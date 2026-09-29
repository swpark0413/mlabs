// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

using namespace Rcpp;

double log_odd_anchor_prior(double anchor, const arma::vec& x,
                            const arma::vec& bd);
arma::vec odd_anchor_bin_probs(const arma::vec& x, const arma::vec& bd,
                               const arma::vec& resid, bool use_residual,
                               double temperature, double uniform_mix);
double sample_odd_anchor(const arma::vec& x, const arma::vec& bd,
                         const arma::vec& bin_prob, double uniform_mix);
double log_odd_anchor_density(double a, const arma::vec& x,
                              const arma::vec& bd,
                              const arma::vec& bin_prob, double uniform_mix);
arma::vec generate_odd_knots(double a, const arma::vec& bd,
                             unsigned int degree);
double odd_knot_log_sides(const arma::vec& knots, unsigned int degree,
                          const arma::vec& bd);

arma::vec B_basis_eval(const arma::vec& x, const unsigned int degree, const arma::vec& knots);
arma::vec tensor_basis(const arma::mat& X, const List& knots_list,
                       const arma::ivec& degrees, const arma::ivec& var_idx,
                       const unsigned int K);
arma::mat build_design_matrix(const arma::mat& X, const List& xi_list,
                              const List& C_list, const List& nu_list,
                              const arma::ivec& K_vec, const unsigned int J);
arma::vec generate_knots_weighted(const arma::vec& x, const arma::vec& bd,
                                  const arma::vec& resid,
                                  const unsigned int degree,
                                  const bool use_residual,
                                  const double temperature,
                                  const double knot_window_frac,
                                  const double odd_anchor_uniform);
double log_knot_proposal_density(const arma::vec& knots, const arma::vec& x,
                                 const arma::vec& bd, const arma::vec& resid,
                                 const unsigned int degree,
                                 const bool use_residual,
                                 const double temperature,
                                 const double knot_window_frac,
                                 const double odd_anchor_uniform);
arma::vec compute_relevance(const arma::mat& X, const arma::vec& resid);
arma::vec compute_relevance_ema(const arma::mat& X, const arma::vec& resid,
                                const arma::vec& rho_cache,
                                const unsigned int iter,
                                const unsigned int update_interval,
                                const double decay);
arma::vec effective_var_probs(const arma::vec& pi_var, const arma::vec& rho,
                              const double boost);
arma::ivec sample_variables(const arma::vec& probs, const unsigned int K,
                            const bool replace);
double pwor_logset_ivec(const arma::ivec& S, const unsigned int K,
                        const arma::vec& pv);
arma::vec knot_anchor_weights(const arma::vec& x, const arma::vec& resid,
                              const double temperature);

inline double log_dnorm(double x, double mean, double sd) {
  return R::dnorm(x, mean, sd, true);
}

inline double g_birth_prob(unsigned int J, unsigned int NB_max, double pb) {
  if (J <= 1)        return 1.0;
  if (J >= NB_max)   return 0.0;
  return pb;
}
inline double g_death_prob(unsigned int J, unsigned int NB_max, double pd) {
  if (J <= 1)        return 0.0;
  if (J >= NB_max)   return 1.0;
  return pd;
}

static double log_marginal_basis(const arma::vec& b,
                                 const arma::vec& r,
                                 double sigma2,
                                 double phi2) {
  double bTb = arma::dot(b, b);
  double bTr = arma::dot(b, r);
  double v2  = 1.0 / (bTb / sigma2 + 1.0 / phi2);
  double m   = v2 * bTr / sigma2;
  return 0.5 * std::log(v2 / phi2) + 0.5 * m * m / v2;
}

static double log_knot_prior(const arma::vec& knots,
                             const unsigned int degree,
                             const arma::vec& bd_k,
                             const arma::vec& x_k) {
  double b1 = bd_k[0] - bd_k[2];
  double b2 = bd_k[1] + bd_k[2];
  double range = b2 - b1;
  if (range <= 0.0) return -std::numeric_limits<double>::infinity();

  if (degree % 2 == 0) {
    unsigned int nk = degree + 2;
    return std::lgamma((double) nk + 1.0) - (double) nk * std::log(range);
  }

  double log_sides = odd_knot_log_sides(knots, degree, bd_k);
  if (!std::isfinite(log_sides)) return log_sides;
  return log_odd_anchor_prior(knots[(degree + 1) / 2], x_k, bd_k)
    + log_sides;
}

// [[Rcpp::export]]
List backfitting_update(const arma::vec& y, const arma::mat& X,
                        arma::vec beta, arma::mat Bmat,
                        List C_list, List xi_list, List nu_list,
                        arma::ivec K_vec, unsigned int J,
                        const double beta0,
                        const double sigmay, const double sigmab,
                        const arma::mat& bd, const unsigned int n,
                        const arma::ivec& allowed_degrees,
                        const double xi_step_size = 0.05,
                        const double knot_temperature = 1.0,
                        const unsigned int n_refresh = 0,
                        const double knot_window_frac = 0.001,
                        const double odd_anchor_uniform = 0.05) {
  double sigmay2 = sigmay * sigmay;
  double sigmab2 = sigmab * sigmab;
  unsigned int nS = allowed_degrees.n_elem;
  int d_max = (nS > 0) ? (int) allowed_degrees[nS - 1] : 2;

  arma::vec fit = beta0 + Bmat * beta;

  for (unsigned int j = 0; j < J; ++j) {

    arma::vec b_j = Bmat.col(j);
    arma::vec partial_fit = fit - beta[j] * b_j;
    arma::vec resid_j = y - partial_fit;

    double bTb = arma::dot(b_j, b_j);
    double bTr = arma::dot(b_j, resid_j);
    double v2 = 1.0 / (bTb / sigmay2 + 1.0 / sigmab2);
    double m = v2 * (bTr / sigmay2);
    if (std::isfinite(m) && std::isfinite(v2) && v2 > 0.0) {
      beta[j] = R::rnorm(m, std::sqrt(v2));
    } else {
      double safe_sd = std::isfinite(sigmab) ? sigmab : 1.0;
      beta[j] = R::rnorm(0.0, safe_sd);
    }
    if (!std::isfinite(beta[j]) ||
        std::abs(beta[j]) > 100.0 * sigmab) {
      beta[j] = R::rnorm(0.0, sigmab);
    }

    fit = partial_fit + beta[j] * b_j;

    List xi_j = xi_list[j];
    arma::ivec C_j = as<arma::ivec>(C_list[j]);
    arma::ivec nu_j = as<arma::ivec>(nu_list[j]);
    unsigned int K_j = (unsigned int) K_vec[j];

    for (unsigned int k = 0; k < K_j; ++k) {
      arma::vec knots = as<arma::vec>(xi_j[k]);
      arma::vec x_k = X.col(nu_j[k]);
      arma::vec bd_k = bd.col(nu_j[k]);
      double bd1 = bd_k[0] - bd_k[2];
      double bd2 = bd_k[1] + bd_k[2];
      double range = bd2 - bd1;

      int c_k = (int) C_j[k];
      int central_idx = (c_k % 2 == 1) ? (c_k + 1) / 2 : -1;

      arma::vec knots_prop = knots;
      for (unsigned int ii = 0; ii < knots.n_elem; ++ii) {
        if ((int) ii == central_idx) continue;
        knots_prop[ii] += R::rnorm(0.0, xi_step_size * range);
      }
      knots_prop = arma::sort(knots_prop);

      bool ok = true;

      if (knots_prop.min() < bd1 || knots_prop.max() > bd2) ok = false;

      if (ok && central_idx >= 0 &&
          knots_prop[central_idx] != knots[central_idx]) {
        ok = false;
      }
      if (ok && central_idx >= 0 &&
          !std::isfinite(odd_knot_log_sides(knots_prop, c_k, bd_k))) {
          ok = false;
      }

      if (ok) {
        List xi_j_prop = clone(xi_j);
        xi_j_prop[k] = knots_prop;
        arma::vec b_j_prop = tensor_basis(X, xi_j_prop, C_j, nu_j, K_j);

        double bTb_prop = arma::dot(b_j_prop, b_j_prop);
        double bTr_prop = arma::dot(b_j_prop, resid_j);
        double log_lik_old = (beta[j] / sigmay2) * bTr
        - 0.5 * beta[j] * beta[j] * bTb / sigmay2;
        double log_lik_new = (beta[j] / sigmay2) * bTr_prop
        - 0.5 * beta[j] * beta[j] * bTb_prop / sigmay2;
        double logA = log_lik_new - log_lik_old;

        if (std::log(R::runif(0.0, 1.0)) < logA) {
          xi_j = xi_j_prop;
          b_j = b_j_prop;
          bTb = bTb_prop;
          bTr = bTr_prop;
        }
      }

      if (central_idx >= 0) {
        arma::vec knots_cur = as<arma::vec>(xi_j[k]);
        arma::vec bin_prob = odd_anchor_bin_probs(
          x_k, bd_k, resid_j, true, knot_temperature, odd_anchor_uniform);
        double anchor_prop = sample_odd_anchor(
          x_k, bd_k, bin_prob, odd_anchor_uniform);
        arma::vec knots_prop2 = generate_odd_knots(anchor_prop, bd_k, c_k);
        double log_g_cur = log_odd_anchor_density(
          knots_cur[central_idx], x_k, bd_k, bin_prob, odd_anchor_uniform);
        double log_g_prop = log_odd_anchor_density(
          knots_prop2[central_idx], x_k, bd_k, bin_prob, odd_anchor_uniform);
        if (std::isfinite(log_g_cur) && std::isfinite(log_g_prop) &&
            std::isfinite(odd_knot_log_sides(knots_prop2, c_k, bd_k))) {
          List xi_j_prop2 = clone(xi_j);
          xi_j_prop2[k] = knots_prop2;
          arma::vec b_j_prop2 = tensor_basis(X, xi_j_prop2, C_j, nu_j, K_j);

          double bTb_p2 = arma::dot(b_j_prop2, b_j_prop2);
          double bTr_p2 = arma::dot(b_j_prop2, resid_j);
          double log_lik_old = (beta[j] / sigmay2) * bTr
          - 0.5 * beta[j] * beta[j] * bTb / sigmay2;
          double log_lik_new = (beta[j] / sigmay2) * bTr_p2
          - 0.5 * beta[j] * beta[j] * bTb_p2 / sigmay2;

          double log_h_cur = log_odd_anchor_prior(
            knots_cur[central_idx], x_k, bd_k);
          double log_h_prop = log_odd_anchor_prior(
            knots_prop2[central_idx], x_k, bd_k);
          double log_anchor_ratio = (log_h_prop - log_h_cur)
            + (log_g_cur - log_g_prop);

          double logA2 = (log_lik_new - log_lik_old) + log_anchor_ratio;

          if (std::log(R::runif(0.0, 1.0)) < logA2) {
            xi_j = xi_j_prop2;
            b_j = b_j_prop2;
            bTb = bTb_p2;
            bTr = bTr_p2;
          }
        }
      }

    }

    xi_list[j] = xi_j;
    C_list[j] = C_j;
    Bmat.col(j) = b_j;
    fit = partial_fit + beta[j] * b_j;
  }

  if (n_refresh > 0 && J > 0) {
    unsigned int n_do = std::min(n_refresh, J);

    arma::uvec perm = arma::randperm(J);

    for (unsigned int ii = 0; ii < n_do; ++ii) {
      unsigned int j = perm[ii];

      arma::vec b_j = Bmat.col(j);
      arma::vec partial_fit_j = fit - beta[j] * b_j;
      arma::vec resid_j = y - partial_fit_j;

      List xi_j = xi_list[j];
      arma::ivec C_j = as<arma::ivec>(C_list[j]);
      arma::ivec nu_j = as<arma::ivec>(nu_list[j]);
      unsigned int K_j = (unsigned int) K_vec[j];

      double log_m_cur = log_marginal_basis(b_j, resid_j, sigmay2, sigmab2);

      double log_q_xi_cur = 0.0;
      double log_prior_xi_cur = 0.0;
      for (unsigned int k = 0; k < K_j; ++k) {
        arma::vec x_k = X.col(nu_j[k]);
        arma::vec bd_k = bd.col(nu_j[k]);
        arma::vec knots_k = as<arma::vec>(xi_j[k]);

        log_q_xi_cur += log_knot_proposal_density(
          knots_k, x_k, bd_k, resid_j, (unsigned int) C_j[k],
                                                         true, knot_temperature, knot_window_frac, odd_anchor_uniform);

        log_prior_xi_cur += log_knot_prior(
          knots_k, (unsigned int) C_j[k], bd_k, x_k);
      }

      arma::ivec C_prop(K_j);
      for (unsigned int k = 0; k < K_j; ++k) {
        unsigned int idx = (unsigned int) std::floor(
          R::runif(0.0, (double) nS));
        if (idx >= nS) idx = nS - 1;
        C_prop[k] = allowed_degrees[idx];
      }

      List xi_prop(K_j);
      double log_q_xi_prop = 0.0;
      double log_prior_xi_prop = 0.0;
      for (unsigned int k = 0; k < K_j; ++k) {
        arma::vec x_k = X.col(nu_j[k]);
        arma::vec bd_k = bd.col(nu_j[k]);

        arma::vec knots_prop = generate_knots_weighted(
          x_k, bd_k, resid_j, (unsigned int) C_prop[k],
                                                   true, knot_temperature, knot_window_frac, odd_anchor_uniform);
        xi_prop[k] = knots_prop;

        log_q_xi_prop += log_knot_proposal_density(
          knots_prop, x_k, bd_k, resid_j, (unsigned int) C_prop[k],
                                                               true, knot_temperature, knot_window_frac, odd_anchor_uniform);

        log_prior_xi_prop += log_knot_prior(
          knots_prop, (unsigned int) C_prop[k], bd_k, x_k);
      }

      arma::vec b_prop = tensor_basis(X, xi_prop, C_prop, nu_j, K_j);

      double log_m_prop = log_marginal_basis(b_prop, resid_j, sigmay2, sigmab2);

      double logA_refresh =
        (log_m_prop - log_m_cur)
        + (log_prior_xi_prop - log_prior_xi_cur)
        + (log_q_xi_cur - log_q_xi_prop);

        if (!std::isfinite(logA_refresh))
          logA_refresh = -std::numeric_limits<double>::infinity();

        if (std::log(R::runif(0.0, 1.0)) < logA_refresh) {

          C_list[j] = C_prop;
          xi_list[j] = xi_prop;

          Bmat.col(j) = b_prop;

          double bTb_new = arma::dot(b_prop, b_prop);
          double bTr_new = arma::dot(b_prop, resid_j);
          double v2_new = 1.0 / (bTb_new / sigmay2 + 1.0 / sigmab2);
          double m_new  = v2_new * (bTr_new / sigmay2);
          if (std::isfinite(m_new) && std::isfinite(v2_new) && v2_new > 0.0) {
            beta[j] = R::rnorm(m_new, std::sqrt(v2_new));
          } else {
            beta[j] = R::rnorm(0.0, sigmab);
          }
          if (!std::isfinite(beta[j])) beta[j] = 0.0;

          fit = partial_fit_j + beta[j] * b_prop;
        }

    }
  }

  for (unsigned int j = 0; j < J; ++j) {
    arma::ivec C_j_check = as<arma::ivec>(C_list[j]);
    List xi_j_check = xi_list[j];
    arma::ivec nu_j_check = as<arma::ivec>(nu_list[j]);
    unsigned int K_j_check = (unsigned int) K_vec[j];
    bool fixed = false;

    for (unsigned int k = 0; k < K_j_check; ++k) {
      arma::vec knots_k = as<arma::vec>(xi_j_check[k]);
      unsigned int expected = (unsigned int) C_j_check[k] + 2;
      if (knots_k.n_elem != expected) {
        arma::vec bd_k = bd.col(nu_j_check[k]);
        double lo = bd_k[0] - bd_k[2];
        double hi = bd_k[1] + bd_k[2];
        arma::vec new_knots(expected);
        for (unsigned int kk = 0; kk < expected; ++kk)
          new_knots[kk] = R::runif(lo, hi);
        new_knots = arma::sort(new_knots);
        xi_j_check[k] = new_knots;
        fixed = true;
      }
    }
    if (fixed) {
      xi_list[j] = xi_j_check;

      arma::vec b_fix = tensor_basis(X, xi_j_check, C_j_check, nu_j_check, K_j_check);
      Bmat.col(j) = b_fix;
      arma::vec fit_fix = beta0 + Bmat * beta - beta[j] * b_fix;
      arma::vec resid_fix = y - fit_fix;
      double bTb_fix = arma::dot(b_fix, b_fix);
      double bTr_fix = arma::dot(b_fix, resid_fix);
      double v2_fix = 1.0 / (bTb_fix / sigmay2 + 1.0 / sigmab2);
      double m_fix  = v2_fix * (bTr_fix / sigmay2);
      if (std::isfinite(m_fix) && v2_fix > 0.0)
        beta[j] = R::rnorm(m_fix, std::sqrt(v2_fix));
      else
        beta[j] = 0.0;
    }
  }

  return List::create(
    Named("beta") = beta,
    Named("Bmat") = Bmat,
    Named("C") = C_list,
    Named("xi") = xi_list);
}

// [[Rcpp::export]]
double update_sigmay(const arma::vec& y, const arma::vec& beta,
                     const arma::mat& Bmat, const double beta0,
                     const double sigr, const double sigR) {
  unsigned int n = y.n_elem;
  arma::vec fit = beta0 + Bmat * beta;
  arma::vec resid = y - fit;
  double rss = arma::dot(resid, resid);
  double new_r = sigr + (double) n;
  double new_R = sigr * sigR + rss;
  double inv_sigma2 = R::rgamma(0.5 * new_r, 2.0 / new_R);
  return 1.0 / std::sqrt(inv_sigma2);
}

// [[Rcpp::export]]
List birth_step_collapsed(const arma::vec& y, const arma::mat& X,
                          arma::vec beta, arma::mat Bmat,
                          List C_list, List xi_list, List nu_list,
                          arma::ivec K_vec, unsigned int J,
                          const double M, const double beta0,
                          const double sigmay, const double sigmab,
                          const arma::mat& bd,
                          const unsigned int maxInter,
                          const unsigned int n, const unsigned int p,
                          const arma::vec& jump_prob,
                          const arma::vec& pi_var,
                          const arma::ivec& allowed_degrees,
                          const double knot_temperature = 1.0,
                          const arma::vec& rho_cache = arma::vec(),
                          const unsigned int iter = 0,
                          const unsigned int rho_update_interval = 10,
                          const double boost = 5.0,
                          const double ref_sigma = 1.0,
                          const unsigned int NB_max = 100,
                          const double lambda_J = 0.0,
                          const arma::vec& k_weights = arma::vec(),
                          const double knot_window_frac = 0.001,
                          const double odd_anchor_uniform = 0.05) {
  double sigmay2 = sigmay * sigmay;
  double sigmab2 = sigmab * sigmab;
  unsigned int nS = allowed_degrees.n_elem;
  int d_max = (nS > 0) ? (int) allowed_degrees[nS - 1] : 2;

  arma::vec fit_current = beta0 + Bmat * beta;
  arma::vec resid = y - fit_current;

  unsigned int K_new;
  if (k_weights.n_elem == maxInter && arma::sum(k_weights) > 0.0) {
    double u = R::runif(0.0, 1.0) * arma::sum(k_weights), c = 0.0; K_new = maxInter;
    for (unsigned int kk = 0; kk < maxInter; ++kk) { c += k_weights[kk]; if (u <= c) { K_new = kk + 1; break; } }
  } else {
    K_new = (unsigned int) std::floor(R::runif(0.0, (double) maxInter)) + 1;
  }
  if (K_new > maxInter) K_new = maxInter;

  arma::vec rho = compute_relevance_ema(X, resid, rho_cache,
                                        iter, rho_update_interval, 0.9);
  arma::vec tilde_pi = effective_var_probs(pi_var, rho, boost);
  arma::ivec cand_nu = sample_variables(tilde_pi, K_new, false);

  arma::ivec cand_C(K_new);
  for (unsigned int k = 0; k < K_new; ++k) {
    unsigned int idx = (unsigned int) std::floor(R::runif(0.0, (double) nS));
    if (idx >= nS) idx = nS - 1;
    cand_C[k] = allowed_degrees[idx];
  }

  List cand_xi(K_new);
  double log_q_xi = 0.0;
  double log_prior_xi = 0.0;
  for (unsigned int k = 0; k < K_new; ++k) {
    arma::vec x_k = X.col(cand_nu[k]);
    arma::vec bd_k = bd.col(cand_nu[k]);
    arma::vec knots_k = generate_knots_weighted(x_k, bd_k, resid, cand_C[k],
                                                true, knot_temperature, knot_window_frac, odd_anchor_uniform);
    cand_xi[k] = knots_k;

    log_q_xi += log_knot_proposal_density(knots_k, x_k, bd_k, resid, cand_C[k],
                                          true, knot_temperature, knot_window_frac, odd_anchor_uniform);

    log_prior_xi += log_knot_prior(knots_k, cand_C[k], bd_k, x_k);
  }

  arma::vec b_new = tensor_basis(X, cand_xi, cand_C, cand_nu, K_new);

  double bTb = arma::dot(b_new, b_new);
  double bTr = arma::dot(b_new, resid);

  if (!std::isfinite(bTb) || !std::isfinite(bTr)) {
    return List::create(
      Named("beta") = beta, Named("Bmat") = Bmat,
      Named("C") = C_list, Named("xi") = xi_list,
      Named("nu") = nu_list, Named("K") = K_vec,
      Named("J") = J, Named("accepted") = false,
      Named("rho") = rho);
  }

  double v2 = 1.0 / (bTb / sigmay2 + 1.0 / sigmab2);
  double m  = v2 * (bTr / sigmay2);
  double log_marg = 0.5 * std::log(v2 / sigmab2) + 0.5 * m * m / v2;

  double log_prior_J = std::log(M) - lambda_J - std::log((double)(J + 1));

  double log_prior_nu    = pwor_logset_ivec(cand_nu, K_new, pi_var);
  double log_proposal_nu = pwor_logset_ivec(cand_nu, K_new, tilde_pi);

  double gb_fwd = g_birth_prob(J,     NB_max, jump_prob[0]);
  double gd_rev = g_death_prob(J + 1, NB_max, jump_prob[1]);
  double log_move_ratio = std::log(gd_rev) - std::log(gb_fwd);

  double logA = log_marg
  + log_prior_J
  + (log_prior_nu - log_proposal_nu)
  + (log_prior_xi - log_q_xi)
  + log_move_ratio;

  if (!std::isfinite(logA)) logA = -std::numeric_limits<double>::infinity();

  double u = R::runif(0.0, 1.0);
  bool accept = (std::log(u) < logA);

  if (accept) {

    double beta_new;
    if (std::isfinite(m) && std::isfinite(v2) && v2 > 0.0) {
      beta_new = R::rnorm(m, std::sqrt(v2));
    } else {
      beta_new = R::rnorm(0.0, sigmab);
    }
    if (!std::isfinite(beta_new)) beta_new = 0.0;

    arma::vec new_beta(beta.n_elem + 1);
    new_beta.head(beta.n_elem) = beta;
    new_beta[beta.n_elem] = beta_new;
    beta = new_beta;

    arma::mat new_Bmat(Bmat.n_rows, Bmat.n_cols + 1);
    new_Bmat.head_cols(Bmat.n_cols) = Bmat;
    new_Bmat.col(Bmat.n_cols) = b_new;
    Bmat = new_Bmat;

    arma::ivec new_K(K_vec.n_elem + 1);
    new_K.head(K_vec.n_elem) = K_vec;
    new_K[K_vec.n_elem] = (int) K_new;
    K_vec = new_K;

    xi_list.push_back(cand_xi);
    C_list.push_back(cand_C);
    nu_list.push_back(cand_nu);
    J = J + 1;
  }

  return List::create(
    Named("beta") = beta, Named("Bmat") = Bmat,
    Named("C") = C_list, Named("xi") = xi_list,
    Named("nu") = nu_list, Named("K") = K_vec,
    Named("J") = J, Named("accepted") = accept,
    Named("rho") = rho);
}

// [[Rcpp::export]]
List death_step_collapsed(const arma::vec& y, const arma::mat& X,
                          arma::vec beta, arma::mat Bmat,
                          List C_list, List xi_list, List nu_list,
                          arma::ivec K_vec, unsigned int J,
                          const double M, const double beta0,
                          const double sigmay, const double sigmab,
                          const arma::mat& bd,
                          const unsigned int n, const unsigned int p,
                          const arma::vec& jump_prob,
                          const arma::vec& pi_var,
                          const arma::ivec& allowed_degrees,
                          const double knot_temperature = 1.0,
                          const arma::vec& rho_cache = arma::vec(),
                          const unsigned int iter = 0,
                          const unsigned int rho_update_interval = 10,
                          const double boost = 5.0,
                          const double ref_sigma = 1.0,
                          const unsigned int NB_max = 100,
                          const double lambda_J = 0.0,
                          const double knot_window_frac = 0.001,
                          const double odd_anchor_uniform = 0.05) {
  double sigmay2 = sigmay * sigmay;
  double sigmab2 = sigmab * sigmab;
  unsigned int nS_death = allowed_degrees.n_elem;
  int d_max_death = (nS_death > 0) ? (int) allowed_degrees[nS_death - 1] : 2;

  unsigned int j_star = (unsigned int) std::floor(R::runif(0.0, (double) J));
  if (j_star >= J) j_star = J - 1;

  double beta_rm = beta[j_star];
  arma::vec b_rm = Bmat.col(j_star);
  List xi_rm = xi_list[j_star];
  arma::ivec C_rm = as<arma::ivec>(C_list[j_star]);
  arma::ivec nu_rm = as<arma::ivec>(nu_list[j_star]);
  unsigned int K_rm = (unsigned int) K_vec[j_star];

  arma::vec partial_fit = beta0 + Bmat * beta - beta_rm * b_rm;
  arma::vec resid_partial = y - partial_fit;

  double bTb = arma::dot(b_rm, b_rm);
  double bTr = arma::dot(b_rm, resid_partial);
  double v2 = 1.0 / (bTb / sigmay2 + 1.0 / sigmab2);
  double m  = v2 * (bTr / sigmay2);
  double log_marg = 0.5 * std::log(v2 / sigmab2) + 0.5 * m * m / v2;

  arma::vec rho_rev = compute_relevance_ema(X, resid_partial, rho_cache,
                                            iter, rho_update_interval, 0.9);
  arma::vec tilde_pi_rev = effective_var_probs(pi_var, rho_rev, boost);

  double log_prior_nu    = pwor_logset_ivec(nu_rm, K_rm, pi_var);
  double log_proposal_nu = pwor_logset_ivec(nu_rm, K_rm, tilde_pi_rev);

  double log_prior_xi = 0.0, log_q_xi = 0.0;
  for (unsigned int k = 0; k < K_rm; ++k) {
    arma::vec x_k = X.col(nu_rm[k]);
    arma::vec bd_k = bd.col(nu_rm[k]);
    arma::vec knots_k = as<arma::vec>(xi_rm[k]);

    log_q_xi += log_knot_proposal_density(knots_k, x_k, bd_k, resid_partial,
                                          C_rm[k], true, knot_temperature, knot_window_frac, odd_anchor_uniform);

    log_prior_xi += log_knot_prior(knots_k, C_rm[k], bd_k, x_k);
  }

  double log_prior_J = std::log((double) J) - std::log(M) + lambda_J;
  double gd_fwd = g_death_prob(J,     NB_max, jump_prob[1]);
  double gb_rev = g_birth_prob(J - 1, NB_max, jump_prob[0]);
  double log_move_ratio = std::log(gb_rev) - std::log(gd_fwd);

  double logA = -log_marg
  + log_prior_J
  - (log_prior_nu - log_proposal_nu)
  - (log_prior_xi - log_q_xi)
  + log_move_ratio;

  if (!std::isfinite(logA)) logA = -std::numeric_limits<double>::infinity();

  double u = R::runif(0.0, 1.0);
  bool accept = (std::log(u) < logA);

  if (accept) {
    arma::vec new_beta(beta.n_elem - 1);
    unsigned int idx = 0;
    for (unsigned int i = 0; i < beta.n_elem; ++i) {
      if (i != j_star) { new_beta[idx++] = beta[i]; }
    }
    beta = new_beta;

    arma::mat new_Bmat(Bmat.n_rows, Bmat.n_cols - 1);
    idx = 0;
    for (unsigned int i = 0; i < Bmat.n_cols; ++i) {
      if (i != j_star) { new_Bmat.col(idx++) = Bmat.col(i); }
    }
    Bmat = new_Bmat;

    arma::ivec new_K(K_vec.n_elem - 1);
    idx = 0;
    for (unsigned int i = 0; i < K_vec.n_elem; ++i) {
      if (i != j_star) { new_K[idx++] = K_vec[i]; }
    }
    K_vec = new_K;

    xi_list.erase(j_star);
    C_list.erase(j_star);
    nu_list.erase(j_star);
    J = J - 1;
  }

  return List::create(
    Named("beta") = beta, Named("Bmat") = Bmat,
    Named("C") = C_list, Named("xi") = xi_list,
    Named("nu") = nu_list, Named("K") = K_vec,
    Named("J") = J, Named("accepted") = accept,
    Named("rho") = rho_rev);
}
