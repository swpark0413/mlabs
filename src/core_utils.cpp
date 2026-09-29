// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

using namespace Rcpp;

static const double ODD_ANCHOR_PRIOR_UNIFORM = 0.05;

// [[Rcpp::export]]
arma::vec B_basis_eval(const arma::vec& x, const unsigned int degree,
                       const arma::vec& knots) {
  unsigned int n = x.n_elem;
  arma::vec res(n, arma::fill::zeros);

  if (degree == 0) {
    arma::uvec cond = arma::find((knots[0] <= x) && (x < knots[1]));
    res.elem(cond).ones();
    return res;
  }

  if (degree == 1) {
    double k1 = knots[0], k2 = knots[1], k3 = knots[2];
    arma::uvec cond1 = arma::find((k1 <= x) && (x < k2));
    arma::uvec cond2 = arma::find((k2 <= x) && (x < k3));
    arma::vec res1(n, arma::fill::zeros), res2(n, arma::fill::zeros);
    res1.elem(cond1).ones();
    res2.elem(cond2).ones();
    res = (x - k1) / (k2 - k1) % res1 + (k3 - x) / (k3 - k2) % res2;
    return res;
  }

  if (degree == 2) {
    double k1 = knots[0], k2 = knots[1], k3 = knots[2], k4 = knots[3];
    double k21 = k2 - k1, k31 = k3 - k1, k32 = k3 - k2, k42 = k4 - k2, k43 = k4 - k3;

    arma::uvec cond1 = arma::find((k1 <= x) && (x < k2));
    arma::uvec cond2 = arma::find((k2 <= x) && (x < k3));
    arma::uvec cond3 = arma::find((k3 <= x) && (x < k4));
    arma::vec r1(n, arma::fill::zeros), r2(n, arma::fill::zeros), r3(n, arma::fill::zeros);
    r1.elem(cond1).ones();
    r2.elem(cond2).ones();
    r3.elem(cond3).ones();

    res = ((x - k1) % (x - k1) / (k31 * k21)) % r1
    + (((x - k1) % (k3 - x)) / (k31 * k32) + ((k4 - x) % (x - k2)) / (k42 * k32)) % r2
    + (((k4 - x) % (k4 - x)) / (k42 * k43)) % r3;
    return res;
  }

  double k1 = knots[0], k2 = knots[1], k3 = knots[2], k4 = knots[3], k5 = knots[4];
  double k21 = k2 - k1, k31 = k3 - k1, k41 = k4 - k1;
  double k32 = k3 - k2, k42 = k4 - k2, k52 = k5 - k2;
  double k43 = k4 - k3, k53 = k5 - k3, k54 = k5 - k4;

  arma::uvec cond1 = arma::find((k1 <= x) && (x < k2));
  arma::uvec cond2 = arma::find((k2 <= x) && (x < k3));
  arma::uvec cond3 = arma::find((k3 <= x) && (x < k4));
  arma::uvec cond4 = arma::find((k4 <= x) && (x < k5));
  arma::vec r1(n, arma::fill::zeros), r2(n, arma::fill::zeros),
  r3(n, arma::fill::zeros), r4(n, arma::fill::zeros);
  r1.elem(cond1).ones();
  r2.elem(cond2).ones();
  r3.elem(cond3).ones();
  r4.elem(cond4).ones();

  res = ((x - k1) % (x - k1) % (x - k1)) / (k41 * k31 * k21) % r1
  + ((x - k1) % (x - k1) % (k3 - x) / (k41 * k31 * k32)
       + (x - k1) % (k4 - x) % (x - k2) / (k41 * k42 * k32)
       + (k5 - x) % (x - k2) % (x - k2) / (k52 * k42 * k32)) % r2
       + ((x - k1) % (k4 - x) % (k4 - x) / (k41 * k42 * k43)
       + (x - k2) % (k5 - x) % (k4 - x) / (k52 * k43 * k42)
       + (k5 - x) % (k5 - x) % (x - k3) / (k52 * k53 * k43)) % r3
       + ((k5 - x) % (k5 - x) % (k5 - x)) / (k52 * k53 * k54) % r4;
       return res;
}

// [[Rcpp::export]]
arma::vec tensor_basis(const arma::mat& X, const List& knots_list,
                       const arma::ivec& degrees, const arma::ivec& var_idx,
                       const unsigned int K) {
  unsigned int n = X.n_rows;
  arma::vec res = arma::ones<arma::vec>(n);
  for (unsigned int k = 0; k < K; ++k) {
    arma::vec kn = as<arma::vec>(knots_list[k]);
    res = res % B_basis_eval(X.col(var_idx[k]), degrees[k], kn);
  }
  return res;
}

// [[Rcpp::export]]
arma::mat build_design_matrix(const arma::mat& X, const List& xi_list,
                              const List& C_list, const List& nu_list,
                              const arma::ivec& K_vec, const unsigned int J) {
  unsigned int n = X.n_rows;
  arma::mat B(n, J);
  for (unsigned int j = 0; j < J; ++j) {
    List kn_j = xi_list[j];
    arma::ivec Cj = as<arma::ivec>(C_list[j]);
    arma::ivec nuj = as<arma::ivec>(nu_list[j]);
    B.col(j) = tensor_basis(X, kn_j, Cj, nuj, K_vec[j]);
  }
  return B;
}

// [[Rcpp::export]]
arma::mat compute_boundary(const arma::mat& X, const double prop) {
  unsigned int p = X.n_cols;
  arma::mat bd(3, p);
  for (unsigned int j = 0; j < p; ++j) {
    double mn = X.col(j).min();
    double mx = X.col(j).max();
    bd(0, j) = mn;
    bd(1, j) = mx;
    bd(2, j) = (mx - mn) * prop;
  }
  return bd;
}

inline arma::vec sort_two(arma::vec v) {
  if (v[0] > v[1]) std::swap(v[0], v[1]);
  return v;
}

arma::vec knot_anchor_weights(const arma::vec& x,
                              const arma::vec& resid,
                              const double T) {
  unsigned int n = x.n_elem;

  if (n == 0) return arma::vec();

  arma::vec abs_r = arma::abs(resid);
  double mean_abs_r = arma::mean(abs_r);

  double eps_floor = 1e-3 * mean_abs_r + 1e-12;

  unsigned int n_bins = (unsigned int) std::floor(std::sqrt((double) n));
  if (n_bins < 5)  n_bins = 5;
  if (n_bins > 20) n_bins = 20;
  if (n_bins > n)  n_bins = n;

  double xmin = x.min();
  double xmax = x.max();
  double range = xmax - xmin;

  if (range < 1e-12) {
    arma::vec w(n, arma::fill::value(1.0 / (double) n));
    return w;
  }

  arma::uvec bin_of(n);
  arma::vec  bin_sum(n_bins, arma::fill::zeros);
  arma::uvec bin_cnt(n_bins, arma::fill::zeros);
  for (unsigned int i = 0; i < n; ++i) {
    double u = (x[i] - xmin) / range;
    unsigned int b = (unsigned int) std::floor(u * (double) n_bins);
    if (b >= n_bins) b = n_bins - 1;
    bin_of[i] = b;
    bin_sum[b] += abs_r[i];
    bin_cnt[b] += 1;
  }

  arma::vec bin_mean(n_bins, arma::fill::zeros);
  for (unsigned int b = 0; b < n_bins; ++b) {

    if (bin_cnt[b] > 0) bin_mean[b] = bin_sum[b] / (double) bin_cnt[b];
  }

  arma::vec w(n);
  for (unsigned int i = 0; i < n; ++i) {
    w[i] = std::pow(bin_mean[bin_of[i]] + eps_floor, 1.0 / T);
  }

  double s = arma::sum(w);
  if (s < 1e-300) {
    w.fill(1.0 / (double) n);
  } else {
    w /= s;
  }
  return w;
}

static unsigned int odd_anchor_bin_index(double a, double xmin,
                                         double xmax, unsigned int nb) {
  unsigned int b = (unsigned int) std::floor(
    ((a - xmin) / (xmax - xmin)) * (double) nb);
  return std::min(b, nb - 1);
}

arma::vec odd_anchor_bin_probs(const arma::vec& x, const arma::vec& bd,
                               const arma::vec& resid, bool use_residual,
                               double temperature, double uniform_mix) {
  if (bd.n_elem != 3 || !bd.is_finite())
    Rcpp::stop("Odd-degree knots require a finite three-element boundary.");
  double lo = bd[0] - bd[2];
  double hi = bd[1] + bd[2];
  if (!std::isfinite(hi - lo) || !(hi > lo))
    Rcpp::stop("Odd-degree knots require a positive finite knot-domain width.");
  if (!std::isfinite(uniform_mix) || uniform_mix <= 0.0 || uniform_mix > 1.0)
    Rcpp::stop("odd_anchor_uniform must be in (0, 1].");

  arma::vec bin_prob;
  if (!use_residual || uniform_mix == 1.0) return bin_prob;
  if (x.n_elem == 0 || !x.is_finite() || resid.n_elem != x.n_elem ||
      !resid.is_finite())
      Rcpp::stop("Odd-degree anchor proposals require finite, equally sized x and residuals.");
  if (!std::isfinite(temperature))
    Rcpp::stop("knot temperature must be finite.");
  double xmin = x.min();
  double xmax = x.max();
  if (xmin < lo || xmax > hi)
    Rcpp::stop("Observed predictors must lie inside the knot domain.");
  if (xmax - xmin < 1e-12) return bin_prob;

  unsigned int nb = (unsigned int) std::floor(std::sqrt((double) x.n_elem));
  nb = std::min((unsigned int) x.n_elem, std::max(5u, std::min(20u, nb)));
  bin_prob.zeros(nb);
  arma::vec w = knot_anchor_weights(x, resid, std::max(1.0, temperature));
  for (unsigned int i = 0; i < x.n_elem; ++i)
    bin_prob[odd_anchor_bin_index(x[i], xmin, xmax, nb)] += w[i];
  bin_prob /= arma::sum(bin_prob);
  return bin_prob;
}

double sample_odd_anchor(const arma::vec& x, const arma::vec& bd,
                         const arma::vec& bin_prob, double uniform_mix) {
  double lo = bd[0] - bd[2];
  double hi = bd[1] + bd[2];
  if (bin_prob.n_elem == 0 || R::runif(0.0, 1.0) < uniform_mix)
    return R::runif(lo, hi);

  double xmin = x.min();
  double xmax = x.max();
  unsigned int nb = bin_prob.n_elem;
  double width = (xmax - xmin) / (double) nb;
  double u = R::runif(0.0, 1.0), cum = 0.0;
  unsigned int b = nb - 1;
  for (unsigned int j = 0; j < nb; ++j) {
    cum += bin_prob[j];
    if (u < cum) { b = j; break; }
  }
  double left = xmin + b * width;
  double right = (b + 1 == nb) ? xmax : xmin + (b + 1) * width;
  return R::runif(left, right);
}

double log_odd_anchor_density(double a, const arma::vec& x,
                              const arma::vec& bd,
                              const arma::vec& bin_prob,
                              double uniform_mix) {
  double lo = bd[0] - bd[2];
  double hi = bd[1] + bd[2];
  if (!std::isfinite(a) || a <= lo || a >= hi)
    return -std::numeric_limits<double>::infinity();
  if (bin_prob.n_elem == 0) return -std::log(hi - lo);

  double density = uniform_mix / (hi - lo);
  double xmin = x.min();
  double xmax = x.max();
  if (a >= xmin && a <= xmax) {
    unsigned int nb = bin_prob.n_elem;
    double width = (xmax - xmin) / (double) nb;
    unsigned int b = odd_anchor_bin_index(a, xmin, xmax, nb);
    density += (1.0 - uniform_mix) * bin_prob[b] / width;
  }
  return std::log(density);
}

double log_odd_anchor_prior(double anchor, const arma::vec& x,
                            const arma::vec& bd) {
  arma::vec zero_resid(x.n_elem, arma::fill::zeros);
  arma::vec prior_bin_prob = odd_anchor_bin_probs(
    x, bd, zero_resid, true, 1.0, ODD_ANCHOR_PRIOR_UNIFORM);
  return log_odd_anchor_density(
    anchor, x, bd, prior_bin_prob, ODD_ANCHOR_PRIOR_UNIFORM);
}

arma::vec generate_odd_knots(double a, const arma::vec& bd,
                             unsigned int degree) {
  unsigned int k = (degree + 1) / 2;
  double lo = bd[0] - bd[2];
  double hi = bd[1] + bd[2];
  arma::vec knots(degree + 2);
  for (unsigned int i = 0; i < k; ++i)
    knots[i] = R::runif(lo, a);
  knots[k] = a;
  for (unsigned int i = 0; i < k; ++i)
    knots[k + 1 + i] = R::runif(a, hi);
  knots.head(k) = arma::sort(arma::vec(knots.head(k)));
  knots.tail(k) = arma::sort(arma::vec(knots.tail(k)));
  return knots;
}

double odd_knot_log_sides(const arma::vec& knots, unsigned int degree,
                          const arma::vec& bd) {
  double neg_inf = -std::numeric_limits<double>::infinity();
  if (bd.n_elem != 3 || !bd.is_finite() || knots.n_elem != degree + 2 ||
      !knots.is_finite()) return neg_inf;
      double lo = bd[0] - bd[2];
      double hi = bd[1] + bd[2];
      if (!(knots[0] > lo && knots[knots.n_elem - 1] < hi)) return neg_inf;
      for (unsigned int i = 1; i < knots.n_elem; ++i)
        if (!(knots[i] > knots[i - 1])) return neg_inf;
        unsigned int k = (degree + 1) / 2;
        double a = knots[k];
        return 2.0 * std::lgamma((double) k + 1.0)
          - k * std::log(a - lo) - k * std::log(hi - a);
}

// [[Rcpp::export]]
arma::vec generate_knots_weighted(const arma::vec& x, const arma::vec& bd,
                                  const arma::vec& resid,
                                  const unsigned int degree,
                                  const bool use_residual,
                                  const double temperature = 1.0,
                                  const double knot_window_frac = 0.01,
                                  const double odd_anchor_uniform = 0.05) {
  double bd1 = bd[0] - bd[2];
  double bd2 = bd[1] + bd[2];
  unsigned int n = x.n_elem;

  double T = std::max(1.0, temperature);

  if (degree % 2 == 1) {

    double anchor_mix = use_residual ? odd_anchor_uniform
    : ODD_ANCHOR_PRIOR_UNIFORM;
    arma::vec anchor_resid = resid;
    if (!use_residual) anchor_resid.zeros(x.n_elem);
    arma::vec bin_prob = odd_anchor_bin_probs(
      x, bd, anchor_resid, true, T, anchor_mix);
    double anchor = sample_odd_anchor(x, bd, bin_prob, anchor_mix);
    return generate_odd_knots(anchor, bd, degree);
  }

  double pts;
  if (use_residual) {
    arma::vec w = knot_anchor_weights(x, resid, T);
    double u = R::runif(0.0, 1.0);
    double cum = 0.0;
    unsigned int idx = n - 1;
    for (unsigned int i = 0; i < n; ++i) {
      cum += w[i];
      if (u < cum) { idx = i; break; }
    }
    pts = x[idx];
  } else {
    unsigned int idx = (unsigned int) std::floor(R::runif(0.0, (double) n));
    if (idx >= n) idx = n - 1;
    pts = x[idx];
  }

  unsigned int n_c = degree + 2;
  arma::vec knots(n_c);

  if (degree == 0) {
    knots[0] = R::runif(bd1, pts);
    knots[1] = R::runif(pts, bd2);
  } else if (degree == 2) {
    arma::vec k1(2), k2(2);
    k1[0] = R::runif(bd1, pts);
    k1[1] = R::runif(bd1, pts);
    k1 = sort_two(k1);
    k2[0] = R::runif(pts, bd2);
    k2[1] = R::runif(pts, bd2);
    k2 = sort_two(k2);
    knots[0] = k1[0]; knots[1] = k1[1]; knots[2] = k2[0]; knots[3] = k2[1];
  }

  return knots;
}

// [[Rcpp::export]]
double log_knot_proposal_density(const arma::vec& knots, const arma::vec& x,
                                 const arma::vec& bd, const arma::vec& resid,
                                 const unsigned int degree,
                                 const bool use_residual,
                                 const double temperature = 1.0,
                                 const double knot_window_frac = 0.01,
                                 const double odd_anchor_uniform = 0.05) {
  double bd1 = bd[0] - bd[2];
  double bd2 = bd[1] + bd[2];
  unsigned int n = x.n_elem;

  double T = std::max(1.0, temperature);

  if (degree % 2 == 1) {
    double log_sides = odd_knot_log_sides(knots, degree, bd);
    if (!std::isfinite(log_sides)) return log_sides;
    if (!use_residual) {
      return log_odd_anchor_prior(knots[(degree + 1) / 2], x, bd)
      + log_sides;
    }
    arma::vec bin_prob = odd_anchor_bin_probs(
      x, bd, resid, use_residual, T, odd_anchor_uniform);
    return log_odd_anchor_density(
      knots[(degree + 1) / 2], x, bd, bin_prob, odd_anchor_uniform) + log_sides;
  }

  arma::vec w(n);
  if (use_residual) {
    w = knot_anchor_weights(x, resid, T);
  } else {
    w.fill(1.0 / (double) n);
  }

  double log_dens;

  if (degree == 0) {

    double sum_dens = 0.0;
    for (unsigned int i = 0; i < n; ++i) {
      if (x[i] > knots[0] && x[i] < knots[1]) {
        double left_range  = x[i] - bd1;
        double right_range = bd2 - x[i];
        if (left_range > 0.0 && right_range > 0.0) {
          sum_dens += w[i] / (left_range * right_range);
        }
      }
    }
    if (sum_dens <= 0.0) return -std::numeric_limits<double>::infinity();
    log_dens = std::log(sum_dens);

  } else if (degree == 2) {

    double sum_dens = 0.0;
    for (unsigned int i = 0; i < n; ++i) {
      if (x[i] > knots[1] && x[i] < knots[2]) {
        double left_range  = x[i] - bd1;
        double right_range = bd2 - x[i];
        if (left_range > 0.0 && right_range > 0.0) {
          double cond = 4.0 / (left_range * left_range * right_range * right_range);
          sum_dens += w[i] * cond;
        }
      }
    }
    if (sum_dens <= 0.0) return -std::numeric_limits<double>::infinity();
    log_dens = std::log(sum_dens);

  }

  return log_dens;
}

// [[Rcpp::export]]
arma::vec compute_relevance(const arma::mat& X, const arma::vec& resid) {
  unsigned int p = X.n_cols;
  unsigned int n = X.n_rows;
  arma::vec rho(p, arma::fill::zeros);

  double rbar = arma::mean(resid);
  arma::vec r_centered = resid - rbar;
  double ss_total = arma::dot(r_centered, r_centered);

  if (ss_total < 1e-12) {

    return rho;
  }

  unsigned int n_bins = (unsigned int) std::floor(std::sqrt((double) n));
  if (n_bins < 5) n_bins = 5;
  if (n_bins > 20) n_bins = 20;

  for (unsigned int l = 0; l < p; ++l) {
    arma::vec xl = X.col(l);
    double xmin = xl.min();
    double xmax = xl.max();
    double range = xmax - xmin;

    if (range < 1e-12) {
      rho[l] = 0.0;
      continue;
    }

    arma::vec bin_sum(n_bins, arma::fill::zeros);
    arma::uvec bin_count(n_bins, arma::fill::zeros);
    arma::uvec bin_of(n);

    for (unsigned int i = 0; i < n; ++i) {
      double u = (xl[i] - xmin) / range;
      unsigned int b = (unsigned int) std::floor(u * (double) n_bins);
      if (b >= n_bins) b = n_bins - 1;
      bin_of[i] = b;
      bin_sum[b] += resid[i];
      bin_count[b] += 1;
    }

    arma::vec bin_mean(n_bins, arma::fill::zeros);
    for (unsigned int b = 0; b < n_bins; ++b) {
      if (bin_count[b] > 0) bin_mean[b] = bin_sum[b] / (double) bin_count[b];
    }

    double ss_within = 0.0;
    for (unsigned int i = 0; i < n; ++i) {
      double d = resid[i] - bin_mean[bin_of[i]];
      ss_within += d * d;
    }

    double r2 = 1.0 - ss_within / ss_total;
    if (r2 < 0.0) r2 = 0.0;
    if (r2 > 1.0) r2 = 1.0;
    rho[l] = r2;
  }

  return rho;
}

// [[Rcpp::export]]
arma::vec compute_relevance_ema(const arma::mat& X, const arma::vec& resid,
                                const arma::vec& rho_cache,
                                const unsigned int iter,
                                const unsigned int update_interval = 10,
                                const double decay = 0.9) {
  unsigned int p = X.n_cols;
  bool do_full = (rho_cache.n_elem != p) || (iter % update_interval == 0);

  if (!do_full) return rho_cache;

  arma::vec rho_new = compute_relevance(X, resid);
  if (rho_cache.n_elem != p) return rho_new;
  return decay * rho_cache + (1.0 - decay) * rho_new;
}

// [[Rcpp::export]]
arma::vec effective_var_probs(const arma::vec& pi_var, const arma::vec& rho,
                              const double boost = 5.0) {

  if (!rho.is_finite() || !pi_var.is_finite()) {
    return arma::vec(pi_var.n_elem, arma::fill::value(1.0 / (double) pi_var.n_elem));
  }
  arma::vec w = pi_var % (1.0 + boost * rho);
  double s = arma::sum(w);
  if (s < 1e-300 || !std::isfinite(s)) {

    return arma::vec(pi_var.n_elem, arma::fill::value(1.0 / (double) pi_var.n_elem));
  }
  return w / s;
}

// [[Rcpp::export]]
arma::ivec sample_variables(const arma::vec& probs, const unsigned int K,
                            const bool replace) {
  unsigned int p = probs.n_elem;
  arma::vec cand = arma::regspace<arma::vec>(0, p - 1);
  arma::vec sampled = Rcpp::RcppArmadillo::sample(cand, K, replace, probs);
  return arma::conv_to<arma::ivec>::from(sampled);
}

double pwor_rec(const std::vector<unsigned int>& rem,
                double used, const arma::vec& pv) {
  if (rem.empty()) return 1.0;
  double tot = 0.0;
  for (size_t i = 0; i < rem.size(); ++i) {
    unsigned int s = rem[i];
    std::vector<unsigned int> nxt;
    nxt.reserve(rem.size() - 1);
    for (size_t t = 0; t < rem.size(); ++t) if (t != i) nxt.push_back(rem[t]);
    tot += (pv[s] / (1.0 - used)) * pwor_rec(nxt, used + pv[s], pv);
  }
  return tot;
}

double pwor_logset(const std::vector<unsigned int>& S,
                   const arma::vec& pv) {
  if (S.empty()) return 0.0;
  if (S.size() == 1) return std::log(pv[S[0]] + 1e-300);
  return std::log(pwor_rec(S, 0.0, pv) + 1e-300);
}

double pwor_logset_ivec(const arma::ivec& S, const unsigned int K,
                        const arma::vec& pv) {
  std::vector<unsigned int> s;
  s.reserve(K);
  for (unsigned int k = 0; k < K; ++k) s.push_back((unsigned int) S[k]);
  return pwor_logset(s, pv);
}

// [[Rcpp::export]]
arma::vec update_pi_mh_wor(const arma::vec& pi_cur, const List& nu_list,
                           const arma::ivec& K_vec, const unsigned int J,
                           const unsigned int p, const double alpha) {

  arma::vec counts(p, arma::fill::zeros);
  std::vector< std::vector<unsigned int> > S;
  S.reserve(J);
  for (unsigned int j = 0; j < J; ++j) {
    arma::ivec nuj = as<arma::ivec>(nu_list[j]);
    unsigned int K = (unsigned int) K_vec[j];
    std::vector<unsigned int> sj; sj.reserve(K);
    for (unsigned int m = 0; m < K; ++m) {
      unsigned int idx = (unsigned int) nuj[m];
      counts[idx] += 1.0;
      sj.push_back(idx);
    }
    if (K > 0) S.push_back(sj);
  }

  arma::vec ap = (alpha / (double) p) + counts;
  arma::vec g(p);
  for (unsigned int l = 0; l < p; ++l) g[l] = R::rgamma(ap[l], 1.0);
  double gs = arma::sum(g);
  if (!(gs > 0.0) || !std::isfinite(gs)) return pi_cur;
  arma::vec pstar = g / gs;

  double a1 = alpha / (double) p - 1.0;
  double logh_star = 0.0, logh_cur = 0.0, logq_star = 0.0, logq_cur = 0.0;
  for (unsigned int l = 0; l < p; ++l) {
    logh_star += a1 * std::log(pstar[l]  + 1e-300);
    logh_cur  += a1 * std::log(pi_cur[l] + 1e-300);
    logq_star += (ap[l] - 1.0) * std::log(pstar[l]  + 1e-300);
    logq_cur  += (ap[l] - 1.0) * std::log(pi_cur[l] + 1e-300);
  }
  for (size_t j = 0; j < S.size(); ++j) {
    logh_star += pwor_logset(S[j], pstar);
    logh_cur  += pwor_logset(S[j], pi_cur);
  }

  double logA = (logh_star - logh_cur) + (logq_cur - logq_star);
  if (std::isfinite(logA) && std::log(R::runif(0.0, 1.0)) < logA) return pstar;
  return pi_cur;
}
