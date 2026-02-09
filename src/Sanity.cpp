#include <Rcpp.h>
#include <cmath>
#include <vector>
#include "PsiGamma.h"
#include "FukushimaLambertW.h"

using namespace Rcpp;

// Forward declaration
double fitfrac(double* f, double q, const double* counts, int C, double gene_size, double v, const double* cell_size, const double b);
double normalization(const std::vector<double>& Q, double q);
double get_epsilon_2(double d, double v, double gene_size, double f);

// [[Rcpp::export]]
List get_gene_expression_level(NumericVector counts, NumericVector cell_size,
                               double vmin, double vmax, int numbin,
                               double a, double b) {
    // Equations mentioned are from supplementary of Sanity paper:
    // https://doi.org/10.1038/s41587-021-00875-x
    //  Nomenclature (subscript g is omitted since this function runs for 1 gene):
    //  - counts = n_{gc}
    //  - cell_size = N_c
    //  - gene_size = n_g
    //  - var_g = v_g
    //  - delta_v = \delta_{gc} as a function of v_g

    /***********************/
    /*** Setup variables ***/
    /***********************/

    // Extract summary data
    const int C = cell_size.length();  // total number of cells
    const double gene_size = sum(counts) + a;  // n_g + a: total number of UMI of gene
    const double lib_size = sum(cell_size);

    // Initialize output variables
    double mu = 0.0;             // mean gene expression (log-scale)
    double var_mu = 0.0;         // variance of mean gene expression
    double var = 0.0;            // variance of gene expression
    NumericVector delta(C);      // logfold of gene expression
    NumericVector var_delta(C);  // variance of logFC
    NumericVector lik(numbin);   // likelihood of each variance level
    NumericVector var_g(numbin); // gene variance

    // Temporary variables
    double q_exp = 0.0;               // expected q_g value to determine mean expression
    std::vector<double> q_v(numbin);  // fitted q value for each bin
    std::vector<double> f(C);
    std::vector<std::vector<double>> delta_v(numbin, std::vector<double>(C));  // logFC as function of variance
    std::vector<std::vector<double>> var_delta_v(numbin, std::vector<double>(C));

    // Initialize Auxiliary Variables and Parameters
    double Lmax = -1e+100;  // to normalize likelihood
    const double dv = log(vmax / vmin) / (numbin - 1.0);  // bin variance log-step
    double q = log(lib_size + b); // initial guess for q

    /*********************************************************/
    /*** Main Loop Compute mu and delta as function of var ***/
    /*********************************************************/
    for (int k = 0; k < numbin; ++k) {
        var_g[k] = vmin * exp(dv * k);  // current gene variance level
        double inv_v = 1.0 / var_g[k];
        double beta = inv_v / gene_size;  // denominator of f_gc (eq 26)

        // fit q(g) = mean number of reads deviating from mean (log-scale)
        // and f(g,c) = ratio of deviating reads attributes to cell c
        // using eq 27
        q[k] = fitfrac(f.data(), q, counts.begin(), C, gene_size, var_g[k], cell_size.begin(), b);

        /*** Compute optimal delta given v_g ***/
        for (int c = 0; c < C; ++c) {
            delta_v[k][c] = v * (counts[c] - f[c] * gene_size); // eq 28 in supp rearranged to avoid logs
        }

        /*** Laplace approx of likelihood at optimal delta ***/
        // Compute  loglik of P(n_c, delta | var_g) (eq 20) at optimal delta (L^*)
        // Terms:
        //   1. C*log(var_g)/2: Gaussian prior normalization term (log-det of cov)
        //   2. sum(delta^2/2/var_g): Gaussian prior kernel (mahanalobis dist)
        //   3. sum(n * delta): Poisson kernel
        //   4. (n + a) * q: Poisson normalization
        double L = 0.5 * C * log(inv_v);  // 1st term (minus obsorved by inv_v)
        double sum_delta_v2 = 0;  // 2nd term
        for (int c = 0; c < C; ++c) {
            sum_delta_v2 += delta_v[k][c] * delta_v[k][c];
            L += counts[c] * delta_v[k][c];  // 3rd term
        }
        L -= 0.5 * sum_delta_v2 * inv_v;  // 2nd term
        L -= gene_size * q[k];      // 4th term

        // Equation 33 of Supp: Compute log-determinant of the Hessian at L*
        // Notes:
        //   - The Hessian is a combination of diagonal + C rank 1 correction terms
        //   - The term 1 - sum((n+a)f^2 / (n+a)f + 1/var_g) = 1 - sum(f^2/(f+beta))
        //     corresponds to the correction terms. We store it (det_corr)
        //     because it's used in var(delta) calculation
        //   - (n+a)f + 1/var_g = (n+a)(f + beta), (n+a) does not depend on cell
        double ldet = 0.0;
        double sum_f2_over_f_plus_beta = 0.0;
        for (int c = 0; c < C; ++c) {
            double denom = f[c] + beta;
            var_delta_v[k][c] = f[c] * f[c] / denom;  // f^2 / (f + beta)

            ldet += log(denom);  // 2nd term (product over c) without (n+a)
            sum_f2_over_f_plus_beta += var_delta_v[k][c];  // sum of 1st term
        }
        double det_corr = 1.0 - sum_f2_over_f_plus_beta; // rank 1 corrections
        ldet += log(det_corr);  // 1st term
        // ldet += C * log(gene_size); // missing (n+a) from 2nd term (not necessary b/c same across cells)

        // P(n_c | var_g) = exp(L^*) / sqrt(det(H)) {eq 32}
        L -= 0.5 * ldet;
        lik[k] = L;
        // update Lmax
        if (L > Lmax) {
            Lmax = L;
        }

        /*** Compute variance of delta given v_g ***/
        // Equation 37 of Supplementary
        for (int c = 0; c < C; ++c) {
            if (counts[c] <= 0.5) {
                // For zero counts compute symmetric CI
                var_delta_v[k][c] = get_epsilon_2(delta_v[k][c], var_g[k], gene_size, f[c]);
            } else {
                // 1 - sum(f_c^2/(f_c + beta) for c != c') == det_corr + f^2_c'/(f_c' + beta)
                var_delta_v[k][c] += det_corr;
                var_delta_v[k][c] /= det_corr * (gene_size * f[c] + inv_v);
            }
        }

        /*** Compute mean expression given v_g ***/
        //mu_v[k] = psi::Psi_0(gene_size) - q; // Equation 45 in Supplementary

    }

    // Normalize likelihood
    double sum_L = 0.0;
    for (int k = 0; k < numbin; ++k) {
        lik[k] -= Lmax; // for numerical stability
        lik[k] = exp(lik[k]);
        sum_L += lik[k];
    }
    for (int k = 0; k < numbin; ++k) {
        lik[k] /= sum_L;
    }

    // Compute average gene expression (mu) and variance
    for (int k = 0; k < numbin; ++k) {
        q_exp += lik[k] * q_v[k];
        //mu += lik[k] * mu_v[k];
        var += lik[k] * var_g[k];
    }
    mu = psi::Psi_0(gene_size) + q_exp;

    // Compute mu variance
    var_mu = psi::Psi_1(gene_size);
    for (int k = 0; k < numbin; ++k) {
        double dev_mu = q_v[k] - q_exp;
        var_mu += lik[k] * dev_mu * dev_mu;
    }

    // Compute delta and its variance
    for (int c = 0; c < C; ++c) {
        delta[c] = 0.0;
        for (int k = 0; k < numbin; ++k) {
            delta[c] += lik[k] * delta_v[k][c];
        }

        var_delta[c] = 0.0;
        for (int k = 0; k < numbin; ++k) {
            double dev_delta = delta_v[k][c] - delta[c];
            var_delta[c] += lik[k] * (dev_delta * dev_delta + var_delta_v[k][c]);
        }
    }

    // Return results as a list
    return List::create(
        Named("mu") = mu,
        Named("var_mu") = var_mu,
        Named("var") = var,
        Named("delta") = delta,
        Named("var_delta") = var_delta,
        Named("lik") = lik,
        Named("var_k") = var_g
    );
}

double fitfrac(double* f, double q, const double* counts, int C, double gene_size, double v,
               const double* cell_size, const double b) {

    double qmin, qmax, nor;
    double dq = log(2.0);
    std::vector<double> Q(C);
    double beta = 1.0 / (gene_size * v);
    double logbeta = log(beta);

    for (int c = 0; c < C; ++c) {
        Q[c] = log(cell_size[c]) + counts[c] * v - logbeta;
    }

    // Get initial normalization
    nor = beta * normalization(Q, q) + b * exp(-q);

    // Bracket the solution
    if (nor < 1) {
        while (nor < 1) {
            qmax = q;
            q = q - dq;
            nor = beta * normalization(Q, q) + b * exp(-q);
        }
        qmin = q;
    } else {
        while (nor > 1) {
            qmin = q;
            q = q + dq;
            nor = beta * normalization(Q, q) + b * exp(-q);
        }
        qmax = q;
    }

    // Bisection method
    const double tol = 1e-7;
    double diff = 1.0;
    while (diff > tol) {
        q = (qmax + qmin) / 2.0;
        nor = beta * normalization(Q, q) + b * exp(-q);
        if (nor > 1) {
            qmin = q;
        } else {
            qmax = q;
        }
        diff = fabs(nor - 1.0);
    }

    q = (qmax + qmin) / 2.0;
    for (int c = 0; c < C; ++c) {
        double x = Q[c] - q;
        if (x > 50) {
            f[c] = beta * (x - log(x) + log(x)/x);
        } else {
            f[c] = beta * Fukushima::LambertW0(exp(x));
        }
    }

    return q;
}

double normalization(const std::vector<double>& Q, double q) {
    double nor = 0;
    for (size_t c = 0; c < Q.size(); ++c) {
        double x = Q[c] - q;
        if (std::isnan(x)) {
            Rcpp::warning("x = nan in normalization!");
            continue;
        }
        if (x > 50) {
            nor += (x - log(x) + log(x)/x);
        } else {
            nor += Fukushima::LambertW0(exp(x));
        }
    }
    return nor;
}

double get_epsilon_2(double d, double v, double gene_size, double f){

    // Solve equation 38 of supplementary
    double sigma_c;
    double dL;
    double s_low = 0.0;
    double vnf = v * gene_size * f;
    // approximate exp(x) - 1 = x (valid for small x) and solve the resulting
    // quadratic:  A = (1 + vnf) / 2, B = d + vnf, C = -v
    double s_high = ( -(d+vnf) + sqrt( (d+vnf)*(d+vnf) + v*(1.0 + vnf) ) )/(1.0 + vnf);

    // bisection method: dL == 1/2
    double tol = 0.0000001;
    double diff = 1.0;
    while(diff > tol){
        sigma_c = 0.5 * (s_high + s_low);
        // approximate log(1 - f*(exp(e) - 1)) = f * (exp(e) - 1)
        dL = sigma_c * (d + 0.5*sigma_c) / v + gene_size * f*expm1(sigma_c);
        if (dL < 0.5) {
            s_low = sigma_c;
        } else {
            s_high = sigma_c;
        }
        diff = fabs(dL-0.5);
    }
    return sigma_c * sigma_c;
}
