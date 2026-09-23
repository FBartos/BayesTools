#ifndef BTLKJCORE_H_
#define BTLKJCORE_H_

namespace bayestools {
namespace lkj {

unsigned int n_cpc(unsigned int K);
unsigned int cpc_index(unsigned int i, unsigned int j);
unsigned int flat_index(unsigned int row, unsigned int column, unsigned int K);

bool valid_u(double const *u, unsigned int length);
bool valid_u_strided(double const *u, unsigned int length, unsigned int stride);
bool valid_alpha(double const *alpha, unsigned int length);

double cpc_from_u(double u);
double pair_alpha(unsigned int i, unsigned int K, double eta);

void fill_alpha(double *alpha, unsigned int K, double eta);
void fill_cholesky_from_u(double *value, double const *u, unsigned int K);
void fill_cholesky_from_u_strided(double *value, double const *u, unsigned int K, unsigned int stride);
void fill_corr_from_cholesky(double *value, double const *L, unsigned int K);
void fill_corr_from_u(double *value, double const *u, unsigned int K);

double log_density_u(double const *u, double const *alpha, unsigned int length);
double log_density_u_strided(double const *u, double const *alpha, unsigned int length, unsigned int stride);
double log_density_u_strided_valid_alpha(double const *u, double const *alpha, unsigned int length, unsigned int stride);

}
}

#endif
