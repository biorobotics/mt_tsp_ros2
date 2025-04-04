#pragma once
#include <Eigen/Dense>
#include <omp.h>

using namespace Eigen;

void parallel_squareform(Ref<Matrix<double, -1, -1, RowMajor>> M, const Ref<const VectorXd>& v, int n, int num_openmp_threads) {
  double *Mptr = M.data();
  const double *vptr = v.data();

  omp_set_num_threads(num_openmp_threads);

  #pragma omp parallel for
  for (size_t i = 1; i < n; ++i) {
    double *it1 = Mptr + 1 + (n + 1)*(i - 1);
    // Populate M(i - 1, i) and everything to the right
    const double *vptr1 = vptr + n*(i - 1) - (i - 1)*(i)/2;
    memcpy(it1, vptr1, (n - i) * sizeof(double));

    // Populate M(i, i - 1) and everything below
    #pragma omp parallel for
    for (size_t j = i; j < n; ++j) {
      double *it2 = Mptr + i * (n + 1) - 1 + n*(j - i);
      const double *vptr2 = vptr1 + j - i;
      *it2 = *vptr2;
    }
  }
}
