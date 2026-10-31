#ifndef GMS_BLAS3_UNROLLED_BLOCKED_DP_HPP
#define GMS_BLAS3_UNROLLED_BLOCKED_DP_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace gms {
namespace blas3 {

inline std::size_t f2d_idx(std::size_t row, std::size_t col, std::size_t ld) noexcept {
    return row + col * ld;
}

inline bool lsame(char ca, char cb) noexcept {
    unsigned char a = static_cast<unsigned char>(ca);
    unsigned char b = static_cast<unsigned char>(cb);
    if (a == b) {
        return true;
    }
    if (a >= 'A' && a <= 'Z') a = static_cast<unsigned char>(a + ('a' - 'A'));
    if (b >= 'A' && b <= 'Z') b = static_cast<unsigned char>(b + ('a' - 'A'));
    return a == b;
}

inline void dgemm_core(bool nota, bool notb,
    int m, int n, int k,
    double alpha,
    const double* A, int lda,
    const double* B, int ldb,
    double beta,
    double* C, int ldc,
    int mb, int nb, int nbt, int kb) {
    const double zero = 0.0;
    const double one = 1.0;
    const std::size_t t1_size = static_cast<std::size_t>(kb) * static_cast<std::size_t>(mb);
    const std::size_t t2_size = static_cast<std::size_t>(kb) * static_cast<std::size_t>(nbt);
    std::unique_ptr<double[]> t1(new double[t1_size]);
    std::unique_ptr<double[]> t2(new double[t2_size]);
    std::fill_n(t1.get(), t1_size, zero);
    std::fill_n(t2.get(), t2_size, zero);

    if (notb) {
        for (int jj = 1; jj <= n; jj += nb) {
            const int jsec = std::min(nb, n - jj + 1);
            const int ujsec = jsec - (jsec % 4);
            for (int ll = 1; ll <= k; ll += kb) {
                const int lsec = std::min(kb, k - ll + 1);
                const int ulsec = lsec - (lsec % 2);
                double delta = one;
                if (ll == 1) delta = beta;
                for (int ii = 1; ii <= m; ii += mb) {
                    const int isec = std::min(mb, m - ii + 1);
                    const int uisec = isec - (isec % 2);
                    if (nota) {
                        for (int l = ll; l <= ll + ulsec - 1; l += 2) {
                            for (int i = ii; i <= ii + uisec - 1; i += 2) {
                                const std::size_t row0 = static_cast<std::size_t>(i - 1);
                                const std::size_t row1 = static_cast<std::size_t>(i);
                                const std::size_t col0 = static_cast<std::size_t>(l - 1);
                                const std::size_t col1 = static_cast<std::size_t>(l);
                                const std::size_t t1row0 = static_cast<std::size_t>(l - ll);
                                const std::size_t t1row1 = static_cast<std::size_t>(l - ll + 1);
                                const std::size_t ti0 = static_cast<std::size_t>(i - ii);
                                const std::size_t ti1 = static_cast<std::size_t>(i - ii + 1);
                                const std::size_t idx00 = f2d_idx(t1row0, ti0, static_cast<std::size_t>(kb));
                                const std::size_t idx10 = f2d_idx(t1row1, ti0, static_cast<std::size_t>(kb));
                                const std::size_t idx01 = f2d_idx(t1row0, ti1, static_cast<std::size_t>(kb));
                                const std::size_t idx11 = f2d_idx(t1row1, ti1, static_cast<std::size_t>(kb));
                                t1[idx00] = alpha * A[f2d_idx(row0, col0, static_cast<std::size_t>(lda))];
                                t1[idx10] = alpha * A[f2d_idx(row0, col1, static_cast<std::size_t>(lda))];
                                t1[idx01] = alpha * A[f2d_idx(row1, col0, static_cast<std::size_t>(lda))];
                                t1[idx11] = alpha * A[f2d_idx(row1, col1, static_cast<std::size_t>(lda))];
                            }
                            if (uisec < isec) {
                                const std::size_t row = static_cast<std::size_t>(ii + isec - 2);
                                const std::size_t t1row0 = static_cast<std::size_t>(l - ll);
                                const std::size_t t1row1 = static_cast<std::size_t>(l - ll + 1);
                                const std::size_t ti = static_cast<std::size_t>(isec - 1);
                                t1[f2d_idx(t1row0, ti, static_cast<std::size_t>(kb))] = alpha * A[f2d_idx(row, col0, static_cast<std::size_t>(lda))];
                                t1[f2d_idx(t1row1, ti, static_cast<std::size_t>(kb))] = alpha * A[f2d_idx(row, col1, static_cast<std::size_t>(lda))];
                            }
                        }
                        if (ulsec < lsec) {
                            for (int i = ii; i <= ii + isec - 1; ++i) {
                                const std::size_t row = static_cast<std::size_t>(i - 1);
                                const std::size_t col = static_cast<std::size_t>(ll + lsec - 2);
                                t1[f2d_idx(static_cast<std::size_t>(lsec - 1), static_cast<std::size_t>(i - ii), static_cast<std::size_t>(kb))] = alpha * A[f2d_idx(row, col, static_cast<std::size_t>(lda))];
                            }
                        }
                    } else {
                        for (int i = ii; i <= ii + uisec - 1; i += 2) {
#pragma omp simd
                            for (int l = ll; l <= ll + ulsec - 1; l += 2) {
                                const std::size_t row0 = static_cast<std::size_t>(l - 1);
                                const std::size_t row1 = static_cast<std::size_t>(l);
                                const std::size_t col0 = static_cast<std::size_t>(i - 1);
                                const std::size_t col1 = static_cast<std::size_t>(i);
                                const std::size_t t1row0 = static_cast<std::size_t>(l - ll);
                                const std::size_t t1row1 = static_cast<std::size_t>(l - ll + 1);
                                const std::size_t ti0 = static_cast<std::size_t>(i - ii);
                                const std::size_t ti1 = static_cast<std::size_t>(i - ii + 1);
                                t1[f2d_idx(t1row0, ti0, static_cast<std::size_t>(kb))] = alpha * A[f2d_idx(row0, col0, static_cast<std::size_t>(lda))];
                                t1[f2d_idx(t1row0, ti1, static_cast<std::size_t>(kb))] = alpha * A[f2d_idx(row0, col1, static_cast<std::size_t>(lda))];
                                t1[f2d_idx(t1row1, ti0, static_cast<std::size_t>(kb))] = alpha * A[f2d_idx(row1, col0, static_cast<std::size_t>(lda))];
                                t1[f2d_idx(t1row1, ti1, static_cast<std::size_t>(kb))] = alpha * A[f2d_idx(row1, col1, static_cast<std::size_t>(lda))];
                            }
                            if (ulsec < lsec) {
                                const std::size_t row = static_cast<std::size_t>(ll + lsec - 2);
                                const std::size_t ti0 = static_cast<std::size_t>(i - ii);
                                const std::size_t ti1 = static_cast<std::size_t>(i - ii + 1);
                                t1[f2d_idx(static_cast<std::size_t>(lsec - 1), ti0, static_cast<std::size_t>(kb))] = alpha * A[f2d_idx(row, static_cast<std::size_t>(i - 1), static_cast<std::size_t>(lda))];
                                t1[f2d_idx(static_cast<std::size_t>(lsec - 1), ti1, static_cast<std::size_t>(kb))] = alpha * A[f2d_idx(row, static_cast<std::size_t>(i), static_cast<std::size_t>(lda))];
                            }
                        }
                        if (uisec < isec) {
#pragma omp simd
                            for (int l = ll; l <= ll + lsec - 1; ++l) {
                                const std::size_t row = static_cast<std::size_t>(l - 1);
                                const std::size_t col = static_cast<std::size_t>(ii + isec - 2);
                                t1[f2d_idx(static_cast<std::size_t>(l - ll), static_cast<std::size_t>(isec - 1), static_cast<std::size_t>(kb))] = alpha * A[f2d_idx(row, col, static_cast<std::size_t>(lda))];
                            }
                        }
                    }

                    const int ui4sec = isec - (isec % 4);
                    for (int j = jj; j <= jj + ujsec - 1; j += 4) {
#pragma omp simd
                        for (int i = ii; i <= ii + ui4sec - 1; i += 4) {
                            const std::size_t ci00 = static_cast<std::size_t>(i - 1);
                            const std::size_t ci10 = static_cast<std::size_t>(i);
                            const std::size_t ci20 = static_cast<std::size_t>(i + 1);
                            const std::size_t ci30 = static_cast<std::size_t>(i + 2);
                            const std::size_t cj00 = static_cast<std::size_t>(j - 1);
                            const std::size_t cj10 = static_cast<std::size_t>(j);
                            const std::size_t cj20 = static_cast<std::size_t>(j + 1);
                            const std::size_t cj30 = static_cast<std::size_t>(j + 2);
                            double f11 = delta * C[f2d_idx(ci00, cj00, static_cast<std::size_t>(ldc))];
                            double f21 = delta * C[f2d_idx(ci10, cj00, static_cast<std::size_t>(ldc))];
                            double f12 = delta * C[f2d_idx(ci00, cj10, static_cast<std::size_t>(ldc))];
                            double f22 = delta * C[f2d_idx(ci10, cj10, static_cast<std::size_t>(ldc))];
                            double f13 = delta * C[f2d_idx(ci00, cj20, static_cast<std::size_t>(ldc))];
                            double f23 = delta * C[f2d_idx(ci10, cj20, static_cast<std::size_t>(ldc))];
                            double f14 = delta * C[f2d_idx(ci00, cj30, static_cast<std::size_t>(ldc))];
                            double f24 = delta * C[f2d_idx(ci10, cj30, static_cast<std::size_t>(ldc))];
                            double f31 = delta * C[f2d_idx(ci20, cj00, static_cast<std::size_t>(ldc))];
                            double f41 = delta * C[f2d_idx(ci30, cj00, static_cast<std::size_t>(ldc))];
                            double f32 = delta * C[f2d_idx(ci20, cj10, static_cast<std::size_t>(ldc))];
                            double f42 = delta * C[f2d_idx(ci30, cj10, static_cast<std::size_t>(ldc))];
                            double f33 = delta * C[f2d_idx(ci20, cj20, static_cast<std::size_t>(ldc))];
                            double f43 = delta * C[f2d_idx(ci30, cj20, static_cast<std::size_t>(ldc))];
                            double f34 = delta * C[f2d_idx(ci20, cj30, static_cast<std::size_t>(ldc))];
                            double f44 = delta * C[f2d_idx(ci30, cj30, static_cast<std::size_t>(ldc))];
                            for (int l = ll; l <= ll + lsec - 1; ++l) {
                                const std::size_t tl = static_cast<std::size_t>(l - ll);
                                const double b00 = B[f2d_idx(static_cast<std::size_t>(l - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldb))];
                                const double b01 = B[f2d_idx(static_cast<std::size_t>(l - 1), static_cast<std::size_t>(j), static_cast<std::size_t>(ldb))];
                                const double b02 = B[f2d_idx(static_cast<std::size_t>(l - 1), static_cast<std::size_t>(j + 1), static_cast<std::size_t>(ldb))];
                                const double b03 = B[f2d_idx(static_cast<std::size_t>(l - 1), static_cast<std::size_t>(j + 2), static_cast<std::size_t>(ldb))];
                                const double t11 = t1[f2d_idx(tl, static_cast<std::size_t>(i - ii), static_cast<std::size_t>(kb))];
                                const double t21 = t1[f2d_idx(tl, static_cast<std::size_t>(i - ii + 1), static_cast<std::size_t>(kb))];
                                const double t31 = t1[f2d_idx(tl, static_cast<std::size_t>(i - ii + 2), static_cast<std::size_t>(kb))];
                                const double t41 = t1[f2d_idx(tl, static_cast<std::size_t>(i - ii + 3), static_cast<std::size_t>(kb))];
                                f11 += t11 * b00;
                                f21 += t21 * b00;
                                f12 += t11 * b01;
                                f22 += t21 * b01;
                                f13 += t11 * b02;
                                f23 += t21 * b02;
                                f14 += t11 * b03;
                                f24 += t21 * b03;
                                f31 += t31 * b00;
                                f41 += t41 * b00;
                                f32 += t31 * b01;
                                f42 += t41 * b01;
                                f33 += t31 * b02;
                                f43 += t41 * b02;
                                f34 += t31 * b03;
                                f44 += t41 * b03;
                            }
                            C[f2d_idx(ci00, cj00, static_cast<std::size_t>(ldc))] = f11;
                            C[f2d_idx(ci10, cj00, static_cast<std::size_t>(ldc))] = f21;
                            C[f2d_idx(ci00, cj10, static_cast<std::size_t>(ldc))] = f12;
                            C[f2d_idx(ci10, cj10, static_cast<std::size_t>(ldc))] = f22;
                            C[f2d_idx(ci00, cj20, static_cast<std::size_t>(ldc))] = f13;
                            C[f2d_idx(ci10, cj20, static_cast<std::size_t>(ldc))] = f23;
                            C[f2d_idx(ci00, cj30, static_cast<std::size_t>(ldc))] = f14;
                            C[f2d_idx(ci10, cj30, static_cast<std::size_t>(ldc))] = f24;
                            C[f2d_idx(ci20, cj00, static_cast<std::size_t>(ldc))] = f31;
                            C[f2d_idx(ci30, cj00, static_cast<std::size_t>(ldc))] = f41;
                            C[f2d_idx(ci20, cj10, static_cast<std::size_t>(ldc))] = f32;
                            C[f2d_idx(ci30, cj10, static_cast<std::size_t>(ldc))] = f42;
                            C[f2d_idx(ci20, cj20, static_cast<std::size_t>(ldc))] = f33;
                            C[f2d_idx(ci30, cj20, static_cast<std::size_t>(ldc))] = f43;
                            C[f2d_idx(ci20, cj30, static_cast<std::size_t>(ldc))] = f34;
                            C[f2d_idx(ci30, cj30, static_cast<std::size_t>(ldc))] = f44;
                        }
                    }

                    if (ujsec < jsec) {
                        for (int j = jj + ujsec; j <= jj + jsec - 1; ++j) {
                            for (int i = ii; i <= ii + isec - 1; ++i) {
                                double f11 = delta * C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))];
                                for (int l = ll; l <= ll + lsec - 1; ++l) {
                                    const std::size_t tl = static_cast<std::size_t>(l - ll);
                                    const std::size_t ti = static_cast<std::size_t>(i - ii);
                                    f11 += t1[f2d_idx(tl, ti, static_cast<std::size_t>(kb))] * B[f2d_idx(static_cast<std::size_t>(l - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldb))];
                                }
                                C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] = f11;
                            }
                        }
                    }
                }
            }
        }
    } else {
        for (int jj = 1; jj <= n; jj += nbt) {
            const int jsec = std::min(nbt, n - jj + 1);
            for (int ll = 1; ll <= k; ll += kb) {
                const int lsec = std::min(kb, k - ll + 1);
                double delta = one;
                if (ll == 1) delta = beta;
                const int ulsec = lsec - (lsec % 2);
                for (int l = ll; l <= ll + ulsec - 1; l += 2) {
                    for (int j = jj; j <= jj + jsec - 1; j += 2) {
                        const std::size_t t2row0 = static_cast<std::size_t>(l - ll);
                        const std::size_t t2row1 = static_cast<std::size_t>(l - ll + 1);
                        const std::size_t t2col0 = static_cast<std::size_t>(j - jj);
                        const std::size_t t2col1 = static_cast<std::size_t>(j - jj + 1);
                        t2[f2d_idx(t2row0, t2col0, static_cast<std::size_t>(kb))] = alpha * B[f2d_idx(static_cast<std::size_t>(j - 1), static_cast<std::size_t>(l - 1), static_cast<std::size_t>(ldb))];
                        t2[f2d_idx(t2row1, t2col0, static_cast<std::size_t>(kb))] = alpha * B[f2d_idx(static_cast<std::size_t>(j - 1), static_cast<std::size_t>(l), static_cast<std::size_t>(ldb))];
                        t2[f2d_idx(t2row0, t2col1, static_cast<std::size_t>(kb))] = alpha * B[f2d_idx(static_cast<std::size_t>(j), static_cast<std::size_t>(l - 1), static_cast<std::size_t>(ldb))];
                        t2[f2d_idx(t2row1, t2col1, static_cast<std::size_t>(kb))] = alpha * B[f2d_idx(static_cast<std::size_t>(j), static_cast<std::size_t>(l), static_cast<std::size_t>(ldb))];
                    }
                    if (jsec % 2 != 0) {
                        const std::size_t jlast = static_cast<std::size_t>(jsec - 1);
                        t2[f2d_idx(static_cast<std::size_t>(l - ll), jlast, static_cast<std::size_t>(kb))] = alpha * B[f2d_idx(static_cast<std::size_t>(jj + jsec - 2), static_cast<std::size_t>(l - 1), static_cast<std::size_t>(ldb))];
                        t2[f2d_idx(static_cast<std::size_t>(l - ll + 1), jlast, static_cast<std::size_t>(kb))] = alpha * B[f2d_idx(static_cast<std::size_t>(jj + jsec - 2), static_cast<std::size_t>(l), static_cast<std::size_t>(ldb))];
                    }
                }
                if (ulsec < lsec) {
                    for (int j = jj; j <= jj + jsec - 1; ++j) {
                        const std::size_t tl = static_cast<std::size_t>(lsec - 1);
                        const std::size_t tj = static_cast<std::size_t>(j - jj);
                        t2[f2d_idx(tl, tj, static_cast<std::size_t>(kb))] = alpha * B[f2d_idx(static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ll + lsec - 2), static_cast<std::size_t>(ldb))];
                    }
                }
                const int ujsec = jsec - (jsec % 4);
                for (int ii = 1; ii <= m; ii += mb) {
                    const int isec = std::min(mb, m - ii + 1);
                    const int ui2sec = isec - (isec % 2);
                    if (nota) {
                        for (int l = ll; l <= ll + ulsec - 1; l += 2) {
                            for (int i = ii; i <= ii + ui2sec - 1; i += 2) {
                                const std::size_t row0 = static_cast<std::size_t>(i - 1);
                                const std::size_t row1 = static_cast<std::size_t>(i);
                                const std::size_t col0 = static_cast<std::size_t>(l - 1);
                                const std::size_t col1 = static_cast<std::size_t>(l);
                                const std::size_t ti0 = static_cast<std::size_t>(i - ii);
                                const std::size_t ti1 = static_cast<std::size_t>(i - ii + 1);
                                const std::size_t t1row0 = static_cast<std::size_t>(l - ll);
                                const std::size_t t1row1 = static_cast<std::size_t>(l - ll + 1);
                                t1[f2d_idx(t1row0, ti0, static_cast<std::size_t>(kb))] = A[f2d_idx(row0, col0, static_cast<std::size_t>(lda))];
                                t1[f2d_idx(t1row1, ti0, static_cast<std::size_t>(kb))] = A[f2d_idx(row0, col1, static_cast<std::size_t>(lda))];
                                t1[f2d_idx(t1row0, ti1, static_cast<std::size_t>(kb))] = A[f2d_idx(row1, col0, static_cast<std::size_t>(lda))];
                                t1[f2d_idx(t1row1, ti1, static_cast<std::size_t>(kb))] = A[f2d_idx(row1, col1, static_cast<std::size_t>(lda))];
                            }
                            if (ui2sec < isec) {
                                const std::size_t row = static_cast<std::size_t>(ii + isec - 2);
                                const std::size_t t1row0 = static_cast<std::size_t>(l - ll);
                                const std::size_t t1row1 = static_cast<std::size_t>(l - ll + 1);
                                const std::size_t ti = static_cast<std::size_t>(isec - 1);
                                t1[f2d_idx(t1row0, ti, static_cast<std::size_t>(kb))] = A[f2d_idx(row, static_cast<std::size_t>(l - 1), static_cast<std::size_t>(lda))];
                                t1[f2d_idx(t1row1, ti, static_cast<std::size_t>(kb))] = A[f2d_idx(row, static_cast<std::size_t>(l), static_cast<std::size_t>(lda))];
                            }
                        }
                        if (ulsec < lsec) {
                            for (int i = ii; i <= ii + isec - 1; ++i) {
                                const std::size_t row = static_cast<std::size_t>(i - 1);
                                const std::size_t col = static_cast<std::size_t>(ll + lsec - 2);
                                t1[f2d_idx(static_cast<std::size_t>(lsec - 1), static_cast<std::size_t>(i - ii), static_cast<std::size_t>(kb))] = A[f2d_idx(row, col, static_cast<std::size_t>(lda))];
                            }
                        }
                    } else {
                        for (int i = ii; i <= ii + ui2sec - 1; i += 2) {
#pragma omp simd
                            for (int l = ll; l <= ll + ulsec - 1; l += 2) {
                                const std::size_t row0 = static_cast<std::size_t>(l - 1);
                                const std::size_t row1 = static_cast<std::size_t>(l);
                                const std::size_t col0 = static_cast<std::size_t>(i - 1);
                                const std::size_t col1 = static_cast<std::size_t>(i);
                                const std::size_t ti0 = static_cast<std::size_t>(i - ii);
                                const std::size_t ti1 = static_cast<std::size_t>(i - ii + 1);
                                const std::size_t t1row0 = static_cast<std::size_t>(l - ll);
                                const std::size_t t1row1 = static_cast<std::size_t>(l - ll + 1);
                                t1[f2d_idx(t1row0, ti0, static_cast<std::size_t>(kb))] = A[f2d_idx(row0, col0, static_cast<std::size_t>(lda))];
                                t1[f2d_idx(t1row0, ti1, static_cast<std::size_t>(kb))] = A[f2d_idx(row0, col1, static_cast<std::size_t>(lda))];
                                t1[f2d_idx(t1row1, ti0, static_cast<std::size_t>(kb))] = A[f2d_idx(row1, col0, static_cast<std::size_t>(lda))];
                                t1[f2d_idx(t1row1, ti1, static_cast<std::size_t>(kb))] = A[f2d_idx(row1, col1, static_cast<std::size_t>(lda))];
                            }
                            if (ulsec < lsec) {
                                const std::size_t col = static_cast<std::size_t>(ll + lsec - 2);
                                const std::size_t ti0 = static_cast<std::size_t>(i - ii);
                                const std::size_t ti1 = static_cast<std::size_t>(i - ii + 1);
                                t1[f2d_idx(static_cast<std::size_t>(lsec - 1), ti0, static_cast<std::size_t>(kb))] = A[f2d_idx(col, static_cast<std::size_t>(i - 1), static_cast<std::size_t>(lda))];
                                t1[f2d_idx(static_cast<std::size_t>(lsec - 1), ti1, static_cast<std::size_t>(kb))] = A[f2d_idx(col, static_cast<std::size_t>(i), static_cast<std::size_t>(lda))];
                            }
                        }
                        if (ui2sec < isec) {
#pragma omp simd
                            for (int l = ll; l <= ll + lsec - 1; ++l) {
                                const std::size_t row = static_cast<std::size_t>(l - 1);
                                const std::size_t col = static_cast<std::size_t>(ii + isec - 2);
                                t1[f2d_idx(static_cast<std::size_t>(l - ll), static_cast<std::size_t>(isec - 1), static_cast<std::size_t>(kb))] = A[f2d_idx(row, col, static_cast<std::size_t>(lda))];
                            }
                        }
                    }

                    const int ui4sec = isec - (isec % 4);
                    for (int j = jj; j <= jj + ujsec - 1; j += 4) {
#pragma omp simd
                        for (int i = ii; i <= ii + ui4sec - 1; i += 4) {
                            const std::size_t ci00 = static_cast<std::size_t>(i - 1);
                            const std::size_t ci10 = static_cast<std::size_t>(i);
                            const std::size_t ci20 = static_cast<std::size_t>(i + 1);
                            const std::size_t ci30 = static_cast<std::size_t>(i + 2);
                            const std::size_t cj00 = static_cast<std::size_t>(j - 1);
                            const std::size_t cj10 = static_cast<std::size_t>(j);
                            const std::size_t cj20 = static_cast<std::size_t>(j + 1);
                            const std::size_t cj30 = static_cast<std::size_t>(j + 2);
                            double f11 = delta * C[f2d_idx(ci00, cj00, static_cast<std::size_t>(ldc))];
                            double f21 = delta * C[f2d_idx(ci10, cj00, static_cast<std::size_t>(ldc))];
                            double f12 = delta * C[f2d_idx(ci00, cj10, static_cast<std::size_t>(ldc))];
                            double f22 = delta * C[f2d_idx(ci10, cj10, static_cast<std::size_t>(ldc))];
                            double f13 = delta * C[f2d_idx(ci00, cj20, static_cast<std::size_t>(ldc))];
                            double f23 = delta * C[f2d_idx(ci10, cj20, static_cast<std::size_t>(ldc))];
                            double f14 = delta * C[f2d_idx(ci00, cj30, static_cast<std::size_t>(ldc))];
                            double f24 = delta * C[f2d_idx(ci10, cj30, static_cast<std::size_t>(ldc))];
                            double f31 = delta * C[f2d_idx(ci20, cj00, static_cast<std::size_t>(ldc))];
                            double f41 = delta * C[f2d_idx(ci30, cj00, static_cast<std::size_t>(ldc))];
                            double f32 = delta * C[f2d_idx(ci20, cj10, static_cast<std::size_t>(ldc))];
                            double f42 = delta * C[f2d_idx(ci30, cj10, static_cast<std::size_t>(ldc))];
                            double f33 = delta * C[f2d_idx(ci20, cj20, static_cast<std::size_t>(ldc))];
                            double f43 = delta * C[f2d_idx(ci30, cj20, static_cast<std::size_t>(ldc))];
                            double f34 = delta * C[f2d_idx(ci20, cj30, static_cast<std::size_t>(ldc))];
                            double f44 = delta * C[f2d_idx(ci30, cj30, static_cast<std::size_t>(ldc))];
                            for (int l = ll; l <= ll + lsec - 1; ++l) {
                                const std::size_t tl = static_cast<std::size_t>(l - ll);
                                const std::size_t ti0 = static_cast<std::size_t>(i - ii);
                                const std::size_t ti1 = static_cast<std::size_t>(i - ii + 1);
                                const std::size_t ti2 = static_cast<std::size_t>(i - ii + 2);
                                const std::size_t ti3 = static_cast<std::size_t>(i - ii + 3);
                                const double c0 = t1[f2d_idx(tl, ti0, static_cast<std::size_t>(kb))];
                                const double c1 = t2[f2d_idx(tl, static_cast<std::size_t>(j - jj), static_cast<std::size_t>(kb))];
                                f11 += c0 * c1;
                                const double c2 = t1[f2d_idx(tl, ti1, static_cast<std::size_t>(kb))];
                                f21 += c2 * c1;
                                const double c3 = t2[f2d_idx(tl, static_cast<std::size_t>(j - jj + 1), static_cast<std::size_t>(kb))];
                                f12 += c0 * c3;
                                f22 += c2 * c3;
                                const double c4 = t1[f2d_idx(tl, ti2, static_cast<std::size_t>(kb))];
                                const double c5 = t2[f2d_idx(tl, static_cast<std::size_t>(j - jj + 2), static_cast<std::size_t>(kb))];
                                f13 += c0 * c4;
                                f23 += c2 * c4;
                                const double c6 = t1[f2d_idx(tl, ti3, static_cast<std::size_t>(kb))];
                                const double c7 = t2[f2d_idx(tl, static_cast<std::size_t>(j - jj + 3), static_cast<std::size_t>(kb))];
                                f14 += c0 * c7;
                                f24 += c2 * c7;
                                f31 += c4 * c1;
                                f41 += c6 * c1;
                                f32 += c4 * c3;
                                f42 += c6 * c3;
                                f33 += c4 * c5;
                                f43 += c6 * c5;
                                f34 += c4 * c7;
                                f44 += c6 * c7;
                            }
                            C[f2d_idx(ci00, cj00, static_cast<std::size_t>(ldc))] = f11;
                            C[f2d_idx(ci10, cj00, static_cast<std::size_t>(ldc))] = f21;
                            C[f2d_idx(ci00, cj10, static_cast<std::size_t>(ldc))] = f12;
                            C[f2d_idx(ci10, cj10, static_cast<std::size_t>(ldc))] = f22;
                            C[f2d_idx(ci00, cj20, static_cast<std::size_t>(ldc))] = f13;
                            C[f2d_idx(ci10, cj20, static_cast<std::size_t>(ldc))] = f23;
                            C[f2d_idx(ci00, cj30, static_cast<std::size_t>(ldc))] = f14;
                            C[f2d_idx(ci10, cj30, static_cast<std::size_t>(ldc))] = f24;
                            C[f2d_idx(ci20, cj00, static_cast<std::size_t>(ldc))] = f31;
                            C[f2d_idx(ci30, cj00, static_cast<std::size_t>(ldc))] = f41;
                            C[f2d_idx(ci20, cj10, static_cast<std::size_t>(ldc))] = f32;
                            C[f2d_idx(ci30, cj10, static_cast<std::size_t>(ldc))] = f42;
                            C[f2d_idx(ci20, cj20, static_cast<std::size_t>(ldc))] = f33;
                            C[f2d_idx(ci30, cj20, static_cast<std::size_t>(ldc))] = f43;
                            C[f2d_idx(ci20, cj30, static_cast<std::size_t>(ldc))] = f34;
                            C[f2d_idx(ci30, cj30, static_cast<std::size_t>(ldc))] = f44;
                        }
                    }
                    if (ujsec < jsec) {
                        for (int j = jj + ujsec; j <= jj + jsec - 1; ++j) {
                            for (int i = ii; i <= ii + isec - 1; ++i) {
                                double f11 = delta * C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))];
                                for (int l = ll; l <= ll + lsec - 1; ++l) {
                                    const std::size_t tl = static_cast<std::size_t>(l - ll);
                                    const std::size_t ti = static_cast<std::size_t>(i - ii);
                                    f11 += t1[f2d_idx(tl, ti, static_cast<std::size_t>(kb))] * t2[f2d_idx(tl, static_cast<std::size_t>(j - jj), static_cast<std::size_t>(kb))];
                                }
                                C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] = f11;
                            }
                        }
                    }
                }
            }
        }
    }
}

inline void dgemm(char transa, char transb,
    int m, int n, int k,
    double alpha,
    const double* A, int lda,
    const double* B, int ldb,
    double beta,
    double* C, int ldc,
    int mb, int nb, int nbt, int kb) {
    const bool nota = lsame(transa, 'N');
    const bool notb = lsame(transb, 'N');
    const int nrowa = nota ? m : k;
    const int nrowb = notb ? k : n;
    if ((!lsame(transa, 'N') && !lsame(transa, 'C') && !lsame(transa, 'T')) ||
        (!lsame(transb, 'N') && !lsame(transb, 'C') && !lsame(transb, 'T')) ||
        m < 0 || n < 0 || k < 0 ||
        lda < std::max(1, nrowa) || ldb < std::max(1, nrowb) || ldc < std::max(1, m)) {
        return;
    }
    if ((m == 0) || (n == 0) || (((alpha == 0.0) || (k == 0)) && (beta == 1.0))) {
        return;
    }
    if ((alpha == 0.0) || (k == 0)) {
        if (beta == 0.0) {
            for (int j = 1; j <= n; ++j) {
                for (int i = 1; i <= m; ++i) {
                    C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] = 0.0;
                }
            }
        } else {
            for (int j = 1; j <= n; ++j) {
                for (int i = 1; i <= m; ++i) {
                    C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] =
                        beta * C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))];
                }
            }
        }
        return;
    }
    dgemm_core(nota, notb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, mb, nb, nbt, kb);
}

inline void dgemm_omp(char transa, char transb,
    int m, int n, int k,
    double alpha,
    const double* A, int lda,
    const double* B, int ldb,
    double beta,
    double* C, int ldc,
    int mb, int nb, int nbt, int kb) {
    const bool nota = lsame(transa, 'N');
    const bool notb = lsame(transb, 'N');
    const int nrowa = nota ? m : k;
    const int nrowb = notb ? k : n;
    if ((!lsame(transa, 'N') && !lsame(transa, 'C') && !lsame(transa, 'T')) ||
        (!lsame(transb, 'N') && !lsame(transb, 'C') && !lsame(transb, 'T')) ||
        m < 0 || n < 0 || k < 0 ||
        lda < std::max(1, nrowa) || ldb < std::max(1, nrowb) || ldc < std::max(1, m)) {
        return;
    }
    if ((m == 0) || (n == 0) || (((alpha == 0.0) || (k == 0)) && (beta == 1.0))) {
        return;
    }
    if ((alpha == 0.0) || (k == 0)) {
        if (beta == 0.0) {
#pragma omp parallel for
            for (int j = 1; j <= n; ++j) {
                for (int i = 1; i <= m; ++i) {
                    C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] = 0.0;
                }
            }
        } else {
#pragma omp parallel for
            for (int j = 1; j <= n; ++j) {
                for (int i = 1; i <= m; ++i) {
                    C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] =
                        beta * C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))];
                }
            }
        }
        return;
    }

    const double zero = 0.0;
    const double one = 1.0;
    const std::size_t t1_size = static_cast<std::size_t>(kb) * static_cast<std::size_t>(mb);
    const std::size_t t2_size = static_cast<std::size_t>(kb) * static_cast<std::size_t>(nbt);
    std::unique_ptr<double[]> t1(new double[t1_size]);
    std::unique_ptr<double[]> t2(new double[t2_size]);
    std::fill_n(t1.get(), t1_size, zero);
    std::fill_n(t2.get(), t2_size, zero);

#pragma omp parallel for default(none) schedule(dynamic) shared(A, B, C, t1, t2, m, n, k, lda, ldb, ldc, mb, nb, nbt, kb, alpha, beta, nota, notb)
    for (int jj = 1; jj <= n; jj += nb) {
        const int jsec = std::min(nb, n - jj + 1);
        const int ujsec = jsec - (jsec % 4);
        for (int ll = 1; ll <= k; ll += kb) {
            const int lsec = std::min(kb, k - ll + 1);
            const int ulsec = lsec - (lsec % 2);
            double delta = one;
            if (ll == 1) delta = beta;
            for (int ii = 1; ii <= m; ii += mb) {
                const int isec = std::min(mb, m - ii + 1);
                const int uisec = isec - (isec % 2);
                if (notb) {
                    if (nota) {
                        for (int l = ll; l <= ll + ulsec - 1; l += 2) {
                            for (int i = ii; i <= ii + uisec - 1; i += 2) {
                                const std::size_t row0 = static_cast<std::size_t>(i - 1);
                                const std::size_t row1 = static_cast<std::size_t>(i);
                                const std::size_t col0 = static_cast<std::size_t>(l - 1);
                                const std::size_t col1 = static_cast<std::size_t>(l);
                                const std::size_t t1row0 = static_cast<std::size_t>(l - ll);
                                const std::size_t t1row1 = static_cast<std::size_t>(l - ll + 1);
                                const std::size_t ti0 = static_cast<std::size_t>(i - ii);
                                const std::size_t ti1 = static_cast<std::size_t>(i - ii + 1);
                                t1[f2d_idx(t1row0, ti0, static_cast<std::size_t>(kb))] = alpha * A[f2d_idx(row0, col0, static_cast<std::size_t>(lda))];
                                t1[f2d_idx(t1row1, ti0, static_cast<std::size_t>(kb))] = alpha * A[f2d_idx(row0, col1, static_cast<std::size_t>(lda))];
                                t1[f2d_idx(t1row0, ti1, static_cast<std::size_t>(kb))] = alpha * A[f2d_idx(row1, col0, static_cast<std::size_t>(lda))];
                                t1[f2d_idx(t1row1, ti1, static_cast<std::size_t>(kb))] = alpha * A[f2d_idx(row1, col1, static_cast<std::size_t>(lda))];
                            }
                        }
                    }
                    const int ui4sec = isec - (isec % 4);
                    for (int j = jj; j <= jj + ujsec - 1; j += 4) {
#pragma omp simd
                        for (int i = ii; i <= ii + ui4sec - 1; i += 4) {
                            double f11 = delta * C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))];
                            double f21 = delta * C[f2d_idx(static_cast<std::size_t>(i), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))];
                            double f12 = delta * C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j), static_cast<std::size_t>(ldc))];
                            double f22 = delta * C[f2d_idx(static_cast<std::size_t>(i), static_cast<std::size_t>(j), static_cast<std::size_t>(ldc))];
                            double f13 = delta * C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j + 1), static_cast<std::size_t>(ldc))];
                            double f23 = delta * C[f2d_idx(static_cast<std::size_t>(i), static_cast<std::size_t>(j + 1), static_cast<std::size_t>(ldc))];
                            double f14 = delta * C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j + 2), static_cast<std::size_t>(ldc))];
                            double f24 = delta * C[f2d_idx(static_cast<std::size_t>(i), static_cast<std::size_t>(j + 2), static_cast<std::size_t>(ldc))];
                            double f31 = delta * C[f2d_idx(static_cast<std::size_t>(i + 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))];
                            double f41 = delta * C[f2d_idx(static_cast<std::size_t>(i + 2), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))];
                            double f32 = delta * C[f2d_idx(static_cast<std::size_t>(i + 1), static_cast<std::size_t>(j), static_cast<std::size_t>(ldc))];
                            double f42 = delta * C[f2d_idx(static_cast<std::size_t>(i + 2), static_cast<std::size_t>(j), static_cast<std::size_t>(ldc))];
                            double f33 = delta * C[f2d_idx(static_cast<std::size_t>(i + 1), static_cast<std::size_t>(j + 1), static_cast<std::size_t>(ldc))];
                            double f43 = delta * C[f2d_idx(static_cast<std::size_t>(i + 2), static_cast<std::size_t>(j + 1), static_cast<std::size_t>(ldc))];
                            double f34 = delta * C[f2d_idx(static_cast<std::size_t>(i + 1), static_cast<std::size_t>(j + 2), static_cast<std::size_t>(ldc))];
                            double f44 = delta * C[f2d_idx(static_cast<std::size_t>(i + 2), static_cast<std::size_t>(j + 2), static_cast<std::size_t>(ldc))];
                            for (int l = ll; l <= ll + lsec - 1; ++l) {
                                const std::size_t tl = static_cast<std::size_t>(l - ll);
                                const double b0 = B[f2d_idx(static_cast<std::size_t>(l - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldb))];
                                const double b1 = B[f2d_idx(static_cast<std::size_t>(l - 1), static_cast<std::size_t>(j), static_cast<std::size_t>(ldb))];
                                const double b2 = B[f2d_idx(static_cast<std::size_t>(l - 1), static_cast<std::size_t>(j + 1), static_cast<std::size_t>(ldb))];
                                const double b3 = B[f2d_idx(static_cast<std::size_t>(l - 1), static_cast<std::size_t>(j + 2), static_cast<std::size_t>(ldb))];
                                const double t11 = t1[f2d_idx(tl, static_cast<std::size_t>(i - ii), static_cast<std::size_t>(kb))];
                                const double t21 = t1[f2d_idx(tl, static_cast<std::size_t>(i - ii + 1), static_cast<std::size_t>(kb))];
                                const double t31 = t1[f2d_idx(tl, static_cast<std::size_t>(i - ii + 2), static_cast<std::size_t>(kb))];
                                const double t41 = t1[f2d_idx(tl, static_cast<std::size_t>(i - ii + 3), static_cast<std::size_t>(kb))];
                                f11 += t11 * b0;
                                f21 += t21 * b0;
                                f12 += t11 * b1;
                                f22 += t21 * b1;
                                f13 += t11 * b2;
                                f23 += t21 * b2;
                                f14 += t11 * b3;
                                f24 += t21 * b3;
                                f31 += t31 * b0;
                                f41 += t41 * b0;
                                f32 += t31 * b1;
                                f42 += t41 * b1;
                                f33 += t31 * b2;
                                f43 += t41 * b2;
                                f34 += t31 * b3;
                                f44 += t41 * b3;
                            }
                            C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] = f11;
                            C[f2d_idx(static_cast<std::size_t>(i), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] = f21;
                            C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j), static_cast<std::size_t>(ldc))] = f12;
                            C[f2d_idx(static_cast<std::size_t>(i), static_cast<std::size_t>(j), static_cast<std::size_t>(ldc))] = f22;
                            C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j + 1), static_cast<std::size_t>(ldc))] = f13;
                            C[f2d_idx(static_cast<std::size_t>(i), static_cast<std::size_t>(j + 1), static_cast<std::size_t>(ldc))] = f23;
                            C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j + 2), static_cast<std::size_t>(ldc))] = f14;
                            C[f2d_idx(static_cast<std::size_t>(i), static_cast<std::size_t>(j + 2), static_cast<std::size_t>(ldc))] = f24;
                            C[f2d_idx(static_cast<std::size_t>(i + 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] = f31;
                            C[f2d_idx(static_cast<std::size_t>(i + 2), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] = f41;
                            C[f2d_idx(static_cast<std::size_t>(i + 1), static_cast<std::size_t>(j), static_cast<std::size_t>(ldc))] = f32;
                            C[f2d_idx(static_cast<std::size_t>(i + 2), static_cast<std::size_t>(j), static_cast<std::size_t>(ldc))] = f42;
                            C[f2d_idx(static_cast<std::size_t>(i + 1), static_cast<std::size_t>(j + 1), static_cast<std::size_t>(ldc))] = f33;
                            C[f2d_idx(static_cast<std::size_t>(i + 2), static_cast<std::size_t>(j + 1), static_cast<std::size_t>(ldc))] = f43;
                            C[f2d_idx(static_cast<std::size_t>(i + 1), static_cast<std::size_t>(j + 2), static_cast<std::size_t>(ldc))] = f34;
                            C[f2d_idx(static_cast<std::size_t>(i + 2), static_cast<std::size_t>(j + 2), static_cast<std::size_t>(ldc))] = f44;
                        }
                    }
                }
            }
        }
    }
}

inline void dsymm(char side, char uplo,
    int m, int n,
    double alpha,
    const double* A, int lda,
    const double* B, int ldb,
    double beta,
    double* C, int ldc,
    int rcb, int mb, int nb, int nbt, int kb) {
    const bool lside = lsame(side, 'L');
    const bool upper = lsame(uplo, 'U');
    const int nrowa = lsame(side, 'L') ? m : n;
    if ((!lside && !lsame(side, 'R')) || (!upper && !lsame(uplo, 'L')) ||
        m < 0 || n < 0 || lda < std::max(1, nrowa) ||
        ldb < std::max(1, m) || ldc < std::max(1, m)) {
        return;
    }
    if ((m == 0) || (n == 0) || ((alpha == 0.0) && (beta == 1.0))) {
        return;
    }
    if (alpha == 0.0) {
        if (beta == 0.0) {
            for (int j = 1; j <= n; ++j) {
                for (int i = 1; i <= m; ++i) {
                    C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] = 0.0;
                }
            }
        } else {
            for (int j = 1; j <= n; ++j) {
                for (int i = 1; i <= m; ++i) {
                    C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] *= beta;
                }
            }
        }
        return;
    }
    if (lside) {
        if (upper) {
            for (int ii = 1; ii <= m; ii += rcb) {
                const int isec = std::min(rcb, m - ii + 1);
                std::unique_ptr<double[]> t1(new double[static_cast<std::size_t>(rcb) * static_cast<std::size_t>(rcb)]);
                std::fill_n(t1.get(), static_cast<std::size_t>(rcb) * static_cast<std::size_t>(rcb), 0.0);
                for (int j = ii; j <= ii + isec - 1; ++j) {
                    for (int i = ii; i <= ii + isec - 1; ++i) {
                        const std::size_t row = static_cast<std::size_t>(i - 1);
                        const std::size_t col = static_cast<std::size_t>(j - 1);
                        t1[f2d_idx(static_cast<std::size_t>(i - ii), static_cast<std::size_t>(j - ii), static_cast<std::size_t>(rcb))] = alpha * A[f2d_idx(row, col, static_cast<std::size_t>(lda))];
                    }
                }
                dgemm('T', 'N', isec, n, isec, 1.0, t1.get(), rcb, B + (ii - 1), ldb, beta, C + (ii - 1), ldc, mb, nb, nbt, kb);
            }
        } else {
            for (int ii = m; ii >= 1; ii -= rcb) {
                const int isec = std::min(rcb, ii);
                std::unique_ptr<double[]> t1(new double[static_cast<std::size_t>(rcb) * static_cast<std::size_t>(rcb)]);
                std::fill_n(t1.get(), static_cast<std::size_t>(rcb) * static_cast<std::size_t>(rcb), 0.0);
                for (int j = ii - isec + 1; j <= ii; ++j) {
                    for (int i = ii - isec + 1; i <= ii; ++i) {
                        const std::size_t row = static_cast<std::size_t>(i - 1);
                        const std::size_t col = static_cast<std::size_t>(j - 1);
                        t1[f2d_idx(static_cast<std::size_t>(i - (ii - isec + 1)), static_cast<std::size_t>(j - (ii - isec + 1)), static_cast<std::size_t>(rcb))] = alpha * A[f2d_idx(row, col, static_cast<std::size_t>(lda))];
                    }
                }
                dgemm('T', 'N', isec, n, isec, 1.0, t1.get(), rcb, B + (ii - isec), ldb, beta, C + (ii - isec), ldc, mb, nb, nbt, kb);
            }
        }
    } else {
        if (upper) {
            for (int ii = 1; ii <= n; ii += rcb) {
                const int isec = std::min(rcb, n - ii + 1);
                std::unique_ptr<double[]> t1(new double[static_cast<std::size_t>(rcb) * static_cast<std::size_t>(rcb)]);
                std::fill_n(t1.get(), static_cast<std::size_t>(rcb) * static_cast<std::size_t>(rcb), 0.0);
                for (int j = ii; j <= ii + isec - 1; ++j) {
                    for (int i = ii; i <= ii + isec - 1; ++i) {
                        const std::size_t row = static_cast<std::size_t>(i - 1);
                        const std::size_t col = static_cast<std::size_t>(j - 1);
                        t1[f2d_idx(static_cast<std::size_t>(i - ii), static_cast<std::size_t>(j - ii), static_cast<std::size_t>(rcb))] = alpha * A[f2d_idx(row, col, static_cast<std::size_t>(lda))];
                    }
                }
                dgemm('N', 'N', m, isec, isec, 1.0, B + (ii - 1), ldb, t1.get(), rcb, beta, C + (ii - 1), ldc, mb, nb, nbt, kb);
            }
        } else {
            for (int ii = n; ii >= 1; ii -= rcb) {
                const int isec = std::min(rcb, ii);
                std::unique_ptr<double[]> t1(new double[static_cast<std::size_t>(rcb) * static_cast<std::size_t>(rcb)]);
                std::fill_n(t1.get(), static_cast<std::size_t>(rcb) * static_cast<std::size_t>(rcb), 0.0);
                for (int j = ii - isec + 1; j <= ii; ++j) {
                    for (int i = ii - isec + 1; i <= ii; ++i) {
                        const std::size_t row = static_cast<std::size_t>(i - 1);
                        const std::size_t col = static_cast<std::size_t>(j - 1);
                        t1[f2d_idx(static_cast<std::size_t>(i - (ii - isec + 1)), static_cast<std::size_t>(j - (ii - isec + 1)), static_cast<std::size_t>(rcb))] = alpha * A[f2d_idx(row, col, static_cast<std::size_t>(lda))];
                    }
                }
                dgemm('N', 'N', m, isec, isec, 1.0, B + (ii - isec), ldb, t1.get(), rcb, beta, C + (ii - isec), ldc, mb, nb, nbt, kb);
            }
        }
    }
}

inline void dsyr2k(char uplo, char trans,
    int n, int k,
    double alpha,
    const double* A, int lda,
    const double* B, int ldb,
    double beta,
    double* C, int ldc,
    int rcb, int mb, int nb, int nbt, int kb) {
    const bool upper = lsame(uplo, 'U');
    const bool notr = lsame(trans, 'N');
    const int nrowa = notr ? n : k;
    if ((!upper && !lsame(uplo, 'L')) || ((!notr && !lsame(trans, 'T') && !lsame(trans, 'C'))) ||
        n < 0 || k < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, nrowa) || ldc < std::max(1, n)) {
        return;
    }
    if ((n == 0) || (((alpha == 0.0) || (k == 0)) && (beta == 1.0))) {
        return;
    }
    if ((alpha == 0.0) || (k == 0)) {
        if (upper) {
            if (beta == 0.0) {
                for (int j = 1; j <= n; ++j) {
                    for (int i = 1; i <= j; ++i) {
                        C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] = 0.0;
                    }
                }
            } else {
                for (int j = 1; j <= n; ++j) {
                    for (int i = 1; i <= j; ++i) {
                        C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] *= beta;
                    }
                }
            }
        } else {
            if (beta == 0.0) {
                for (int j = 1; j <= n; ++j) {
                    for (int i = j; i <= n; ++i) {
                        C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] = 0.0;
                    }
                }
            } else {
                for (int j = 1; j <= n; ++j) {
                    for (int i = j; i <= n; ++i) {
                        C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] *= beta;
                    }
                }
            }
        }
        return;
    }
    std::unique_ptr<double[]> t1(new double[static_cast<std::size_t>(rcb) * static_cast<std::size_t>(rcb)]);
    std::fill_n(t1.get(), static_cast<std::size_t>(rcb) * static_cast<std::size_t>(rcb), 0.0);
    for (int ii = 1; ii <= n; ii += rcb) {
        const int isec = std::min(rcb, n - ii + 1);
        if (notr) {
            dgemm('N', 'T', isec, isec, k, alpha, A + (ii - 1), lda, B + (ii - 1), ldb, 0.0, t1.get(), rcb, mb, nb, nbt, kb);
            for (int j = ii + isec - 2; j >= ii - 1; j -= 2) {
                const int ui = j - ii + 1 - ((j - ii + 1) % 2);
                for (int i = ii; i <= ii + ui - 1; i += 2) {
                    const std::size_t row0 = static_cast<std::size_t>(i - ii);
                    const std::size_t row1 = static_cast<std::size_t>(i - ii + 1);
                    const std::size_t col0 = static_cast<std::size_t>(j - ii);
                    const std::size_t col1 = static_cast<std::size_t>(j - ii + 1);
                    C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] =
                        beta * C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] +
                        t1[f2d_idx(row0, col0, static_cast<std::size_t>(rcb))] + t1[f2d_idx(col0, row0, static_cast<std::size_t>(rcb))];
                    C[f2d_idx(static_cast<std::size_t>(i), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] =
                        beta * C[f2d_idx(static_cast<std::size_t>(i), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] +
                        t1[f2d_idx(row1, col0, static_cast<std::size_t>(rcb))] + t1[f2d_idx(col0, row1, static_cast<std::size_t>(rcb))];
                }
            }
        } else {
            dgemm('T', 'N', isec, isec, k, alpha, A + (ii - 1), lda, B + (ii - 1), ldb, 0.0, t1.get(), rcb, mb, nb, nbt, kb);
            for (int j = ii + isec - 2; j >= ii - 1; j -= 2) {
                const int ui = j - ii + 1 - ((j - ii + 1) % 2);
                for (int i = ii; i <= ii + ui - 1; i += 2) {
                    const std::size_t row0 = static_cast<std::size_t>(i - ii);
                    const std::size_t row1 = static_cast<std::size_t>(i - ii + 1);
                    const std::size_t col0 = static_cast<std::size_t>(j - ii);
                    const std::size_t col1 = static_cast<std::size_t>(j - ii + 1);
                    C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] =
                        beta * C[f2d_idx(static_cast<std::size_t>(i - 1), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] +
                        t1[f2d_idx(row0, col0, static_cast<std::size_t>(rcb))] + t1[f2d_idx(col0, row0, static_cast<std::size_t>(rcb))];
                    C[f2d_idx(static_cast<std::size_t>(i), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] =
                        beta * C[f2d_idx(static_cast<std::size_t>(i), static_cast<std::size_t>(j - 1), static_cast<std::size_t>(ldc))] +
                        t1[f2d_idx(row1, col0, static_cast<std::size_t>(rcb))] + t1[f2d_idx(col0, row1, static_cast<std::size_t>(rcb))];
                }
            }
        }
    }
}

} // namespace blas3
} // namespace gms

#endif // GMS_BLAS3_UNROLLED_BLOCKED_DP_HPP
