// #include <bits/stdc++.h>
// using namespace std;
#include <vector>
#include <cstring>
#include <iostream>
#include "matrix.h"
using bf16 = __bf16;
using fp16 = _Float16;
using fp32 = float;
using fp64 = double;
const int RESTART = 20;
using std::cerr;

template<const int N>
struct GMRES {
    Mtx<bf16, N> A_bf16;
    Mtx<fp16, N> A_fp16;
    Mtx<fp32, N> A_fp32;
    Mtx<fp64, N> A_fp64;
    Vec<fp64, N> b;
    GMRES(Mtx<fp64, N> _A, const fp64* _b) :
    A_bf16(_A), A_fp16(_A), A_fp32(_A), A_fp64(_A), b(_b) {}
    template<typename T>
    Vec<fp64, N> work(Vec<fp64, N> x, fp64 eps, int IterMax) {
        Mtx<T, N> A;
        if constexpr (std::is_same_v<T, bf16>) {
            A = A_bf16;
        } else if constexpr (std::is_same_v<T, fp16>) {
            A = A_fp16;
        } else if constexpr (std::is_same_v<T, fp32>) {
            A = A_fp32;
        } else if constexpr (std::is_same_v<T, fp64>) {
            A = A_fp64;
        } else {
            static_assert(false, "type of T error!");
        }

        Vec<fp64, N> r = b - A * x;
        fp64 beta = r.nrm2(), tol = eps * b.nrm2();
        if (beta < tol) return x;

        Vec<T, N> V[RESTART + 1];
        Vec<T, RESTART> H[RESTART + 1];
        for (int iter = 0; iter * RESTART < IterMax; ++iter) {
            memset(V, 0, sizeof V);
            memset(H, 0, sizeof H);
            V[0] = 1 / beta * r;
            
            Vec<T, RESTART + 1> s;
            s[0] = beta;

            std::vector<std::pair<T, T>> givens(RESTART);
            int j = 0;
            for (; j < RESTART; ++j) {
                Vec<T, N> w = A * V[j];
                for (int k = 0; k <= j; ++k) {
                    H[k][j] = dot(w, V[k]); //?
                    w = w - H[k][j] * V[k];
                }
                H[j + 1][j] = w.nrm2();
                if (H[j + 1][j] < tol) {
                    break;
                }
                V[j + 1] = 1 / H[j + 1][j] * w;

                givens_rotation(H, givens, j);

                s[j + 1] = -givens[j].second * s[j];
                s[j] = givens[j].first * s[j];

                T res = std::abs(s[j + 1]);
                if (res < tol) break;
            }

            std::vector<T> y(j + 1, 0);
            for (int i = j; i >= 0; --i) {
                y[i] = s[i];
                for (int k = i + 1; k <= j; ++k) {
                    y[i] -= H[i][k] * y[k];
                }
                y[i] /= H[i][i];
            }

            for (int k = 0; k <= j; ++k) {
                x = x + y[k] * V[k];
            }

            r = b - A * x;
            beta = r.nrm2();
            if (beta < tol) {
                return x;
            }
        }
        return x;
    }
};
int main() {
    const int N = 6;
    fp64 A[N][N] = {
        {4, -1, 0, -1, 0, 0},
        {-1, 4, -1, 0, -1, 0},
        {0, -1, 4, 0, 0, -1},
        {-1, 0, 0, 4, -1, 0},
        {0, -1, 0, -1, 4, -1},
        {0, 0, -1, 0, -1, 4}
    };
    fp64 b[N] = {0, 5, 0, 6, -2, 6};
    fp64 x[N] = {0, 0, 0, 0, 0, 0};
    GMRES<N> proc(Mtx<fp64, N>(A), b);
    auto ans = proc.work<fp64>(Vec<fp64, N>(x), 1e-6, 1000);
    std::cout << ans << '\n';
}