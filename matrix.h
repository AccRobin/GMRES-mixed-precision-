#include <iostream>
#include <vector>
#include <cmath>
template<typename T, const int N>
struct Mtx {
    T a[N][N];
    Mtx() : a() {}
    template<typename P>
    Mtx(P _a[][N]) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) a[i][j] = _a[i][j];
        }
    }
    template<typename P>
    Mtx(const Mtx<P, N>& _a) {
        *this = Mtx(_a.a);
    }
    T* operator[](int x) {
        return a[x];
    }
    const T* operator[](int x) const {
        return a[x];
    }
};

template<typename T, const int N>
struct Vec {
    T a[N];
    Vec() : a() {}
    template<typename P>
    Vec(P* _a) {
        for (int i = 0; i < N; ++i) a[i] = _a[i];
    }
    template<typename P>
    Vec(const Vec<P, N>& _a) : a(_a.a) {}
    T& operator[](int x) {
        return a[x];
    }
    const T& operator[](int x) const {
        return a[x];
    }
    Vec<T, N> operator-() const {
        Vec<T, N> r;
        for (int i = 0; i < N; ++i) r[i] = -a[i];
        return r;
    }
    Vec<T, N> operator+(const Vec<T, N>& b) const {
        Vec<T, N> r;
        for (int i = 0; i < N; ++i) r[i] = a[i] + b[i];
        return r;
    }
    Vec<T, N> operator-(const Vec<T, N>& b) const {
        return *this + -b;
    }
    T nrm2() const ;
};
template<typename T, const int N>
T dot(const Vec<T, N>& a, const Vec<T, N>& b){
    T res = 0;
    for (int i = 0; i < N; ++i) res += a[i] * b[i];
    return res;
}
template<typename T, const int N>
T Vec<T, N>::nrm2() const {
    return sqrt(dot(*this, *this));
}
template<typename P, typename T, const int N>
Vec<T, N> operator*(const P& a, const Vec<T, N>& x) {
    Vec<T, N> r = x;
    for (int i = 0; i < N; ++i) r[i] = a * r[i];
    return r;
}
template<typename T, const int N>
Vec<T, N> operator*(const Mtx<T, N>& A, const Vec<T, N>& x) {
    Vec<T, N> r;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            r[i] += A[i][j] * x[j];
        }
    }
    return r;
}

template<typename T, const int N>
void givens_rotation(Vec<T, N>* H, std::vector<std::pair<T, T>>& givens, int i) {
    for (int k = 0; k < i; ++k) {
        T temp = givens[k].first * H[k][i] + givens[k].second * H[k+1][i];
        H[k+1][i] = -givens[k].second * H[k][i] + givens[k].first * H[k+1][i];
        H[k][i] = temp;
    }
    
    T h1 = H[i][i];
    T h2 = H[i+1][i];
    T r = sqrt(h1 * h1 + h2 * h2);
    if (r == 0) {
        throw "Zero divisor in Givens rotation";
    }
    givens[i].first = h1 / r;
    givens[i].second = h2 / r;
    
    H[i][i] = r;
    H[i+1][i] = 0;
}

template<typename T, const int N>
std::ostream& operator<<(std::ostream& out, const Vec<T, N>& x) {
    out << '[';
    for (int i = 0; i < N; ++i) out << x[i] << ",]"[i + 1 == N];
    return out;
}