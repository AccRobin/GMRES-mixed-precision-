#include <cmath>

using bf16 = __bf16;
using fp16 = _Float16;
using fp32 = float;
using fp64 = double;

fp16 abs(const fp16& x) {
    return x < 0 ? -x : x;
}
bf16 abs(const bf16& x) {
    return x < 0 ? -x : x;
}
template<typename T>
T abs(const T& x) {
    return std::abs(x);
}

fp16 sqrt(const fp16& x) {
    return (fp16)std::sqrt((fp32)x);
}
bf16 sqrt(const bf16& x) {
    return (bf16)std::sqrt((fp32)x);
}