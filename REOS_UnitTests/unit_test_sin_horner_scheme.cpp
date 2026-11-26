#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <type_traits>

template <typename T>
constexpr T sin_horner_scheme(T x) noexcept {
    static_assert(std::is_floating_point_v<T>, "fast_sin requires a floating-point type");

    constexpr T two_over_pi = static_cast<T>(0.6366197723675813430755350534900574);
    constexpr T pi_over_2_hi = static_cast<T>(1.57079637050628662109375);
    constexpr T pi_over_2_lo = static_cast<T>(-4.37113900018624283e-8);

    T scaled = x * two_over_pi;
    auto k = static_cast<std::int64_t>(scaled + (scaled >= T(0) ? T(0.5) : T(-0.5)));
    T kf = static_cast<T>(k);

    T r = (x - kf * pi_over_2_hi) - kf * pi_over_2_lo;
    T r2 = r * r;

    T result{};
    if ((k & 1) == 0) {
        constexpr T c3 = static_cast<T>(-0.166666666666666324348);
        constexpr T c5 = static_cast<T>( 0.00833333333332248946124);
        constexpr T c7 = static_cast<T>(-0.000198412698298579493134);
        result = r * (T(1) + r2 * (c3 + r2 * (c5 + r2 * c7)));
    } else {
        constexpr T c2 = static_cast<T>(-0.499999999999999944489);
        constexpr T c4 = static_cast<T>( 0.0416666666666666019037);
        constexpr T c6 = static_cast<T>(-0.00138888888888741095749);
        result = T(1) + r2 * (c2 + r2 * (c4 + r2 * c6));
    }

    return (k & 2) ? -result : result;
}

// Example of unit tests covering major implementation areas.

// Compile-time verification
constexpr double ct_zero = sin_horner_scheme(0.0);
constexpr float ct_half_pi = sin_horner_scheme(1.5707963f);
static_assert(ct_zero == 0.0, "Compile-time fast_sin(0.0) failed");
static_assert(ct_half_pi > 0.999999f && ct_half_pi <= 1.000001f, "Compile-time fast_sin(pi/2) failed");

template <typename T>
void run_tests(T max_allowed_error, const char* type_name) {
    constexpr T pi = static_cast<T>(3.1415926535897932384626433832795);

    //Exact & Cardinal Angles
    assert(sin_horner_scheme(T(0)) == T(0));
    assert(sin_horner_scheme(T(-0.0)) == T(0));
    assert(std::abs(sin_horner_scheme(pi / T(2)) - T(1)) < max_allowed_error);
    assert(std::abs(sin_horner_scheme(-pi / T(2)) - T(-1)) < max_allowed_error);
    assert(std::abs(sin_horner_scheme(pi) - T(0)) < max_allowed_error);
    assert(std::abs(sin_horner_scheme(T(2) * pi) - T(0)) < max_allowed_error);

    //Odd Symmetry: sin(-x) == -sin(x)
    for (T x = T(0.1); x < T(100.0); x += T(1.37)) {
        assert(std::abs(sin_horner_scheme(-x) + sin_horner_scheme(x)) < std::numeric_limits<T>::epsilon() * T(4));
    }

    //Primary Range Sweep [-2*pi, 2*pi]
    T max_err_primary = T(0);
    for (T x = -T(2) * pi; x <= T(2) * pi; x += T(0.001)) {
        T err = std::abs(sin_horner_scheme(x) - std::sin(x));
        if (err > max_err_primary) max_err_primary = err;
        assert(err < max_allowed_error);
    }

    //Large Range Sweep [-10,000, 10,000] radians
    T max_err_large = T(0);
    for (T x = T(-10000); x <= T(10000); x += T(0.73)) {
        T err = std::abs(sin_horner_scheme(x) - std::sin(x));
        if (err > max_err_large) max_err_large = err;
        assert(err < max_allowed_error * T(10)); // Slight tolerance growth for large k
    }

    std::cout << "[" << type_name << "] Passed! "
              << "Max Primary Err: " << max_err_primary
              << " | Max Large-Range Err: " << max_err_large << '\n';
}

int main() 
{
    run_tests<float>(1e-6f, "float");
    run_tests<double>(1e-10, "double");
    std::cout << "All unit tests passed successfully.\n";
    return 0;
}