#ifndef BTGAMMARANGE_H_
#define BTGAMMARANGE_H_

#include <cfloat>
#include <cmath>
#include <limits>

#if defined(__MINGW32__)
#include <_mingw.h>
#endif
#if defined(__GNUC__) && defined(__x86_64__)
#include <cpuid.h>
#include <xmmintrin.h>
#endif

namespace bayestools {
namespace gamma {
namespace range {

// This certificate covers the inspected MinGW-w64 11.x static x87 log route.
// Other backends keep the existing unavailable FAR result. Installed proving
// must also check the actual DLL binding; the UCRT label alone is insufficient.
struct environment {
  bool compiled_profile;
  bool documented_cpu;
  bool hypervisor;
  unsigned short control_word;
  unsigned int mxcsr;
  long double ln2;
  bool available;
};

inline environment current_environment()
{
  environment out = {false, false, false, 0, 0, 0.0L, false};
#if defined(__GNUC__) && !defined(__clang__) && defined(__MINGW32__) && \
    defined(__x86_64__) && defined(__SSE2__) && \
    defined(__MINGW64_VERSION_MAJOR) && __MINGW64_VERSION_MAJOR == 11 && \
    !defined(__FAST_MATH__) && __FINITE_MATH_ONLY__ == 0 && FLT_EVAL_METHOD == 0
  out.compiled_profile = DBL_MANT_DIG == 53 && DBL_MAX_EXP == 1024 &&
    LDBL_MANT_DIG == 64 && LDBL_MAX_EXP == 16384 &&
    std::numeric_limits<double>::is_iec559;
  unsigned int eax, ebx, ecx, edx;
  if(__get_cpuid(0, &eax, &ebx, &ecx, &edx)){
    out.documented_cpu =
      (ebx == 0x756e6547 && edx == 0x49656e69 && ecx == 0x6c65746e) ||
      (ebx == 0x68747541 && edx == 0x69746e65 && ecx == 0x444d4163);
  }
  if(__get_cpuid(1, &eax, &ebx, &ecx, &edx)) out.hypervisor = (ecx & (1U << 31)) != 0;
  __asm__ __volatile__("fnstcw %0" : "=m" (out.control_word));
  out.mxcsr = _mm_getcsr();
  __asm__ __volatile__("fldln2; fstpt %0" : "=m" (out.ln2));
  const unsigned int precision = out.control_word & 0x0300;
  // Exact extended limits were checked against Arb ln(2): every value in this
  // interval differs from ln(2) by less than 2^-60. Compare in extended precision.
  const long double lower = 0xb.17217f7d1cf79a4p-4L;
  const long double upper = 0xb.17217f7d1cf79b4p-4L;
  out.available = out.compiled_profile && out.documented_cpu &&
    (out.control_word & 0x0c00) == 0 && (precision == 0x0200 || precision == 0x0300) &&
    (out.mxcsr & (0x6000 | 0x8000 | 0x0040)) == 0 &&
    out.ln2 >= lower && out.ln2 <= upper &&
    std::log(2.0) == 0x1.62e42fefa39efp-1;
#endif
  return out;
}

const double far_threshold = 0x1.6447141f9342ap+9;
const double log_error = 0x1p-42;
const double twice_unit_roundoff = 0x1p-52;

inline double unavailable()
{
  return std::numeric_limits<double>::quiet_NaN();
}

inline double add(double x, double y)
{
  if(!(std::isfinite(x) && x >= 0.0 && std::isfinite(y) && y >= 0.0)) return unavailable();
  double value = x + y;
  return std::isfinite(value) ?
    std::nextafter(value, std::numeric_limits<double>::infinity()) : unavailable();
}

inline double multiply(double x, double y)
{
  if(!(std::isfinite(x) && x >= 0.0 && std::isfinite(y) && y >= 0.0)) return unavailable();
  double value = x * y;
  return std::isfinite(value) ?
    std::nextafter(value, std::numeric_limits<double>::infinity()) : unavailable();
}

inline double rounding_error(double value)
{
  return add(multiply(twice_unit_roundoff, std::fabs(value)),
    std::numeric_limits<double>::denorm_min());
}

inline double invgamma_error(double h)
{
  return add(add(log_error, log_error), rounding_error(h));
}

inline double nonlocal_error(double order, double m, double h, bool inverse)
{
  double distance_error = add(log_error, twice_unit_roundoff);
  double operand_error = add(log_error, multiply(2.0, distance_error));
  if(!inverse) operand_error = add(operand_error, log_error);
  double inner = add(operand_error, rounding_error(m));
  return add(inverse ? multiply(order, inner) : inner, rounding_error(h));
}

inline bool certified(double shape, double h, double error)
{
  // A<2a (A=a for inverse-Gamma), E<1/4, h>T and finite RN(a*h)
  // imply A*ell<3M and r>8M. Prefix<-5M+3 and logQ<-5M+4;
  // finite supported Jacobians/normalizers cannot recover finite log range.
  return std::isfinite(shape) && shape > 0.0 && std::isfinite(h) &&
    h > far_threshold && std::isfinite(error) && error >= 0.0 && error < 0.25 &&
    std::isfinite(shape * h);
}

}
}
}
#endif
