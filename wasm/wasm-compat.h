// Adds whatever libc/libc++ functionality is needed by either reimplenting it,
// or importing it from javascript.
#pragma once

#define COMPAT_INLINE inline __attribute__((always_inline))

//-----------------------------------------------------------------------------
// stdint.h
typedef unsigned long size_t;
typedef unsigned long uintptr_t;
typedef long ptrdiff_t;
typedef long ssize_t;
typedef long intptr_t;

typedef signed char        int8_t;
typedef short              int16_t;
typedef int                int32_t;
typedef long long          int64_t;
typedef long long          intmax_t;
typedef unsigned char      uint8_t;
typedef unsigned short     uint16_t;
typedef unsigned int       uint32_t;
typedef unsigned long long uint64_t;
typedef unsigned long long u_int64_t;
typedef unsigned long long uintmax_t;

//-----------------------------------------------------------------------------
// bits/alltypes.h
typedef float float_t;
typedef double double_t;

//-----------------------------------------------------------------------------
// float.h
const float FLT_MIN = 1.17549435082228750797e-38F;
const float FLT_MAX = 3.40282346638528859812e+38F;
const double DBL_MIN = 2.22507385850720138309e-308;
const double DBL_MAX = 1.79769313486231570815e+308;

//-----------------------------------------------------------------------------
// stdarg.h
typedef __builtin_va_list va_list;
#define va_start(v,l)   __builtin_va_start(v,l)
#define va_end(v)       __builtin_va_end(v)
#define va_arg(v,l)     __builtin_va_arg(v,l)
#define va_copy(d,s)    __builtin_va_copy(d,s)

//-----------------------------------------------------------------------------
// stdlib.h
extern "C" void abort(); // Javascript must define this

//-----------------------------------------------------------------------------
// memcpy.c
extern "C" {
void *memset(void *dest, int c, size_t n);
void *memcpy(void *__restrict dest, const void *__restrict src, size_t n);
}

//-----------------------------------------------------------------------------
// strnlen.c
void *memchr(const void *src, int c, size_t n);
COMPAT_INLINE size_t strnlen(const char *s, size_t n) {
	const char *p = (const char *)memchr(s, 0, n);
	return p ? p-s : n;
}

//-----------------------------------------------------------------------------
// assert.h
void assert_impl(const char* expression, const char* file, int line);
#define assert(expression) (void)((!!(expression)) || (assert_impl((#expression), __FILE__, __LINE__), 0))

//-----------------------------------------------------------------------------
// stdio.h
extern "C" void puts(const char*); // Javascript must define this

//-----------------------------------------------------------------------------
#include "printf.h"

//-----------------------------------------------------------------------------
// math.h
#define M_E             2.7182818284590452354   /* e */
#define M_LOG2E         1.4426950408889634074   /* log_2 e */
#define M_LOG10E        0.43429448190325182765  /* log_10 e */
#define M_LN2           0.69314718055994530942  /* log_e 2 */
#define M_LN10          2.30258509299404568402  /* log_e 10 */
#define M_PI            3.14159265358979323846  /* pi */
#define M_PI_2          1.57079632679489661923  /* pi/2 */
#define M_PI_4          0.78539816339744830962  /* pi/4 */
#define M_1_PI          0.31830988618379067154  /* 1/pi */
#define M_2_PI          0.63661977236758134308  /* 2/pi */
#define M_2_SQRTPI      1.12837916709551257390  /* 2/sqrt(pi) */
#define M_SQRT2         1.41421356237309504880  /* sqrt(2) */
#define M_SQRT1_2       0.70710678118654752440  /* 1/sqrt(2) */

#define isnan(x) ( \
	sizeof(x) == sizeof(float) ? (__FLOAT_BITS(x) & 0x7fffffff) > 0x7f800000 : (__DOUBLE_BITS(x) & -1ULL>>1) > 0x7ffULL<<52)

// Use builtins where we can:
COMPAT_INLINE int abs(int x) { return __builtin_abs(x); }
COMPAT_INLINE float fabsf(float x) { return __builtin_fabsf(x); }
COMPAT_INLINE double fabs(double x) { return __builtin_fabs(x); }
COMPAT_INLINE float ceilf(float x) { return __builtin_ceilf(x); }
COMPAT_INLINE double ceil(double x) { return __builtin_ceil(x); }
COMPAT_INLINE float floorf(float x) { return __builtin_floorf(x); }
COMPAT_INLINE double floor(double x) { return __builtin_floor(x); }
COMPAT_INLINE float fmaxf(float x, float y) { return x > y ? x : y; } // fast-math-style (ignores NaN and signed 0)
COMPAT_INLINE double fmax(double x, double y) { return x > y ? x : y; }
COMPAT_INLINE float fminf(float x, float y) { return x < y ? x : y; }
COMPAT_INLINE double fmin(double x, double y) { return x < y ? x : y; }
COMPAT_INLINE float sqrtf(float x) { return __builtin_sqrtf(x); }
COMPAT_INLINE double sqrt(double x) { return __builtin_sqrt(x); }
COMPAT_INLINE float truncf(float x) { return __builtin_truncf(x); }
COMPAT_INLINE double trunc(double x) { return __builtin_trunc(x); }

// Expect javascript to provide these
extern "C" double acos(double x);
extern "C" double acosh(double x);
extern "C" double asin(double x);
extern "C" double asinh(double x);
extern "C" double atan(double x);
extern "C" double atanh(double x);
extern "C" double atan2(double y, double x);
extern "C" double cos(double x);
extern "C" double cosh(double x);
extern "C" double exp(double x);
extern "C" double exp2(double x);
extern "C" double log(double x);
extern "C" double log2(double x);
extern "C" double pow(double x, double y);
extern "C" double round(double x);
extern "C" double sin(double x);
extern "C" double sinh(double x);
extern "C" double tan(double x);
extern "C" double tanh(double x);

extern "C" {
// Provide float versions/overloads for everything
COMPAT_INLINE float acosf(float x) { return (float)acos((double)x); }
COMPAT_INLINE float acoshf(float x) { return (float)acosh((double)x); }
COMPAT_INLINE float asinf(float x) { return (float)asin((double)x); }
COMPAT_INLINE float asinhf(float x) { return (float)asinh((double)x); }
COMPAT_INLINE float atanf(float x) { return (float)atan((double)x); }
COMPAT_INLINE float atanhf(float x) { return (float)atanh((double)x); }
COMPAT_INLINE float atan2f(float y, float x) { return (float)atan2((double)x, (double)y); }
COMPAT_INLINE float cosf(float x) { return (float)cos((double)x); }
COMPAT_INLINE float coshf(float x) { return (float)cosh((double)x); }
COMPAT_INLINE float expf(float x) { return (float)exp((double)x); }
COMPAT_INLINE float exp2f(float x) { return (float)exp2((double)x); }
COMPAT_INLINE float fmaf(float x, float y, float z) { return (float)((double)x * (double)y + (double)z); } // Close enough?
COMPAT_INLINE float logf(float x) { return (float)log((double)x); }
COMPAT_INLINE float log2f(float x) { return (float)log2((double)x); }
float modff(float x, float *iptr); // Javscript doesn't provide this, so we use the musl version
COMPAT_INLINE float powf(float x, float y) { return (float)pow((double)x, (double)y); }
COMPAT_INLINE float roundf(float x) { return (float)round((double)x); }
COMPAT_INLINE float sinf(float x) { return (float)sin((double)x); }
COMPAT_INLINE float sinhf(float x) { return (float)sinh((double)x); }
COMPAT_INLINE float tanf(float x) { return (float)tan((double)x); }
COMPAT_INLINE float tanhf(float x) { return (float)tanh((double)x); }
}

COMPAT_INLINE float abs(float x) { return fabsf(x); }
COMPAT_INLINE double abs(double x) { return fabs(x); }
COMPAT_INLINE float acos(float x) { return acosf(x); }
COMPAT_INLINE float acosh(float x) { return acoshf(x); }
COMPAT_INLINE float asin(float x) { return asinf(x); }
COMPAT_INLINE float asinh(float x) { return asinhf(x); }
COMPAT_INLINE float atan(float x) { return atanf(x); }
COMPAT_INLINE float atanh(float x) { return atanhf(x); }
COMPAT_INLINE float atan2(float y, float x) { return atan2f(y, x); }
COMPAT_INLINE float ceil(float x) { return ceilf(x); }
COMPAT_INLINE float cos(float x) { return cosf(x); }
COMPAT_INLINE float cosh(float x) { return coshf(x); }
COMPAT_INLINE float exp(float x) { return expf(x); }
COMPAT_INLINE float exp2(float x) { return exp2f(x); }
COMPAT_INLINE float floor(float x) { return floorf(x); }
COMPAT_INLINE float fma(float a, float b, float c) { return fmaf(a, b, c); }
COMPAT_INLINE float log(float x) { return logf(x); }
COMPAT_INLINE float log2(float x) { return log2f(x); }
//COMPAT_INLINE float max(float x, float y) { return fmaxf(x, y); }
//COMPAT_INLINE double max(double x, double y) { return fmax(x, y); }
//COMPAT_INLINE float min(float x, float y) { return fminf(x, y); }
//COMPAT_INLINE double min(double x, double y) { return fmin(x, y); }
COMPAT_INLINE float modf(float x, float* i) { return modff(x, i); }
COMPAT_INLINE float pow(float x, float y) { return powf(x, y); }
COMPAT_INLINE float round(float x) { return roundf(x); }
COMPAT_INLINE float sin(float x) { return sinf(x); }
COMPAT_INLINE float sinh(float x) { return sinhf(x); }
COMPAT_INLINE float sqrt(float x) { return sqrtf(x); }
COMPAT_INLINE float tan(float x) { return tanf(x); }
COMPAT_INLINE float tanh(float x) { return tanhf(x); }
COMPAT_INLINE float trunc(float x) { return truncf(x); }

//-----------------------------------------------------------------------------
// <new>
inline void* operator new  (size_t, void* __p) { return __p; }

//-----------------------------------------------------------------------------
// Custom add-ons
extern "C" void debugln(const char* cstr, size_t length);
extern "C" double performance_now();
