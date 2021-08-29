//-----------------------------------------------------------------------------
// Highly GLSL-compatible vectormath, including swizzling.
//-----------------------------------------------------------------------------
#pragma once

#include <stddef.h>
#include <stdint.h>
#include <cmath>
inline double max(double a, double b) { return fmax(a, b); }
inline double min(double a, double b) { return fmin(a, b); }
inline float max(float a, float b) { return fmaxf(a, b); }
inline float min(float a, float b) { return fminf(a, b); }

class vec2;
class vec3;
class vec4;
class mat2;
class mat3;
class mat4;

//-----------------------------------------------------------------------------
// swiz struct templates (allow vector fields to be accessed out of order)
//-----------------------------------------------------------------------------
template <uint32_t i, uint32_t j> struct swiz2 {
	swiz2<i, j>& operator = (const swiz2<i, j>& other) { ((float*)this)[i] = ((float*)&other)[i]; ((float*)this)[j] = ((float*)&other)[j]; return *this; }
	swiz2<i, j>& operator = (const vec2& other) { ((float*)this)[i] = ((float*)&other)[0]; ((float*)this)[j] = ((float*)&other)[1]; return *this; }
};
// Big reduction in PCH size by keeping these operators external
template <uint32_t i, uint32_t j> swiz2<i, j>& operator += (swiz2<i, j>& a, const vec2& b) { ((float*)&a)[i] += ((float*)&b)[0]; ((float*)&a)[j] += ((float*)&b)[1]; return a; }
template <uint32_t i, uint32_t j> swiz2<i, j>& operator -= (swiz2<i, j>& a, const vec2& b) { ((float*)&a)[i] -= ((float*)&b)[0]; ((float*)&a)[j] -= ((float*)&b)[1]; return a; }
template <uint32_t i, uint32_t j> swiz2<i, j>& operator *= (swiz2<i, j>& a, const vec2& b) { ((float*)&a)[i] *= ((float*)&b)[0]; ((float*)&a)[j] *= ((float*)&b)[1]; return a; }
template <uint32_t i, uint32_t j> swiz2<i, j>& operator /= (swiz2<i, j>& a, const vec2& b) { ((float*)&a)[i] /= ((float*)&b)[0]; ((float*)&a)[j] /= ((float*)&b)[1]; return a; }
template <uint32_t i, uint32_t j> swiz2<i, j>& operator += (swiz2<i, j>& a, float b) { ((float*)&a)[i] += b; ((float*)&a)[j] += b; return a; }
template <uint32_t i, uint32_t j> swiz2<i, j>& operator -= (swiz2<i, j>& a, float b) { ((float*)&a)[i] -= b; ((float*)&a)[j] -= b; return a; }
template <uint32_t i, uint32_t j> swiz2<i, j>& operator *= (swiz2<i, j>& a, float b) { ((float*)&a)[i] *= b; ((float*)&a)[j] *= b; return a; }
template <uint32_t i, uint32_t j> swiz2<i, j>& operator /= (swiz2<i, j>& a, float b) { ((float*)&a)[i] /= b; ((float*)&a)[j] /= b; return a; }

template <uint32_t i, uint32_t j, uint32_t k> struct swiz3 {
	swiz3<i, j, k>& operator = (const swiz3<i, j, k>& other) { ((float*)this)[i] = ((float*)&other)[i]; ((float*)this)[j] = ((float*)&other)[j]; ((float*)this)[k] = ((float*)&other)[k]; return *this; }
	swiz3<i, j, k>& operator = (const vec3& other) { ((float*)this)[i] = ((float*)&other)[0]; ((float*)this)[j] = ((float*)&other)[1]; ((float*)this)[k] = ((float*)&other)[2]; return *this; }
};
template <uint32_t i, uint32_t j, uint32_t k> swiz3<i, j, k>& operator += (swiz3<i, j, k>& a, const vec3& b) { ((float*)&a)[i] += ((float*)&b)[0]; ((float*)&a)[j] += ((float*)&b)[1]; ((float*)&a)[k] += ((float*)&b)[2]; return a; }
template <uint32_t i, uint32_t j, uint32_t k> swiz3<i, j, k>& operator -= (swiz3<i, j, k>& a, const vec3& b) { ((float*)&a)[i] -= ((float*)&b)[0]; ((float*)&a)[j] -= ((float*)&b)[1]; ((float*)&a)[k] -= ((float*)&b)[2]; return a; }
template <uint32_t i, uint32_t j, uint32_t k> swiz3<i, j, k>& operator *= (swiz3<i, j, k>& a, const vec3& b) { ((float*)&a)[i] *= ((float*)&b)[0]; ((float*)&a)[j] *= ((float*)&b)[1]; ((float*)&a)[k] *= ((float*)&b)[2]; return a; }
template <uint32_t i, uint32_t j, uint32_t k> swiz3<i, j, k>& operator /= (swiz3<i, j, k>& a, const vec3& b) { ((float*)&a)[i] /= ((float*)&b)[0]; ((float*)&a)[j] /= ((float*)&b)[1]; ((float*)&a)[k] /= ((float*)&b)[2]; return a; }
template <uint32_t i, uint32_t j, uint32_t k> swiz3<i, j, k>& operator += (swiz3<i, j, k>& a, float b) { ((float*)&a)[i] += b; ((float*)&a)[j] += b; ((float*)&a)[k] += b; return a; }
template <uint32_t i, uint32_t j, uint32_t k> swiz3<i, j, k>& operator -= (swiz3<i, j, k>& a, float b) { ((float*)&a)[i] -= b; ((float*)&a)[j] -= b; ((float*)&a)[k] -= b; return a; }
template <uint32_t i, uint32_t j, uint32_t k> swiz3<i, j, k>& operator *= (swiz3<i, j, k>& a, float b) { ((float*)&a)[i] *= b; ((float*)&a)[j] *= b; ((float*)&a)[k] *= b; return a; }
template <uint32_t i, uint32_t j, uint32_t k> swiz3<i, j, k>& operator /= (swiz3<i, j, k>& a, float b) { ((float*)&a)[i] /= b; ((float*)&a)[j] /= b; ((float*)&a)[k] /= b; return a; }

template <uint32_t i, uint32_t j, uint32_t k, uint32_t l> struct swiz4 {
	swiz4<i, j, k, l>& operator = (const swiz4<i, j, k, l>& other) { ((float*)this)[i] = ((float*)&other)[i]; ((float*)this)[j] = ((float*)&other)[j]; ((float*)this)[k] = ((float*)&other)[k]; ((float*)this)[l] = ((float*)&other)[l]; return *this; }
	swiz4<i, j, k, l>& operator = (const vec4& other) { ((float*)this)[i] = ((float*)&other)[0]; ((float*)this)[j] = ((float*)&other)[1]; ((float*)this)[k] = ((float*)&other)[2]; ((float*)this)[l] = ((float*)&other)[3]; return *this; }
};
template <uint32_t i, uint32_t j, uint32_t k, uint32_t l> swiz4<i, j, k, l>& operator += (swiz4<i, j, k, l>& a, const vec4& b) { ((float*)&a)[i] += ((float*)&b)[0]; ((float*)&a)[j] += ((float*)&b)[1]; ((float*)&a)[k] += ((float*)&b)[2]; ((float*)&a)[l] += ((float*)&b)[3]; return a; }
template <uint32_t i, uint32_t j, uint32_t k, uint32_t l> swiz4<i, j, k, l>& operator -= (swiz4<i, j, k, l>& a, const vec4& b) { ((float*)&a)[i] -= ((float*)&b)[0]; ((float*)&a)[j] -= ((float*)&b)[1]; ((float*)&a)[k] -= ((float*)&b)[2]; ((float*)&a)[l] -= ((float*)&b)[3]; return a; }
template <uint32_t i, uint32_t j, uint32_t k, uint32_t l> swiz4<i, j, k, l>& operator *= (swiz4<i, j, k, l>& a, const vec4& b) { ((float*)&a)[i] *= ((float*)&b)[0]; ((float*)&a)[j] *= ((float*)&b)[1]; ((float*)&a)[k] *= ((float*)&b)[2]; ((float*)&a)[l] *= ((float*)&b)[3]; return a; }
template <uint32_t i, uint32_t j, uint32_t k, uint32_t l> swiz4<i, j, k, l>& operator /= (swiz4<i, j, k, l>& a, const vec4& b) { ((float*)&a)[i] /= ((float*)&b)[0]; ((float*)&a)[j] /= ((float*)&b)[1]; ((float*)&a)[k] /= ((float*)&b)[2]; ((float*)&a)[l] /= ((float*)&b)[3]; return a; }
template <uint32_t i, uint32_t j, uint32_t k, uint32_t l> swiz4<i, j, k, l>& operator += (swiz4<i, j, k, l>& a, float b) { ((float*)&a)[i] += b; ((float*)&a)[j] += b; ((float*)&a)[k] += b; ((float*)&a)[l] += b; return a; }
template <uint32_t i, uint32_t j, uint32_t k, uint32_t l> swiz4<i, j, k, l>& operator -= (swiz4<i, j, k, l>& a, float b) { ((float*)&a)[i] -= b; ((float*)&a)[j] -= b; ((float*)&a)[k] -= b; ((float*)&a)[l] -= b; return a; }
template <uint32_t i, uint32_t j, uint32_t k, uint32_t l> swiz4<i, j, k, l>& operator *= (swiz4<i, j, k, l>& a, float b) { ((float*)&a)[i] *= b; ((float*)&a)[j] *= b; ((float*)&a)[k] *= b; ((float*)&a)[l] *= b; return a; }
template <uint32_t i, uint32_t j, uint32_t k, uint32_t l> swiz4<i, j, k, l>& operator /= (swiz4<i, j, k, l>& a, float b) { ((float*)&a)[i] /= b; ((float*)&a)[j] /= b; ((float*)&a)[k] /= b; ((float*)&a)[l] /= b; return a; }

// Macros for generating all the different permutations of swizzles
#define LEQ12(a) a
#define LEQ13(a) a
#define LEQ14(a) a
#define LEQ22(a) a
#define LEQ23(a) a
#define LEQ24(a) a
#define LEQ32(a) 
#define LEQ33(a) a
#define LEQ34(a) a
#define LEQ42(a) 
#define LEQ43(a) 
#define LEQ44(a) a

#define SWIZ2_H2(n, I,X,R,S, J,Y,G,T, ...) swiz2<I, J> X##Y, R##G, S##T;
#define SWIZ2_H1(n, ...) LEQ1##n(SWIZ2_H2(n, 0,x,r,s, __VA_ARGS__)) LEQ2##n(SWIZ2_H2(n, 1,y,g,t, __VA_ARGS__)) LEQ3##n(SWIZ2_H2(n, 2,z,b,p, __VA_ARGS__)) LEQ4##n(SWIZ2_H2(n, 3,w,a,q, __VA_ARGS__))
#define SWIZ2(n, ...)    LEQ1##n(SWIZ2_H1(n, 0,x,r,s, __VA_ARGS__)) LEQ2##n(SWIZ2_H1(n, 1,y,g,t, __VA_ARGS__)) LEQ3##n(SWIZ2_H1(n, 2,z,b,p, __VA_ARGS__)) LEQ4##n(SWIZ2_H1(n, 3,w,a,q, __VA_ARGS__))

#define SWIZ3_H3(n, I,X,R,S, J,Y,G,T, K,Z,B,P, ...) swiz3<I, J, K> X##Y##Z, R##G##B, S##T##P;
#define SWIZ3_H2(n, ...) LEQ1##n(SWIZ3_H3(n, 0,x,r,s, __VA_ARGS__)) LEQ2##n(SWIZ3_H3(n, 1,y,g,t, __VA_ARGS__)) LEQ3##n(SWIZ3_H3(n, 2,z,b,p, __VA_ARGS__)) LEQ4##n(SWIZ3_H3(n, 3,w,a,q, __VA_ARGS__))
#define SWIZ3_H1(n, ...) LEQ1##n(SWIZ3_H2(n, 0,x,r,s, __VA_ARGS__)) LEQ2##n(SWIZ3_H2(n, 1,y,g,t, __VA_ARGS__)) LEQ3##n(SWIZ3_H2(n, 2,z,b,p, __VA_ARGS__)) LEQ4##n(SWIZ3_H2(n, 3,w,a,q, __VA_ARGS__))
#define SWIZ3(n, ...)    LEQ1##n(SWIZ3_H1(n, 0,x,r,s, __VA_ARGS__)) LEQ2##n(SWIZ3_H1(n, 1,y,g,t, __VA_ARGS__)) LEQ3##n(SWIZ3_H1(n, 2,z,b,p, __VA_ARGS__)) LEQ4##n(SWIZ3_H1(n, 3,w,a,q, __VA_ARGS__))

#define SWIZ4_H4(n, I,X,R,S, J,Y,G,T, K,Z,B,P, L,W,A,Q, ...) swiz4<I, J, K, L> X##Y##Z##W, R##G##B##A, S##T##P##Q;
#define SWIZ4_H3(n, ...) LEQ1##n(SWIZ4_H4(n, 0,x,r,s, __VA_ARGS__)) LEQ2##n(SWIZ4_H4(n, 1,y,g,t, __VA_ARGS__)) LEQ3##n(SWIZ4_H4(n, 2,z,b,p, __VA_ARGS__)) LEQ4##n(SWIZ4_H4(n, 3,w,a,q, __VA_ARGS__))
#define SWIZ4_H2(n, ...) LEQ1##n(SWIZ4_H3(n, 0,x,r,s, __VA_ARGS__)) LEQ2##n(SWIZ4_H3(n, 1,y,g,t, __VA_ARGS__)) LEQ3##n(SWIZ4_H3(n, 2,z,b,p, __VA_ARGS__)) LEQ4##n(SWIZ4_H3(n, 3,w,a,q, __VA_ARGS__))
#define SWIZ4_H1(n, ...) LEQ1##n(SWIZ4_H2(n, 0,x,r,s, __VA_ARGS__)) LEQ2##n(SWIZ4_H2(n, 1,y,g,t, __VA_ARGS__)) LEQ3##n(SWIZ4_H2(n, 2,z,b,p, __VA_ARGS__)) LEQ4##n(SWIZ4_H2(n, 3,w,a,q, __VA_ARGS__))
#define SWIZ4(n, ...)    LEQ1##n(SWIZ4_H1(n, 0,x,r,s, __VA_ARGS__)) LEQ2##n(SWIZ4_H1(n, 1,y,g,t, __VA_ARGS__)) LEQ3##n(SWIZ4_H1(n, 2,z,b,p, __VA_ARGS__)) LEQ4##n(SWIZ4_H1(n, 3,w,a,q, __VA_ARGS__))

//-----------------------------------------------------------------------------
// vec classes
//-----------------------------------------------------------------------------
class alignas(8) vec2 {
public:
	union {
		float data[2];
		struct { float x, y; };
		struct { float r, g; };
		struct { float s, t; };
		SWIZ2(2);
		SWIZ3(2);
		SWIZ4(2);
	};

	inline explicit vec2() {}
	inline explicit vec2(float s) { data[0] = s; data[1] = s; }
	//inline explicit vec2(float x_, float y_) { data[0] = x_; data[1] = y_; }
	inline explicit constexpr vec2(float x_, float y_) : x(x_), y(y_) {}
	template <uint32_t i, uint32_t j>
	inline vec2(const swiz2<i, j>& xy) { data[0] = ((float*)&xy)[i]; data[1] = ((float*)&xy)[j]; }

	inline vec2& operator = (const vec2& other) { x = other.x; y = other.y; return *this; }
	inline float& operator [] (size_t i) { return data[i]; }
	inline float operator [] (size_t i) const { return data[i]; }
};

class alignas(16) vec3 {
public:
	union {
		float data[3];
		struct { float x, y, z; };
		struct { float r, g, b; };
		struct { float s, t, p; };
		SWIZ2(3);
		SWIZ3(3);
		SWIZ4(3);
	};

	inline explicit vec3() {}
	inline explicit vec3(float s) { data[0] = s; data[1] = s; data[2] = s; }
	//inline explicit vec3(float x_, float y_, float z_) { data[0] = x_; data[1] = y_;  data[2] = z_; }
	inline explicit constexpr vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}
	inline explicit vec3(const vec2& xy, float z_) { data[0] = xy.x; data[1] = xy.y;  data[2] = z_; }
	inline explicit vec3(float x_, const vec2& yz) { data[0] = x_; data[1] = yz.x;  data[2] = yz.y; }
	template <uint32_t i, uint32_t j, uint32_t k>
	inline vec3(const swiz3<i, j, k>& xyz) { data[0] = ((float*)&xyz)[i]; data[1] = ((float*)&xyz)[j]; data[2] = ((float*)&xyz)[k]; }

	inline vec3& operator = (const vec3& other) { x = other.x; y = other.y; z = other.z; return *this; }
	inline float& operator [] (size_t i) { return data[i]; }
	inline float operator [] (size_t i) const { return data[i]; }
};

class alignas(16) vec4 {
public:
	union {
		float data[4];
		struct { float x, y, z, w; };
		struct { float r, g, b, a; };
		struct { float s, t, p, q; };
		SWIZ2(4);
		SWIZ3(4);
		SWIZ4(4);
	};

	inline explicit vec4() {}
	inline explicit vec4(float s) { data[0] = s; data[1] = s; data[2] = s; data[3] = s; }
	inline explicit vec4(float x_, float y_, float z_, float w_) { data[0] = x_; data[1] = y_;  data[2] = z_; data[3] = w_; }
	inline explicit vec4(const vec2& xy, float z_, float w_) { data[0] = xy.x; data[1] = xy.y;  data[2] = z_; data[3] = w_; }
	inline explicit vec4(float x_, const vec2& yz, float w_) { data[0] = x_; data[1] = yz.x;  data[2] = yz.y; data[3] = w_; }
	inline explicit vec4(float x_, float y_, const vec2& zw) { data[0] = x_; data[1] = y_;  data[2] = zw.x; data[3] = zw.y; }
	inline explicit vec4(const vec2& xy, const vec2& zw) { data[0] = xy.x; data[1] = xy.y;  data[2] = zw.x; data[3] = zw.y; }
	inline explicit vec4(const vec3& xyz, float w_) { data[0] = xyz.x; data[1] = xyz.y;  data[2] = xyz.z; data[3] = w_; }
	inline explicit vec4(float x_, const vec3& yzw) { data[0] = x_; data[1] = yzw.x;  data[2] = yzw.y; data[3] = yzw.z; }
	template <uint32_t i, uint32_t j, uint32_t k, uint32_t l>
	inline vec4(const swiz4<i, j, k, l>& xyzw) { data[0] = ((float*)&xyzw)[i]; data[1] = ((float*)&xyzw)[j]; data[2] = ((float*)&xyzw)[k]; data[3] = ((float*)&xyzw)[l]; }

	inline vec4& operator = (const vec4& other) { x = other.x; y = other.y; z = other.z; w = other.w; return *this; }
	inline float& operator [] (size_t i) { return data[i]; }
	inline float operator [] (size_t i) const { return data[i]; }
};

//-----------------------------------------------------------------------------
// mat classes
//-----------------------------------------------------------------------------
class mat2 {
public:
	vec2 cols[2];

	inline mat2() {}
	inline mat2(float s) { cols[0] = vec2(s, 0.0f); cols[1] = vec2(0.0f, s); }
	inline mat2(const vec2& c0, const vec2& c1) { cols[0] = c0; cols[1] = c1; }

	inline vec2& operator [] (size_t i) { return cols[i]; }
	inline vec2 operator [] (size_t i) const { return cols[i]; }
};

class mat3 {
public:
	vec3 cols[3];

	inline mat3() {}
	inline mat3(float s) { cols[0] = vec3(s, 0.0f, 0.0f); cols[1] = vec3(0.0f, s, 0.0f); cols[2] = vec3(0.0f, 0.0f, s); }
	inline mat3(const vec3& c0, const vec3& c1, const vec3& c2) { cols[0] = c0; cols[1] = c1; cols[2] = c2; }

	inline vec3& operator [] (size_t i) { return cols[i]; }
	inline vec3 operator [] (size_t i) const { return cols[i]; }
};

class mat4 {
public:
	vec4 cols[4];

	inline mat4() {}
	inline mat4(float s) { cols[0] = vec4(s, 0.0f, 0.0f, 0.0f); cols[1] = vec4(0.0f, s, 0.0f, 0.0f); cols[2] = vec4(0.0f, 0.0f, s, 0.0f); cols[3] = vec4(0.0f, 0.0f, 0.0f, s); }
	inline mat4(const vec4& c0, const vec4& c1, const vec4& c2, const vec4& c3) { cols[0] = c0; cols[1] = c1; cols[2] = c2; cols[3] = c3; }

	inline vec4& operator [] (size_t i) { return cols[i]; }
	inline vec4 operator [] (size_t i) const { return cols[i]; }
};

//-----------------------------------------------------------------------------
// Helpers for defining vec interfaces
//-----------------------------------------------------------------------------
#if _MSC_VER && !__INTEL_COMPILER
#define VEC_FORCEINLINE __forceinline
#else
#define VEC_FORCEINLINE __attribute__((always_inline))
#endif

inline VEC_FORCEINLINE vec2 vec_map(float(*f)(float), vec2 x) { return vec2(f(x[0]), f(x[1])); }
inline VEC_FORCEINLINE vec2 vec_map(float(*f)(float, float), vec2 x, vec2 y) { return vec2(f(x[0], y[0]), f(x[1], y[1])); }
inline VEC_FORCEINLINE vec2 vec_map(float(*f)(float, float, float), vec2 x, vec2 y, vec2 z) { return vec2(f(x[0], y[0], z[0]), f(x[1], y[1], z[1])); }
inline VEC_FORCEINLINE vec3 vec_map(float(*f)(float), vec3 x) { return vec3(f(x[0]), f(x[1]), f(x[2])); }
inline VEC_FORCEINLINE vec3 vec_map(float(*f)(float, float), vec3 x, vec3 y) { return vec3(f(x[0], y[0]), f(x[1], y[1]), f(x[2], y[2])); }
inline VEC_FORCEINLINE vec3 vec_map(float(*f)(float, float, float), vec3 x, vec3 y, vec3 z) { return vec3(f(x[0], y[0], z[0]), f(x[1], y[1], z[1]), f(x[2], y[2], z[2])); }
inline VEC_FORCEINLINE vec4 vec_map(float(*f)(float), vec4 x) { return vec4(f(x[0]), f(x[1]), f(x[2]), f(x[3])); }
inline VEC_FORCEINLINE vec4 vec_map(float(*f)(float, float), vec4 x, vec4 y) { return vec4(f(x[0], y[0]), f(x[1], y[1]), f(x[2], y[2]), f(x[3], y[3])); }
inline VEC_FORCEINLINE vec4 vec_map(float(*f)(float, float, float), vec4 x, vec4 y, vec4 z) { return vec4(f(x[0], y[0], z[0]), f(x[1], y[1], z[1]), f(x[2], y[2], z[2]), f(x[3], y[3], z[3])); }

//-----------------------------------------------------------------------------
// float GLSL functions missing from <cmath>
//-----------------------------------------------------------------------------
inline VEC_FORCEINLINE float clamp(float x, float vmin, float vmax) { return min(max(vmin, x), vmax); }
inline VEC_FORCEINLINE float degrees(float radians) { return 180.0f * radians / 3.14159265358979323846f; }
inline VEC_FORCEINLINE float distance(float x, float y) { return abs(x - y); }
inline VEC_FORCEINLINE float dot(float x, float y) { return x * y; }
inline VEC_FORCEINLINE float fract(float x) { return x - floor(x); }
inline VEC_FORCEINLINE float inversesqrt(float x) { return sqrt(1.0f / x); }
inline VEC_FORCEINLINE float length(float x) { return abs(x); }
inline VEC_FORCEINLINE float mix(float x, float y, float a) { return x * (1.0f - a) + y * a; }
inline VEC_FORCEINLINE float mod(float x, float y) { return x - y * floor(x / y); }
inline VEC_FORCEINLINE float modf(float x, float& i) { return modff(x, &i); }
inline VEC_FORCEINLINE float normalize(float x) { return x * (1.0f / length(x)); }
inline VEC_FORCEINLINE float safeNormalize(float x) { float len = max(0.0000001f, length(x)); return x * (1.0f / len); }
inline VEC_FORCEINLINE float radians(float degrees) { return 3.14159265358979323846f * degrees / 180.0f; }
inline VEC_FORCEINLINE float sign(float x) { return (x > 0.0f) ? 1.0f : ((x < 0.0f) ? -1.0f : 0.0f); }
inline VEC_FORCEINLINE float smoothstep(float edge0, float edge1, float x) { float t = clamp((x - edge0) / (edge1 - edge0), 0.0f, 1.0f); return t * t * (3.0f - 2.0f * t); }
inline VEC_FORCEINLINE float step(float edge, float x) { return x < edge ? 0.0f : 1.0f; }

//-----------------------------------------------------------------------------
// vec2 operators and functions
//-----------------------------------------------------------------------------
inline vec2 operator + (const vec2& x) { return x; }
inline vec2 operator - (const vec2& x) { return vec2(-x.data[0], -x.data[1]); }
inline vec2 operator + (const vec2& x, const vec2& y) { return vec2(x.data[0] + y.data[0], x.data[1] + y.data[1]); }
inline vec2 operator - (const vec2& x, const vec2& y) { return vec2(x.data[0] - y.data[0], x.data[1] - y.data[1]); }
inline vec2 operator * (const vec2& x, const vec2& y) { return vec2(x.data[0] * y.data[0], x.data[1] * y.data[1]); }
inline vec2 operator / (const vec2& x, const vec2& y) { return vec2(x.data[0] / y.data[0], x.data[1] / y.data[1]); }
inline vec2 operator + (const vec2& x, float s) { return x + vec2(s); }
inline vec2 operator + (float s, const vec2& x) { return vec2(s) + x; }
inline vec2 operator - (const vec2& x, float s) { return x - vec2(s); }
inline vec2 operator - (float s, const vec2& x) { return vec2(s) - x; }
inline vec2 operator * (const vec2& x, float s) { return x * vec2(s); }
inline vec2 operator * (float s, const vec2& x) { return vec2(s) * x; }
inline vec2 operator / (const vec2& x, float s) { return x / vec2(s); }
inline vec2 operator / (float s, const vec2& x) { return vec2(s) / x; }
inline vec2& operator += (vec2& x, const vec2& y) { x = x + y; return x; }
inline vec2& operator -= (vec2& x, const vec2& y) { x = x - y; return x; }
inline vec2& operator *= (vec2& x, const vec2& y) { x = x * y; return x; }
inline vec2& operator /= (vec2& x, const vec2& y) { x = x / y; return x; }
inline vec2& operator += (vec2& x, float s) { x = x + s; return x; }
inline vec2& operator -= (vec2& x, float s) { x = x - s; return x; }
inline vec2& operator *= (vec2& x, float s) { x = x * s; return x; }
inline vec2& operator /= (vec2& x, float s) { x = x / s; return x; }

inline vec2 abs(vec2 x) { return vec_map(abs, x); }
inline vec2 acos(vec2 x) { return vec_map(acos, x); }
inline vec2 acosh(vec2 x) { return vec_map(acosh, x); }
inline vec2 asin(vec2 x) { return vec_map(asin, x); }
inline vec2 asinh(vec2 x) { return vec_map(asinh, x); }
inline vec2 atan(vec2 x) { return vec_map(atan, x); }
inline vec2 atanh(vec2 x) { return vec_map(atanh, x); }
inline vec2 ceil(vec2 x) { return vec_map(ceil, x); }
inline vec2 clamp(vec2 x, vec2 vmin, vec2 vmax) { return vec_map(clamp, x, vmin, vmax); }
inline vec2 clamp(vec2 x, float smin, float smax) { return clamp(x, vec2(smin), vec2(smax)); }
inline vec2 cos(vec2 x) { return vec_map(cos, x); }
inline vec2 cosh(vec2 x) { return vec_map(cosh, x); }
inline vec2 degrees(vec2 radians) { return vec_map(degrees, radians); }
inline float distance(vec2 x, vec2 y) { vec2 a = x - y; return sqrt(a[0] * a[0] + a[1] * a[1]); }
inline float dot(vec2 x, vec2 y) { return x[0] * y[0] + x[1] * y[1]; }
inline vec2 exp(vec2 x) { return vec_map(exp, x); }
inline vec2 exp2(vec2 x) { return vec_map(exp2, x); }
inline vec2 floor(vec2 x) { return vec_map(floor, x); }
inline vec2 fma(vec2 a, vec2 b, vec2 c) { return vec_map(fma, a, b, c); }
inline vec2 fract(vec2 x) { return vec_map(fract, x); }
inline vec2 inversesqrt(vec2 x) { return vec_map(inversesqrt, x); }
inline float length(vec2 x) { return sqrt(dot(x, x)); }
inline vec2 log(vec2 x) { return vec_map(log, x); }
inline vec2 log2(vec2 x) { return vec_map(log2, x); }
inline vec2 max(vec2 x, vec2 y) { return vec_map(max, x, y); }
inline vec2 min(vec2 x, vec2 y) { return vec_map(min, x, y); }
inline vec2 mix(vec2 x, vec2 y, vec2 a) { return x * (1.0f - a) + y * a; }
inline vec2 mix(vec2 x, vec2 y, float a) { return x * (1.0f - a) + y * a; }
inline vec2 mod(vec2 x, vec2 y) { return vec_map(mod, x, y); }
inline vec2 mod(vec2 x, float y) { return mod(x, vec2(y)); }
inline vec2 modf(vec2 x, vec2& i) { return vec2(modf(x[0], &i[0]), modf(x[1], &i[1])); }
inline vec2 normalize(vec2 x) { return x * (1.0f / length(x)); }
inline vec2 safeNormalize(vec2 x) { float len = max(0.0000001f, length(x)); return x * (1.0f / len); }
inline vec2 pow(vec2 x, vec2 y) { return vec_map(pow, x, y); }
inline vec2 radians(vec2 degrees) { return vec_map(radians, degrees); }
inline vec2 round(vec2 x) { return vec_map(round, x); }
inline vec2 sign(vec2 x) { return vec_map(sign, x); }
inline vec2 sin(vec2 x) { return vec_map(sin, x); }
inline vec2 sinh(vec2 x) { return vec_map(sinh, x); }
inline vec2 smoothstep(vec2 edge0, vec2 edge1, vec2 x) { return vec_map(smoothstep, edge0, edge1, x); }
inline vec2 sqrt(vec2 x) { return vec_map(sqrt, x); }
inline vec2 step(vec2 edge, vec2 x) { return vec_map(step, edge, x); }
inline vec2 step(float edge, vec2 x) { return step(vec2(edge), x); }
inline vec2 tan(vec2 x) { return vec_map(tan, x); }
inline vec2 tanh(vec2 x) { return vec_map(tanh, x); }
inline vec2 trunc(vec2 x) { return vec_map(trunc, x); }

//-----------------------------------------------------------------------------
// vec3 operators and functions
//-----------------------------------------------------------------------------
inline vec3 operator + (const vec3& x) { return x; }
inline vec3 operator - (const vec3& x) { return vec3(-x.data[0], -x.data[1], -x.data[2]); }
inline vec3 operator + (const vec3& x, const vec3& y) { return vec3(x.data[0] + y.data[0], x.data[1] + y.data[1], x.data[2] + y.data[2]); }
inline vec3 operator - (const vec3& x, const vec3& y) { return vec3(x.data[0] - y.data[0], x.data[1] - y.data[1], x.data[2] - y.data[2]); }
inline vec3 operator * (const vec3& x, const vec3& y) { return vec3(x.data[0] * y.data[0], x.data[1] * y.data[1], x.data[2] * y.data[2]); }
inline vec3 operator / (const vec3& x, const vec3& y) { return vec3(x.data[0] / y.data[0], x.data[1] / y.data[1], x.data[2] / y.data[2]); }
inline vec3 operator + (const vec3& x, float s) { return x + vec3(s); }
inline vec3 operator + (float s, const vec3& x) { return vec3(s) + x; }
inline vec3 operator - (const vec3& x, float s) { return x - vec3(s); }
inline vec3 operator - (float s, const vec3& x) { return vec3(s) - x; }
inline vec3 operator * (const vec3& x, float s) { return x * vec3(s); }
inline vec3 operator * (float s, const vec3& x) { return vec3(s) * x; }
inline vec3 operator / (const vec3& x, float s) { return x / vec3(s); }
inline vec3 operator / (float s, const vec3& x) { return vec3(s) / x; }
inline vec3& operator += (vec3& x, const vec3& y) { x = x + y; return x; }
inline vec3& operator -= (vec3& x, const vec3& y) { x = x - y; return x; }
inline vec3& operator *= (vec3& x, const vec3& y) { x = x * y; return x; }
inline vec3& operator /= (vec3& x, const vec3& y) { x = x / y; return x; }
inline vec3& operator += (vec3& x, float s) { x = x + s; return x; }
inline vec3& operator -= (vec3& x, float s) { x = x - s; return x; }
inline vec3& operator *= (vec3& x, float s) { x = x * s; return x; }
inline vec3& operator /= (vec3& x, float s) { x = x / s; return x; }

inline vec3 abs(vec3 x) { return vec_map(abs, x); }
inline vec3 acos(vec3 x) { return vec_map(acos, x); }
inline vec3 acosh(vec3 x) { return vec_map(acosh, x); }
inline vec3 asin(vec3 x) { return vec_map(asin, x); }
inline vec3 asinh(vec3 x) { return vec_map(asinh, x); }
inline vec3 atan(vec3 x) { return vec_map(atan, x); }
inline vec3 atanh(vec3 x) { return vec_map(atanh, x); }
inline vec3 ceil(vec3 x) { return vec_map(ceil, x); }
inline vec3 clamp(vec3 x, vec3 vmin, vec3 vmax) { return vec_map(clamp, x, vmin, vmax); }
inline vec3 clamp(vec3 x, float smin, float smax) { return clamp(x, vec3(smin), vec3(smax)); }
inline vec3 cos(vec3 x) { return vec_map(cos, x); }
inline vec3 cosh(vec3 x) { return vec_map(cosh, x); }
inline vec3 cross(vec3 x, vec3 y) { return vec3(x[1] * y[2] - y[1] * x[2], x[2] * y[0] - y[2] * x[0], x[0] * y[1] - y[0] * x[1]); }
inline vec3 degrees(vec3 radians) { return vec_map(degrees, radians); }
inline float distance(vec3 x, vec3 y) { vec3 a = x - y; return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }
inline float dot(vec3 x, vec3 y) { return x[0] * y[0] + x[1] * y[1] + x[2] * y[2]; }
inline vec3 exp(vec3 x) { return vec_map(exp, x); }
inline vec3 exp2(vec3 x) { return vec_map(exp2, x); }
inline vec3 floor(vec3 x) { return vec_map(floor, x); }
inline vec3 fma(vec3 a, vec3 b, vec3 c) { return vec_map(fma, a, b, c); }
inline vec3 fract(vec3 x) { return vec_map(fract, x); }
inline vec3 inversesqrt(vec3 x) { return vec_map(inversesqrt, x); }
inline float length(vec3 x) { return sqrt(dot(x, x)); }
inline vec3 log(vec3 x) { return vec_map(log, x); }
inline vec3 log2(vec3 x) { return vec_map(log2, x); }
inline vec3 max(vec3 x, vec3 y) { return vec_map(max, x, y); }
inline vec3 min(vec3 x, vec3 y) { return vec_map(min, x, y); }
inline vec3 mix(vec3 x, vec3 y, vec3 a) { return x * (1.0f - a) + y * a; }
inline vec3 mix(vec3 x, vec3 y, float a) { return x * (1.0f - a) + y * a; }
inline vec3 mod(vec3 x, vec3 y) { return vec_map(mod, x, y); }
inline vec3 mod(vec3 x, float y) { return mod(x, vec3(y)); }
inline vec3 modf(vec3 x, vec3& i) { return vec3(modf(x[0], &i[0]), modf(x[1], &i[1]), modf(x[2], &i[2])); }
inline vec3 normalize(vec3 x) { return x * (1.0f / length(x)); }
inline vec3 safeNormalize(vec3 x) { float len = max(0.0000001f, length(x)); return x * (1.0f / len); }
inline vec3 pow(vec3 x, vec3 y) { return vec_map(pow, x, y); }
inline vec3 radians(vec3 degrees) { return vec_map(radians, degrees); }
inline vec3 round(vec3 x) { return vec_map(round, x); }
inline vec3 sign(vec3 x) { return vec_map(sign, x); }
inline vec3 sin(vec3 x) { return vec_map(sin, x); }
inline vec3 sinh(vec3 x) { return vec_map(sinh, x); }
inline vec3 smoothstep(vec3 edge0, vec3 edge1, vec3 x) { return vec_map(smoothstep, edge0, edge1, x); }
inline vec3 sqrt(vec3 x) { return vec_map(sqrt, x); }
inline vec3 step(vec3 edge, vec3 x) { return vec_map(step, edge, x); }
inline vec3 step(float edge, vec3 x) { return step(vec3(edge), x); }
inline vec3 tan(vec3 x) { return vec_map(tan, x); }
inline vec3 tanh(vec3 x) { return vec_map(tanh, x); }
inline vec3 trunc(vec3 x) { return vec_map(trunc, x); }

//-----------------------------------------------------------------------------
// vec4 operators and functions
//-----------------------------------------------------------------------------
inline vec4 operator + (const vec4& x) { return x; }
inline vec4 operator - (const vec4& x) { return vec4(-x.data[0], -x.data[1], -x.data[2], -x.data[3]); }
inline vec4 operator + (const vec4& x, const vec4& y) { return vec4(x.data[0] + y.data[0], x.data[1] + y.data[1], x.data[2] + y.data[2], x.data[3] + y.data[3]); }
inline vec4 operator - (const vec4& x, const vec4& y) { return vec4(x.data[0] - y.data[0], x.data[1] - y.data[1], x.data[2] - y.data[2], x.data[3] - y.data[3]); }
inline vec4 operator * (const vec4& x, const vec4& y) { return vec4(x.data[0] * y.data[0], x.data[1] * y.data[1], x.data[2] * y.data[2], x.data[3] * y.data[3]); }
inline vec4 operator / (const vec4& x, const vec4& y) { return vec4(x.data[0] / y.data[0], x.data[1] / y.data[1], x.data[2] / y.data[2], x.data[3] / y.data[3]); }
inline vec4 operator + (const vec4& x, float s) { return x + vec4(s); }
inline vec4 operator + (float s, const vec4& x) { return vec4(s) + x; }
inline vec4 operator - (const vec4& x, float s) { return x - vec4(s); }
inline vec4 operator - (float s, const vec4& x) { return vec4(s) - x; }
inline vec4 operator * (const vec4& x, float s) { return x * vec4(s); }
inline vec4 operator * (float s, const vec4& x) { return vec4(s) * x; }
inline vec4 operator / (const vec4& x, float s) { return x / vec4(s); }
inline vec4 operator / (float s, const vec4& x) { return vec4(s) / x; }
inline vec4& operator += (vec4& x, const vec4& y) { x = x + y; return x; }
inline vec4& operator -= (vec4& x, const vec4& y) { x = x - y; return x; }
inline vec4& operator *= (vec4& x, const vec4& y) { x = x * y; return x; }
inline vec4& operator /= (vec4& x, const vec4& y) { x = x / y; return x; }
inline vec4& operator += (vec4& x, float s) { x = x + s; return x; }
inline vec4& operator -= (vec4& x, float s) { x = x - s; return x; }
inline vec4& operator *= (vec4& x, float s) { x = x * s; return x; }
inline vec4& operator /= (vec4& x, float s) { x = x / s; return x; }

inline vec4 abs(vec4 x) { return vec_map(abs, x); }
inline vec4 acos(vec4 x) { return vec_map(acos, x); }
inline vec4 acosh(vec4 x) { return vec_map(acosh, x); }
inline vec4 asin(vec4 x) { return vec_map(asin, x); }
inline vec4 asinh(vec4 x) { return vec_map(asinh, x); }
inline vec4 atan(vec4 x) { return vec_map(atan, x); }
inline vec4 atanh(vec4 x) { return vec_map(atanh, x); }
inline vec4 ceil(vec4 x) { return vec_map(ceil, x); }
inline vec4 clamp(vec4 x, vec4 vmin, vec4 vmax) { return vec_map(clamp, x, vmin, vmax); }
inline vec4 clamp(vec4 x, float smin, float smax) { return clamp(x, vec4(smin), vec4(smax)); }
inline vec4 cos(vec4 x) { return vec_map(cos, x); }
inline vec4 cosh(vec4 x) { return vec_map(cosh, x); }
inline vec4 degrees(vec4 radians) { return vec_map(degrees, radians); }
inline float distance(vec4 x, vec4 y) { vec4 a = x - y; return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2] + a[3] * a[3]); }
inline float dot(vec4 x, vec4 y) { return x[0] * y[0] + x[1] * y[1] + x[2] * y[2] + x[3] * y[3]; }
inline vec4 exp(vec4 x) { return vec_map(exp, x); }
inline vec4 exp2(vec4 x) { return vec_map(exp2, x); }
inline vec4 floor(vec4 x) { return vec_map(floor, x); }
inline vec4 fma(vec4 a, vec4 b, vec4 c) { return vec_map(fma, a, b, c); }
inline vec4 fract(vec4 x) { return vec_map(fract, x); }
inline vec4 inversesqrt(vec4 x) { return vec_map(inversesqrt, x); }
inline float length(vec4 x) { return sqrt(dot(x, x)); }
inline vec4 log(vec4 x) { return vec_map(log, x); }
inline vec4 log2(vec4 x) { return vec_map(log2, x); }
inline vec4 max(vec4 x, vec4 y) { return vec_map(max, x, y); }
inline vec4 min(vec4 x, vec4 y) { return vec_map(min, x, y); }
inline vec4 mix(vec4 x, vec4 y, vec4 a) { return vec_map(mix, x, y, a); }
inline vec4 mix(vec4 x, vec4 y, float a) { return mix(x, y, vec4(a)); }
inline vec4 mod(vec4 x, vec4 y) { return vec_map(mod, x, y); }
inline vec4 mod(vec4 x, float y) { return mod(x, vec4(y)); }
inline vec4 modf(vec4 x, vec4& i) { return vec4(modf(x[0], &i[0]), modf(x[1], &i[1]), modf(x[2], &i[2]), modf(x[3], &i[3])); }
inline vec4 normalize(vec4 x) { return x * (1.0f / length(x)); }
inline vec4 safeNormalize(vec4 x) { float len = max(0.0000001f, length(x)); return x * (1.0f / len); }
inline vec4 pow(vec4 x, vec4 y) { return vec_map(pow, x, y); }
inline vec4 radians(vec4 degrees) { return vec_map(radians, degrees); }
inline vec4 round(vec4 x) { return vec_map(round, x); }
inline vec4 sign(vec4 x) { return vec_map(sign, x); }
inline vec4 sin(vec4 x) { return vec_map(sin, x); }
inline vec4 sinh(vec4 x) { return vec_map(sinh, x); }
inline vec4 smoothstep(vec4 edge0, vec4 edge1, vec4 x) { return vec_map(smoothstep, edge0, edge1, x); }
inline vec4 sqrt(vec4 x) { return vec_map(sqrt, x); }
inline vec4 step(vec4 edge, vec4 x) { return vec_map(step, edge, x); }
inline vec4 step(float edge, vec4 x) { return step(vec4(edge), x); }
inline vec4 tan(vec4 x) { return vec_map(tan, x); }
inline vec4 tanh(vec4 x) { return vec_map(tanh, x); }
inline vec4 trunc(vec4 x) { return vec_map(trunc, x); }

//-----------------------------------------------------------------------------
// mat2 operators and functions
//-----------------------------------------------------------------------------
inline mat2 operator + (const mat2& x, const mat2& y) { return mat2(x[0] + y[0], x[1] + y[1]); }
inline mat2 operator - (const mat2& x, const mat2& y) { return mat2(x[0] - y[0], x[1] - y[1]); }
inline mat2 operator * (const mat2& m, float s) { return mat2(m[0] * s, m[1] * s); }
inline mat2 operator * (float s, const mat2& m) { return m * s; }
inline mat2 operator / (const mat2& m, float s) { return m * (1.0f / s); }
inline vec2 operator * (const mat2& m, const vec2& v) {
	vec2 r0 = vec2(m[0][0], m[1][0]);
	vec2 r1 = vec2(m[0][1], m[1][1]);
	return vec2(dot(r0, v), dot(r1, v));
}
inline mat2 operator * (const mat2& x, const mat2& y) { return mat2(x * y[0], x * y[1]); }
inline mat2& operator += (mat2& x, const mat2& y) { x = x + y; return x; }
inline mat2& operator -= (mat2& x, const mat2& y) { x = x - y; return x; }
inline mat2& operator *= (mat2& m, float s) { m = m * s; return m; }
inline mat2& operator /= (mat2& m, float s) { m = m / s; return m; }
inline mat2& operator *= (mat2& x, const mat2& y) { x = x * y; return x; }

inline mat2 matrixCompMult(mat2 x, mat2 y) { return mat2(x[0] * y[0], x[1] * y[1]); }
inline mat2 outerProduct(vec2 c, vec2 r) { return mat2(c * r[0], c * r[1]); }

mat2 transpose(const mat2& m);
float determinant(const mat2& m);
mat2 inverse(const mat2& m);

//-----------------------------------------------------------------------------
// mat3 operators and functions
//-----------------------------------------------------------------------------
inline mat3 operator + (const mat3& x, const mat3& y) { return mat3(x[0] + y[0], x[1] + y[1], x[2] + y[2]); }
inline mat3 operator - (const mat3& x, const mat3& y) { return mat3(x[0] - y[0], x[1] - y[1], x[2] - y[2]); }
inline mat3 operator * (const mat3& m, float s) { return mat3(m[0] * s, m[1] * s, m[2] * s); }
inline mat3 operator * (float s, const mat3& m) { return m * s; }
inline mat3 operator / (const mat3& m, float s) { return m * (1.0f / s); }
inline vec3 operator * (const mat3& m, const vec3& v) {
	vec3 r0 = vec3(m[0][0], m[1][0], m[2][0]);
	vec3 r1 = vec3(m[0][1], m[1][1], m[2][1]);
	vec3 r2 = vec3(m[0][2], m[1][2], m[2][2]);
	return vec3(dot(r0, v), dot(r1, v), dot(r2, v));
}
inline mat3 operator * (const mat3& x, const mat3& y) { return mat3(x * y[0], x * y[1], x * y[2]); }
inline mat3& operator += (mat3& x, const mat3& y) { x = x + y; return x; }
inline mat3& operator -= (mat3& x, const mat3& y) { x = x - y; return x; }
inline mat3& operator *= (mat3& m, float s) { m = m * s; return m; }
inline mat3& operator /= (mat3& m, float s) { m = m / s; return m; }
inline mat3& operator *= (mat3& x, const mat3& y) { x = x * y; return x; }

inline mat3 matrixCompMult(mat3 x, mat3 y) { return mat3(x[0] * y[0], x[1] * y[1], x[2] * y[2]); }
inline mat3 outerProduct(vec3 c, vec3 r) { return mat3(c * r[0], c * r[1], c * r[2]); }

mat3 transpose(const mat3& m);
float determinant(const mat3& m);
mat3 inverse(const mat3& m);

//-----------------------------------------------------------------------------
// mat4 operators and functions
//-----------------------------------------------------------------------------
inline mat4 operator + (const mat4& x, const mat4& y) { return mat4(x[0] + y[0], x[1] + y[1], x[2] + y[2], x[3] + y[3]); }
inline mat4 operator - (const mat4& x, const mat4& y) { return mat4(x[0] - y[0], x[1] - y[1], x[2] - y[2], x[3] - y[3]); }
inline mat4 operator * (const mat4& m, float s) { return mat4(m[0] * s, m[1] * s, m[2] * s, m[3] * s); }
inline mat4 operator * (float s, const mat4& m) { return m * s; }
inline mat4 operator / (const mat4& m, float s) { return m * (1.0f / s); }
inline vec4 operator * (const mat4& m, const vec4& v) {
	vec4 r0 = vec4(m[0][0], m[1][0], m[2][0], m[3][0]);
	vec4 r1 = vec4(m[0][1], m[1][1], m[2][1], m[3][1]);
	vec4 r2 = vec4(m[0][2], m[1][2], m[2][2], m[3][2]);
	vec4 r3 = vec4(m[0][3], m[1][3], m[2][3], m[3][3]);
	return vec4(dot(r0, v), dot(r1, v), dot(r2, v), dot(r3, v));
}
inline mat4 operator * (const mat4& x, const mat4& y) { return mat4(x * y[0], x * y[1], x * y[2], x * y[3]); }
inline mat4& operator += (mat4& x, const mat4& y) { x = x + y; return x; }
inline mat4& operator -= (mat4& x, const mat4& y) { x = x - y; return x; }
inline mat4& operator *= (mat4& m, float s) { m = m * s; return m; }
inline mat4& operator /= (mat4& m, float s) { m = m / s; return m; }
inline mat4& operator *= (mat4& x, const mat4& y) { x = x * y; return x; }

inline mat4 matrixCompMult(mat4 x, mat4 y) { return mat4(x[0] * y[0], x[1] * y[1], x[2] * y[2], x[3] * y[3]); }
inline mat4 outerProduct(vec4 c, vec4 r) { return mat4(c * r[0], c * r[1], c * r[2], c * r[3]); }

mat4 transpose(const mat4& m);
float determinant(const mat4& m);
mat4 inverse(const mat4& m);

//-----------------------------------------------------------------------------
// Partial dvec implementation
//-----------------------------------------------------------------------------
class alignas(16) dvec2 {
public:
	using vec = vec2;

	union {
		double data[2];
		struct { double x, y; };
		struct { double r, g; };
		struct { double s, t; };
	};

	inline explicit dvec2() {}
	inline explicit dvec2(double s) { data[0] = s; data[1] = s; }
	inline explicit dvec2(double x_, double y_) { data[0] = x_; data[1] = y_; }
	inline explicit dvec2(const vec2& v) { data[0] = (double)v.x; data[1] = (double)v.y; }

	inline explicit operator vec2() { return vec2((float)data[0], (float)data[1]); }

	inline dvec2& operator = (const dvec2& other) { x = other.x; y = other.y; return *this; }
	inline double& operator [] (size_t i) { return data[i]; }
	inline double operator [] (size_t i) const { return data[i]; }
};

inline dvec2 operator + (const dvec2& x) { return x; }
inline dvec2 operator - (const dvec2& x) { return dvec2(-x.data[0], -x.data[1]); }
inline dvec2 operator + (const dvec2& x, const dvec2& y) { return dvec2(x.data[0] + y.data[0], x.data[1] + y.data[1]); }
inline dvec2 operator - (const dvec2& x, const dvec2& y) { return dvec2(x.data[0] - y.data[0], x.data[1] - y.data[1]); }
inline dvec2 operator * (const dvec2& x, const dvec2& y) { return dvec2(x.data[0] * y.data[0], x.data[1] * y.data[1]); }
inline dvec2 operator / (const dvec2& x, const dvec2& y) { return dvec2(x.data[0] / y.data[0], x.data[1] / y.data[1]); }
inline dvec2 operator + (const dvec2& x, double s) { return x + dvec2(s); }
inline dvec2 operator + (double s, const dvec2& x) { return dvec2(s) + x; }
inline dvec2 operator - (const dvec2& x, double s) { return x - dvec2(s); }
inline dvec2 operator - (double s, const dvec2& x) { return dvec2(s) - x; }
inline dvec2 operator * (const dvec2& x, double s) { return x * dvec2(s); }
inline dvec2 operator * (double s, const dvec2& x) { return dvec2(s) * x; }
inline dvec2 operator / (const dvec2& x, double s) { return x / dvec2(s); }
inline dvec2 operator / (double s, const dvec2& x) { return dvec2(s) / x; }
inline dvec2& operator += (dvec2& x, const dvec2& y) { x = x + y; return x; }
inline dvec2& operator -= (dvec2& x, const dvec2& y) { x = x - y; return x; }
inline dvec2& operator *= (dvec2& x, const dvec2& y) { x = x * y; return x; }
inline dvec2& operator /= (dvec2& x, const dvec2& y) { x = x / y; return x; }
inline dvec2& operator += (dvec2& x, double s) { x = x + s; return x; }
inline dvec2& operator -= (dvec2& x, double s) { x = x - s; return x; }
inline dvec2& operator *= (dvec2& x, double s) { x = x * s; return x; }
inline dvec2& operator /= (dvec2& x, double s) { x = x / s; return x; }

class alignas(32) dvec3 {
public:
	using vec = vec3;

	union {
		double data[3];
		struct { double x, y, z; };
		struct { double r, g, b; };
		struct { double s, t, p; };
	};

	inline explicit dvec3() {}
	inline explicit dvec3(double s) { data[0] = s; data[1] = s; data[2] = s; }
	inline explicit dvec3(double x_, double y_, double z_) { data[0] = x_; data[1] = y_; data[2] = z_; }
	inline explicit dvec3(const vec3& v) { data[0] = (double)v.x; data[1] = (double)v.y; data[2] = (double)v.z; }

	inline explicit operator vec3() { return vec3((float)data[0], (float)data[1], (float)data[2]); }

	inline dvec3& operator = (const dvec3& other) { x = other.x; y = other.y; z = other.z; return *this; }
	inline double& operator [] (size_t i) { return data[i]; }
	inline double operator [] (size_t i) const { return data[i]; }
};

inline dvec3 operator + (const dvec3& x) { return x; }
inline dvec3 operator - (const dvec3& x) { return dvec3(-x.data[0], -x.data[1], -x.data[2]); }
inline dvec3 operator + (const dvec3& x, const dvec3& y) { return dvec3(x.data[0] + y.data[0], x.data[1] + y.data[1], x.data[2] + y.data[2]); }
inline dvec3 operator - (const dvec3& x, const dvec3& y) { return dvec3(x.data[0] - y.data[0], x.data[1] - y.data[1], x.data[2] - y.data[2]); }
inline dvec3 operator * (const dvec3& x, const dvec3& y) { return dvec3(x.data[0] * y.data[0], x.data[1] * y.data[1], x.data[2] * y.data[2]); }
inline dvec3 operator / (const dvec3& x, const dvec3& y) { return dvec3(x.data[0] / y.data[0], x.data[1] / y.data[1], x.data[2] / y.data[2]); }
inline dvec3 operator + (const dvec3& x, double s) { return x + dvec3(s); }
inline dvec3 operator + (double s, const dvec3& x) { return dvec3(s) + x; }
inline dvec3 operator - (const dvec3& x, double s) { return x - dvec3(s); }
inline dvec3 operator - (double s, const dvec3& x) { return dvec3(s) - x; }
inline dvec3 operator * (const dvec3& x, double s) { return x * dvec3(s); }
inline dvec3 operator * (double s, const dvec3& x) { return dvec3(s) * x; }
inline dvec3 operator / (const dvec3& x, double s) { return x / dvec3(s); }
inline dvec3 operator / (double s, const dvec3& x) { return dvec3(s) / x; }
inline dvec3& operator += (dvec3& x, const dvec3& y) { x = x + y; return x; }
inline dvec3& operator -= (dvec3& x, const dvec3& y) { x = x - y; return x; }
inline dvec3& operator *= (dvec3& x, const dvec3& y) { x = x * y; return x; }
inline dvec3& operator /= (dvec3& x, const dvec3& y) { x = x / y; return x; }
inline dvec3& operator += (dvec3& x, double s) { x = x + s; return x; }
inline dvec3& operator -= (dvec3& x, double s) { x = x - s; return x; }
inline dvec3& operator *= (dvec3& x, double s) { x = x * s; return x; }
inline dvec3& operator /= (dvec3& x, double s) { x = x / s; return x; }

//-----------------------------------------------------------------------------
// Additional non-standard helpers
//-----------------------------------------------------------------------------
inline vec3 vec3_xAxis() { return vec3(1.0f, 0.0f, 0.0f); }
inline vec3 vec3_yAxis() { return vec3(0.0f, 1.0f, 0.0f); }
inline vec3 vec3_zAxis() { return vec3(0.0f, 0.0f, 1.0f); }
inline vec4 vec4_xAxis() { return vec4(1.0f, 0.0f, 0.0f, 0.0f); }
inline vec4 vec4_yAxis() { return vec4(0.0f, 1.0f, 0.0f, 0.0f); }
inline vec4 vec4_zAxis() { return vec4(0.0f, 0.0f, 1.0f, 0.0f); }
inline vec4 vec4_wAxis() { return vec4(0.0f, 0.0f, 0.0f, 1.0f); }

mat2 rotation(float angleRadians);

mat3 axisAngleRotation(const vec3& axis, float angleRadians);
mat3 axisCrossAxisRotation(const vec3& trueAxis, const vec3& hintAxis, uint32_t i, uint32_t j, uint32_t k);

inline mat3 xyRotation(const vec3& x, const vec3& y) { return axisCrossAxisRotation(x, y, 0, 1, 2); }
inline mat3 xzRotation(const vec3& x, const vec3& z) { return axisCrossAxisRotation(x, z, 0, 2, 1); }
inline mat3 yzRotation(const vec3& y, const vec3& z) { return axisCrossAxisRotation(y, z, 1, 2, 0); }
inline mat3 yxRotation(const vec3& y, const vec3& x) { return axisCrossAxisRotation(y, x, 1, 0, 2); }
inline mat3 zxRotation(const vec3& z, const vec3& x) { return axisCrossAxisRotation(z, x, 2, 0, 1); }
inline mat3 zyRotation(const vec3& z, const vec3& y) { return axisCrossAxisRotation(z, y, 2, 1, 0); }

mat3 transform(const mat2& upper2x2, const vec2& translation);
mat4 transform(const mat3& upper3x3, const vec3& translation);
mat4 projection(float near, float verticalFovDegrees, float width, float height);
mat4 ortho(float left, float right, float bottom, float top, float near, float far);

//-----------------------------------------------------------------------------
// Debug printing
//-----------------------------------------------------------------------------
#if DEBUG_OUT
namespace DVAL_ADL {
	inline int Stringify(char* buffer, size_t count, const vec2& v) { return snprintf(buffer, count, "(%f, %f)", v[0], v[1]); }
	inline int Stringify(char* buffer, size_t count, const vec3& v) { return snprintf(buffer, count, "(%f, %f, %f)", v[0], v[1], v[2]); }
	inline int Stringify(char* buffer, size_t count, const vec4& v) { return snprintf(buffer, count, "(%f, %f, %f, %f)", v[0], v[1], v[2], v[3]); }
	inline int Stringify(char* buffer, size_t count, const mat2& m) { return snprintf(buffer, count, "[%f, %f]\n[%f, %f]\n", m[0][0], m[1][0], m[0][1], m[1][1]); }
	inline int Stringify(char* buffer, size_t count, const mat3& m) { return snprintf(buffer, count, "[%f, %f, %f]\n[%f, %f, %f]\n[%f, %f, %f]\n", m[0][0], m[1][0], m[2][0], m[0][1], m[1][1], m[2][1], m[0][2], m[1][2], m[2][2]); }
	inline int Stringify(char* buffer, size_t count, const mat4& m) { return snprintf(buffer, count, "[%f, %f, %f, %f]\n[%f, %f, %f, %f]\n[%f, %f, %f, %f]\n[%f, %f, %f, %f]\n", m[0][0], m[1][0], m[2][0], m[3][0], m[0][1], m[1][1], m[2][1], m[3][1], m[0][2], m[1][2], m[2][2], m[3][2], m[0][3], m[1][3], m[2][3], m[3][3]); }
}
#endif

//-----------------------------------------------------------------------------
// Testing
//-----------------------------------------------------------------------------
inline void GlslVectormathTest() {
#if 0
	vec2 v2a(0.1f, 0.2f), v2b(1.1f, 1.2f);
	vec3 v3a(2.1f, 2.2f, 2.3f), v3b(3.1f, 3.2f, 3.3f);
	vec4 v4a(4.1f, 4.2f, 4.3f, 4.4f), v4b(5.1f, 5.2f, 5.3f, 5.4f);
	vec4 v4c = vec4(v2a, v2b);
	v2a.xx = v2b.yy;
	v3a.rbg += v2b.xyx;
	v3a.rbg *= 2.0f;
	printf("v2a: %f, %f\n", v2a.x, v2a.y);
	printf("v3a: %f, %f, %f\n", v3a.x, v3a.y, v3a.z);
	printf("v4c: %f, %f, %f, %f\n", v4c.x, v4c.y, v4c.z, v4c.w);

	vec2 v2c(1.0f, 2.0f);
	vec2 v2d = -v2c;
	v2c += 3.0f * v2d;
	printf("v2c: %f, %f\n", v2c.x, v2c.y);

	mat3 m3(vec3(7, 0, -3), vec3(2, 3, 4), vec3(1, -1, -2));
	mat3 m3i = inverse(m3) * m3;
	printf("m3i: %.2f %.2f %.2f\n     %.2f %.2f %.2f\n     %.2f %.2f %.2f\n",
		m3i[0][0], m3i[1][0], m3i[2][0],
		m3i[0][1], m3i[1][1], m3i[2][1],
		m3i[0][2], m3i[1][2], m3i[2][2]);

	mat4 m4(vec4(1, 0, 2, 1), vec4(1, 3, 3, 0), vec4(1, 1, 1, 2), vec4(0, 2, 0, 1));
	mat4 m4i = inverse(m4) * m4;
	printf("m4i: %.2f %.2f %.2f %.2f\n     %.2f %.2f %.2f %.2f\n     %.2f %.2f %.2f %.2f\n     %.2f %.2f %.2f %.2f\n",
		m4i[0][0], m4i[1][0], m4i[2][0], m4i[3][0],
		m4i[0][1], m4i[1][1], m4i[2][1], m4i[3][1],
		m4i[0][2], m4i[1][2], m4i[2][2], m4i[3][2],
		m4i[0][3], m4i[1][3], m4i[2][3], m4i[3][3]);
#endif
}
