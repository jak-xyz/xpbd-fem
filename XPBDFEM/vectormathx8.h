//-----------------------------------------------------------------------------
// Partial AVX2 structure of arrays implementation of vectormath.
//-----------------------------------------------------------------------------
#pragma once

#include <immintrin.h>
#include <stddef.h>
#include <stdint.h>
#include <cmath>
#include "vectormath.h"

class floatx8;
class vec2x8;
class vec3x8;
class vec4x8;
class mat2x8;
class mat3x8;
class mat4x8;

#ifndef VECTORCALL
#define VECTORCALL __vectorcall
#endif

//-----------------------------------------------------------------------------
// scalar classes
//-----------------------------------------------------------------------------
class floatx8 {
public:
	__m256 v;

	inline explicit floatx8() {}
	inline explicit floatx8(float s) { v = _mm256_set1_ps(s); }
	inline explicit floatx8(float f0, float f1, float f2, float f3, float f4, float f5, float f6, float f7) { v = _mm256_setr_ps(f0, f1, f2, f3, f4, f5, f6, f7); }
	inline explicit floatx8(const __m256& x) { v = x; }

	inline floatx8& operator = (const floatx8& other) { v = other.v; return *this; }
	inline floatx8& operator = (const __m256& x) { v = x; return *this; }

	floatx8& load(const float* p) { v = _mm256_loadu_ps(p); return *this; }
	void store(float* p) const { _mm256_storeu_ps(p, v); }
};

class fboolx8 {
public:
	__m256 v;

	inline explicit fboolx8() {}
	inline explicit fboolx8(const __m256& x) { v = x; }

	inline fboolx8& operator = (const fboolx8& other) { v = other.v; return *this; }
	inline fboolx8& operator = (const __m256& x) { v = x; return *this; }
};

class doublex8 {
public:
	__m256d v0;
	__m256d v1;

	inline explicit doublex8() {}
	inline explicit doublex8(double s) { v0 = _mm256_set1_pd(s); v1 = _mm256_set1_pd(s); }
	inline explicit doublex8(double d0, double d1, double d2, double d3, double d4, double d5, double d6, double d7) { v0 = _mm256_setr_pd(d0, d1, d2, d3); v1 = _mm256_setr_pd(d4, d5, d6, d7); }
	inline explicit doublex8(const __m256d& x, const __m256d& y) { v0 = x; v1 = y; }
	inline explicit doublex8(const floatx8& x) {
		v0 = _mm256_cvtps_pd(_mm256_extractf128_ps(x.v, 0));
		v1 = _mm256_cvtps_pd(_mm256_extractf128_ps(x.v, 1));
	}
	inline explicit operator floatx8() {
		return floatx8(_mm256_set_m128(_mm256_cvtpd_ps(v1), _mm256_cvtpd_ps(v0)));
	}

	doublex8& load(const double* p) { v0 = _mm256_loadu_pd(p); v1 = _mm256_loadu_pd(p + 4); return *this; }
	void store(double* p) const { _mm256_storeu_pd(p, v0); _mm256_storeu_pd(p + 4, v1); }
};

//-----------------------------------------------------------------------------
// vec classes
//-----------------------------------------------------------------------------
class vec2x8 {
public:
	union {
		floatx8 data[2];
		struct { floatx8 x, y; };
		struct { floatx8 r, g; };
		struct { floatx8 s, t; };
	};

	inline explicit vec2x8() {}
	inline explicit vec2x8(floatx8 s) { data[0] = s; data[1] = s; }
	inline explicit vec2x8(floatx8 x_, floatx8 y_) { data[0] = x_; data[1] = y_; }

	inline vec2x8& operator = (const vec2x8& other) { x = other.x; y = other.y; return *this; }
	inline floatx8& operator [] (size_t i) { return data[i]; }
	inline floatx8 operator [] (size_t i) const { return data[i]; }

	vec2x8& load(const vec2* v0, const vec2* v1, const vec2* v2, const vec2* v3, const vec2* v4, const vec2* v5, const vec2* v6, const vec2* v7) {
		x = floatx8(v0->x, v1->x, v2->x, v3->x, v4->x, v5->x, v6->x, v7->x);
		y = floatx8(v0->y, v1->y, v2->y, v3->y, v4->y, v5->y, v6->y, v7->y);
		return *this;
	}
	void store(vec2* v0, vec2* v1, vec2* v2, vec2* v3, vec2* v4, vec2* v5, vec2* v6, vec2* v7) const {
		float x[8], y[8];
		data[0].store(x);
		data[1].store(y);
		v0->x = x[0]; v1->x = x[1]; v2->x = x[2]; v3->x = x[3]; v4->x = x[4]; v5->x = x[5]; v6->x = x[6]; v7->x = x[7];
		v0->y = y[0]; v1->y = y[1]; v2->y = y[2]; v3->y = y[3]; v4->y = y[4]; v5->y = y[5]; v6->y = y[6]; v7->y = y[7];
	}
};

class vec3x8 {
public:
	union {
		floatx8 data[3];
		struct { floatx8 x, y, z; };
		struct { floatx8 r, g, b; };
		struct { floatx8 s, t, p; };
	};

	inline explicit vec3x8() {}
	inline explicit vec3x8(floatx8 s) { data[0] = s; data[1] = s; data[2] = s; }
	inline explicit vec3x8(floatx8 x_, floatx8 y_, floatx8 z_) { data[0] = x_; data[1] = y_;  data[2] = z_; }
	inline explicit vec3x8(const vec2x8& xy, floatx8 z_) { data[0] = xy.x; data[1] = xy.y;  data[2] = z_; }
	inline explicit vec3x8(floatx8 x_, const vec2x8& yz) { data[0] = x_; data[1] = yz.x;  data[2] = yz.y; }

	inline vec3x8& operator = (const vec3x8& other) { x = other.x; y = other.y; z = other.z; return *this; }
	inline floatx8& operator [] (size_t i) { return data[i]; }
	inline floatx8 operator [] (size_t i) const { return data[i]; }

	vec3x8& load(const vec3* v0, const vec3* v1, const vec3* v2, const vec3* v3, const vec3* v4, const vec3* v5, const vec3* v6, const vec3* v7) {
		x = floatx8(v0->x, v1->x, v2->x, v3->x, v4->x, v5->x, v6->x, v7->x);
		y = floatx8(v0->y, v1->y, v2->y, v3->y, v4->y, v5->y, v6->y, v7->y);
		z = floatx8(v0->z, v1->z, v2->z, v3->z, v4->z, v5->z, v6->z, v7->z);
		return *this;
	}
	void store(vec3* v0, vec3* v1, vec3* v2, vec3* v3, vec3* v4, vec3* v5, vec3* v6, vec3* v7) const {
		float x[8], y[8], z[8];
		data[0].store(x);
		data[1].store(y);
		data[2].store(z);
		v0->x = x[0]; v1->x = x[1]; v2->x = x[2]; v3->x = x[3]; v4->x = x[4]; v5->x = x[5]; v6->x = x[6]; v7->x = x[7];
		v0->y = y[0]; v1->y = y[1]; v2->y = y[2]; v3->y = y[3]; v4->y = y[4]; v5->y = y[5]; v6->y = y[6]; v7->y = y[7];
		v0->z = z[0]; v1->z = z[1]; v2->z = z[2]; v3->z = z[3]; v4->z = z[4]; v5->z = z[5]; v6->z = z[6]; v7->z = z[7];
	}
};

class vec4x8 {
public:
	union {
		floatx8 data[4];
		struct { floatx8 x, y, z, w; };
		struct { floatx8 r, g, b, a; };
		struct { floatx8 s, t, p, q; };
	};

	inline explicit vec4x8() {}
	inline explicit vec4x8(floatx8 s) { data[0] = s; data[1] = s; data[2] = s; data[3] = s; }
	inline explicit vec4x8(floatx8 x_, floatx8 y_, floatx8 z_, floatx8 w_) { data[0] = x_; data[1] = y_;  data[2] = z_; data[3] = w_; }
	inline explicit vec4x8(const vec2x8& xy, floatx8 z_, floatx8 w_) { data[0] = xy.x; data[1] = xy.y;  data[2] = z_; data[3] = w_; }
	inline explicit vec4x8(floatx8 x_, const vec2x8& yz, floatx8 w_) { data[0] = x_; data[1] = yz.x;  data[2] = yz.y; data[3] = w_; }
	inline explicit vec4x8(floatx8 x_, floatx8 y_, const vec2x8& zw) { data[0] = x_; data[1] = y_;  data[2] = zw.x; data[3] = zw.y; }
	inline explicit vec4x8(const vec2x8& xy, const vec2x8& zw) { data[0] = xy.x; data[1] = xy.y;  data[2] = zw.x; data[3] = zw.y; }
	inline explicit vec4x8(const vec3x8& xyz, floatx8 w_) { data[0] = xyz.x; data[1] = xyz.y;  data[2] = xyz.z; data[3] = w_; }
	inline explicit vec4x8(floatx8 x_, const vec3x8& yzw) { data[0] = x_; data[1] = yzw.x;  data[2] = yzw.y; data[3] = yzw.z; }

	inline vec4x8& operator = (const vec4x8& other) { x = other.x; y = other.y; z = other.z; w = other.w; return *this; }
	inline floatx8& operator [] (size_t i) { return data[i]; }
	inline floatx8 operator [] (size_t i) const { return data[i]; }

	vec4x8& load(const vec4* v0, const vec4* v1, const vec4* v2, const vec4* v3, const vec4* v4, const vec4* v5, const vec4* v6, const vec4* v7) {
		x = floatx8(v0->x, v1->x, v2->x, v3->x, v4->x, v5->x, v6->x, v7->x);
		y = floatx8(v0->y, v1->y, v2->y, v3->y, v4->y, v5->y, v6->y, v7->y);
		z = floatx8(v0->z, v1->z, v2->z, v3->z, v4->z, v5->z, v6->z, v7->z);
		w = floatx8(v0->w, v1->w, v2->w, v3->w, v4->w, v5->w, v6->w, v7->w);
		return *this;
	}
	void store(vec4* v0, vec4* v1, vec4* v2, vec4* v3, vec4* v4, vec4* v5, vec4* v6, vec4* v7) const {
		float x[8], y[8], z[8], w[8];
		data[0].store(x);
		data[1].store(y);
		data[2].store(z);
		data[2].store(w);
		v0->x = x[0]; v1->x = x[1]; v2->x = x[2]; v3->x = x[3]; v4->x = x[4]; v5->x = x[5]; v6->x = x[6]; v7->x = x[7];
		v0->y = y[0]; v1->y = y[1]; v2->y = y[2]; v3->y = y[3]; v4->y = y[4]; v5->y = y[5]; v6->y = y[6]; v7->y = y[7];
		v0->z = z[0]; v1->z = z[1]; v2->z = z[2]; v3->z = z[3]; v4->z = z[4]; v5->z = z[5]; v6->z = z[6]; v7->z = z[7];
		v0->w = w[0]; v1->w = w[1]; v2->w = w[2]; v3->w = w[3]; v4->w = w[4]; v5->w = w[5]; v6->w = w[6]; v7->w = w[7];
	}
};

//-----------------------------------------------------------------------------
// mat classes
//-----------------------------------------------------------------------------
class mat2x8 {
public:
	vec2x8 cols[2];

	inline mat2x8() {}
	inline explicit mat2x8(const vec2x8& c0, const vec2x8& c1) { cols[0] = c0; cols[1] = c1; }

	inline vec2x8& operator [] (size_t i) { return cols[i]; }
	inline vec2x8 operator [] (size_t i) const { return cols[i]; }

	mat2x8& load(const mat2* m0, const mat2* m1, const mat2* m2, const mat2* m3, const mat2* m4, const mat2* m5, const mat2* m6, const mat2* m7) {
		cols[0].load(&(*m0)[0], &(*m1)[0], &(*m2)[0], &(*m3)[0], &(*m4)[0], &(*m5)[0], &(*m6)[0], &(*m7)[0]);
		cols[1].load(&(*m0)[1], &(*m1)[1], &(*m2)[1], &(*m3)[1], &(*m4)[1], &(*m5)[1], &(*m6)[1], &(*m7)[1]);
		return *this;
	}
	void store(mat2* m0, mat2* m1, mat2* m2, mat2* m3, mat2* m4, mat2* m5, mat2* m6, mat2* m7) const {
		cols[0].store(&(*m0)[0], &(*m1)[0], &(*m2)[0], &(*m3)[0], &(*m4)[0], &(*m5)[0], &(*m6)[0], &(*m7)[0]);
		cols[1].store(&(*m0)[1], &(*m1)[1], &(*m2)[1], &(*m3)[1], &(*m4)[1], &(*m5)[1], &(*m6)[1], &(*m7)[1]);
	}
};

class mat3x8 {
public:
	vec3x8 cols[3];

	inline mat3x8() {}
	inline explicit mat3x8(const vec3x8& c0, const vec3x8& c1, const vec3x8& c2) { cols[0] = c0; cols[1] = c1; cols[2] = c2; }

	inline vec3x8& operator [] (size_t i) { return cols[i]; }
	inline vec3x8 operator [] (size_t i) const { return cols[i]; }

	mat3x8& load(const mat3* m0, const mat3* m1, const mat3* m2, const mat3* m3, const mat3* m4, const mat3* m5, const mat3* m6, const mat3* m7) {
		cols[0].load(&(*m0)[0], &(*m1)[0], &(*m2)[0], &(*m3)[0], &(*m4)[0], &(*m5)[0], &(*m6)[0], &(*m7)[0]);
		cols[1].load(&(*m0)[1], &(*m1)[1], &(*m2)[1], &(*m3)[1], &(*m4)[1], &(*m5)[1], &(*m6)[1], &(*m7)[1]);
		cols[2].load(&(*m0)[2], &(*m1)[2], &(*m2)[2], &(*m3)[2], &(*m4)[2], &(*m5)[2], &(*m6)[2], &(*m7)[2]);
		return *this;
	}
	void store(mat3* m0, mat3* m1, mat3* m2, mat3* m3, mat3* m4, mat3* m5, mat3* m6, mat3* m7) const {
		cols[0].store(&(*m0)[0], &(*m1)[0], &(*m2)[0], &(*m3)[0], &(*m4)[0], &(*m5)[0], &(*m6)[0], &(*m7)[0]);
		cols[1].store(&(*m0)[1], &(*m1)[1], &(*m2)[1], &(*m3)[1], &(*m4)[1], &(*m5)[1], &(*m6)[1], &(*m7)[1]);
		cols[2].store(&(*m0)[2], &(*m1)[2], &(*m2)[2], &(*m3)[2], &(*m4)[2], &(*m5)[2], &(*m6)[2], &(*m7)[2]);
	}
};

class mat4x8 {
public:
	vec4x8 cols[4];

	inline mat4x8() {}
	inline explicit mat4x8(const vec4x8& c0, const vec4x8& c1, const vec4x8& c2, const vec4x8& c3) { cols[0] = c0; cols[1] = c1; cols[2] = c2; cols[3] = c3; }

	inline vec4x8& operator [] (size_t i) { return cols[i]; }
	inline vec4x8 operator [] (size_t i) const { return cols[i]; }

	mat4x8& load(const mat4* m0, const mat4* m1, const mat4* m2, const mat4* m3, const mat4* m4, const mat4* m5, const mat4* m6, const mat4* m7) {
		cols[0].load(&(*m0)[0], &(*m1)[0], &(*m2)[0], &(*m3)[0], &(*m4)[0], &(*m5)[0], &(*m6)[0], &(*m7)[0]);
		cols[1].load(&(*m0)[1], &(*m1)[1], &(*m2)[1], &(*m3)[1], &(*m4)[1], &(*m5)[1], &(*m6)[1], &(*m7)[1]);
		cols[2].load(&(*m0)[2], &(*m1)[2], &(*m2)[2], &(*m3)[2], &(*m4)[2], &(*m5)[2], &(*m6)[2], &(*m7)[2]);
		cols[3].load(&(*m0)[3], &(*m1)[3], &(*m2)[3], &(*m3)[3], &(*m4)[3], &(*m5)[3], &(*m6)[3], &(*m7)[3]);
		return *this;
	}
	void store(mat4* m0, mat4* m1, mat4* m2, mat4* m3, mat4* m4, mat4* m5, mat4* m6, mat4* m7) const {
		cols[0].store(&(*m0)[0], &(*m1)[0], &(*m2)[0], &(*m3)[0], &(*m4)[0], &(*m5)[0], &(*m6)[0], &(*m7)[0]);
		cols[1].store(&(*m0)[1], &(*m1)[1], &(*m2)[1], &(*m3)[1], &(*m4)[1], &(*m5)[1], &(*m6)[1], &(*m7)[1]);
		cols[2].store(&(*m0)[2], &(*m1)[2], &(*m2)[2], &(*m3)[2], &(*m4)[2], &(*m5)[2], &(*m6)[2], &(*m7)[2]);
		cols[3].store(&(*m0)[3], &(*m1)[3], &(*m2)[3], &(*m3)[3], &(*m4)[3], &(*m5)[3], &(*m6)[3], &(*m7)[3]);
	}
};

//-----------------------------------------------------------------------------
// Helpers for defining vec interfaces
//-----------------------------------------------------------------------------
inline VEC_FORCEINLINE vec2x8 VECTORCALL vec_map(floatx8(VECTORCALL*f)(floatx8), vec2x8 x) { return vec2x8(f(x[0]), f(x[1])); }
inline VEC_FORCEINLINE vec2x8 VECTORCALL vec_map(floatx8(VECTORCALL*f)(floatx8, floatx8), vec2x8 x, vec2x8 y) { return vec2x8(f(x[0], y[0]), f(x[1], y[1])); }
inline VEC_FORCEINLINE vec2x8 VECTORCALL vec_map(floatx8(VECTORCALL*f)(floatx8, floatx8, floatx8), vec2x8 x, vec2x8 y, vec2x8 z) { return vec2x8(f(x[0], y[0], z[0]), f(x[1], y[1], z[1])); }
inline VEC_FORCEINLINE vec3x8 VECTORCALL vec_map(floatx8(VECTORCALL*f)(floatx8), vec3x8 x) { return vec3x8(f(x[0]), f(x[1]), f(x[2])); }
inline VEC_FORCEINLINE vec3x8 VECTORCALL vec_map(floatx8(VECTORCALL*f)(floatx8, floatx8), vec3x8 x, vec3x8 y) { return vec3x8(f(x[0], y[0]), f(x[1], y[1]), f(x[2], y[2])); }
inline VEC_FORCEINLINE vec3x8 VECTORCALL vec_map(floatx8(VECTORCALL*f)(floatx8, floatx8, floatx8), vec3x8 x, vec3x8 y, vec3x8 z) { return vec3x8(f(x[0], y[0], z[0]), f(x[1], y[1], z[1]), f(x[2], y[2], z[2])); }
inline VEC_FORCEINLINE vec4x8 VECTORCALL vec_map(floatx8(VECTORCALL*f)(floatx8), vec4x8 x) { return vec4x8(f(x[0]), f(x[1]), f(x[2]), f(x[3])); }
inline VEC_FORCEINLINE vec4x8 VECTORCALL vec_map(floatx8(VECTORCALL*f)(floatx8, floatx8), vec4x8 x, vec4x8 y) { return vec4x8(f(x[0], y[0]), f(x[1], y[1]), f(x[2], y[2]), f(x[3], y[3])); }
inline VEC_FORCEINLINE vec4x8 VECTORCALL vec_map(floatx8(VECTORCALL*f)(floatx8, floatx8, floatx8), vec4x8 x, vec4x8 y, vec4x8 z) { return vec4x8(f(x[0], y[0], z[0]), f(x[1], y[1], z[1]), f(x[2], y[2], z[2]), f(x[3], y[3], z[3])); }

//-----------------------------------------------------------------------------
// floatx8 operators and functions
//-----------------------------------------------------------------------------
inline floatx8 VECTORCALL operator + (floatx8 x) { return x; }
inline floatx8 VECTORCALL operator - (floatx8 x) { return floatx8(_mm256_sub_ps(_mm256_set1_ps(0.0f), x.v)); }
inline floatx8 VECTORCALL operator + (floatx8 x, floatx8 y) { return floatx8(_mm256_add_ps(x.v, y.v)); }
inline floatx8 VECTORCALL operator - (floatx8 x, floatx8 y) { return floatx8(_mm256_sub_ps(x.v, y.v)); }
inline floatx8 VECTORCALL operator * (floatx8 x, floatx8 y) { return floatx8(_mm256_mul_ps(x.v, y.v)); }
inline floatx8 VECTORCALL operator / (floatx8 x, floatx8 y) { return floatx8(_mm256_div_ps(x.v, y.v)); }
inline floatx8 VECTORCALL operator + (floatx8 x, float s) { return x + floatx8(s); }
inline floatx8 VECTORCALL operator + (float s, floatx8 x) { return floatx8(s) + x; }
inline floatx8 VECTORCALL operator - (floatx8 x, float s) { return x - floatx8(s); }
inline floatx8 VECTORCALL operator - (float s, floatx8 x) { return floatx8(s) - x; }
inline floatx8 VECTORCALL operator * (floatx8 x, float s) { return x * floatx8(s); }
inline floatx8 VECTORCALL operator * (float s, floatx8 x) { return floatx8(s) * x; }
inline floatx8 VECTORCALL operator / (floatx8 x, float s) { return x / floatx8(s); }
inline floatx8 VECTORCALL operator / (float s, floatx8 x) { return floatx8(s) / x; }
inline floatx8& VECTORCALL operator += (floatx8& x, floatx8 y) { x = x + y; return x; }
inline floatx8& VECTORCALL operator -= (floatx8& x, floatx8 y) { x = x - y; return x; }
inline floatx8& VECTORCALL operator *= (floatx8& x, floatx8 y) { x = x * y; return x; }
inline floatx8& VECTORCALL operator /= (floatx8& x, floatx8 y) { x = x / y; return x; }
inline floatx8& VECTORCALL operator += (floatx8& x, float s) { x = x + s; return x; }
inline floatx8& VECTORCALL operator -= (floatx8& x, float s) { x = x - s; return x; }
inline floatx8& VECTORCALL operator *= (floatx8& x, float s) { x = x * s; return x; }
inline floatx8& VECTORCALL operator /= (floatx8& x, float s) { x = x / s; return x; }

inline fboolx8 VECTORCALL operator == (floatx8 x, floatx8 y) { return fboolx8(_mm256_cmp_ps(x.v, y.v, _CMP_EQ_OQ)); }
inline fboolx8 VECTORCALL operator != (floatx8 x, floatx8 y) { return fboolx8(_mm256_cmp_ps(x.v, y.v, _CMP_NEQ_OQ)); }
inline fboolx8 VECTORCALL operator < (floatx8 x, floatx8 y) { return fboolx8(_mm256_cmp_ps(x.v, y.v, _CMP_LT_OQ)); }
inline fboolx8 VECTORCALL operator <= (floatx8 x, floatx8 y) { return fboolx8(_mm256_cmp_ps(x.v, y.v, _CMP_LE_OQ)); }
inline fboolx8 VECTORCALL operator > (floatx8 x, floatx8 y) { return fboolx8(_mm256_cmp_ps(x.v, y.v, _CMP_GT_OQ)); }
inline fboolx8 VECTORCALL operator >= (floatx8 x, floatx8 y) { return fboolx8(_mm256_cmp_ps(x.v, y.v, _CMP_GE_OQ)); }
inline fboolx8 VECTORCALL operator == (fboolx8 x, fboolx8 y) { return fboolx8(_mm256_cmp_ps(x.v, y.v, _CMP_EQ_OQ)); }
inline fboolx8 VECTORCALL operator != (fboolx8 x, fboolx8 y) { return fboolx8(_mm256_cmp_ps(x.v, y.v, _CMP_NEQ_OQ)); }
inline fboolx8 VECTORCALL operator && (fboolx8 x, fboolx8 y) { return fboolx8(_mm256_and_ps(x.v, y.v)); }
inline fboolx8 VECTORCALL operator || (fboolx8 x, fboolx8 y) { return fboolx8(_mm256_or_ps(x.v, y.v)); }
inline bool VECTORCALL all(fboolx8 x) { return 0xff == _mm256_movemask_ps(x.v); }
inline bool VECTORCALL any(fboolx8 x) { return _mm256_movemask_ps(x.v); }
inline floatx8 VECTORCALL select(fboolx8 cond, floatx8 onTrue, floatx8 onFalse) { return floatx8(_mm256_blendv_ps(onFalse.v, onTrue.v, cond.v)); };

inline floatx8 VECTORCALL abs(floatx8 x) { return floatx8(_mm256_and_ps(_mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff)), x.v)); }
inline floatx8 VECTORCALL ceil(floatx8 x) { return floatx8(_mm256_ceil_ps(x.v)); }
inline floatx8 VECTORCALL clamp(floatx8 x, floatx8 fmin, floatx8 fmax) { return floatx8(_mm256_min_ps(_mm256_max_ps(fmin.v, x.v), fmax.v)); }
inline floatx8 VECTORCALL floor(floatx8 x) { return floatx8(_mm256_floor_ps(x.v)); }
inline floatx8 VECTORCALL fma(floatx8 a, floatx8 b, floatx8 c) { return floatx8(_mm256_fmadd_ps(a.v, b.v, c.v)); }
inline floatx8 VECTORCALL fract(floatx8 x) { return x - floor(x); }
inline floatx8 VECTORCALL log(floatx8 x) { return floatx8(_mm256_log_ps(x.v)); }
inline floatx8 VECTORCALL max(floatx8 x, floatx8 y) { return floatx8(_mm256_max_ps(x.v, y.v)); }
inline floatx8 VECTORCALL min(floatx8 x, floatx8 y) { return floatx8(_mm256_min_ps(x.v, y.v)); }
inline floatx8 VECTORCALL mix(floatx8 x, floatx8 y, const floatx8& a) { return x * (1.0f - a) + y * a; }
inline floatx8 VECTORCALL rcp_est(floatx8 x) { return floatx8(_mm256_rcp_ps(x.v)); }
inline floatx8 VECTORCALL rsqrt_est(floatx8 x) { return floatx8(_mm256_rsqrt_ps(x.v)); }
inline floatx8 VECTORCALL round(floatx8 x) { return floatx8(_mm256_round_ps(x.v, (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC))); }
inline floatx8 VECTORCALL smoothstep(floatx8 edge0, floatx8 edge1, floatx8 x) { floatx8 t = clamp((x - edge0) / (edge1 - edge0), floatx8(0.0f), floatx8(1.0f)); return t * t * (3.0f - 2.0f * t); }
inline floatx8 VECTORCALL sqrt(floatx8 x) { return floatx8(_mm256_sqrt_ps(x.v)); }
inline floatx8 VECTORCALL inversesqrt(floatx8 x) { return 1.0f / sqrt(x); }
inline floatx8 VECTORCALL step(floatx8 edge, floatx8 x) { return select(x < edge, floatx8(0.0f), floatx8(1.0f)); }
inline floatx8 VECTORCALL trunc(floatx8 x) { return floatx8(_mm256_round_ps(x.v, (_MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC))); }

//-----------------------------------------------------------------------------
// doublex8 operators and functions
//-----------------------------------------------------------------------------
inline doublex8 VECTORCALL operator + (doublex8 x) { return x; }
inline doublex8 VECTORCALL operator - (doublex8 x) { return doublex8(_mm256_sub_pd(_mm256_set1_pd(0.0), x.v0), _mm256_sub_pd(_mm256_set1_pd(0.0), x.v1)); }
inline doublex8 VECTORCALL operator + (doublex8 x, doublex8 y) { return doublex8(_mm256_add_pd(x.v0, y.v0), _mm256_add_pd(x.v1, y.v1)); }
inline doublex8 VECTORCALL operator - (doublex8 x, doublex8 y) { return doublex8(_mm256_sub_pd(x.v0, y.v0), _mm256_sub_pd(x.v1, y.v1)); }
inline doublex8 VECTORCALL operator * (doublex8 x, doublex8 y) { return doublex8(_mm256_mul_pd(x.v0, y.v0), _mm256_mul_pd(x.v1, y.v1)); }
inline doublex8 VECTORCALL operator / (doublex8 x, doublex8 y) { return doublex8(_mm256_div_pd(x.v0, y.v0), _mm256_div_pd(x.v1, y.v1)); }
inline doublex8 VECTORCALL operator + (doublex8 x, double s) { return x + doublex8(s); }
inline doublex8 VECTORCALL operator + (double s, doublex8 x) { return doublex8(s) + x; }
inline doublex8 VECTORCALL operator - (doublex8 x, double s) { return x - doublex8(s); }
inline doublex8 VECTORCALL operator - (double s, doublex8 x) { return doublex8(s) - x; }
inline doublex8 VECTORCALL operator * (doublex8 x, double s) { return x * doublex8(s); }
inline doublex8 VECTORCALL operator * (double s, doublex8 x) { return doublex8(s) * x; }
inline doublex8 VECTORCALL operator / (doublex8 x, double s) { return x / doublex8(s); }
inline doublex8 VECTORCALL operator / (double s, doublex8 x) { return doublex8(s) / x; }
inline doublex8& VECTORCALL operator += (doublex8& x, doublex8 y) { x = x + y; return x; }
inline doublex8& VECTORCALL operator -= (doublex8& x, doublex8 y) { x = x - y; return x; }
inline doublex8& VECTORCALL operator *= (doublex8& x, doublex8 y) { x = x * y; return x; }
inline doublex8& VECTORCALL operator /= (doublex8& x, doublex8 y) { x = x / y; return x; }
inline doublex8& VECTORCALL operator += (doublex8& x, double s) { x = x + s; return x; }
inline doublex8& VECTORCALL operator -= (doublex8& x, double s) { x = x - s; return x; }
inline doublex8& VECTORCALL operator *= (doublex8& x, double s) { x = x * s; return x; }
inline doublex8& VECTORCALL operator /= (doublex8& x, double s) { x = x / s; return x; }

inline doublex8 VECTORCALL max(doublex8 x, doublex8 y) { return doublex8(_mm256_max_pd(x.v0, y.v0), _mm256_max_pd(x.v1, y.v1)); }
inline doublex8 VECTORCALL sqrt(doublex8 x) { return doublex8(_mm256_sqrt_pd(x.v0), _mm256_sqrt_pd(x.v1)); }

//-----------------------------------------------------------------------------
// vec2x8 operators and functions
//-----------------------------------------------------------------------------
inline vec2x8 VECTORCALL operator + (vec2x8 x) { return x; }
inline vec2x8 VECTORCALL operator - (vec2x8 x) { return vec2x8(-x.data[0], -x.data[1]); }
inline vec2x8 VECTORCALL operator + (vec2x8 x, vec2x8 y) { return vec2x8(x.data[0] + y.data[0], x.data[1] + y.data[1]); }
inline vec2x8 VECTORCALL operator - (vec2x8 x, vec2x8 y) { return vec2x8(x.data[0] - y.data[0], x.data[1] - y.data[1]); }
inline vec2x8 VECTORCALL operator * (vec2x8 x, vec2x8 y) { return vec2x8(x.data[0] * y.data[0], x.data[1] * y.data[1]); }
inline vec2x8 VECTORCALL operator / (vec2x8 x, vec2x8 y) { return vec2x8(x.data[0] / y.data[0], x.data[1] / y.data[1]); }
inline vec2x8 VECTORCALL operator + (vec2x8 x, floatx8 s) { return x + vec2x8(s); }
inline vec2x8 VECTORCALL operator + (floatx8 s, vec2x8 x) { return vec2x8(s) + x; }
inline vec2x8 VECTORCALL operator - (vec2x8 x, floatx8 s) { return x - vec2x8(s); }
inline vec2x8 VECTORCALL operator - (floatx8 s, vec2x8 x) { return vec2x8(s) - x; }
inline vec2x8 VECTORCALL operator * (vec2x8 x, floatx8 s) { return x * vec2x8(s); }
inline vec2x8 VECTORCALL operator * (floatx8 s, vec2x8 x) { return vec2x8(s) * x; }
inline vec2x8 VECTORCALL operator / (vec2x8 x, floatx8 s) { return x / vec2x8(s); }
inline vec2x8 VECTORCALL operator / (floatx8 s, vec2x8 x) { return vec2x8(s) / x; }
inline vec2x8& VECTORCALL operator += (vec2x8& x, vec2x8 y) { x = x + y; return x; }
inline vec2x8& VECTORCALL operator -= (vec2x8& x, vec2x8 y) { x = x - y; return x; }
inline vec2x8& VECTORCALL operator *= (vec2x8& x, vec2x8 y) { x = x * y; return x; }
inline vec2x8& VECTORCALL operator /= (vec2x8& x, vec2x8 y) { x = x / y; return x; }
inline vec2x8& VECTORCALL operator += (vec2x8& x, floatx8 s) { x = x + s; return x; }
inline vec2x8& VECTORCALL operator -= (vec2x8& x, floatx8 s) { x = x - s; return x; }
inline vec2x8& VECTORCALL operator *= (vec2x8& x, floatx8 s) { x = x * s; return x; }
inline vec2x8& VECTORCALL operator /= (vec2x8& x, floatx8 s) { x = x / s; return x; }

inline vec2x8 VECTORCALL abs(vec2x8 x) { return vec_map(abs, x); }
inline vec2x8 VECTORCALL ceil(vec2x8 x) { return vec_map(ceil, x); }
inline vec2x8 VECTORCALL clamp(vec2x8 x, vec2x8 vmin, vec2x8 vmax) { return vec_map(clamp, x, vmin, vmax); }
inline vec2x8 VECTORCALL clamp(vec2x8 x, floatx8 smin, floatx8 smax) { return clamp(x, vec2x8(smin), vec2x8(smax)); }
inline floatx8 VECTORCALL distance(vec2x8 x, vec2x8 y) { vec2x8 a = x - y; return sqrt(a[0] * a[0] + a[1] * a[1]); }
inline floatx8 VECTORCALL dot(vec2x8 x, vec2x8 y) { return x[0] * y[0] + x[1] * y[1]; }
inline vec2x8 VECTORCALL floor(vec2x8 x) { return vec_map(floor, x); }
inline vec2x8 VECTORCALL fma(vec2x8 a, vec2x8 b, vec2x8 c) { return vec_map(fma, a, b, c); }
inline vec2x8 VECTORCALL fract(vec2x8 x) { return vec_map(fract, x); }
inline vec2x8 VECTORCALL inversesqrt(vec2x8 x) { return vec_map(inversesqrt, x); }
inline floatx8 VECTORCALL length(vec2x8 x) { return sqrt(dot(x, x)); }
inline vec2x8 VECTORCALL max(vec2x8 x, vec2x8 y) { return vec_map(max, x, y); }
inline vec2x8 VECTORCALL min(vec2x8 x, vec2x8 y) { return vec_map(min, x, y); }
inline vec2x8 VECTORCALL mix(vec2x8 x, vec2x8 y, vec2x8 a) { return x * (vec2x8(floatx8(1.0f)) - a) + y * a; }
inline vec2x8 VECTORCALL mix(vec2x8 x, vec2x8 y, floatx8 a) { return x * (floatx8(1.0f) - a) + y * a; }
inline vec2x8 VECTORCALL normalize(vec2x8 x) { return x * (1.0f / length(x)); }
inline vec2x8 VECTORCALL safeNormalize(vec2x8 x) { floatx8 len = max(floatx8(1.0e-11f), length(x)); return x * (floatx8(1.0f) / len); }
inline vec2x8 VECTORCALL round(vec2x8 x) { return vec_map(round, x); }
inline vec2x8 VECTORCALL smoothstep(vec2x8 edge0, vec2x8 edge1, vec2x8 x) { return vec_map(smoothstep, edge0, edge1, x); }
inline vec2x8 VECTORCALL sqrt(vec2x8 x) { return vec_map(sqrt, x); }
inline vec2x8 VECTORCALL step(vec2x8 edge, vec2x8 x) { return vec_map(step, edge, x); }
inline vec2x8 VECTORCALL step(floatx8 edge, vec2x8 x) { return step(vec2x8(edge), x); }
inline vec2x8 VECTORCALL trunc(vec2x8 x) { return vec_map(trunc, x); }

//-----------------------------------------------------------------------------
// vec3x8 operators and functions
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// vec4x8 operators and functions
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// mat2x8 operators and functions
//-----------------------------------------------------------------------------
inline mat2x8 VECTORCALL operator + (const mat2x8& x, const mat2x8& y) { return mat2x8(x[0] + y[0], x[1] + y[1]); }
inline mat2x8 VECTORCALL operator - (const mat2x8& x, const mat2x8& y) { return mat2x8(x[0] - y[0], x[1] - y[1]); }
inline mat2x8 VECTORCALL operator * (const mat2x8& m, floatx8 s) { return mat2x8(m[0] * s, m[1] * s); }
inline mat2x8 VECTORCALL operator * (floatx8 s, const mat2x8& m) { return m * s; }
inline mat2x8 VECTORCALL operator / (const mat2x8& m, floatx8 s) { return m * (1.0f / s); }
inline vec2x8 VECTORCALL operator * (const mat2x8& m, vec2x8 v) {
	return vec2x8(
		m[0][0] * v[0] + m[1][0] * v[1],
		m[0][1] * v[0] + m[1][1] * v[1]);
}
inline mat2x8 VECTORCALL operator * (const mat2x8& x, const mat2x8& y) { return mat2x8(x * y[0], x * y[1]); }
inline mat2x8& VECTORCALL operator += (mat2x8& x, const mat2x8& y) { x = x + y; return x; }
inline mat2x8& VECTORCALL operator -= (mat2x8& x, const mat2x8& y) { x = x - y; return x; }
inline mat2x8& VECTORCALL operator *= (mat2x8& m, floatx8 s) { m = m * s; return m; }
inline mat2x8& VECTORCALL operator /= (mat2x8& m, floatx8 s) { m = m / s; return m; }
inline mat2x8& VECTORCALL operator *= (mat2x8& x, const mat2x8& y) { x = x * y; return x; }

inline mat2x8 VECTORCALL transpose(const mat2x8& m) {
	return mat2x8(
		vec2x8(m[0][0], m[1][0]),
		vec2x8(m[0][1], m[1][1]));
}

inline floatx8 VECTORCALL determinant(const mat2x8& m) {
	return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

inline mat2x8 VECTORCALL inverse(const mat2x8& m) {
	return 1.0f / determinant(m) * mat2x8(
		vec2x8(m[1][1], -m[0][1]),
		vec2x8(-m[1][0], m[0][0]));
}

//-----------------------------------------------------------------------------
// mat3x8 operators and functions
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// mat4x8 operators and functions
//-----------------------------------------------------------------------------
