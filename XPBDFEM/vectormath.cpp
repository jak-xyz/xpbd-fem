#include "vectormath.h"
#define _USE_MATH_DEFINES
#include <math.h>

//-----------------------------------------------------------------------------
// mat2 operators and functions
//-----------------------------------------------------------------------------
mat2 transpose(const mat2& m) {
	return mat2(
		vec2(m[0][0], m[1][0]),
		vec2(m[0][1], m[1][1]));
}

float determinant(const mat2& m) {
	return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

mat2 inverse(const mat2& m) {
	return 1.0f / determinant(m) * mat2(
		vec2(m[1][1], -m[0][1]),
		vec2(-m[1][0], m[0][0]));
}

//-----------------------------------------------------------------------------
// mat3 operators and functions
//-----------------------------------------------------------------------------
mat3 transpose(const mat3& m) {
	return mat3(
		vec3(m[0][0], m[1][0], m[2][0]),
		vec3(m[0][1], m[1][1], m[2][1]),
		vec3(m[0][2], m[1][2], m[2][2]));
}

float determinant(const mat3& m) {
	float a = m[0][0], b = m[1][0], c = m[2][0];
	float d = m[0][1], e = m[1][1], f = m[2][1];
	float g = m[0][2], h = m[1][2], i = m[2][2];
	return a * (e*i - f * h) - b * (d*i - f * g) + c * (d*h - e * g);
}

mat3 inverse(const mat3& m) {
	double m00 = m[0][0], m01 = m[0][1], m02 = m[0][2];
	double m10 = m[1][0], m11 = m[1][1], m12 = m[1][2];
	double m20 = m[2][0], m21 = m[2][1], m22 = m[2][2];

	mat3 adj;
	adj[0][0] = (float)+(m11 * m22 - m12 * m21);
	adj[0][1] = (float)-(m01 * m22 - m02 * m21);
	adj[0][2] = (float)+(m01 * m12 - m02 * m11);
	adj[1][0] = (float)-(m10 * m22 - m12 * m20);
	adj[1][1] = (float)+(m00 * m22 - m02 * m20);
	adj[1][2] = (float)-(m00 * m12 - m02 * m10);
	adj[2][0] = (float)+(m10 * m21 - m11 * m20);
	adj[2][1] = (float)-(m00 * m21 - m01 * m20);
	adj[2][2] = (float)+(m00 * m11 - m01 * m10);

	float det = dot(m[0], vec3(adj[0][0], adj[1][0], adj[2][0]));

	return adj / det;
}

//-----------------------------------------------------------------------------
// mat4 operators and functions
//-----------------------------------------------------------------------------
mat4 transpose(const mat4& m) {
	return mat4(
		vec4(m[0][0], m[1][0], m[2][0], m[3][0]),
		vec4(m[0][1], m[1][1], m[2][1], m[3][1]),
		vec4(m[0][2], m[1][2], m[2][2], m[3][2]),
		vec4(m[0][3], m[1][3], m[2][3], m[3][3]));
}

float determinant(const mat4& m) {
	float a = m[0][0], b = m[1][0], c = m[2][0], d = m[3][0];
	mat3 ma = mat3(m[1].yzw, m[2].yzw, m[3].yzw);
	mat3 mb = mat3(m[0].yzw, m[2].yzw, m[3].yzw);
	mat3 mc = mat3(m[0].yzw, m[1].yzw, m[3].yzw);
	mat3 md = mat3(m[0].yzw, m[1].yzw, m[2].yzw);
	return a * determinant(ma) - b * determinant(mb) + c * determinant(mc) - d * determinant(md);
}

mat4 inverse(const mat4& m) {
	double m00 = m[0][0], m01 = m[0][1], m02 = m[0][2], m03 = m[0][3];
	double m10 = m[1][0], m11 = m[1][1], m12 = m[1][2], m13 = m[1][3];
	double m20 = m[2][0], m21 = m[2][1], m22 = m[2][2], m23 = m[2][3];
	double m30 = m[3][0], m31 = m[3][1], m32 = m[3][2], m33 = m[3][3];

	double s00 = m22 * m33 - m32 * m23;
	double s01 = m21 * m33 - m31 * m23;
	double s02 = m21 * m32 - m31 * m22;
	double s03 = m20 * m33 - m30 * m23;
	double s04 = m20 * m32 - m30 * m22;
	double s05 = m20 * m31 - m30 * m21;
	double s06 = m12 * m33 - m32 * m13;
	double s07 = m11 * m33 - m31 * m13;
	double s08 = m11 * m32 - m31 * m12;
	double s09 = m10 * m33 - m30 * m13;
	double s10 = m10 * m32 - m30 * m12;
	double s11 = m11 * m33 - m31 * m13;
	double s12 = m10 * m31 - m30 * m11;
	double s13 = m12 * m23 - m22 * m13;
	double s14 = m11 * m23 - m21 * m13;
	double s15 = m11 * m22 - m21 * m12;
	double s16 = m10 * m23 - m20 * m13;
	double s17 = m10 * m22 - m20 * m12;
	double s18 = m10 * m21 - m20 * m11;

	mat4 adj;
	adj[0][0] = (float)+(m11 * s00 - m12 * s01 + m13 * s02);
	adj[0][1] = (float)-(m01 * s00 - m02 * s01 + m03 * s02);
	adj[0][2] = (float)+(m01 * s06 - m02 * s07 + m03 * s08);
	adj[0][3] = (float)-(m01 * s13 - m02 * s14 + m03 * s15);
	adj[1][0] = (float)-(m10 * s00 - m12 * s03 + m13 * s04);
	adj[1][1] = (float)+(m00 * s00 - m02 * s03 + m03 * s04);
	adj[1][2] = (float)-(m00 * s06 - m02 * s09 + m03 * s10);
	adj[1][3] = (float)+(m00 * s13 - m02 * s16 + m03 * s17);
	adj[2][0] = (float)+(m10 * s01 - m11 * s03 + m13 * s05);
	adj[2][1] = (float)-(m00 * s01 - m01 * s03 + m03 * s05);
	adj[2][2] = (float)+(m00 * s11 - m01 * s09 + m03 * s12);
	adj[2][3] = (float)-(m00 * s14 - m01 * s16 + m03 * s18);
	adj[3][0] = (float)-(m10 * s02 - m11 * s04 + m12 * s05);
	adj[3][1] = (float)+(m00 * s02 - m01 * s04 + m02 * s05);
	adj[3][2] = (float)-(m00 * s08 - m01 * s10 + m02 * s12);
	adj[3][3] = (float)+(m00 * s15 - m01 * s17 + m02 * s18);

	float det = dot(m[0], vec4(adj[0][0], adj[1][0], adj[2][0], adj[3][0]));

	return adj / det;
}

//------------------------------------------------------------------------------
// Additional non-standard helpers
//------------------------------------------------------------------------------
mat2 rotation(float angleRadians) {
	float c = cos(angleRadians);
	float s = sin(angleRadians);
	return mat2(vec2(c, s), vec2(-s, c));
}

mat3 axisAngleRotation(const vec3& axis, float angleRadians) {
	vec3 u = safeNormalize(axis);
	float x = u.x, y = u.y, z = u.z;
	float c = cos(angleRadians);
	float s = sin(angleRadians);

	return mat3(
		vec3(c + x * x*(1.0f - c), y*x*(1.0f - c) + z * s, z*x*(1.0f - c) - y * s),
		vec3(x*y*(1.0f - c) - z * s, c + y * y*(1.0f - c), z*y*(1.0f - c) + x * s),
		vec3(x*z*(1.0f - c) + y * s, y*z*(1.0f - c) - x * s, c + z * z*(1.0f - c)));
}

mat3 axisCrossAxisRotation(const vec3& trueAxis, const vec3& hintAxis, uint32_t i, uint32_t j, uint32_t k) {
	bool naturalOrder = (i + 1 == j) || (i == 2 && j == 0);
	mat3 m;
	m[i] = safeNormalize(trueAxis);
	if (naturalOrder) {
		m[k] = safeNormalize(cross(m[i], hintAxis));
		m[j] = safeNormalize(cross(m[k], m[i]));
	}
	else {
		m[k] = safeNormalize(cross(hintAxis, m[i]));
		m[j] = safeNormalize(cross(m[i], m[k]));
	}
	return m;
}

mat3 transform(const mat2& upper2x2, const vec2& translation) {
	return mat3(
		vec3(upper2x2[0], 0.0f),
		vec3(upper2x2[1], 0.0f),
		vec3(translation, 1.0f));
}

mat4 transform(const mat3& upper3x3, const vec3& translation) {
	return mat4(
		vec4(upper3x3[0], 0.0f),
		vec4(upper3x3[1], 0.0f),
		vec4(upper3x3[2], 0.0f),
		vec4(translation, 1.0f));
}

mat4 projection(float near, float verticalFovDegrees, float width, float height) {
	float t = tan(0.5f * verticalFovDegrees * ((float)M_PI / 180.0f));
	float h = width / height * near * t;
	float v = near * t;
	return mat4( // Remember, mat4 is column major, so this is visually transposed ...
		vec4(near / h, 0.0f, 0.0f, 0.0f),
		vec4(0.0f, near / v, 0.0f, 0.0f),
		vec4(0.0f, 0.0f, -1.0f, -1.0f),
		vec4(0.0f, 0.0f, -2.0f * near, 0.0f));
}

mat4 ortho(float left, float right, float bottom, float top, float near, float far) {
	mat4 m = mat4(vec4(0.0f), vec4(0.0f), vec4(0.0f), vec4(0.0f));
	m[0][0] = 2.0f / (right - left);
	m[1][1] = 2.0f / (top - bottom);
	m[2][2] = 2.0f / (far - near);
	m[3][0] = -(right + left) / (right - left);
	m[3][1] = -(top + bottom) / (top - bottom);
	m[3][2] = -(far + near) / (far - near);
	m[3][3] = 1.0f;
	return m;
}
