//-----------------------------------------------------------------------------
// Function templates for running the XPBD integration algorithm, both in its
// regular purely Gauss-Seidel form, and partially simultaneous.
//
// (And as a bonus, contains the damping algorithm described in the original
//  PDB paper.)
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include "vectormath.h"
#include "Ldlt.h"

//-----------------------------------------------------------------------------
template<typename Dvec, typename Vec, uint32_t Nodes>
inline void XpbdConstrain(Dvec* X, Dvec* O, float* w, const uint32_t(&is)[Nodes], float C, Vec(&g)[Nodes], float dt, float compliance, float damping, float overRelaxation) {
	float alpha = compliance / (dt * dt);
	float beta = damping * (dt * dt);
	float gamma = alpha * beta / dt;
	float ldiv = 0.0f;
	float damp = 0.0f;
	for (uint32_t i = 0; i < Nodes; i++) {
		ldiv += w[is[i]] * dot(g[i], g[i]);
		damp += dot(Vec(X[is[i]] - O[is[i]]), g[i]);
	}
	ldiv = ldiv * (1.0f + gamma) + alpha;
	float lambda = (-C - gamma * damp) / ldiv;
	for (uint32_t i = 0; i < Nodes; i++) {
		X[is[i]] += Dvec(overRelaxation * w[is[i]] * lambda * g[i]);
	}
}

template<typename Dvec, typename Vec, uint32_t Constraints, uint32_t Nodes>
inline void XpbdConstrainSimultaneous(
	Dvec* X, Dvec* O, float* w, const uint32_t(&is)[Nodes],
	float(&C)[Constraints], Vec(&g)[Constraints][Nodes],
	float dt, float(&compliance)[Constraints], float(&damping)[Constraints], float overRelaxation
) {
	const uint32_t N = Constraints;

	float alpha[N], gamma[N];
	for (uint32_t i = 0; i < N; i++) {
		alpha[i] = compliance[i] / (dt * dt);
		float beta = damping[i] * (dt * dt);
		gamma[i] = alpha[i] * beta / dt;
	}

	auto Dot = [&](uint32_t i, uint32_t j) {
		float r = 0.0f;
		for (uint32_t k = 0; k < Nodes; k++) {
			r += w[is[k]] * dot(g[i][k], g[j][k]);
		}
		return r;
	};
	float A[((N + 1) * N) / 2];
	float b[N];
	uint32_t k = 0;
	for (uint32_t i = 0; i < N; i++) {
		for (uint32_t j = 0; j <= i; j++) {
			A[k++] = Dot(j, i);
		}
		if (gamma[i] > 0.0f) {
			A[k - 1] += alpha[i] / (1.0f + gamma[i]);
			b[i] = -C[i];
			for (uint32_t k = 0; k < Nodes; k++) {
				b[i] -= gamma[i] * dot(g[i][k], Vec(X[is[k]] - O[is[k]]));
			}
			b[i] /= (1.0f + gamma[i]);
		} else {
			A[k - 1] += alpha[i];
			b[i] = -C[i];
		}
	}
	float L[N * N];
	uint32_t r[N];
	float invD[N];
	LDLTDecomposition(N, A, L, r, invD);
	float lambda[N];
	LDLTSolve(N, L, r, invD, b, lambda);

	for (uint32_t i = 0; i < Nodes; i++) {
		Vec lambdagAccum = Vec(0.0f);
		for (uint32_t j = 0; j < N; j++) { lambdagAccum += lambda[j] * g[j][i]; }
		X[is[i]] += Dvec(overRelaxation * w[is[i]] * lambdagAccum);
	}
}

// Automatically fall back to simple XpbdConstrain when we just have one constraint
template<typename Dvec, typename Vec, uint32_t Nodes>
inline void XpbdConstrainSimultaneous(
	Dvec* X, Dvec* O, float* w, const uint32_t(&is)[Nodes],
	float(&C)[1], Vec(&g)[1][Nodes],
	float dt, float(&compliance)[1], float(&damping)[1], float overRelaxation
) {
	XpbdConstrain(X, O, w, is, C[0], g[0], dt, compliance[0], damping[0], overRelaxation);
}

template<uint32_t N>
inline void PbdDamp(const dvec2* Xd, dvec2* Vd, float* W, const uint32_t(&is)[N], float damping) {
	vec2 X[N], V[N];
	float M[N];
	vec2 Xcm = vec2(0.0f), Vcm = vec2(0.0f);
	float Msum = 0.0f;
	for (uint32_t i = 0; i < N; i++) {
		X[i] = vec2(Xd[is[i]] - Xd[is[N -1]]);
		V[i] = vec2(Vd[is[i]]);
		M[i] = 1.0f / max(1.0e-12f, W[is[i]]);
		Xcm += X[i] * M[i];
		Vcm += V[i] * M[i];
		Msum += M[i];
	}
	Xcm /= Msum;
	Vcm /= Msum;

	vec2 r[N];
	float L = 0.0f;
	float I = 0.0f;
	for (uint32_t i = 0; i < N; i++) {
		r[i] = X[i] - Xcm;
		auto cross = [](vec2 a, vec2 b) { return a.x * b.y - a.y * b.x; };
		L += cross(r[i], M[i] * V[i]);
		I += dot(r[i], r[i]) * M[i];
	}

	float w = L / I;
	for (uint32_t i = 0; i < N; i++) {
		vec2 dV = Vcm + w * vec2(-r[i].y, r[i].x) - V[i];
		Vd[is[i]] += dvec2(damping * dV);
	}
}

template<uint32_t N>
inline void PbdDamp(const dvec3* Xd, dvec3* Vd, float* W, const uint32_t(&is)[N], float damping) {
	vec3 X[N], V[N];
	float M[N];
	vec3 Xcm = vec3(0.0f), Vcm = vec3(0.0f);
	float Msum = 0.0f;
	for (uint32_t i = 0; i < N; i++) {
		X[i] = vec3(Xd[is[i]] - Xd[is[N - 1]]);
		V[i] = vec3(Vd[is[i]]);
		M[i] = 1.0f / max(1.0e-12f, W[is[i]]);
		Xcm += X[i] * M[i];
		Vcm += V[i] * M[i];
		Msum += M[i];
	}
	Xcm /= Msum;
	Vcm /= Msum;

	vec3 r[N], L = vec3(0.0f);
	mat3 I = mat3(0.0f);
	for (uint32_t i = 0; i < N; i++) {
		r[i] = X[i] - Xcm;
		L += cross(r[i], M[i] * V[i]);
		vec3 rr = r[i] * r[i];
		vec3 rp = r[i] * r[i].yzx;
		I += M[i] * mat3(
			vec3(rr.z + rr.y, -rp[0], -rp[2]),
			vec3(-rp[0], rr.z + rr.x, -rp[1]),
			vec3(-rp[2], -rp[1], rr.y + rr.x));
	}

	vec3 w = inverse(I) * L;
	for (uint32_t i = 0; i < N; i++) {
		vec3 dV = Vcm + cross(w, r[i]) - V[i];
		Vd[is[i]] += dvec3(damping * dV);
	}
}
