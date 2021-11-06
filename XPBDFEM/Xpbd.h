//-----------------------------------------------------------------------------
// Function templates for running the XPBD integration algorithm, both in its
// regular purely Gauss-Seidel form, and a partially simultaneous form. Also
// includes a modifed EnergyXpbdConstrain version which operates on potential
// energies (with their compliances factored out) and their gradients.
//
// As a bonus, contains two damping algorithms. The Rayleigh damping from the
// XPBD paper isolated into its own function, as well as the non-rigid-body
// damping described in the original PDB paper.
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include "vectormath.h"
#include "Ldlt.h"
//@HACK: For Trying different damping types
#include "Settings.h"

//-----------------------------------------------------------------------------
template<typename Dvec, typename Vec>
inline void XpbdConstrain(uint32_t nodeCount, Dvec* X, const Dvec* O, const float* w, const uint32_t* is, float C, const Vec* g, float dt, float compliance, float damping) {
	float alpha = compliance / (dt * dt);
	float beta = damping * (dt * dt);
	float gamma = alpha * beta / dt;
	float wgg = 1.0e-22f; // Prevent divide by zero
	float gV = 0.0f;
	for (uint32_t i = 0; i < nodeCount; i++) {
		wgg += w[is[i]] * dot(g[i], g[i]);
		gV += dot(g[i], Vec(X[is[i]] - O[is[i]]));
	}
	float lambda = (-C - gamma * gV) / (wgg * (1.0f + gamma) + alpha);
	for (uint32_t i = 0; i < nodeCount; i++) {
		X[is[i]] += Dvec(w[is[i]] * lambda * g[i]);
	}
}

template<typename Dvec, typename Vec, uint32_t Constraints>
inline void XpbdConstrainSimultaneous(
	uint32_t nodeCount,
	Dvec* X, const Dvec* O, const float* w, const uint32_t* is, // expected: is[nodeCount]
	const float(&C)[Constraints], const Vec* g, // expected: g[Constraints * nodeCount]
	float dt, const float(&compliance)[Constraints], const float(&damping)[Constraints]
) {
	const uint32_t N = Constraints;

	float alpha[N], gamma[N];
	for (uint32_t i = 0; i < N; i++) {
		alpha[i] = compliance[i] / (dt * dt);
		float beta = damping[i] * (dt * dt);
		gamma[i] = alpha[i] * beta / dt;
	}

	float A[((N + 1) * N) / 2];
	float b[N];
	uint32_t k = 0;
	for (uint32_t i = 0; i < N; i++) {
		for (uint32_t j = 0; j <= i; j++) {
			float wgg = 1.0e-22f; // Prevent divide by zero
			for (uint32_t n = 0; n < nodeCount; n++) {
				wgg += w[is[n]] * dot(g[i * nodeCount + n], g[j * nodeCount + n]);
			}
			A[k++] = wgg;
		}
		if (gamma[i] > 0.0f) {
			A[k - 1] += alpha[i] / (1.0f + gamma[i]);
			b[i] = -C[i];
			for (uint32_t k = 0; k < nodeCount; k++) {
				b[i] -= gamma[i] * dot(g[i * nodeCount + k], Vec(X[is[k]] - O[is[k]]));
			}
			b[i] /= (1.0f + gamma[i]);
		} else {
			A[k - 1] += alpha[i];
			b[i] = -C[i];
		}
	}
	float lambda[N];
	LDLTDecomposeAndSolve(N, A, b, lambda);

	for (uint32_t i = 0; i < nodeCount; i++) {
		Vec lambdagAccum = Vec(0.0f);
		for (uint32_t j = 0; j < N; j++) { lambdagAccum += lambda[j] * g[j * nodeCount + i]; }
		X[is[i]] += Dvec(w[is[i]] * lambdagAccum);
	}
}

template<typename Dvec, typename Vec>
inline void EnergyXpbdConstrain(uint32_t nodeCount, Dvec* X, const Dvec* O, const float* w, const uint32_t* is, float U, const Vec* g, float dt, float compliance, float dampingGamma, uint32_t dampingType) {
	//IMPORTANT: U must be >= 0!
	float alpha = compliance / (dt * dt);
	float gamma = dampingGamma / dt;
	float lambdaPrime;
	if (dampingGamma > 0.0f && dampingType < Rayleigh_Post) {
		float wgg = 1.0e-22f; // Prevent divide by zero, except in the very unlikely case where compliance and gradients are 0 and U is large
		float gV = 0.0f;
		for (uint32_t i = 0; i < nodeCount; i++) {
			wgg += w[is[i]] * dot(g[i], g[i]);
			gV += dot(Vec(X[is[i]] - O[is[i]]), g[i]);
		}
		if (dampingType == Rayleigh_Paper) {
			float A = wgg * (1.0f + gamma) + (2.0f * U * alpha);
			lambdaPrime = (-2.0f * U - gamma * gV) / A;
		} else {
			lambdaPrime = (-2.0f * U) / (wgg + (2.0f * U * alpha));
			float invBeta = compliance / (dt * dampingGamma);
			float A = wgg + 2.0f * U * invBeta;
			float b = -gV - wgg * lambdaPrime; // (wgg * lambdaPrime) is the contribution to gV added by the main constraint
			b *= A / max(A, 4.0f * wgg);
			lambdaPrime += b / A;
		}
	} else {
		float wgg = 1.0e-22f;
		for (uint32_t i = 0; i < nodeCount; i++) {
			wgg += w[is[i]] * dot(g[i], g[i]);
		}
		lambdaPrime = (-2.0f * U) / (wgg + (2.0f * U * alpha));
	}
	for (uint32_t i = 0; i < nodeCount; i++) {
		X[is[i]] += Dvec(w[is[i]] * lambdaPrime * g[i]);
	}
}

template<typename Dvec, typename Vec, uint32_t Constraints>
inline void EnergyXpbdConstrainSimultaneous(
	uint32_t nodeCount,
	Dvec* X, const Dvec* O, const float* w, const uint32_t* is, // expected: is[nodeCount]
	const float(&U)[Constraints], const Vec* g, // expected: g[Constraints * nodeCount]
	float dt, const float(&compliance)[Constraints], const float(&dampingGamma)[Constraints],
	uint32_t dampingType
) {
	//IMPORTANT: All U must be >= 0!
	const uint32_t N = Constraints;

	float alpha[N], gamma[N];
	bool anyDamping = false;
	for (uint32_t i = 0; i < N; i++) {
		alpha[i] = compliance[i] / (dt * dt);
		gamma[i] = dampingGamma[i] / dt;
		anyDamping = anyDamping || dampingGamma[i] > 0.0f;
	}

	float gV[N] = { 0.0f };
	if (anyDamping && dampingType < Rayleigh_Post) {
		for (uint32_t n = 0; n < nodeCount; n++) {
			Vec v = Vec(X[is[n]] - O[is[n]]);
			for (uint32_t i = 0; i < N; i++) {
				gV[i] += dot(g[i * nodeCount + n], v);
			}
		}
	}
	float wgg[((N + 1) * N) / 2];
	float A[((N + 1) * N) / 2];
	float b[N];
	uint32_t k = 0;
	for (uint32_t i = 0; i < N; i++) {
		for (uint32_t j = 0; j <= i; j++) {
			wgg[k] = 1.0e-22f; // Prevent divide by zero, except in the very unlikely case where compliance and gradients are 0 and U is large
			for (uint32_t n = 0; n < nodeCount; n++) {
				wgg[k] += w[is[n]] * dot(g[i * nodeCount + n], g[j * nodeCount + n]);
			}
			A[k] = wgg[k];
			++k;
		}
		if (dampingType == Rayleigh_Paper && gamma[i] > 0.0f) {
			A[k - 1] += (2.0f * U[i] * alpha[i]) / (1.0f + gamma[i]);
			b[i] = (-2.0f * U[i] - gamma[i] * gV[i]) / (1.0f + gamma[i]);
		} else {
			A[k - 1] += 2.0f * U[i] * alpha[i];
			b[i] = -2.0f * U[i];
		}
	}

	float lambdaPrime[N];
	if constexpr (N == 2) {
		// Cramer's rule (massaged to prevent floating point over/underflow)
		float invA00 = 1.0f / A[0];
		float invA11 = 1.0f / A[2];
		float invDet = 1.0f / max(0.00000001f, 1.0f - (A[1] * invA00) * (A[1] * invA11));
		lambdaPrime[0] = invDet * ((b[0] * invA00) - (A[1] * invA00) * (b[1] * invA11));
		lambdaPrime[1] = invDet * ((b[1] * invA11) - (b[0] * invA00) * (A[1] * invA11));
	} else {
		LDLTDecomposeAndSolve(N, A, b, lambdaPrime);
	}

	if (dampingType == Rayleigh_Limit && gamma[0] > 0.0f) {
		uint32_t k = 0;
		for (uint32_t i = 0; i < N; i++) {
			k += i + 1;
			float invBeta = compliance[i] / (dt * dampingGamma[i]);
			A[k - 1] = wgg[k - 1] + 2.0f * U[i] * invBeta;
			b[i] *= (-1.0f - alpha[i] * lambdaPrime[i]); // This is equivalent to setting b = -(wgg * lambdaPrime)
			b[i] += -gV[i];                              // (where (wgg * lambdaPrime) is the contribution to gV added by the main constraint)
			// Apply the limit
			b[i] *= A[k - 1] / max(A[k - 1], 8.0f * wgg[k - 1]);
		}
		if constexpr (N == 2) {
			// Cramer's rule (massaged to prevent floating point over/underflow)
			float invA00 = 1.0f / A[0];
			float invA11 = 1.0f / A[2];
			float invDet = 1.0f / max(0.00000001f, 1.0f - (A[1] * invA00) * (A[1] * invA11));
			lambdaPrime[0] += invDet * ((b[0] * invA00) - (A[1] * invA00) * (b[1] * invA11));
			lambdaPrime[1] += invDet * ((b[1] * invA11) - (b[0] * invA00) * (A[1] * invA11));
		} else {
			float lambdaDamp[N];
			LDLTDecomposeAndSolve(N, A, b, lambdaDamp);
			for (uint32_t n = 0; n < N; n++) { lambdaPrime[n] += lambdaDamp[n]; }
		}
	}

	for (uint32_t i = 0; i < nodeCount; i++) {
		Vec lambdagAccum = Vec(0.0f);
		for (uint32_t j = 0; j < N; j++) { lambdagAccum += lambdaPrime[j] * g[j * nodeCount + i]; }
		X[is[i]] += Dvec(w[is[i]] * lambdagAccum);
	}
}

template<typename Dvec, typename Vec>
inline void RayleighDamp(uint32_t nodeCount, Dvec* V, const float* w, const uint32_t* is,  float U, const Vec* g, float dt, float compliance, float dampingGamma) {
	float invBeta = compliance / (dt * dampingGamma);
	float wgg = 1.0e-22f;
	float gV = 0.0f;
	for (uint32_t n = 0; n < nodeCount; n++) {
		wgg += w[is[n]] * dot(g[n], g[n]);
		gV += dot(g[n], Vec(V[is[n]]));
	}
	float lambda = -gV / (2.0f * U * invBeta + wgg);
	for (uint32_t n = 0; n < nodeCount; n++) {
		V[is[n]] += Dvec(w[is[n]] * lambda * g[n]);
	}
}

template<typename Dvec, typename Vec>
inline void RayleighDamp(
	uint32_t nodeCount,
	Dvec* V, const float* w, const uint32_t* is, // expected: is[nodeCount]
	const float(&U)[2], const Vec* g, // expected: g[2 * nodeCount]
	float dt, const float(&compliance)[2], const float(&dampingGamma)[2]
) {
	float invBeta[2] = { compliance[0] / (dt * dampingGamma[0]), compliance[1] / (dt * dampingGamma[1]) };
	float A[3] = { 1.0e-22f, 1.0e-22f, 1.0e-22f };
	float b[2] = { 0.0f };
	for (uint32_t n = 0; n < nodeCount; n++) {
		A[0] += w[is[n]] * dot(g[0 * nodeCount + n], g[0 * nodeCount + n]);
		A[1] += w[is[n]] * dot(g[0 * nodeCount + n], g[1 * nodeCount + n]);
		A[2] += w[is[n]] * dot(g[1 * nodeCount + n], g[1 * nodeCount + n]);
		Vec v = Vec(V[is[n]]);
		b[0] -= dot(g[0 * nodeCount + n], v);
		b[1] -= dot(g[1 * nodeCount + n], v);
	}
	A[0] += 2.0f * U[0] * invBeta[0];
	A[2] += 2.0f * U[1] * invBeta[1];

	float lambdaPrime[2];
	float invA00 = 1.0f / A[0];
	float invA11 = 1.0f / A[2];
	float invDet = 1.0f / max(0.00000001f, 1.0f - (A[1] * invA00) * (A[1] * invA11));
	lambdaPrime[0] = invDet * ((b[0] * invA00) - (A[1] * invA00) * (b[1] * invA11));
	lambdaPrime[1] = invDet * ((b[1] * invA11) - (b[0] * invA00) * (A[1] * invA11));

	for (uint32_t n = 0; n < nodeCount; n++) {
		Vec lambdagAccum = lambdaPrime[0] * g[0 * nodeCount + n] + lambdaPrime[1] * g[1 * nodeCount + n];
		V[is[n]] += Dvec(w[is[n]] * lambdagAccum);
	}
}

template<uint32_t N>
inline void PbdDamp(const dvec2* Xd, dvec2* Vd, const float* W, const uint32_t(&is)[N], float damping) {
	vec2 X[N], V[N];
	float M[N];
	vec2 Xcm = vec2(0.0f), Vcm = vec2(0.0f);
	float Msum = 0.0f;

	// Numerically, we get better results if mass is around 1, and we limit the heaviness of the heaviest particles
	float Wsum = 1.0e-8f;
	for (uint32_t i = 0; i < N; i++) { Wsum += W[is[i]]; }
	float Waverage = Wsum * (1.0f / (float)N);
	float Wmin = 0.0001f * Waverage;

	for (uint32_t i = 0; i < N; i++) {
		X[i] = vec2(Xd[is[i]] - Xd[is[N -1]]);
		V[i] = vec2(Vd[is[i]]);
		M[i] = Waverage / max(Wmin, W[is[i]]);
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
		L += M[i] * cross(r[i], V[i]);
		I += M[i] * dot(r[i], r[i]);
	}

	float w = L / I;
	for (uint32_t i = 0; i < N; i++) {
		vec2 dV = Vcm + w * vec2(-r[i].y, r[i].x) - V[i];
		Vd[is[i]] += dvec2(damping * dV);
	}
}

template<uint32_t N>
inline void PbdDamp(const dvec3* Xd, dvec3* Vd, const float* W, const uint32_t(&is)[N], float damping) {
	vec3 X[N], V[N];
	float M[N];
	vec3 Xcm = vec3(0.0f), Vcm = vec3(0.0f);
	float Msum = 0.0f;

	// Numerically, we get better results if mass is around 1, and we limit the heaviness of the heaviest particles
	float Wsum = 1.0e-8f;
	for (uint32_t i = 0; i < N; i++) { Wsum += W[is[i]]; }
	float Waverage = Wsum * (1.0f / (float)N);
	float Wmin = 0.0001f * Waverage;

	for (uint32_t i = 0; i < N; i++) {
		X[i] = vec3(Xd[is[i]] - Xd[is[N - 1]]);
		V[i] = vec3(Vd[is[i]]);
		M[i] = Waverage / max(Wmin, W[is[i]]);
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
		L += M[i] * cross(r[i], V[i]);
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
