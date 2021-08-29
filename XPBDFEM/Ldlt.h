//-----------------------------------------------------------------------------
// Helpers for LDL^T decomposition and solving of symmetric positive definite
// linear systems.
//
// Usage example:
//
// const uint32_t N = 3;
// float A[] = {
//	1.0f,
//	2.0f, 3.0f,
//	4.0f, 5.0f, 6.0f,
// };
// float b[] = { 7.0f, 8.0f, 9.0f, };
// 
// float L[N * N];
// uint32_t r[N];
// float invD[N];
// LDLTDecomposition(N, A, L, r, invD);
// float X[N];
// LDLTSolve(N, L, r, invD, b, X);
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include <inttypes.h>

namespace LDLTPrivate {
	const uint32_t r[] = { 0, 0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91, 105, 120, 136, 153, 171, 190, 210, 231, 253, 276, 300, };
}

inline void LDLTDecomposition(uint32_t n, float* ASym, float* outL, uint32_t* outr, float* outinvD) {
	using namespace LDLTPrivate;

	// Calculate row offsets
	//uint32_t* r = outr;
	//r[0] = 0;
	//for (uint32_t i = 1; i < n; i++) {
	//	r[i] = r[i - 1] + i - 1;
	//}

	float* L = outL;
	float* D = outinvD;
	for (uint32_t j = 0; j < n; j++) { // each columns
		D[j] = ASym[(r[j] + j) + j];
		for (uint32_t k = 0; k < j; k++) { // sum
			D[j] -= L[r[j] + k] * L[r[j] + k] * D[k];
		}

		for (uint32_t i = j + 1; i < n; i++) { // remaining rows
			L[r[i] + j] = ASym[(r[i] + i) + j];
			for (uint32_t k = 0; k < j; k++) { // sum
				L[r[i] + j] -= L[r[i] + k] * L[r[j] + k] * D[k];
			}
			L[r[i] + j] /= D[j];
		}
	}

	for (uint32_t j = 0; j < n; j++) {
		outinvD[j] = 1.0f / D[j];
	}
}

inline void LDLTSolve(uint32_t n, float* L, uint32_t* r_dummy, float* invD, float* b, float* outX) {
	using namespace LDLTPrivate;

	float* x = outX;

	for (uint32_t i = 0; i < n; i++) {
		x[i] = b[i];
		for (uint32_t j = 0; j < i; j++) {
			x[i] -= L[r[i] + j] * x[j];
		}
	}
	for (uint32_t i = 0; i < n; i++) {
		x[i] *= invD[i];
	}
	for (uint32_t i = 0; i < n; i++) {
		for (uint32_t j = 0; j < i; j++) {
			x[n - 1 - i] -= L[r[n - 1 - j] + n - 1 - i] * x[n - 1 - j];
		}
	}
}
