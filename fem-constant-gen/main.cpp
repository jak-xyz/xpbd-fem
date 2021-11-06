// clear && clang++-12 main.cpp -std=c++17 -O2 && ./a.out

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include "../XPBDFEM/vectormath.h"
#include "../XPBDFEM/vectormath.cpp"

//-----------------------------------------------------------------------------
// Quad/Hex element node layouts
template<typename Vec> struct NodeLayout;
template<> struct NodeLayout<vec2> { uint32_t x, y; };
template<> struct NodeLayout<vec3> { uint32_t x, y, z; };

template<typename Vec, uint32_t Nodes> constexpr uint32_t NodeLayouts = 0;

template<> constexpr NodeLayout<vec2> NodeLayouts<vec2, 4>[4] = {
	{ 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 },
};
template<> constexpr NodeLayout<vec2> NodeLayouts<vec2, 9>[9] = {
	{ 0, 0 }, { 2, 0 }, { 2, 2 }, { 0, 2 }, { 1, 0 }, { 2, 1 }, { 1, 2 }, { 0, 1 }, { 1, 1 },
};
template<> constexpr NodeLayout<vec3> NodeLayouts<vec3, 8>[8] = {
	{ 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 },
	{ 0, 0, 1 }, { 1, 0, 1 }, { 0, 1, 1 }, { 1, 1, 1 },
};
template<> constexpr NodeLayout<vec3> NodeLayouts<vec3, 27>[27] = {
	{ 0, 0, 0 }, { 1, 0, 0 }, { 2, 0, 0 }, { 0, 1, 0 }, { 1, 1, 0 }, { 2, 1, 0 }, { 0, 2, 0 }, { 1, 2, 0 }, { 2, 2, 0 },
	{ 0, 0, 1 }, { 1, 0, 1 }, { 2, 0, 1 }, { 0, 1, 1 }, { 1, 1, 1 }, { 2, 1, 1 }, { 0, 2, 1 }, { 1, 2, 1 }, { 2, 2, 1 },
	{ 0, 0, 2 }, { 1, 0, 2 }, { 2, 0, 2 }, { 0, 1, 2 }, { 1, 1, 2 }, { 2, 1, 2 }, { 0, 2, 2 }, { 1, 2, 2 }, { 2, 2, 2 },
};

//-----------------------------------------------------------------------------
// Shape functions

enum class FemShape { Tri = 0, Tet, Quad, Hex };
const char* ShapeName[] = { "Tri", "Tet", "Quad", "Hex" };

// Tri/Tet:
constexpr float TriLinear_N(uint32_t N, float L1, float L2, float L3) {
	switch (N) {
		case 0: return L1; break;
		case 1: return L2; break;
		case 2: return L3; break;
		default: assert(false); return 0.0f;
	}
}
constexpr vec2 TriLinear_dN(uint32_t N, float L1, float L2, float L3) {
	switch (N) {
		case 0: return vec2(1.0f, 0.0f); break;
		case 1: return vec2(0.0f, 1.0f); break;
		case 2: return vec2(-1.0f, -1.0f); break;
		default: assert(false); return vec2(0.0f);
	}
}
constexpr float TriQuadratic_N(uint32_t N, float L1, float L2, float L3) {
	switch (N) {
		case 0: return (2.0f * L1 - 1.0f) * L1; break;
		case 1: return (2.0f * L2 - 1.0f) * L2; break;
		case 2: return (2.0f * L3 - 1.0f) * L3; break;
		case 3: return 4.0f * L1 * L2; break;
		case 4: return 4.0f * L2 * L3; break;
		case 5: return 4.0f * L3 * L1; break;
		default: assert(false); return 0.0f;
	}
}
constexpr vec2 TriQuadratic_dN(uint32_t N, float L1, float L2, float L3) {
	switch (N) {
		case 0: return vec2((4.0f * L1 - 1.0f), 0.0f); break;
		case 1: return vec2(0.0f, (4.0f * L2 - 1.0f)); break;
		case 2: return vec2(-(4.0f * L3 - 1.0f), -(4.0f * L3 - 1.0f)); break;
		case 3: return vec2((4.0f * L2), (4.0f * L1)); break;
		case 4: return vec2(-(4.0f * L2), (4.0f * L3) - (4.0f * L2)); break;
		case 5: return vec2((4.0f * L3) - (4.0f * L1), -(4.0f * L1)); break;
		default: assert(false); return vec2(0.0f);
	}
}
constexpr float TetLinear_N(uint32_t N, float L1, float L2, float L3, float L4) {
	switch (N) {
		case 0: return L1; break;
		case 1: return L2; break;
		case 2: return L3; break;
		case 3: return L4; break;
		default: assert(false); return 0.0f;
	}
}
constexpr vec3 TetLinear_dN(uint32_t N, float L1, float L2, float L3, float L4) {
	switch (N) {
		case 0: return vec3(1.0f, 0.0f, 0.0f); break;
		case 1: return vec3(0.0f, 1.0f, 0.0f); break;
		case 2: return vec3(0.0f, 0.0f, 1.0f); break;
		case 3: return vec3(-1.0f, -1.0f, -1.0f); break;
		default: assert(false); return vec3(0.0f);
	}
}
constexpr float TetQuadratic_N(uint32_t N, float L1, float L2, float L3, float L4) {
	switch (N) {
		case 0: return (2.0f * L1 - 1.0f) * L1; break;
		case 1: return (2.0f * L2 - 1.0f) * L2; break;
		case 2: return (2.0f * L3 - 1.0f) * L3; break;
		case 3: return (2.0f * L4 - 1.0f) * L4; break;
		case 4: return 4.0f * L1 * L2; break;
		case 5: return 4.0f * L1 * L3; break;
		case 6: return 4.0f * L1 * L4; break;
		case 7: return 4.0f * L2 * L3; break;
		case 8: return 4.0f * L2 * L4; break;
		case 9: return 4.0f * L3 * L4; break;
		default: assert(false); return 0.0f;
	}
}
constexpr vec3 TetQuadratic_dN(uint32_t N, float L1, float L2, float L3, float L4) {
	switch (N) {
		case 0: return vec3((4.0f * L1 - 1.0f), 0.0f, 0.0f); break;
		case 1: return vec3(0.0f, (4.0f * L2 - 1.0f), 0.0f); break;
		case 2: return vec3(0.0f, 0.0f, (4.0f * L3 - 1.0f)); break;
		case 3: return vec3(-(4.0f * L4 - 1.0f), -(4.0f * L4 - 1.0f), -(4.0f * L4 - 1.0f)); break;
		case 4: return vec3((4.0f * L2), (4.0f * L1), 0.0f); break;
		case 5: return vec3((4.0f * L3), 0.0f, (4.0f * L1)); break;
		case 6: return vec3((4.0f * L4) - (4.0f * L1), -(4.0f * L1), -(4.0f * L1)); break;
		case 7: return vec3(0.0f, (4.0f * L3), (4.0f * L2)); break;
		case 8: return vec3(-(4.0f * L2), (4.0f * L4) - (4.0f * L2), -(4.0f * L2)); break;
		case 9: return vec3(-(4.0f * L3), -(4.0f * L3), (4.0f * L4) - (4.0f * L3)); break;
		default: assert(false); return vec3(0.0f);
	}
}

// Quad/Hex:
constexpr float LagrangeLinear_N0(float t) { return (1.0f - t) / 2.0f; }
constexpr float LagrangeLinear_N1(float t) { return (1.0f + t) / 2.0f; }
constexpr float LagrangeLinear_dN0(float t) { return -1.0f / 2.0f; }
constexpr float LagrangeLinear_dN1(float t) { return 1.0f / 2.0f; }
constexpr float LagrangeLinear_N(uint32_t N, float t) { return N == 0 ? LagrangeLinear_N0(t) : LagrangeLinear_N1(t); }
constexpr float LagrangeLinear_dN(uint32_t N, float t) { return N == 0 ? LagrangeLinear_dN0(t) : LagrangeLinear_dN1(t); }

constexpr float LagrangeQuadratic_N0(float t) { return (t * (t - 1.0f)) / 2.0f; }
constexpr float LagrangeQuadratic_N1(float t) { return -(t + 1.0f) * (t - 1.0f); }
constexpr float LagrangeQuadratic_N2(float t) { return ((t + 1.0f) * t) / 2.0f; }
constexpr float LagrangeQuadratic_dN0(float t) { return (2.0f * t - 1.0f) / 2.0f; }
constexpr float LagrangeQuadratic_dN1(float t) { return -2.0f * t; }
constexpr float LagrangeQuadratic_dN2(float t) { return (2.0f * t + 1.0f) / 2.0f; }
constexpr float LagrangeQuadratic_N(uint32_t N, float t) { return N == 0 ? LagrangeQuadratic_N0(t) : N == 1 ? LagrangeQuadratic_N1(t) : LagrangeQuadratic_N2(t); }
constexpr float LagrangeQuadratic_dN(uint32_t N, float t) { return N == 0 ? LagrangeQuadratic_dN0(t) : N == 1 ? LagrangeQuadratic_dN1(t) : LagrangeQuadratic_dN2(t); }

template<FemShape Shape, uint32_t Nodes, typename Uv> float ShapeFunction(uint32_t N, Uv uv);
template<> float ShapeFunction<FemShape::Tri, 3, vec3>(uint32_t N, vec3 L) {
	return TriLinear_N(N, L[0], L[1], L[2]);
}
template<> float ShapeFunction<FemShape::Tri, 6, vec3>(uint32_t N, vec3 L) {
	return TriQuadratic_N(N, L[0], L[1], L[2]);
}
template<> float ShapeFunction<FemShape::Tet, 4, vec4>(uint32_t N, vec4 L) {
	return TetLinear_N(N, L[0], L[1], L[2], L[3]);
}
template<> float ShapeFunction<FemShape::Tet, 10, vec4>(uint32_t N, vec4 L) {
	return TetQuadratic_N(N, L[0], L[1], L[2], L[3]);
}
template<> float ShapeFunction<FemShape::Quad, 4, vec2>(uint32_t N, vec2 uv) {
	const auto& NL = NodeLayouts<vec2, 4>[N];
	return LagrangeLinear_N(NL.x, uv.x) * LagrangeLinear_N(NL.y, uv.y);
}
template<> float ShapeFunction<FemShape::Quad, 9, vec2>(uint32_t N, vec2 uv) {
	const auto& NL = NodeLayouts<vec2, 9>[N];
	return LagrangeQuadratic_N(NL.x, uv.x) * LagrangeQuadratic_N(NL.y, uv.y);
}
template<> float ShapeFunction<FemShape::Hex, 8, vec3>(uint32_t N, vec3 uv) {
	const auto& NL = NodeLayouts<vec3, 8>[N];
	return LagrangeLinear_N(NL.x, uv.x) * LagrangeLinear_N(NL.y, uv.y) * LagrangeLinear_N(NL.z, uv.z);
}
template<> float ShapeFunction<FemShape::Hex, 27, vec3>(uint32_t N, vec3 uv) {
	const auto& NL = NodeLayouts<vec3, 27>[N];
	return LagrangeQuadratic_N(NL.x, uv.x) * LagrangeQuadratic_N(NL.y, uv.y) * LagrangeQuadratic_N(NL.z, uv.z);
}

template<FemShape Shape, uint32_t Nodes, typename Uv> auto ShapeFunctionGradient(uint32_t N, Uv uv);
template<> auto ShapeFunctionGradient<FemShape::Tri, 3, vec3>(uint32_t N, vec3 L) {
	return TriLinear_dN(N, L[0], L[1], L[2]);
}
template<> auto ShapeFunctionGradient<FemShape::Tri, 6, vec3>(uint32_t N, vec3 L) {
	return TriQuadratic_dN(N, L[0], L[1], L[2]);
}
template<> auto ShapeFunctionGradient<FemShape::Tet, 4, vec4>(uint32_t N, vec4 L) {
	return TetLinear_dN(N, L[0], L[1], L[2], L[3]);
}
template<> auto ShapeFunctionGradient<FemShape::Tet, 10, vec4>(uint32_t N, vec4 L) {
	return TetQuadratic_dN(N, L[0], L[1], L[2], L[3]);
}
template<> auto ShapeFunctionGradient<FemShape::Quad, 4, vec2>(uint32_t N, vec2 uv) {
	const auto& NL = NodeLayouts<vec2, 4>[N];
	return vec2(
		LagrangeLinear_dN(NL.x, uv.x) * LagrangeLinear_N(NL.y, uv.y),
		LagrangeLinear_N(NL.x, uv.x) * LagrangeLinear_dN(NL.y, uv.y));
}
template<> auto ShapeFunctionGradient<FemShape::Quad, 9, vec2>(uint32_t N, vec2 uv) {
	const auto& NL = NodeLayouts<vec2, 9>[N];
	return vec2(
		LagrangeQuadratic_dN(NL.x, uv.x) * LagrangeQuadratic_N(NL.y, uv.y),
		LagrangeQuadratic_N(NL.x, uv.x) * LagrangeQuadratic_dN(NL.y, uv.y));
}
template<> auto ShapeFunctionGradient<FemShape::Hex, 8, vec3>(uint32_t N, vec3 uv) {
	const auto& NL = NodeLayouts<vec3, 8>[N];
	return vec3(
		LagrangeLinear_dN(NL.x, uv.x) * LagrangeLinear_N(NL.y, uv.y) * LagrangeLinear_N(NL.z, uv.z),
		LagrangeLinear_N(NL.x, uv.x) * LagrangeLinear_dN(NL.y, uv.y) * LagrangeLinear_N(NL.z, uv.z),
		LagrangeLinear_N(NL.x, uv.x) * LagrangeLinear_N(NL.y, uv.y) * LagrangeLinear_dN(NL.z, uv.z));
}
template<> auto ShapeFunctionGradient<FemShape::Hex, 27, vec3>(uint32_t N, vec3 uv) {
	const auto& NL = NodeLayouts<vec3, 27>[N];
	return vec3(
		LagrangeQuadratic_dN(NL.x, uv.x) * LagrangeQuadratic_N(NL.y, uv.y) * LagrangeQuadratic_N(NL.z, uv.z),
		LagrangeQuadratic_N(NL.x, uv.x) * LagrangeQuadratic_dN(NL.y, uv.y) * LagrangeQuadratic_N(NL.z, uv.z),
		LagrangeQuadratic_N(NL.x, uv.x) * LagrangeQuadratic_N(NL.y, uv.y) * LagrangeQuadratic_dN(NL.z, uv.z));
}

//-----------------------------------------------------------------------------
// Gaussian Quadrature helpers
template<FemShape Shape, uint32_t Points> uint32_t QuadratureUvs = 0;
template<> vec3 QuadratureUvs<FemShape::Tri, 1>[1] = { vec3(1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f), };
template<> vec3 QuadratureUvs<FemShape::Tri, 3>[3] = { vec3(0.5f, 0.5f, 0.0f), vec3(0.5f, 0.0f, 0.5f), vec3(0.0f, 0.5f, 0.5f), };
template<> vec3 QuadratureUvs<FemShape::Tri, 7>[7] = { vec3(1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f), vec3(0.0597158717f, 0.4701420641f, 0.4701420641f), vec3(0.4701420641f, 0.0597158717f, 0.4701420641f), vec3(0.4701420641f, 0.4701420641f, 0.0597158717f), vec3(0.7974269853f, 0.1012865073f, 0.1012865073f), vec3(0.1012865073f, 0.7974269853f, 0.1012865073f), vec3(0.1012865073f, 0.1012865073f, 0.7974269853f), };
template<> vec4 QuadratureUvs<FemShape::Tet, 1>[1] = { vec4(1.0f / 4.0f, 1.0f / 4.0f, 1.0f / 4.0f, 1.0f / 4.0f), };
template<> vec4 QuadratureUvs<FemShape::Tet, 4>[4] = { vec4(0.58541020f, 0.13819660f, 0.13819660f, 0.13819660f), vec4(0.13819660f, 0.58541020f, 0.13819660f, 0.13819660f), vec4(0.13819660f, 0.13819660f, 0.58541020f, 0.13819660f), vec4(0.13819660f, 0.13819660f, 0.13819660f, 0.58541020f), };
template<> vec4 QuadratureUvs<FemShape::Tet, 10>[10] = { vec4(0.77849529f, 0.07383490f, 0.07383490f, 0.07383490f), vec4(0.07383490f, 0.77849529f, 0.07383490f, 0.07383490f), vec4(0.07383490f, 0.07383490f, 0.77849529f, 0.07383490f), vec4(0.07383490f, 0.07383490f, 0.07383490f, 0.77849529f), vec4(0.40624434f, 0.40624434f, 0.09375565f, 0.09375565f), vec4(0.40624434f, 0.09375565f, 0.40624434f, 0.09375565f), vec4(0.40624434f, 0.09375565f, 0.09375565f, 0.40624434f), vec4(0.09375565f, 0.40624434f, 0.40624434f, 0.09375565f), vec4(0.09375565f, 0.40624434f, 0.09375565f, 0.40624434f), vec4(0.09375565f, 0.09375565f, 0.40624434f, 0.40624434f), };
template<> vec2 QuadratureUvs<FemShape::Quad, 1>[1];
template<> vec2 QuadratureUvs<FemShape::Quad, 4>[4];
template<> vec2 QuadratureUvs<FemShape::Quad, 9>[9];
template<> vec3 QuadratureUvs<FemShape::Hex, 1>[1];
template<> vec3 QuadratureUvs<FemShape::Hex, 8>[8];
template<> vec3 QuadratureUvs<FemShape::Hex, 27>[27];

template<FemShape Shape, uint32_t Points> float QuadratureWeights = 0;
template<> float QuadratureWeights<FemShape::Tri, 1>[1] = { 1.0f, };
template<> float QuadratureWeights<FemShape::Tri, 3>[3] = { 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f, };
template<> float QuadratureWeights<FemShape::Tri, 7>[7] = { 0.2250000000f, 0.1323941527f, 0.1323941527f, 0.1323941527f, 0.1259391805f, 0.1259391805f, 0.1259391805f, };
template<> float QuadratureWeights<FemShape::Tet, 1>[1] = { 1.0f, };
template<> float QuadratureWeights<FemShape::Tet, 4>[4] = { 1.0f / 4.0f, 1.0f / 4.0f, 1.0f / 4.0f, 1.0f / 4.0f, };
template<> float QuadratureWeights<FemShape::Tet, 10>[10] = { 0.04763313f, 0.04763313f, 0.04763313f, 0.04763313f, 0.13491124f, 0.13491124f, 0.13491124f, 0.13491124f, 0.13491124f, 0.13491124f, };
template<> float QuadratureWeights<FemShape::Quad, 1>[1];
template<> float QuadratureWeights<FemShape::Quad, 4>[4];
template<> float QuadratureWeights<FemShape::Quad, 9>[9];
template<> float QuadratureWeights<FemShape::Hex, 1>[1];
template<> float QuadratureWeights<FemShape::Hex, 8>[8];
template<> float QuadratureWeights<FemShape::Hex, 27>[27];

// Fill out the quad and hex arrays programmatically
constexpr float LinearQuadratureUvs[6][5] = {
	{  0.0f,       0.0f,       0.0f,       0.0f,       0.0f      },
	{  0.0f,       0.0f,       0.0f,       0.0f,       0.0f      },
	{ -0.577350f,  0.577350f,  0.0f,       0.0f,       0.0f      },
	{ -0.774597f,  0.0f,       0.774597f,  0.0f,       0.0f      },
	{ -0.861136f, -0.339981f,  0.339981f,  0.861136f,  0.0f      },
	{ -0.90618f,  -0.538469f,  0.0f,       0.538469f,  0.90618f  },
};
constexpr float LinearQuadratureWeights[6][5] = {
	{  0.0f,       0.0f,       0.0f,       0.0f,       0.0f      },
	{  2.0f,       0.0f,       0.0f,       0.0f,       0.0f      },
	{  1.0f,       1.0f,       0.0f,       0.0f,       0.0f      },
	{  0.555556f,  0.888889f,  0.555556f,  0.0f,       0.0f      },
	{  0.347855f,  0.652145f,  0.652145f,  0.347855f,  0.0f      },
	{  0.236927f,  0.478629f,  0.568889f,  0.478629f,  0.236927f },
};
template <uint32_t Order>
void InitQuadratureUvsAndWeightsByOrder() {
	for (uint32_t y = 0; y < Order; y++) {
		for (uint32_t x = 0; x < Order; x++) {
			vec2 uv = vec2(LinearQuadratureUvs[Order][x], LinearQuadratureUvs[Order][y]);
			float weight = LinearQuadratureWeights[Order][x] * LinearQuadratureWeights[Order][y];
			QuadratureUvs<FemShape::Quad, Order * Order>[x + y * Order] = uv;
			QuadratureWeights<FemShape::Quad, Order * Order>[x + y * Order] = weight;
		}
	}
	for (uint32_t z = 0; z < Order; z++) {
		for (uint32_t y = 0; y < Order; y++) {
			for (uint32_t x = 0; x < Order; x++) {
				vec3 uv = vec3(LinearQuadratureUvs[Order][x], LinearQuadratureUvs[Order][y], LinearQuadratureUvs[Order][z]);
				float weight = LinearQuadratureWeights[Order][x] * LinearQuadratureWeights[Order][y] * LinearQuadratureWeights[Order][z];
				QuadratureUvs<FemShape::Hex, Order * Order * Order>[x + y * Order + z * Order * Order] = uv;
				QuadratureWeights<FemShape::Hex, Order * Order * Order>[x + y * Order + z * Order * Order] = weight;
			}
		}
	}
}
void InitQuadratureUvsAndWeights() {
	InitQuadratureUvsAndWeightsByOrder<1>();
	InitQuadratureUvsAndWeightsByOrder<2>();
	InitQuadratureUvsAndWeightsByOrder<3>();
}

//-----------------------------------------------------------------------------
// Code printing functions
template<FemShape Shape, uint32_t Points>
void PrintQuadratureWeights() {
	printf("template<> constexpr float QuadratureWeights<FemShape::%s, %d>[%d] = { ", ShapeName[(int)Shape], Points, Points);
	for (uint32_t p = 0; p < Points; ++p) {
		printf("%ff, ", QuadratureWeights<Shape, Points>[p]);
	}
	printf("};\n");
}
template<FemShape Shape, uint32_t Nodes, uint32_t Points>
void PrintShapeFunctions() {
	printf("template<> constexpr float ShapeFunctions<FemShape::%s, %d, %d>[%d][%d] = { ", ShapeName[(int)Shape], Nodes, Points, Points, Nodes);
	for (uint32_t p = 0; p < Points; ++p) {
		if (Nodes == 27 && Points == 27 && p % 9 == 0) { printf("\n	"); }
		printf("{ ");
		for (uint32_t node = 0; node < Nodes; node++) {
			float N = ShapeFunction<Shape, Nodes>(node, QuadratureUvs<Shape, Points>[p]);
			printf("%ff, ", N);
		}
		printf("}, ");
	}
	printf("};\n");
}
template<FemShape Shape, uint32_t Nodes, uint32_t Points>
void PrintShapeFunctionGradients() {
	const char* type = (Shape == FemShape::Tet || Shape == FemShape::Hex) ? "vec3" : "vec2";
	printf("template<> constexpr %s ShapeFunctionGradients<FemShape::%s, %d, %d>[%d][%d] = { ", type, ShapeName[(int)Shape], Nodes, Points, Points, Nodes);
	for (uint32_t p = 0; p < Points; ++p) {
		if (Nodes == 27 && Points == 27 && p % 9 == 0) { printf("\n	"); }
		printf("{ ");
		for (uint32_t node = 0; node < Nodes; node++) {
			auto dN = ShapeFunctionGradient<Shape, Nodes>(node, QuadratureUvs<Shape, Points>[p]);
			if constexpr (sizeof(dN) == sizeof(vec2)) {
				printf("vec2(% -ff,% -ff), ", dN.x, dN.y);
			} else {
				printf("vec3(% -ff,% -ff,% -ff), ", dN.x, dN.y, dN.z);
			}
		}
		printf("}, ");
	}
	printf("};\n");
}

//-----------------------------------------------------------------------------
// Program entry
int main() {
	InitQuadratureUvsAndWeights();

	printf("template<FemShape Shape, uint32_t Points> constexpr uint32_t QuadratureWeights = 0;\n");
	PrintQuadratureWeights<FemShape::Tri, 1>();
	PrintQuadratureWeights<FemShape::Tri, 3>();
	PrintQuadratureWeights<FemShape::Tri, 7>();
	PrintQuadratureWeights<FemShape::Tet, 1>();
	PrintQuadratureWeights<FemShape::Tet, 4>();
	PrintQuadratureWeights<FemShape::Tet, 10>();
	PrintQuadratureWeights<FemShape::Quad, 1>();
	PrintQuadratureWeights<FemShape::Quad, 4>();
	PrintQuadratureWeights<FemShape::Quad, 9>();
	PrintQuadratureWeights<FemShape::Hex, 1>();
	PrintQuadratureWeights<FemShape::Hex, 8>();
	PrintQuadratureWeights<FemShape::Hex, 27>();
	printf("\n");

	printf("template<FemShape Shape, uint32_t Nodes, uint32_t Points> constexpr uint32_t ShapeFunctions = 0;\n");
	PrintShapeFunctions<FemShape::Tri, 3, 1>();
	PrintShapeFunctions<FemShape::Tri, 3, 3>();
	PrintShapeFunctions<FemShape::Tri, 3, 7>();
	PrintShapeFunctions<FemShape::Tri, 6, 1>();
	PrintShapeFunctions<FemShape::Tri, 6, 3>();
	PrintShapeFunctions<FemShape::Tri, 6, 7>();
	PrintShapeFunctions<FemShape::Tet, 4, 1>();
	PrintShapeFunctions<FemShape::Tet, 4, 4>();
	PrintShapeFunctions<FemShape::Tet, 4, 10>();
	PrintShapeFunctions<FemShape::Tet, 10, 1>();
	PrintShapeFunctions<FemShape::Tet, 10, 4>();
	PrintShapeFunctions<FemShape::Tet, 10, 10>();
	PrintShapeFunctions<FemShape::Quad, 4, 1>();
	PrintShapeFunctions<FemShape::Quad, 4, 4>();
	PrintShapeFunctions<FemShape::Quad, 4, 9>();
	PrintShapeFunctions<FemShape::Quad, 9, 1>();
	PrintShapeFunctions<FemShape::Quad, 9, 4>();
	PrintShapeFunctions<FemShape::Quad, 9, 9>();
	PrintShapeFunctions<FemShape::Hex, 8, 1>();
	PrintShapeFunctions<FemShape::Hex, 8, 8>();
	PrintShapeFunctions<FemShape::Hex, 8, 27>();
	PrintShapeFunctions<FemShape::Hex, 27, 1>();
	PrintShapeFunctions<FemShape::Hex, 27, 8>();
	PrintShapeFunctions<FemShape::Hex, 27, 27>();
	printf("\n");

	printf("template<FemShape Shape, uint32_t Nodes, uint32_t Points> constexpr uint32_t ShapeFunctionGradients = 0;\n");
	PrintShapeFunctionGradients<FemShape::Tri, 3, 1>();
	PrintShapeFunctionGradients<FemShape::Tri, 3, 3>();
	PrintShapeFunctionGradients<FemShape::Tri, 3, 7>();
	PrintShapeFunctionGradients<FemShape::Tri, 6, 1>();
	PrintShapeFunctionGradients<FemShape::Tri, 6, 3>();
	PrintShapeFunctionGradients<FemShape::Tri, 6, 7>();
	PrintShapeFunctionGradients<FemShape::Tet, 4, 1>();
	PrintShapeFunctionGradients<FemShape::Tet, 4, 4>();
	PrintShapeFunctionGradients<FemShape::Tet, 4, 10>();
	PrintShapeFunctionGradients<FemShape::Tet, 10, 1>();
	PrintShapeFunctionGradients<FemShape::Tet, 10, 4>();
	PrintShapeFunctionGradients<FemShape::Tet, 10, 10>();
	PrintShapeFunctionGradients<FemShape::Quad, 4, 1>();
	PrintShapeFunctionGradients<FemShape::Quad, 4, 4>();
	PrintShapeFunctionGradients<FemShape::Quad, 4, 9>();
	PrintShapeFunctionGradients<FemShape::Quad, 9, 1>();
	PrintShapeFunctionGradients<FemShape::Quad, 9, 4>();
	PrintShapeFunctionGradients<FemShape::Quad, 9, 9>();
	PrintShapeFunctionGradients<FemShape::Hex, 8, 1>();
	PrintShapeFunctionGradients<FemShape::Hex, 8, 8>();
	PrintShapeFunctionGradients<FemShape::Hex, 8, 27>();
	PrintShapeFunctionGradients<FemShape::Hex, 27, 1>();
	PrintShapeFunctionGradients<FemShape::Hex, 27, 8>();
	PrintShapeFunctionGradients<FemShape::Hex, 27, 27>();
	printf("\n");

	return 0;
}
