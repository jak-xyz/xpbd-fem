// clear && clang++-12 main.cpp -std=c++17 -O2 && ./a.out

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include "../XPBDFEM/vectormath.h"
#include "../XPBDFEM/vectormath.cpp"

template<typename Vec, uint32_t Points> struct PointIterator;
template<uint32_t Points> struct PointIterator<vec2, Points> {
	uint32_t x = 0, y = 0, i = 0;
	PointIterator<vec2, Points>& operator++() {
		i++;
		y = i / Points;
		x = i - y * Points;
		return *this;
	}
};
template<uint32_t Points> struct PointIterator<vec3, Points> {
	uint32_t x = 0, y = 0, z = 0, i = 0;
	PointIterator<vec3, Points>& operator++() {
		i++;
		z = i / (Points * Points);
		y = (i - z * (Points * Points)) / Points;
		x = i - y * Points - z * (Points * Points);
		return* this;
	}
};

//-----------------------------------------------------------------------------
// Element node layouts
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

template<typename Vec, uint32_t Nodes> float Lagrange_N(uint32_t N, float t);
template<> float Lagrange_N<vec2, 4>(uint32_t N, float t) { return LagrangeLinear_N(N, t); }
template<> float Lagrange_N<vec2, 9>(uint32_t N, float t) { return LagrangeQuadratic_N(N, t); }
template<> float Lagrange_N<vec3, 8>(uint32_t N, float t) { return LagrangeLinear_N(N, t); }
template<> float Lagrange_N<vec3, 27>(uint32_t N, float t) { return LagrangeQuadratic_N(N, t); }
template<typename Vec, uint32_t Nodes> float Lagrange_dN(uint32_t N, float t);
template<> float Lagrange_dN<vec2, 4>(uint32_t N, float t) { return LagrangeLinear_dN(N, t); }
template<> float Lagrange_dN<vec2, 9>(uint32_t N, float t) { return LagrangeQuadratic_dN(N, t); }
template<> float Lagrange_dN<vec3, 8>(uint32_t N, float t) { return LagrangeLinear_dN(N, t); }
template<> float Lagrange_dN<vec3, 27>(uint32_t N, float t) { return LagrangeQuadratic_dN(N, t); }


//-----------------------------------------------------------------------------
// Gaussian Quadrature helpers
constexpr float QuadratureUvs[6][5] = {
	{  0.0f,       0.0f,       0.0f,       0.0f,       0.0f      },
	{  0.0f,       0.0f,       0.0f,       0.0f,       0.0f      },
	{ -0.577350f,  0.577350f,  0.0f,       0.0f,       0.0f      },
	{ -0.774597f,  0.0f,       0.774597f,  0.0f,       0.0f      },
	{ -0.861136f, -0.339981f,  0.339981f,  0.861136f,  0.0f      },
	{ -0.90618f,  -0.538469f,  0.0f,       0.538469f,  0.90618f  },
};
constexpr float QuadratureWeights[6][5] = {
	{  0.0f,       0.0f,       0.0f,       0.0f,       0.0f      },
	{  2.0f,       0.0f,       0.0f,       0.0f,       0.0f      },
	{  1.0f,       1.0f,       0.0f,       0.0f,       0.0f      },
	{  0.555556f,  0.888889f,  0.555556f,  0.0f,       0.0f      },
	{  0.347855f,  0.652145f,  0.652145f,  0.347855f,  0.0f      },
	{  0.236927f,  0.478629f,  0.568889f,  0.478629f,  0.236927f },
};

template<uint32_t Points> float CalculateQuadratureWeight(PointIterator<vec2, Points> p) {
	return QuadratureWeights[Points][p.x] * QuadratureWeights[Points][p.y];
}
template<uint32_t Points> float CalculateQuadratureWeight(PointIterator<vec3, Points> p) {
	return QuadratureWeights[Points][p.x] * QuadratureWeights[Points][p.y] * QuadratureWeights[Points][p.z];
}
template<uint32_t Points, uint32_t Nodes> float CalculateShapeFunction(PointIterator<vec2, Points> p, NodeLayout<vec2> n) {
	return Lagrange_N<vec2, Nodes>(n.x, QuadratureUvs[Points][p.x]) * Lagrange_N<vec2, Nodes>(n.y, QuadratureUvs[Points][p.y]);
}
template<uint32_t Points, uint32_t Nodes> float CalculateShapeFunction(PointIterator<vec3, Points> p, NodeLayout<vec3> n) {
	return Lagrange_N<vec3, Nodes>(n.x, QuadratureUvs[Points][p.x]) * Lagrange_N<vec3, Nodes>(n.y, QuadratureUvs[Points][p.y]) * Lagrange_N<vec3, Nodes>(n.z, QuadratureUvs[Points][p.z]);
}
template<uint32_t Points, uint32_t Nodes> vec2 CalculateShapeFunctionGradient(PointIterator<vec2, Points> p, NodeLayout<vec2> n) {
	return vec2(
		Lagrange_dN<vec2, Nodes>(n.x, QuadratureUvs[Points][p.x]) * Lagrange_N<vec2, Nodes>(n.y, QuadratureUvs[Points][p.y]),
		Lagrange_N<vec2, Nodes>(n.x, QuadratureUvs[Points][p.x]) * Lagrange_dN<vec2, Nodes>(n.y, QuadratureUvs[Points][p.y]));
}
template<uint32_t Points, uint32_t Nodes> vec3 CalculateShapeFunctionGradient(PointIterator<vec3, Points> p, NodeLayout<vec3> n) {
	return vec3(
		Lagrange_dN<vec3, Nodes>(n.x, QuadratureUvs[Points][p.x]) * Lagrange_N<vec3, Nodes>(n.y, QuadratureUvs[Points][p.y]) * Lagrange_N<vec3, Nodes>(n.z, QuadratureUvs[Points][p.z]),
		Lagrange_N<vec3, Nodes>(n.x, QuadratureUvs[Points][p.x]) * Lagrange_dN<vec3, Nodes>(n.y, QuadratureUvs[Points][p.y]) * Lagrange_N<vec3, Nodes>(n.z, QuadratureUvs[Points][p.z]),
		Lagrange_N<vec3, Nodes>(n.x, QuadratureUvs[Points][p.x]) * Lagrange_N<vec3, Nodes>(n.y, QuadratureUvs[Points][p.y]) * Lagrange_dN<vec3, Nodes>(n.z, QuadratureUvs[Points][p.z]));
}

template<uint32_t Points>
void PrintQuadratureWeights2d() {
	printf("template<> constexpr float QuadratureWeights<vec2, %d>[%d] = { ", Points * Points, Points * Points);
	for (PointIterator<vec2, Points> p; p.i < Points * Points; ++p) {
		printf("%ff, ", CalculateQuadratureWeight(p));
	}
	printf("};\n");
}
template<uint32_t Points>
void PrintQuadratureWeights3d() {
	printf("template<> constexpr float QuadratureWeights<vec3, %d>[%d] = { ", Points * Points * Points, Points * Points * Points);
	for (PointIterator<vec3, Points> p; p.i < Points * Points * Points; ++p) {
		printf("%ff, ", CalculateQuadratureWeight(p));
	}
	printf("};\n");
}

template<uint32_t Nodes, uint32_t Points>
void PrintShapeFunctions2d() {
	printf("template<> constexpr float ShapeFunctions<vec2, %d, %d>[%d][%d] = { ", Nodes, Points * Points, Points * Points, Nodes);
	for (PointIterator<vec2, Points> p; p.i < Points * Points; ++p) {
		printf("{ ");
		for (uint32_t node = 0; node < Nodes; node++) {
			printf("%ff, ", CalculateShapeFunction<Points, Nodes>(p, NodeLayouts<vec2, Nodes>[node]));
		}
		printf("}, ");
	}
	printf("};\n");
}
template<uint32_t Nodes, uint32_t Points>
void PrintShapeFunctions3d() {
	printf("template<> constexpr float ShapeFunctions<vec3, %d, %d>[%d][%d] = { ", Nodes, Points * Points * Points, Points * Points * Points, Nodes);
	uint32_t py = -1;
	for (PointIterator<vec3, Points> p; p.i < Points * Points * Points; ++p) {
		if (py != p.y) {
			py = p.y;
			printf("\n    ");
		}
		printf("{ ");
		for (uint32_t node = 0; node < Nodes; node++) {
			printf("%ff, ", CalculateShapeFunction<Points, Nodes>(p, NodeLayouts<vec3, Nodes>[node]));
		}
		printf("}, ");
	}
	printf("};\n");
}

template<uint32_t Nodes, uint32_t Points>
void PrintShapeFunctionGradients2d() {
	printf("template<> constexpr vec2 ShapeFunctionGradients<vec2, %d, %d>[%d][%d] = { ", Nodes, Points * Points, Points * Points, Nodes);
	for (PointIterator<vec2, Points> p; p.i < Points * Points; ++p) {
		printf("{ ");
		for (uint32_t node = 0; node < Nodes; node++) {
			vec2 g = CalculateShapeFunctionGradient<Points, Nodes>(p, NodeLayouts<vec2, Nodes>[node]);
			printf("vec2(% -ff,% -ff), ", g.x, g.y);
		}
		printf("}, ");
	}
	printf("};\n");
}
template<uint32_t Nodes, uint32_t Points>
void PrintShapeFunctionGradients3d() {
	printf("template<> constexpr vec3 ShapeFunctionGradients<vec3, %d, %d>[%d][%d] = { ", Nodes, Points * Points * Points, Points * Points * Points, Nodes);
	uint32_t py = -1;
	for (PointIterator<vec3, Points> p; p.i < Points * Points * Points; ++p) {
		if (py != p.y) {
			py = p.y;
			printf("\n    ");
		}
		printf("{ ");
		for (uint32_t node = 0; node < Nodes; node++) {
			vec3 g = CalculateShapeFunctionGradient<Points, Nodes>(p, NodeLayouts<vec3, Nodes>[node]);
			printf("vec3(% -ff,% -ff,% -ff), ", g.x, g.y, g.z);
		}
		printf("}, ");
	}
	printf("};\n");
}

int main() {
	printf("template<typename Vec, uint32_t Points> constexpr uint32_t QuadratureWeights = 0;\n");
	PrintQuadratureWeights2d<1>();
	PrintQuadratureWeights2d<2>();
	PrintQuadratureWeights2d<3>();
	// PrintQuadratureWeights2d<4>();
	// PrintQuadratureWeights2d<5>();
	PrintQuadratureWeights3d<1>();
	PrintQuadratureWeights3d<2>();
	PrintQuadratureWeights3d<3>();
	// PrintQuadratureWeights3d<4>();
	// PrintQuadratureWeights3d<5>();
	printf("\n");

	printf("template<typename Vec, uint32_t Nodes, uint32_t Points> constexpr uint32_t ShapeFunctions = 0;\n");
	PrintShapeFunctions2d<4, 1>();
	PrintShapeFunctions2d<4, 2>();
	PrintShapeFunctions2d<4, 3>();
	// PrintShapeFunctions2d<4, 4>();
	// PrintShapeFunctions2d<4, 5>();
	PrintShapeFunctions2d<9, 1>();
	PrintShapeFunctions2d<9, 2>();
	PrintShapeFunctions2d<9, 3>();
	// PrintShapeFunctions2d<9, 4>();
	// PrintShapeFunctions2d<9, 5>();
	PrintShapeFunctions3d<8, 1>();
	PrintShapeFunctions3d<8, 2>();
	PrintShapeFunctions3d<8, 3>();
	// PrintShapeFunctions3d<8, 4>();
	// PrintShapeFunctions3d<8, 5>();
	PrintShapeFunctions3d<27, 1>();
	PrintShapeFunctions3d<27, 2>();
	PrintShapeFunctions3d<27, 3>();
	// PrintShapeFunctions3d<27, 4>();
	// PrintShapeFunctions3d<27, 5>();
	printf("\n");

	printf("template<typename Vec, uint32_t Nodes, uint32_t Points> constexpr uint32_t ShapeFunctionGradients = 0;\n");
	PrintShapeFunctionGradients2d<4, 1>();
	PrintShapeFunctionGradients2d<4, 2>();
	PrintShapeFunctionGradients2d<4, 3>();
	// PrintShapeFunctionGradients2d<4, 4>();
	// PrintShapeFunctionGradients2d<4, 5>();
	PrintShapeFunctionGradients2d<9, 1>();
	PrintShapeFunctionGradients2d<9, 2>();
	PrintShapeFunctionGradients2d<9, 3>();
	// PrintShapeFunctionGradients2d<9, 4>();
	// PrintShapeFunctionGradients2d<9, 5>();
	PrintShapeFunctionGradients3d<8, 1>();
	PrintShapeFunctionGradients3d<8, 2>();
	PrintShapeFunctionGradients3d<8, 3>();
	// PrintShapeFunctionGradients3d<8, 4>();
	// PrintShapeFunctionGradients3d<8, 5>();
	PrintShapeFunctionGradients3d<27, 1>();
	PrintShapeFunctionGradients3d<27, 2>();
	PrintShapeFunctionGradients3d<27, 3>();
	// PrintShapeFunctionGradients3d<27, 4>();
	// PrintShapeFunctionGradients3d<27, 5>();
	printf("\n");

	return 0;
}
