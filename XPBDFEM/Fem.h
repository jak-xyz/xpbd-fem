//-----------------------------------------------------------------------------
// The guts of the simulation code.
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include "vectormath.h"
#include "Settings.h"

//-----------------------------------------------------------------------------
enum FemShape { Tri, Tet, Quad, Hex };

//-----------------------------------------------------------------------------
template<typename Vec, uint32_t Nodes, uint32_t Points> struct QuadratureCoefficients;
// Quads
template<uint32_t Nodes, uint32_t Points>
struct QuadratureCoefficients<vec2, Nodes, Points> {
	float w[Points]; // Scaled weights for the quadrature points
	mat2 Ji[Points]; // Inverse map jacobian at each quadrature point
};
// Hexahedra
template<uint32_t Nodes, uint32_t Points>
struct QuadratureCoefficients<vec3, Nodes, Points> {
	float w[Points]; // Scaled weights for the quadrature points
	mat3 Ji[Points]; // Inverse map jacobian at each quadrature point
};

template<typename Vec, uint32_t Nodes>
struct PrefactoredIntegrationCoefficients { // Incompressible Neo-Hooken prefactored integration coefficients
	float QQ[Nodes - 1];
	float QR[((Nodes - 1) * (Nodes - 2)) / 2];
};

struct ElementParamsLinearTri {
	static const FemShape Shape = FemShape::Tri;
	static const uint32_t PointsLo = 1;
	static const uint32_t PointsHi = 1;
	mat2 Qi;
	PrefactoredIntegrationCoefficients<vec2, 3> ic;
	float volume;
	float surfaceArea;
};
struct ElementParamsQuadraticTri {
	static const FemShape Shape = FemShape::Tri;
	static const uint32_t PointsLo = 3;
	static const uint32_t PointsHi = 3;
	mat2 Qi;
	PrefactoredIntegrationCoefficients<vec2, 6> ic;
	float volume;
	float surfaceArea;
};
struct ElementParamsLinearTet {
	static const FemShape Shape = FemShape::Tet;
	static const uint32_t PointsLo = 1;
	static const uint32_t PointsHi = 1;
	mat3 Qi;
	PrefactoredIntegrationCoefficients<vec3, 4> ic;
	float volume;
	float surfaceArea;
};
struct ElementParamsQuadraticTet {
	static const FemShape Shape = FemShape::Tet;
	static const uint32_t PointsLo = 4;
	static const uint32_t PointsHi = 4;
	mat3 Qi;
	PrefactoredIntegrationCoefficients<vec3, 10> ic;
	float volume;
	float surfaceArea;
};
struct ElementParamsLinearQuad {
	static const FemShape Shape = FemShape::Quad;
	static const uint32_t PointsLo = 1;
	static const uint32_t PointsHi = 4;
	QuadratureCoefficients<vec2, 4, PointsLo> qcLo;
	QuadratureCoefficients<vec2, 4, PointsHi> qcHi;
	PrefactoredIntegrationCoefficients<vec2, 4> ic;
	float volume;
	float surfaceArea;
};
struct ElementParamsQuadraticQuad {
	static const FemShape Shape = FemShape::Quad;
	static const uint32_t PointsLo = 4;
	static const uint32_t PointsHi = 9;
	QuadratureCoefficients<vec2, 9, PointsLo> qcLo;
	QuadratureCoefficients<vec2, 9, PointsHi> qcHi;
	PrefactoredIntegrationCoefficients<vec2, 9> ic;
	float volume;
	float surfaceArea;
};
struct ElementParamsLinearHex {
	static const FemShape Shape = FemShape::Hex;
	static const uint32_t PointsLo = 1;
	static const uint32_t PointsHi = 8;
	QuadratureCoefficients<vec3, 8, PointsLo> qcLo;
	QuadratureCoefficients<vec3, 8, PointsHi> qcHi;
	PrefactoredIntegrationCoefficients<vec3, 8> ic;
	float volume;
	float surfaceArea;
	float dim;
};
struct ElementParamsQuadraticHex {
	static const FemShape Shape = FemShape::Hex;
	static const uint32_t PointsLo = 8;
	static const uint32_t PointsHi = 27;
	QuadratureCoefficients<vec3, 27, PointsLo> qcLo;
	QuadratureCoefficients<vec3, 27, PointsHi> qcHi;
	PrefactoredIntegrationCoefficients<vec3, 27> ic;
	float volume;
	float surfaceArea;
	float dim;
};

void InitFiniteElement(const dvec2* X, const uint32_t(&is)[3], float density, float* m, ElementParamsLinearTri* ep);
void InitFiniteElement(const dvec2* X, const uint32_t(&is)[6], float density, float* m, ElementParamsQuadraticTri* ep);
void InitFiniteElement(const dvec3* X, const uint32_t(&is)[4], float density, float* m, ElementParamsLinearTet* ep);
void InitFiniteElement(const dvec3* X, const uint32_t(&is)[10], float density, float* m, ElementParamsQuadraticTet* ep);
void InitFiniteElement(const dvec2* X, const uint32_t(&is)[4], float density, float* m, ElementParamsLinearQuad* ep);
void InitFiniteElement(const dvec2* X, const uint32_t(&is)[9], float density, float* m, ElementParamsQuadraticQuad* ep);
void InitFiniteElement(const dvec3* X, const uint32_t(&is)[8], float density, float* m, ElementParamsLinearHex* ep);
void InitFiniteElement(const dvec3* X, const uint32_t(&is)[27], float density, float* m, ElementParamsQuadraticHex* ep);

void SolveElement(float dt, dvec2* X, const dvec2* O, const float* w, const uint32_t(&is)[3], const ElementParamsLinearTri& ep, const Settings& settings);
void SolveElement(float dt, dvec2* X, const dvec2* O, const float* w, const uint32_t(&is)[6], const ElementParamsQuadraticTri& ep, const Settings& settings);
void SolveElement(float dt, dvec3* X, const dvec3* O, const float* w, const uint32_t(&is)[4], const ElementParamsLinearTet& ep, const Settings& settings);
void SolveElement(float dt, dvec3* X, const dvec3* O, const float* w, const uint32_t(&is)[10], const ElementParamsQuadraticTet& ep, const Settings& settings);
void SolveElement(float dt, dvec2* X, const dvec2* O, const float* w, const uint32_t(&is)[4], const ElementParamsLinearQuad& ep, const Settings& settings);
void SolveElement(float dt, dvec2* X, const dvec2* O, const float* w, const uint32_t(&is)[9], const ElementParamsQuadraticQuad& ep, const Settings& settings);
void SolveElement(float dt, dvec3* X, const dvec3* O, const float* w, const uint32_t(&is)[8], const ElementParamsLinearHex& ep, const Settings& settings);
void SolveElement(float dt, dvec3* X, const dvec3* O, const float* w, const uint32_t(&is)[27], const ElementParamsQuadraticHex& ep, const Settings& settings);

void SolveElementVolume(float dt, dvec2* X, const dvec2* O, const float* w, const uint32_t(&is)[3], const ElementParamsLinearTri& ep, const Settings& settings);
void SolveElementVolume(float dt, dvec2* X, const dvec2* O, const float* w, const uint32_t(&is)[6], const ElementParamsQuadraticTri& ep, const Settings& settings);
void SolveElementVolume(float dt, dvec3* X, const dvec3* O, const float* w, const uint32_t(&is)[4], const ElementParamsLinearTet& ep, const Settings& settings);
void SolveElementVolume(float dt, dvec3* X, const dvec3* O, const float* w, const uint32_t(&is)[10], const ElementParamsQuadraticTet& ep, const Settings& settings);
void SolveElementVolume(float dt, dvec2* X, const dvec2* O, const float* w, const uint32_t(&is)[4], const ElementParamsLinearQuad& ep, const Settings& settings);
void SolveElementVolume(float dt, dvec2* X, const dvec2* O, const float* w, const uint32_t(&is)[9], const ElementParamsQuadraticQuad& ep, const Settings& settings);
void SolveElementVolume(float dt, dvec3* X, const dvec3* O, const float* w, const uint32_t(&is)[8], const ElementParamsLinearHex& ep, const Settings& settings);
void SolveElementVolume(float dt, dvec3* X, const dvec3* O, const float* w, const uint32_t(&is)[27], const ElementParamsQuadraticHex& ep, const Settings& settings);

void DampElement(float dt, dvec2* X, dvec2* V, const float* w, const uint32_t(&is)[3], const ElementParamsLinearTri& ep, const Settings& settings);
void DampElement(float dt, dvec2* X, dvec2* V, const float* w, const uint32_t(&is)[6], const ElementParamsQuadraticTri& ep, const Settings& settings);
void DampElement(float dt, dvec3* X, dvec3* V, const float* w, const uint32_t(&is)[4], const ElementParamsLinearTet& ep, const Settings& settings);
void DampElement(float dt, dvec3* X, dvec3* V, const float* w, const uint32_t(&is)[10], const ElementParamsQuadraticTet& ep, const Settings& settings);
void DampElement(float dt, dvec2* X, dvec2* V, const float* w, const uint32_t(&is)[4], const ElementParamsLinearQuad& ep, const Settings& settings);
void DampElement(float dt, dvec2* X, dvec2* V, const float* w, const uint32_t(&is)[9], const ElementParamsQuadraticQuad& ep, const Settings& settings);
void DampElement(float dt, dvec3* X, dvec3* V, const float* w, const uint32_t(&is)[8], const ElementParamsLinearHex& ep, const Settings& settings);
void DampElement(float dt, dvec3* X, dvec3* V, const float* w, const uint32_t(&is)[27], const ElementParamsQuadraticHex& ep, const Settings& settings);

void CalculateElementEnergyAndGradients(float dt, const dvec2* X, const dvec2* O, const float* w, const uint32_t(&is)[3], const ElementParamsLinearTri& ep, uint32_t pressurePoint, float oneRingVolume, const Settings& settings, float(&__restrict U)[2], vec2* (&__restrict g)[2][3]);
void CalculateElementEnergyAndGradients(float dt, const dvec2* X, const dvec2* O, const float* w, const uint32_t(&is)[6], const ElementParamsQuadraticTri& ep, uint32_t pressurePoint, float oneRingVolume, const Settings& settings, float(&__restrict U)[2], vec2* (&__restrict g)[2][6]);
void CalculateElementEnergyAndGradients(float dt, const dvec2* X, const dvec2* O, const float* w, const uint32_t(&is)[4], const ElementParamsLinearQuad& ep, uint32_t pressurePoint, float oneRingVolume, const Settings& settings, float(&__restrict U)[2], vec2* (&__restrict g)[2][4]);
void CalculateElementEnergyAndGradients(float dt, const dvec2* X, const dvec2* O, const float* w, const uint32_t(&is)[9], const ElementParamsQuadraticQuad& ep, uint32_t pressurePoint, float oneRingVolume, const Settings& settings, float(&__restrict U)[2], vec2* (&__restrict g)[2][9]);
void CalculateElementEnergyAndGradients(float dt, const dvec3* X, const dvec3* O, const float* w, const uint32_t(&is)[4], const ElementParamsLinearTet& ep, uint32_t pressurePoint, float oneRingVolume, const Settings& settings, float(&__restrict U)[2], vec3* (&__restrict g)[2][4]);
void CalculateElementEnergyAndGradients(float dt, const dvec3* X, const dvec3* O, const float* w, const uint32_t(&is)[10], const ElementParamsQuadraticTet& ep, uint32_t pressurePoint, float oneRingVolume, const Settings& settings, float(&__restrict U)[2], vec3* (&__restrict g)[2][10]);
void CalculateElementEnergyAndGradients(float dt, const dvec3* X, const dvec3* O, const float* w, const uint32_t(&is)[8], const ElementParamsLinearHex& ep, uint32_t pressurePoint, float oneRingVolume, const Settings& settings, float(&__restrict U)[2], vec3* (&__restrict g)[2][8]);
void CalculateElementEnergyAndGradients(float dt, const dvec3* X, const dvec3* O, const float* w, const uint32_t(&is)[27], const ElementParamsQuadraticHex& ep, uint32_t pressurePoint, float oneRingVolume, const Settings& settings, float(&__restrict U)[2], vec3* (&__restrict g)[2][27]);

float CalculateElementVolume(const dvec2* X, const uint32_t(&is)[3], const ElementParamsLinearTri& ep);
float CalculateElementVolume(const dvec2* X, const uint32_t(&is)[6], const ElementParamsQuadraticTri& ep);
float CalculateElementVolume(const dvec3* X, const uint32_t(&is)[4], const ElementParamsLinearTet& ep);
float CalculateElementVolume(const dvec3* X, const uint32_t(&is)[10], const ElementParamsQuadraticTet& ep);
float CalculateElementVolume(const dvec2* X, const uint32_t(&is)[4], const ElementParamsLinearQuad& ep);
float CalculateElementVolume(const dvec2* X, const uint32_t(&is)[9], const ElementParamsQuadraticQuad& ep);
float CalculateElementVolume(const dvec3* X, const uint32_t(&is)[8], const ElementParamsLinearHex& ep);
float CalculateElementVolume(const dvec3* X, const uint32_t(&is)[27], const ElementParamsQuadraticHex& ep);
