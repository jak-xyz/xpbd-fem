//-----------------------------------------------------------------------------
// The guts of the simulation code.
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include "vectormath.h"
#include "Settings.h"

//-----------------------------------------------------------------------------
template<typename Vec, uint32_t Nodes, uint32_t Points> struct QuadratureCoefficients;

// Triangle 
template<>
struct QuadratureCoefficients<vec2, 3, 1> {
	float w[1];
	mat2 Qi[1];
};
// Tetrahedron
template<>
struct QuadratureCoefficients<vec3, 4, 1> {
	float w[1];
	mat3 Qi[1];
};
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
	float C;
	float QQ[Nodes - 1];
	float QR[((Nodes - 1) * (Nodes - 2)) / 2];
};

struct ElementParamsTriangle {
	QuadratureCoefficients<vec2, 3, 1> qcLo;
	QuadratureCoefficients<vec2, 3, 1> qcHi;
	PrefactoredIntegrationCoefficients<vec2, 3> ic;
	float volume;
	float surfaceArea;
};
struct ElementParamsLinear2d {
	QuadratureCoefficients<vec2, 4, 1> qcLo;
	QuadratureCoefficients<vec2, 4, 4> qcHi;
	PrefactoredIntegrationCoefficients<vec2, 4> ic;
	float volume;
	float surfaceArea;
};
struct ElementParamsQuadratic2d {
	QuadratureCoefficients<vec2, 9, 4> qcLo;
	QuadratureCoefficients<vec2, 9, 9> qcHi;
	PrefactoredIntegrationCoefficients<vec2, 9> ic;
	float volume;
	float surfaceArea;
};
struct ElementParamsLinear3d {
	QuadratureCoefficients<vec3, 8, 1> qcLo;
	QuadratureCoefficients<vec3, 8, 8> qcHi;
	PrefactoredIntegrationCoefficients<vec3, 8> ic;
	float volume;
	float surfaceArea;
};
struct ElementParamsQuadratic3d {
	QuadratureCoefficients<vec3, 27, 8> qcLo;
	QuadratureCoefficients<vec3, 27, 27> qcHi;
	PrefactoredIntegrationCoefficients<vec3, 27> ic;
	float volume;
	float surfaceArea;
};

void GenerateMasses(float(&m)[3], float density, const vec2(&P)[3]);
void GenerateMasses(float(&m)[4], float density, const vec2(&P)[4]);
void GenerateMasses(float(&m)[9], float density, const vec2(&P)[4]);
void GenerateMasses(float(&m)[8], float density, const vec3(&P)[8]);
void GenerateMasses(float(&m)[27], float density, const vec3(&P)[8]);

void GenerateElementParams(ElementParamsTriangle* ep, const vec2(&P)[3]);
void GenerateElementParams(ElementParamsLinear2d* ep, const vec2(&P)[4]);
void GenerateElementParams(ElementParamsQuadratic2d* ep, const vec2(&P)[4]);
void GenerateElementParams(ElementParamsLinear3d* ep, const vec3(&P)[8]);
void GenerateElementParams(ElementParamsQuadratic3d* ep, const vec3(&P)[8]);

void SolveElement(float dt, dvec2* X, dvec2* O, float* w, const uint32_t(&is)[3], const ElementParamsTriangle& ep, const Settings& settings);
void SolveElement(float dt, dvec2* X, dvec2* O, float* w, const uint32_t(&is)[4], const ElementParamsLinear2d& ep, const Settings& settings);
void SolveElement(float dt, dvec2* X, dvec2* O, float* w, const uint32_t(&is)[9], const ElementParamsQuadratic2d& ep, const Settings& settings);
void SolveElement(float dt, dvec3* X, dvec3* O, float* w, const uint32_t(&is)[8], const ElementParamsLinear3d& ep, const Settings& settings);
void SolveElement(float dt, dvec3* X, dvec3* O, float* w, const uint32_t(&is)[27], const ElementParamsQuadratic3d& ep, const Settings& settings);

void SolveElementVolume(float dt, dvec2* X, dvec2* O, float* w, const uint32_t(&is)[3], const ElementParamsTriangle& ep, const Settings& settings);
void SolveElementVolume(float dt, dvec2* X, dvec2* O, float* w, const uint32_t(&is)[4], const ElementParamsLinear2d& ep, const Settings& settings);
void SolveElementVolume(float dt, dvec2* X, dvec2* O, float* w, const uint32_t(&is)[9], const ElementParamsQuadratic2d& ep, const Settings& settings);
void SolveElementVolume(float dt, dvec3* X, dvec3* O, float* w, const uint32_t(&is)[8], const ElementParamsLinear3d& ep, const Settings& settings);
void SolveElementVolume(float dt, dvec3* X, dvec3* O, float* w, const uint32_t(&is)[27], const ElementParamsQuadratic3d& ep, const Settings& settings);

float CalculateElementVolume(const dvec2* X, const uint32_t(&is)[4], const ElementParamsLinear2d& ep);
float CalculateElementVolume(const dvec2* X, const uint32_t(&is)[9], const ElementParamsQuadratic2d& ep);
float CalculateElementVolume(const dvec3* X, const uint32_t(&is)[8], const ElementParamsLinear3d& ep);
float CalculateElementVolume(const dvec3* X, const uint32_t(&is)[27], const ElementParamsQuadratic3d& ep);
