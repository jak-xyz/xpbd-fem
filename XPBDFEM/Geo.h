//-----------------------------------------------------------------------------
// High level orchestration of the creation and simulation of FEM geometry.
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include "Allocator.h"
#include "Fem.h"
#include "Manipulator.h"
#include "Connectivity.h"
#include "Settings.h"
#include "vectormath.h"

//-----------------------------------------------------------------------------
struct Geo {
	static const uint8_t Left = 1 << 0;
	static const uint8_t Right = 1 << 1;
	static const uint8_t Surface = 1 << 2;
	static const uint8_t Edge = 1 << 3;
	static const uint8_t Pickable = 1 << 4;

	virtual void Substep(const Settings& settings, const Manipulator& manip, float dt) = 0;
	virtual void Transform(mat3 transform) = 0;
	virtual void Pick(vec3 rayOrigin, vec3 rayDir, vec3* outNearestPoint, uint32_t* outNearestPointIdx, float* outDistance) = 0;

	virtual void Constrain(const Settings& settings, float dt) = 0;
	virtual void Damp(const Settings& settings, float dt) = 0;
	virtual void Render(const Settings& settings) = 0;
	virtual float CalculateVolume() const = 0;
	virtual uint32_t VertCount() const = 0;
	virtual uint32_t ElementCount() const = 0;

	vec3 color = vec3(1.0f);
	float volume0;
};

//-----------------------------------------------------------------------------
struct Geo2d : public Geo {
	void Substep(const Settings& settings, const Manipulator& manip, float dt) final;
	void Transform(mat3 transform) final;
	void Pick(vec3 rayOrigin, vec3 rayDir, vec3* outNearestPoint, uint32_t* outNearestPointIdx, float* outDistance) final;
	uint32_t VertCount() const final { return vertCount; }

	dvec2* X0;
	dvec2* X;
	dvec2* O;
	dvec2* V;
	float* w;
	uint8_t* flags;
	uint32_t vertCount = 0;

	vec2 origin = vec2(0.0f);
};

struct Geo3d : public Geo {
	void Substep(const Settings& settings, const Manipulator& manip, float dt) final;
	void Transform(mat3 transform) final;
	void Pick(vec3 rayOrigin, vec3 rayDir, vec3* outNearestPoint, uint32_t* outNearestPointIdx, float* outDistance) final;
	uint32_t VertCount() const final { return vertCount; }

	dvec3* X0;
	dvec3* X;
	dvec3* O;
	dvec3* V;
	float* w;
	uint8_t* flags;
	uint32_t vertCount = 0;

	vec3 origin = vec3(0.0f);
};

//-----------------------------------------------------------------------------
struct T3 {
	uint32_t i[3];
	ElementParamsLinearTri ep;
};
struct Q4 {
	uint32_t i[4];
	ElementParamsLinearQuad ep;
};
struct T6 {
	uint32_t i[6];
	ElementParamsQuadraticTri ep;
};
struct Q9 {
	uint32_t i[9];
	ElementParamsQuadraticQuad ep;
};
struct T4 {
	uint32_t i[4];
	ElementParamsLinearTet ep;
};
struct H8 {
	uint32_t i[8];
	ElementParamsLinearHex ep;
};
struct T10 {
	uint32_t i[10];
	ElementParamsQuadraticTet ep;
};
struct H27 {
	uint32_t i[27];
	ElementParamsQuadraticHex ep;
};

struct GeoLinear2d : public Geo2d {
	void Init(Allocator* alloc, float density, const float* nodeData, uint32_t nodeDataCount, const uint32_t* idxData, uint32_t idxDataCount, bool autoResize);

	void Constrain(const Settings& settings, float dt) final;
	void Damp(const Settings& settings, float dt) final;
	void Render(const Settings& settings) final;
	float CalculateVolume() const final;
	uint32_t ElementCount() const final { return con.triCount + con.quadCount; }

	Connectivity2d con;
	T3* t;
	Q4* q;

	OneRing** rings;

	uint32_t* tOrder;
	uint32_t* qOrder;
	uint32_t* vOrder;
};

struct GeoQuadratic2d : public Geo2d {
	void Init(Allocator* alloc, float density, const float* nodeData, uint32_t nodeDataCount, const uint32_t* idxData, uint32_t idxDataCount, bool autoResize);

	void Constrain(const Settings& settings, float dt) final;
	void Damp(const Settings& settings, float dt) final;
	void Render(const Settings& settings) final;
	float CalculateVolume() const final;
	uint32_t ElementCount() const final { return con.triCount + con.quadCount; }

	Connectivity2d con;
	T6* t;
	Q9* q;

	OneRing** rings;

	uint32_t* tOrder;
	uint32_t* qOrder;
	uint32_t* vOrder;
};

struct GeoLinear3d : public Geo3d {
	void Init(Allocator* alloc, float density, const float* nodeData, uint32_t nodeDataCount, const uint32_t* idxData, uint32_t idxDataCount, bool autoResize);

	void Constrain(const Settings& settings, float dt) final;
	void Damp(const Settings& settings, float dt) final;
	void Render(const Settings& settings) final;
	float CalculateVolume() const final;
	uint32_t ElementCount() const final { return con.tetCount + con.hexCount; }

	Connectivity3d con;
	T4* t;
	H8* h;

	OneRing** rings;

	uint32_t* tOrder;
	uint32_t* hOrder;
	uint32_t* vOrder;
};

struct GeoQuadratic3d : public Geo3d {
	void Init(Allocator* alloc, float density, const float* nodeData, uint32_t nodeDataCount, const uint32_t* idxData, uint32_t idxDataCount, bool autoResize);

	void Constrain(const Settings& settings, float dt) final;
	void Damp(const Settings& settings, float dt) final;
	void Render(const Settings& settings) final;
	float CalculateVolume() const final;
	uint32_t ElementCount() const final { return con.tetCount + con.hexCount; }

	Connectivity3d con;
	T10* t;
	H27* h;

	OneRing** rings;

	uint32_t* tOrder;
	uint32_t* hOrder;
	uint32_t* vOrder;
};
