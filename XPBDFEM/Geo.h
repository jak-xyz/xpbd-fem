//-----------------------------------------------------------------------------
// Base geometry class.
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include "Allocator.h"
#include "Manipulator.h"
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
	virtual void Damp(float damping) = 0;
	virtual void Render(const Settings& settings) = 0;
	virtual float CalculateVolume() const = 0;

	vec3 color = vec3(1.0f);
	float volume0;
};

struct Geo2d : public Geo {
	void Substep(const Settings& settings, const Manipulator& manip, float dt) final;
	void Transform(mat3 transform) final;
	void Pick(vec3 rayOrigin, vec3 rayDir, vec3* outNearestPoint, uint32_t* outNearestPointIdx, float* outDistance) final;

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

	dvec3* X0;
	dvec3* X;
	dvec3* O;
	dvec3* V;
	float* w;
	uint8_t* flags;
	uint32_t vertCount = 0;

	vec3 origin = vec3(0.0f);
};
