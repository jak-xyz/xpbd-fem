//-----------------------------------------------------------------------------
// Quadratic FEM rectangular prism.
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include "Fem.h"
#include "Geo.h"

//-----------------------------------------------------------------------------
struct H27 {
	uint32_t i[27];
	ElementParamsQuadratic3d ep;
};

struct H27Block : Geo3d {
	void InitBlock(const Settings& settings, Allocator* alloc, uint32_t width, uint32_t height, uint32_t depth, vec3 elemScale, float mass);

	virtual void Constrain(const Settings& settings, float dt) override;
	virtual void Damp(float damping) override;
	virtual void Render(const Settings& settings) override;
	virtual float CalculateVolume() const override;

	uint32_t width, height, depth;
	H27* e;
	uint32_t elemCount = 0;
};
