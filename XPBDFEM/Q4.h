//-----------------------------------------------------------------------------
// Linear FEM quad.
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include "Fem.h"
#include "Geo.h"

//-----------------------------------------------------------------------------
struct Q4 {
	uint32_t i[4];
	ElementParamsLinear2d ep;
};

struct Q4Block : public Geo2d {
	void InitBlock(const Settings& settings, Allocator* alloc, uint32_t width, uint32_t height, vec2 elemScale, float mass);

	virtual void Constrain(const Settings& settings, float dt) override;
	virtual void Damp(float damping) override;
	virtual void Render(const Settings& settings) override;
	virtual float CalculateVolume() const override;

	uint32_t width, height;
	Q4* e;
	uint32_t elemCount = 0;
};
