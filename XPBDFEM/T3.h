//-----------------------------------------------------------------------------
// Linear FEM triangle.
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include "Fem.h"
#include "Geo.h"

//-----------------------------------------------------------------------------
struct T3 {
	uint32_t i[3];
	ElementParamsTriangle ep;
};

struct T3Block : public Geo2d {
	void InitBlock(const Settings& settings, Allocator* alloc, uint32_t width, uint32_t height, vec2 elemScale, float mass);

	virtual void Constrain(const Settings& settings, float dt) override;
	virtual void Damp(float damping) override;
	virtual void Render(const Settings& settings) override;
	virtual float CalculateVolume() const override;

	uint32_t width, height;
	T3* e;
	uint32_t* eOrderLinear;
	uint32_t* eOrderCheckerboard;
	uint32_t elemCount = 0;
};
