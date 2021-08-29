//------------------------------------------------------------------------------
// Easy access to more-or-less constant device information.
//------------------------------------------------------------------------------
#pragma once

#include "vectormath.h"

inline class Device* g_device = nullptr;
inline const Device& GetDevice() { return *g_device; }

class Device {
public:
	Device() {
		g_device = this;
	}

	vec2 viewportDims = vec2(100.0f);
	float ppi = 96.0f;
};

inline float PixelsToInches(float f) { return f / g_device->ppi; }
inline vec2 PixelsToInches(const vec2& v) { return v / g_device->ppi; }
inline float InchesToPixels(float f) { return f * g_device->ppi; }
inline vec2 InchesToPixels(const vec2& v) { return v * g_device->ppi; }
inline vec2 PixelsToNdc(const vec2& v) {
	return vec2(2.0f, -2.0f) * (v / g_device->viewportDims - vec2(0.5));
}
