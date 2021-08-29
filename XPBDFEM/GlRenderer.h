//-----------------------------------------------------------------------------
// Emulates the wasm rendering setup. Simply creates some big GL buffers you
// can stuff point and line data into, and knows how to draw them.
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include "vectormath.h"

//-----------------------------------------------------------------------------
struct GlSimpleVertex {
	float pos[3];
	uint32_t color;
};

//-----------------------------------------------------------------------------
class GlRenderer {
public:
	void Initialize(uint32_t maxPointVerts, uint32_t maxLineVerts, uint32_t maxLineIdxs);
	void Terminate();

	void Render(
		GlSimpleVertex* pointVerts, uint32_t pointVertCount,
		GlSimpleVertex* lineVerts, uint32_t lineVertCount,
		uint32_t* lineIdxs, uint32_t lineIdxCount,
		const mat4& viewFromWorld, const mat4& projFromView
	);

	uint32_t maxPointVerts;
	uint32_t maxLineVerts;
	uint32_t maxLineIdxs;

	static const uint32_t kBufferCount = 3;
	uint32_t pvbos[kBufferCount];
	uint32_t lvbos[kBufferCount];
	uint32_t libos[kBufferCount];
	uint32_t ubos[kBufferCount];
	uint32_t bufferIdx;

	uint32_t program = 0;
};
