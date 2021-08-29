//-----------------------------------------------------------------------------
// Provides an easy global interface for other systems to throw in debug
// geometry (eg. DLINE(), DSPHERE(), etc.)
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include "vectormath.h"

//-----------------------------------------------------------------------------
inline uint32_t ColorVec3ToUint32(const vec3& v) {
	uint32_t r = (uint32_t)clamp(256.0f * v.r, 0.0f, 255.0f);
	uint32_t g = (uint32_t)clamp(256.0f * v.g, 0.0f, 255.0f);
	uint32_t b = (uint32_t)clamp(256.0f * v.b, 0.0f, 255.0f);
	uint32_t a = 255;
	return r | (g << 8) | (b << 16) | (a << 24);
}

struct DebugVertex {
	float pos[3];
	uint32_t color;
};

inline class DebugGeo* g_debugGeo = nullptr;

class DebugGeo {
public:
	DebugGeo() {
		g_debugGeo = this;
	}

	void Initialize();
	
	void BeginFrame();

	void AddPoint(const vec3& a, uint32_t color) {
		points[pointCount].pos[0] = a[0];
		points[pointCount].pos[1] = a[1];
		points[pointCount].pos[2] = a[2];
		points[pointCount].color = color;
		++pointCount;
	}
	void AddVert(const vec3& a, uint32_t color) {
		verts[vertCount].pos[0] = a[0];
		verts[vertCount].pos[1] = a[1];
		verts[vertCount].pos[2] = a[2];
		verts[vertCount].color = color;
		++vertCount;
	}
	void AddLine(const vec3& a, const vec3& b, const vec3& color) {
		uint32_t ucolor = ColorVec3ToUint32(color);
		AddVert(a, ucolor);
		AddVert(b, ucolor);
		idxs[idxCount++] = vertCount - 2;
		idxs[idxCount++] = vertCount - 1;
	}
	void AddLines(const vec3* a, uint32_t count, const vec3& color);
	void AddSphere(const vec3& o, float r, const vec3& color);
	void AddBox(const vec3& center, const vec3& dims, const mat3& rot, const vec3& color);
	void AddXyz(const vec3& o, const mat3& rot);
	void AddText(const vec3& o, const vec2& dims, const mat3& rot, const vec3& color, const char* str);

	static const uint32_t kMaxVerts = 128 * 1024;
	static const uint32_t kMaxIdxs = 4 * kMaxVerts;
	uint32_t pointCount = 0;
	DebugVertex points[kMaxVerts];
	uint32_t vertCount = 0;
	DebugVertex verts[kMaxVerts];
	uint32_t idxCount = 0;
	uint32_t idxs[kMaxIdxs];

	uint32_t glyphOffsets[256];
};

inline void DPOINT(const vec3& a, const vec3& color) {
	g_debugGeo->AddPoint(a, ColorVec3ToUint32(color));
}

inline void DLINE(const vec3& a, const vec3& b, const vec3& color) {
	g_debugGeo->AddLine(a, b, color);
}

inline void DLINES(const vec3* a, uint32_t count, const vec3& color) {
	g_debugGeo->AddLines(a, count, color);
}

inline void DVEC(const vec3& o, const vec3& e, const vec3& color) {
	g_debugGeo->AddLine(o, o + e, color);
}

inline void DSPHERE(const vec3& o, float r, const vec3& color) {
	g_debugGeo->AddSphere(o, r, color);
}

inline void DBOX(const vec3& center, const vec3& dims, const mat3& rot, const vec3& color) {
	g_debugGeo->AddBox(center, dims, rot, color);
}

inline void DAABB(const vec3& aabbMin, const vec3& aabbMax, const vec3& color) {
	DBOX(0.5f * (aabbMin + aabbMax), 0.5f * (aabbMax - aabbMin), mat3(1.0f), color);
}

inline void DXYZ(const vec3& o, const mat3& rot) {
	g_debugGeo->AddXyz(o, rot);
}

void DTEXT(const vec3& o, const vec2& dims, const mat3& rot, const vec3& color, const char* format, ...);
