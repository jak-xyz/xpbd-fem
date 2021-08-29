#include "DebugGeo.h"

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdarg.h>
#include <stdio.h>

//-----------------------------------------------------------------------------
const int8_t kEndChar = -128;
const int8_t kLineBreak = -127;
const int8_t kFont[] = {
	'\0', 0,0, 50,0, 50,100, 0,100, 0,0, kEndChar,
	' ', 25,0, 25,0, kEndChar,
	'%', 0,72, 8,72, 8,80, 0,80, 0,72, kLineBreak, 0,0, 40,80, kLineBreak, 32,0, 40,0, 40,8, 32,8, 32,0, kEndChar,
	'(', 22,0, 10,10, 0,50, 10,90, 22,100, kEndChar,
	')', 0,0, 12,10, 22,50, 12,90, 0,100, kEndChar,
	'+', 0,40, 40,40, kLineBreak, 20,60, 20,20, kEndChar,
	',', 15,5, 0,-20, kEndChar,
	'-', 0,40, 40,40, kEndChar,
	'.', 0,0, 6,0, 6,6, 0,6, 0,0, kEndChar,
	'/', 0,0, 40,100, kEndChar,
	'0', 25,0, 0,15, 0,85, 25,100, 50,85, 50,15, 25,0, kEndChar,
	'1', 0,0, 30,0, kLineBreak, 15,0, 15,100, 0,80, kEndChar,
	'2', 0,75, 10,100, 33,100, 45,80, 45,65, 25,30, 0,0, 50,0, kEndChar,
	'3', 4,83, 15,100, 38,100, 46,85, 46,70, 30,55, 20,55, 30,55, 50,38, 50,18, 30,0, 20,0, 0,22, kEndChar,
	'4', 40,0, 40,100, 0,30, 50,30, kEndChar,
	'5', 50,100, 10,100, 0,48, 33,55, 50,25, 25,0, 0,10, kEndChar,
	'6', 50,85, 40,100, 20,100, 5,90, 0,70, 0,20, 15,0, 35,0, 50,15, 50,40, 35,57, 15,57, 0,45, kEndChar,
	'7', 0,100, 50,100, 20,0, kEndChar,
	'8', 5,88, 15,100, 35,100, 45,88, 45,67, 3,45, 0,15, 15,0, 35,0, 50,15, 47,45, 5,67, 5,88, kEndChar,
	'9', 2,10, 12,0, 35,0, 45,10, 50,30, 50,85, 38,100, 15,100, 5,85, 5,65, 13,48, 40,48, 50,58, kEndChar,
	':', 0,10, 6,10, 6,16, 0,16, 0,10, kLineBreak, 0,55, 6,55, 6,61, 0,61, 0,55, kEndChar,
	'A', 0,0, 25,100, 50,0, kLineBreak, 10,40, 40,40, kEndChar,
	'B', 0,0, 0,100, 25,100, 50,75, 25,50, 0,50, 25,50, 50,25, 25,0, 0,0, kEndChar,
	'a', 42,50, 35,60, 15,50, 0,25, 12,0, 30,5, 42,50, 44,25, 50,0, kEndChar,
	'b', 0,0, 0,100, kLineBreak, 0,30, 25,60, 45,50, 50,15, 25,0, 0,10, kEndChar,
	'c', 47,45, 36,60, 10,57, 0,40, 0,15, 25,0, 40,0, 50,10, kEndChar,
	'd', 50,0, 50,100, kLineBreak, 50,40, 25,60, 5,50, 0,15, 25,0, 50,20, kEndChar,
	'e', 5,28, 50,28, 40,50, 25,60, 0,40, 0,20, 10,4, 35,0, 50,10, kEndChar,
	'f', -4,40, 22,40, kLineBreak, 10,0, 10,80, 20,100, 30,80, kEndChar,
	'g', 40,60, 43,40, 35,58, 12,60, 0,35, 0,15, 15,0, 35,0, 43,40, 55,-30, 45,-40, 5,-25, 40,0, 55,5, kEndChar,
	'h', 0,0, 0,100, kLineBreak, 0,30, 22,60, 40,60, 50,45, 50,0, kEndChar,
	'i', 0,60, 5,50, 5,0, kLineBreak, 10,80, 10,100, kEndChar,
	'j', 5,60, 5,-20, -5,-40,  -15,-20, kLineBreak, 10,80, 10,100, kEndChar,
	'k', 0,100, 0,0, kLineBreak, 38,65, 0,35, 45,0, kEndChar,
	'l', 0,100, 0,0, 15,15, kEndChar,
	'm', -5,60, 0,40, 0,0, kLineBreak, 0,40, 19,65, 38,35, 38,0, kLineBreak, 38,35, 55,60, 70,35, 70,0, kEndChar,
	'n', -5,60, 0,40, 0,0, kLineBreak, 0,40, 23,60, 40,58, 46,36, 46,0, kEndChar,
	'o', 0,30, 7,51, 25,60, 41,51, 50,30, 41,9, 25,0, 7,9, 0,30, kEndChar,
	'p', 0,60, 0,-40, kLineBreak, 0,40, 25,60, 50,50, 50,20, 25,0, 0,20, kEndChar,
	'q', 50,40, 40,56, 20,60, 5,40, 0,20, 20,0, 50,40, 15,-35, 35,-40, 50,-20, 35,0, kEndChar,
	'r', 0,60, 0,0, kLineBreak, 0,45, 25,60, 40,45, kEndChar,
	's', 38,47, 22,60, 3,44, 38,30, 40,10, 12,0, -10,5, 0,18, kEndChar,
	't', -4,40, 22,40, kLineBreak, 10,100, 10,5, 15,0, 30,15, kEndChar,
	'u', 0,60, 0,15, 15,0, 30,0, 40,20, 40,60, kLineBreak, 40,20, 45,0, 50,10, kEndChar,
	'v', 0,60, 23,0, 46,60, kEndChar,
	'w', 0,60, 15,0, 30,60, 45,0, 60,60, kEndChar,
	'x', 0,0, 50,60, kLineBreak, 0,60, 50,0, kEndChar,
	'y', 0,60, 23,0, kLineBreak, 50,60, 0,-40, kEndChar,
	'z', 0,60, 50,60, 0,0, 50,0, kEndChar,
};

//-----------------------------------------------------------------------------
void DebugGeo::Initialize() {
	for (uint32_t i = 0; i < 256; i++) { glyphOffsets[i] = 1; }
	for (uint32_t i = 0; i < sizeof(kFont); i++) {
		glyphOffsets[kFont[i]] = i + 1;
		while (kFont[i] != kEndChar) { ++i; }
	}
}

void DebugGeo::BeginFrame() {
	// DebugGeo is immediate mode, so reset verts and idxs every frame
	pointCount = 0;
	vertCount = 0;
	idxCount = 0;
}

void DebugGeo::AddLines(const vec3* a, uint32_t count, const vec3& color) {
	uint32_t ucolor = ColorVec3ToUint32(color);
	for (uint32_t i = 0; i < count; i++) {
		AddVert(a[i], ucolor);
		if (i > 0) {
			idxs[idxCount++] = vertCount - 2;
			idxs[idxCount++] = vertCount - 1;
		}
	}
}

void DebugGeo::AddSphere(const vec3& o, float r, const vec3& color) {
	uint32_t ucolor = ColorVec3ToUint32(color);
	uint32_t baseIdx = vertCount;
	const uint32_t kWarp = 7;
	const uint32_t kWeft = 3;
	for (uint32_t j = 0; j < kWeft; j++) {
		float phi = (float)M_PI * ((float)j - 0.5f * (float)(kWeft - 1)) / (float)(kWeft + 1);
		float sy = sin(phi);
		float cy = cos(phi);
		for (uint32_t i = 0; i < kWarp; i++) {
			float theta = 2.0f * (float)M_PI * (float)i / (float)kWarp;
			float sx = sin(theta);
			float cx = cos(theta);

			vec3 v = o + r * vec3(sx * -cy, sy, cx * -cy);
			AddVert(v, ucolor);

			idxs[idxCount++] = baseIdx + i + kWarp * j;
			idxs[idxCount++] = baseIdx + (i + 1) % kWarp + kWarp * j;
			if (j > 0) {
				idxs[idxCount++] = baseIdx + i + kWarp * j;
				idxs[idxCount++] = baseIdx + i + kWarp * (j - 1);
			}
		}
	}
	// Caps
	AddVert(o - r * vec3(0.0f, 1.0f, 0.0f), ucolor);
	AddVert(o + r * vec3(0.0f, 1.0f, 0.0f), ucolor);
	for (uint32_t i = 0; i < kWarp; i++) {
		idxs[idxCount++] = vertCount - 2;
		idxs[idxCount++] = baseIdx + i;
		idxs[idxCount++] = vertCount - 1;
		idxs[idxCount++] = baseIdx + i + kWarp * (kWeft - 1);
	}
}

void DebugGeo::AddBox(const vec3& center, const vec3& dims, const mat3& rot, const vec3& color) {
	uint32_t ucolor = ColorVec3ToUint32(color);
	vec3 x = dims.x * rot[0];
	vec3 y = dims.y * rot[1];
	vec3 z = dims.z * rot[2];
	AddVert(center - x - y - z, ucolor);
	AddVert(center + x - y - z, ucolor);
	AddVert(center - x + y - z, ucolor);
	AddVert(center + x + y - z, ucolor);
	AddVert(center - x - y + z, ucolor);
	AddVert(center + x - y + z, ucolor);
	AddVert(center - x + y + z, ucolor);
	AddVert(center + x + y + z, ucolor);
	idxs[idxCount++] = vertCount - 8;
	idxs[idxCount++] = vertCount - 7;
	idxs[idxCount++] = vertCount - 8;
	idxs[idxCount++] = vertCount - 6;
	idxs[idxCount++] = vertCount - 8;
	idxs[idxCount++] = vertCount - 4;
	idxs[idxCount++] = vertCount - 5;
	idxs[idxCount++] = vertCount - 6;
	idxs[idxCount++] = vertCount - 5;
	idxs[idxCount++] = vertCount - 7;
	idxs[idxCount++] = vertCount - 5;
	idxs[idxCount++] = vertCount - 1;
	idxs[idxCount++] = vertCount - 3;
	idxs[idxCount++] = vertCount - 4;
	idxs[idxCount++] = vertCount - 3;
	idxs[idxCount++] = vertCount - 1;
	idxs[idxCount++] = vertCount - 3;
	idxs[idxCount++] = vertCount - 7;
	idxs[idxCount++] = vertCount - 2;
	idxs[idxCount++] = vertCount - 1;
	idxs[idxCount++] = vertCount - 2;
	idxs[idxCount++] = vertCount - 4;
	idxs[idxCount++] = vertCount - 2;
	idxs[idxCount++] = vertCount - 6;
}

void DebugGeo::AddXyz(const vec3& o, const mat3& rot) {
	AddLine(o, o + safeNormalize(rot[0]), vec3(1.0f, 0.0f, 0.0f));
	AddLine(o, o + safeNormalize(rot[1]), vec3(0.0f, 1.0f, 0.0f));
	AddLine(o, o + safeNormalize(rot[2]), vec3(0.0f, 0.0f, 1.0f));
}

void DebugGeo::AddText(const vec3& o, const vec2& dims, const mat3& rot, const vec3& color, const char* str) {
	uint32_t ucolor = ColorVec3ToUint32(color);

	vec3 cursor = o;
	vec3 lineStart = o;
	for (const char* c = str; *c; c++) {
		if (*c == '\n') {
			lineStart -= 1.5f * dims.y * rot[1];
			cursor = lineStart;
			continue;
		}

		bool prevVert = false;
		float xMax = 0.0f;
		for (uint32_t i = glyphOffsets[*c]; kFont[i] != kEndChar; i++) {
			if (kFont[i] == kLineBreak) {
				prevVert = false;
				continue;
			}
			++i;
			vec2 fv = dims * 0.01f * vec2((float)kFont[i - 1], (float)kFont[i]);
			xMax = max(xMax, fv.x);
			AddVert(cursor + rot * vec3(fv, 0.0f), ucolor);
			if (prevVert) {
				idxs[idxCount++] = vertCount - 2;
				idxs[idxCount++] = vertCount - 1;
			}
			prevVert = true;
		}
		cursor += (xMax + 0.19f * dims.x) * rot[0];
	}
}

void DTEXT(const vec3& o, const vec2& dims, const mat3& rot, const vec3& color, const char* format, ...) {
	char buf[1024];
	va_list args;
	va_start(args, format);
	vsnprintf(buf, sizeof(buf), format, args);
	va_end(args);
	g_debugGeo->AddText(o, dims, rot, color, buf);
}
