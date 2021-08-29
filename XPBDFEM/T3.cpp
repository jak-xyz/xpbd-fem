#include "T3.h"

#include <assert.h>
#include "debugout.h"
#include "DebugGeo.h"
#include "Xpbd.h"

//-----------------------------------------------------------------------------
void T3Block::InitBlock(const Settings& settings, Allocator* alloc, uint32_t width, uint32_t height, vec2 elemScale, float mass) {
	this->width = width;
	this->height = height;
	vertCount = width * height;
	X0 = alloc->Alloc<dvec2>(vertCount);
	X = alloc->Alloc<dvec2>(vertCount);
	O = alloc->Alloc<dvec2>(vertCount);
	V = alloc->Alloc<dvec2>(vertCount);
	w = alloc->Alloc<float>(vertCount);
	flags = alloc->Alloc<uint8_t>(vertCount);
	elemCount = 2 * (width - 1) * (height - 1);
	e = alloc->Alloc<T3>(elemCount);
	eOrderLinear = alloc->Alloc<uint32_t>(elemCount);
	eOrderCheckerboard = alloc->Alloc<uint32_t>(elemCount);

	vec2 bottomLeft = elemScale * vec2(-0.5f * (float)(width - 1), -0.5f * (float)(height - 1));
	float density = mass / ((elemScale.x * elemScale.y) * (float)((width - 1) * (height - 1)));

	// Add verts
	uint32_t randState = 1;
	auto RandF = [&]() {
		randState = ((uint64_t)randState * 48271) % 0x7fffffff;
		return (float)randState * 9.3132258e-10f - 1.0f;
	};
	for (uint32_t y = 0; y < height; y++) {
		uint32_t yEdge = y == 0 || y == height - 1 ? 1 : 0;
		for (uint32_t x = 0; x < width; x++) {
			uint32_t xEdge = x == 0 || x == width - 1 ? 1 : 0;
			uint32_t i = x + y * width;
			vec2 wiggle = settings.wonkiness * vec2(xEdge ? 0.0f : RandF(), yEdge ? 0.0f : RandF());
			O[i] = X[i] = X0[i] = dvec2(bottomLeft + elemScale * (vec2((float)x, (float)y) + 0.5f * wiggle));
			V[i] = dvec2(0.0);
			w[i] = 0.0f; // We will accumulate un-inverted mass in here as we construct elements
			flags[i] = x == 0 ? Geo::Left : x == width - 1 ? Geo::Right : 0;
		}
	}

	// Add elements
	uint32_t ie = 0;
	for (uint32_t y = 0; y < height - 1; y++) {
		for (uint32_t x = 0; x < width - 1; x++) {
			uint32_t i = x + y * width;
			uint32_t i10 = (x + 1) + (y + 0) * width;
			uint32_t i01 = (x + 0) + (y + 1) * width;
			uint32_t i11 = (x + 1) + (y + 1) * width;
			e[ie + 0].i[0] = i10; e[ie + 0].i[1] = i01; e[ie + 0].i[2] = i;
			e[ie + 1].i[0] = i01; e[ie + 1].i[1] = i10; e[ie + 1].i[2] = i11;
			for (uint32_t j = 0; j < 2; j++) {
				uint32_t* is = e[ie + j].i;
				vec2 P[3] = { vec2(X[is[0]] - X[is[2]]), vec2(X[is[1]] - X[is[2]]),  vec2(0.0f) };
				float m[3];
				GenerateMasses(m, density, P);
				for (uint32_t i = 0; i < 3; i++) { w[is[i]] += m[i]; }
				GenerateElementParams(&e[ie + j].ep, P);
			}
			ie += 2;
		}
	}
	assert(ie == elemCount);

	// Flip masses into inverse masses
	for (uint32_t i = 0; i < vertCount; i++) { w[i] = 1.0f / w[i]; }

	// Add alternate constraint orders
	for (uint32_t i = 0; i < elemCount; i++) { eOrderLinear[i] = i; }
	ie = 0;
	for (uint32_t color2 = 0; color2 < 2; color2++) {
		for (uint32_t color = 0; color < 2; color++) {
			for (uint32_t y = 0; y < height - 1; y++) {
				for (uint32_t x = (color + y) % 2; x < width - 1; x += 2) {
					eOrderCheckerboard[ie++] = color2 + 2 * (x + y * (width - 1));
				}
			}
		}
	}
}

void T3Block::Constrain(const Settings& settings, float dt) {
	uint32_t constraintOrder = (settings.flags >> Settings_ConstraintOrderBit) & Settings_ConstraintOrderMask;
	const uint32_t* order = constraintOrder == 1 ? eOrderCheckerboard : eOrderLinear;
	for (uint32_t i = 0; i < elemCount; i++) {
		SolveElement(dt, X, O, w, e[order[i]].i, e[order[i]].ep, settings);
	}
	for (uint32_t itr = 0; itr < settings.volumePasses; itr++) {
		for (uint32_t i = 0; i < elemCount; i++) {
			SolveElementVolume(dt, X, O, w, e[order[i]].i, e[order[i]].ep, settings);
		}
	}

}

void T3Block::Damp(float damping) {
	for (uint32_t i = 0; i < elemCount; i++) {
		PbdDamp(X, V, w, e[i].i, min(1.0f, damping / e[i].ep.volume));
	}
}

void T3Block::Render(const Settings& settings) {
	for (uint32_t y = 0; y < height; y++) {
		for (uint32_t x = 0; x < width; x++) {
			uint32_t i = x + y * width;
			bool locked =
				((settings.flags & Settings_LockLeft) && (flags[i] & Geo::Left)) ||
				((settings.flags & Settings_LockRight) && (flags[i] & Geo::Right));
			if (x < width - 1) { DLINE(vec3(vec2(X[i]), 0.0f), vec3(vec2(X[i + 1]), 0.0f), color); }
			if (y < height - 1) { DLINE(vec3(vec2(X[i]), 0.0f), vec3(vec2(X[i + width]), 0.0f), locked ? vec3(0.7f, 0.0f, 0.0f) : color); }
			if (x < width - 1 && y < height - 1) {
				DLINE(vec3(vec2(X[i]), 0.0f), vec3(vec2(X[i + 1 + width]), 0.0f), color);
			}
		}
	}
	return;


	for (uint32_t ei = 0; ei < elemCount; ei++) {
		uint32_t i = e[ei].i[0], j = e[ei].i[1], k = e[ei].i[2];
		vec2 centroid = vec2((1.0f / 3.0f) * (X[i] + X[j] + X[k]));
		vec3 A = vec3(mix(centroid, vec2(X[i]), 0.8f), 0.0f);
		vec3 B = vec3(mix(centroid, vec2(X[j]), 0.8f), 0.0f);
		vec3 C = vec3(mix(centroid, vec2(X[k]), 0.8f), 0.0f);
		DLINE(A, B, color);
		DLINE(B, C, color);
		DLINE(C, A, color);
	}
}

float T3Block::CalculateVolume() const {
	float area = 0.0f;
	for (uint32_t i = 0; i < elemCount; i++) {
		vec2 A = vec2(X[e[i].i[1]] - X[e[i].i[0]]);
		vec2 B = vec2(X[e[i].i[2]] - X[e[i].i[0]]);
		area += 0.5f * (A.x * B.y - A.y * B.x);
	}
	return area;
}
