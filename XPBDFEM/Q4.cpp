#include "Q4.h"

#include <assert.h>
#include "DebugGeo.h"
#include "Xpbd.h"

//-----------------------------------------------------------------------------
void Q4Block::InitBlock(const Settings& settings, Allocator* alloc, uint32_t width, uint32_t height, vec2 elemScale, float mass) {
	this->width = width;
	this->height = height;
	vertCount = width * height;
	X0 = alloc->Alloc<dvec2>(vertCount);
	X = alloc->Alloc<dvec2>(vertCount);
	O = alloc->Alloc<dvec2>(vertCount);
	V = alloc->Alloc<dvec2>(vertCount);
	w = alloc->Alloc<float>(vertCount);
	flags = alloc->Alloc<uint8_t>(vertCount);
	elemCount = (width - 1) * (height - 1);
	e = alloc->Alloc<Q4>(elemCount);

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
			V[i] = dvec2(0.0f);
			w[i] = 0.0f; // We will accumulate un-inverted mass in here as we construct elements
			flags[i] = x == 0 ? Geo::Left : x == width - 1 ? Geo::Right : 0;
		}
	}

	// Add elements
	uint32_t ie = 0;
	for (uint32_t y = 0; y < height - 1; y++) {
		for (uint32_t x = 0; x < width - 1; x++) {
			e[ie].i[0] = (x + 0) + (y + 0) * width;
			e[ie].i[1] = (x + 1) + (y + 0) * width;
			e[ie].i[2] = (x + 1) + (y + 1) * width;
			e[ie].i[3] = (x + 0) + (y + 1) * width;
			const uint32_t* is = e[ie].i;

			// The 4 corners of the shape we're mapping to
			vec2 P[4] = { vec2(0.0f), vec2(X[is[1]] - X[is[0]]), vec2(X[is[2]] - X[is[0]]), vec2(X[is[3]] - X[is[0]]) };

			float m[4];
			GenerateMasses(m, density, P);
			for (uint32_t i = 0; i < 4; i++) { w[is[i]] += m[i]; }

			GenerateElementParams(&e[ie].ep, P);

			++ie;
		}
	}
	assert(ie == elemCount);

	// Flip masses into inverse masses
	for (uint32_t i = 0; i < vertCount; i++) { w[i] = 1.0f / w[i]; }
}

void Q4Block::Constrain(const Settings& settings, float dt) {
	for (uint32_t i = 0; i < elemCount; i++) {
		SolveElement(dt, X, O, w, e[i].i, e[i].ep, settings);
	}
	for (uint32_t itr = 0; itr < settings.volumePasses; itr++) {
		for (uint32_t i = 0; i < elemCount; i++) {
			SolveElementVolume(dt, X, O, w, e[i].i, e[i].ep, settings);
		}
	}
}

void Q4Block::Damp(float damping) {
	for (uint32_t i = 0; i < elemCount; i++) {
		PbdDamp(X, V, w, e[i].i, min(1.0f, damping / e[i].ep.volume));
	}
}

void Q4Block::Render(const Settings& settings) {
	for (uint32_t y = 0; y < height; y++) {
		for (uint32_t x = 0; x < width; x++) {
			uint32_t i = x + y * width;
			bool locked =
				((settings.flags & Settings_LockLeft) && (flags[i] & Geo::Left)) ||
				((settings.flags & Settings_LockRight) && (flags[i] & Geo::Right));
			if (x < width - 1) { DLINE(vec3(vec2(X[i]), 0.0f), vec3(vec2(X[i + 1]), 0.0f), color); }
			if (y < height - 1) { DLINE(vec3(vec2(X[i]), 0.0f), vec3(vec2(X[i + width]), 0.0f), locked ? vec3(0.7f, 0.0f, 0.0f) : color); }
		}
	}
}

float Q4Block::CalculateVolume() const {
	float volume = 0.0f;
	for (uint32_t i = 0; i < elemCount; i++) {
		volume += CalculateElementVolume(X, e[i].i, e[i].ep);
	}
	return volume;
}
