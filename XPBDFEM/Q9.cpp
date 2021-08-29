#include "Q9.h"

#include <assert.h>
#include "DebugGeo.h"
#include "Xpbd.h"

//-----------------------------------------------------------------------------
void Q9Block::InitBlock(const Settings& settings, Allocator* alloc, uint32_t width, uint32_t height, vec2 elemScale, float mass) {
	assert((width + 1) % 2 == 0);
	assert((height + 1) % 2 == 0);
	this->width = width;
	this->height = height;
	vertCount = width * height;
	X0 = alloc->Alloc<dvec2>(vertCount);
	X = alloc->Alloc<dvec2>(vertCount);
	O = alloc->Alloc<dvec2>(vertCount);
	V = alloc->Alloc<dvec2>(vertCount);
	w = alloc->Alloc<float>(vertCount);
	flags = alloc->Alloc<uint8_t>(vertCount);
	elemCount = (width - 1) * (height - 1) / 4;
	e = alloc->Alloc<Q9>(elemCount);
	eOrderCheckerboard = alloc->Alloc<uint32_t>(elemCount);
	eOrderSimd = alloc->Alloc<uint32_t>(elemCount);

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
			O[i] = X[i] = X0[i] = dvec2(bottomLeft + elemScale * (vec2((float)x, (float)y) + wiggle));
			V[i] = dvec2(0.0f);
			w[i] = 0.0f; // We will accumulate un-inverted mass in here as we construct elements
			flags[i] = x == 0 ? Geo::Left : x == width - 1 ? Geo::Right : 0;
		}
	}
	// Fix up middle verts
	for (uint32_t y = 0; y < height; y += 2) {
		for (uint32_t x = 0; x < width; x += 2) {
			if (x < width - 2) {
				uint32_t i = (x + 1) + (y + 0) * width;
				O[i] = X[i] = X0[i] = 0.5f * (X[(x + 0) + (y + 0) * width] + X[(x + 2) + (y + 0) * width]);
			}
			if (y < height - 2) {
				uint32_t i = (x + 0) + (y + 1) * width;
				O[i] = X[i] = X0[i] = 0.5f * (X[(x + 0) + (y + 0) * width] + X[(x + 0) + (y + 2) * width]);
			}
			if (x < width - 2 && y < height - 2) {
				uint32_t i = (x + 1) + (y + 1) * width;
				O[i] = X[i] = X0[i] = 0.25f * (X[(x + 0) + (y + 0) * width] + X[(x + 2) + (y + 0) * width] + X[(x + 0) + (y + 2) * width] + X[(x + 2) + (y + 2) * width]);
			}
		}
	}

	// Add elements
	uint32_t ie = 0;
	for (uint32_t y = 0; y < height - 1; y += 2) {
		for (uint32_t x = 0; x < width - 1; x += 2) {
			// Element node layout:
			// 3---6---2
			// |       |
			// 7   8   5
			// |       |
			// 0---4---1
			e[ie].i[0] = (x + 0) + (y + 0) * width;
			e[ie].i[1] = (x + 2) + (y + 0) * width;
			e[ie].i[2] = (x + 2) + (y + 2) * width;
			e[ie].i[3] = (x + 0) + (y + 2) * width;
			e[ie].i[4] = (x + 1) + (y + 0) * width;
			e[ie].i[5] = (x + 2) + (y + 1) * width;
			e[ie].i[6] = (x + 1) + (y + 2) * width;
			e[ie].i[7] = (x + 0) + (y + 1) * width;
			e[ie].i[8] = (x + 1) + (y + 1) * width;
			const uint32_t* is = e[ie].i;

			// The 4 corners of the shape we're mapping to
			vec2 P[4] = { vec2(0.0f), vec2(X[is[1]] - X[is[0]]), vec2(X[is[2]] - X[is[0]]), vec2(X[is[3]] - X[is[0]]) };

			float m[9];
			GenerateMasses(m, density, P);
			for (uint32_t i = 0; i < 9; i++) { w[is[i]] += m[i]; }

			GenerateElementParams(&e[ie].ep, P);

			++ie;
		}
	}
	assert(ie == elemCount);

	// Flip masses into inverse masses
	for (uint32_t i = 0; i < vertCount; i++) { w[i] = 1.0f / w[i]; }

	// Add alternate constraint orders
	ie = 0;
	for (uint32_t color = 0; color < 2; color++) {
		for (uint32_t y = 0; y < (height - 1) / 2; y++) {
			for (uint32_t x = (color + y) % 2; x < (width - 1) / 2; x += 2) {
				eOrderCheckerboard[ie++] = x + y * ((width - 1) / 2);
			}
		}
	}
}

void Q9Block::Constrain(const Settings& settings, float dt) {
	uint32_t constraintOrder = (settings.flags >> Settings_ConstraintOrderBit) & Settings_ConstraintOrderMask;
	if (constraintOrder == 0) {
		for (uint32_t i = 0; i < elemCount; i++) {
			SolveElement(dt, X, O, w, e[i].i, e[i].ep, settings);
		}
	} else {
		const uint32_t* order = constraintOrder == 1 ? eOrderCheckerboard : eOrderSimd;
		for (uint32_t i = 0; i < elemCount; i++) {
			SolveElement(dt, X, O, w, e[i].i, e[i].ep, settings);
		}
	}
	for (uint32_t itr = 0; itr < settings.volumePasses; itr++) {
		for (uint32_t i = 0; i < elemCount; i++) {
			SolveElementVolume(dt, X, O, w, e[i].i, e[i].ep, settings);
		}
	}
}

void Q9Block::Damp(float damping) {
	for (uint32_t i = 0; i < elemCount; i++) {
		PbdDamp(X, V, w, e[i].i, min(1.0f, damping / e[i].ep.volume));
	}
}

void Q9Block::Render(const Settings& settings) {
	auto QUADRATICLINE = [&](uint32_t i, uint32_t j, uint32_t k, vec3 color) {
		const uint32_t kCount = 7;
		vec3 P[kCount];
		for (uint32_t x = 0; x < kCount; x++) {
			float t = (float)x * (1.0f / (float)(kCount - 1)) * 2.0f - 1.0f;
			float Ni = (t * (t - 1.0f)) * 0.5f, Nj = -(t + 1.0f) * (t - 1.0f), Nk = ((t + 1.0f) * t) * 0.5f;
			P[x] = vec3(Ni * vec2(X[i]) + Nj * vec2(X[j]) + Nk * vec2(X[k]), 0.0f);
		}
		DLINES(P, kCount, color);
	};
	for (uint32_t y = 0; y < height; y += 2) {
		for (uint32_t x = 0; x < width; x += 2) {
			uint32_t i = x + y * width;
			bool locked =
				((settings.flags & Settings_LockLeft) && (flags[i] & Geo::Left)) ||
				((settings.flags & Settings_LockRight) && (flags[i] & Geo::Right));
			if (x < width - 2) { QUADRATICLINE(i, i + 1, i + 2, color); }
			if (y < height - 2) { QUADRATICLINE(i, i + width, i + 2 * width, locked ? vec3(0.7f, 0.0f, 0.0f) : color); }
			if (x < width - 2 && y < height - 2) { DPOINT(vec3(vec2(X[i + 1 + width]), 0.0f), color); }
		}
	}
}

float Q9Block::CalculateVolume() const {
	float volume = 0.0f;
	for (uint32_t i = 0; i < elemCount; i++) {
		volume += CalculateElementVolume(X, e[i].i, e[i].ep);
	}
	return volume;
}
