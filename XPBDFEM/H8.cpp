#include "H8.h"

#include <assert.h>
#include "DebugGeo.h"
#include "Xpbd.h"

//-----------------------------------------------------------------------------
void H8Block::InitBlock(const Settings& settings, Allocator* alloc, uint32_t width, uint32_t height, uint32_t depth, vec3 elemScale, float mass) {
	this->width = width;
	this->height = height;
	this->depth = depth;
	vertCount = width * height * depth;
	X0 = alloc->Alloc<dvec3>(vertCount);
	X = alloc->Alloc<dvec3>(vertCount);
	O = alloc->Alloc<dvec3>(vertCount);
	V = alloc->Alloc<dvec3>(vertCount);
	w = alloc->Alloc<float>(vertCount);
	flags = alloc->Alloc<uint8_t>(vertCount);
	elemCount = (width - 1) * (height - 1) * (depth - 1);
	e = alloc->Alloc<H8>(elemCount);

	vec3 minX = elemScale * vec3(-0.5f * (float)(width - 1), -0.5f * (float)(height - 1), -0.5f * (float)(depth - 1));
	float density = mass / ((elemScale.x * elemScale.y * elemScale.z) * (float)((width - 1) * (height - 1) * (depth - 1)));

	// Add verts
	uint32_t randState = 1;
	auto RandF = [&]() {
		randState = ((uint64_t)randState * 48271) % 0x7fffffff;
		return (float)randState * 9.3132258e-10f - 1.0f;
	};
	for (uint32_t z = 0; z < depth; z++) {
		uint32_t zEdge = z == 0 || z == depth - 1 ? 1 : 0;
		for (uint32_t y = 0; y < height; y++) {
			uint32_t yEdge = y == 0 || y == height - 1 ? 1 : 0;
			for (uint32_t x = 0; x < width; x++) {
				uint32_t xEdge = x == 0 || x == width - 1 ? 1 : 0;
				uint32_t i = x + y * width + z * width * height;
				vec3 wiggle = settings.wonkiness * vec3(xEdge ? 0.0f : RandF(), yEdge ? 0.0f : RandF(), zEdge ? 0.0f : RandF());
				O[i] = X[i] = X0[i] = dvec3(minX + elemScale * (vec3((float)x, (float)y, (float)z) + 0.5f * wiggle));
				V[i] = dvec3(0.0);
				w[i] = 0.0f; // We will accumulate un-inverted mass in here as we construct elements
				flags[i] = x == 0 ? Geo::Left : x == width - 1 ? Geo::Right : 0;
				flags[i] |= xEdge + yEdge + zEdge >= 1 ? Geo::Surface : 0;
				flags[i] |= xEdge + yEdge + zEdge >= 2 ? Geo::Edge : 0;
				flags[i] |= Geo::Pickable;
			}
		}
	}

	// Add elements
	uint32_t ie = 0;
	for (uint32_t z = 0; z < depth - 1; z++) {
		for (uint32_t y = 0; y < height - 1; y++) {
			for (uint32_t x = 0; x < width - 1; x++) {
				for (uint32_t k = 0; k < 2; k++) {
					for (uint32_t j = 0; j < 2; j++) {
						for (uint32_t i = 0; i < 2; i++) {
							e[ie].i[i + j * 2 + k * 4] = (x + i) + (y + j) * width + (z + k) * width * height;
						}
					}
				}

				vec3 P[8];
				for (uint32_t i = 0; i < 8; i++) {
					P[i] = vec3(X[e[ie].i[i]] - X[e[ie].i[7]]);
				}

				float m[8];
				GenerateMasses(m, density, P);
				for (uint32_t i = 0; i < 8; i++) { w[e[ie].i[i]] += m[i]; }

				GenerateElementParams(&e[ie].ep, P);

				++ie;
			}
		}
	}
	assert(ie == elemCount);

	// Flip masses into inverse masses
	for (uint32_t i = 0; i < vertCount; i++) { w[i] = 1.0f / w[i]; }
}

void H8Block::Constrain(const Settings& settings, float dt) {
	for (uint32_t i = 0; i < elemCount; i++) {
		SolveElement(dt, X, O, w, e[i].i, e[i].ep, settings);
	}
	for (uint32_t itr = 0; itr < settings.volumePasses; itr++) {
		for (uint32_t i = 0; i < elemCount; i++) {
			SolveElementVolume(dt, X, O, w, e[i].i, e[i].ep, settings);
		}
	}
}

void H8Block::Damp(float damping) {
	for (uint32_t i = 0; i < elemCount; i++) {
		PbdDamp(X, V, w, e[i].i, min(1.0f, damping / e[i].ep.surfaceArea));
	}
}

void H8Block::Render(const Settings& settings) {
	auto H8DLINE = [&](uint32_t i0, uint32_t i1) {
		uint8_t f = flags[i0] & flags[i1];
		//if (!(f & Geo::Surface)) { return; }
		vec3 c = (f & Geo::Edge) ? color : (f & Geo::Surface) ? mix(color, vec3(1.0f), 0.2f) : mix(color, vec3(1.0f), 0.7f);
		if ((settings.flags & Settings_LockLeft) && (f & Geo::Left)) { c = vec3(0.7f, 0.0f, 0.0f); }
		if ((settings.flags & Settings_LockRight) && (f & Geo::Right)) { c = vec3(0.7f, 0.0f, 0.0f); }
		DLINE(vec3(X[i0]), vec3(X[i1]), c);
	};
	for (uint32_t z = 0; z < depth; z++) {
		for (uint32_t y = 0; y < height; y++) {
			for (uint32_t x = 0; x < width; x++) {
				uint32_t i = x + y * width + z * width * height;
				if (x < width - 1) { H8DLINE(i, i + 1); }
				if (y < height - 1) { H8DLINE(i, i + width); }
				if (z < depth - 1) { H8DLINE(i, i + width * height); }
			}
		}
	}
}

float H8Block::CalculateVolume() const {
	float volume = 0.0f;
	for (uint32_t i = 0; i < elemCount; i++) {
		volume += CalculateElementVolume(X, e[i].i, e[i].ep);
	}
	return volume;
}
