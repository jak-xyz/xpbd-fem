#include "H27.h"

#include <assert.h>
#include "DebugGeo.h"
#include "Xpbd.h"

//-----------------------------------------------------------------------------
void H27Block::InitBlock(const Settings& settings, Allocator* alloc, uint32_t width, uint32_t height, uint32_t depth, vec3 elemScale, float mass) {
	assert((width + 1) % 2 == 0);
	assert((height + 1) % 2 == 0);
	assert((depth + 1) % 2 == 0);
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
	elemCount = (width - 1) * (height - 1) * (depth - 1) / 8;
	e = alloc->Alloc<H27>(elemCount);

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
				flags[i] |= ((x % 2) == 1 && (y % 2) == 1 && (z % 2) == 1) ? Pickable : 0;
			}
		}
	}
	// Fix up middle verts
	auto Idx = [&](uint32_t x, uint32_t y, uint32_t z) { return x + y * width + z * width * height; };
	for (uint32_t z = 0; z < depth; z += 2) {
		for (uint32_t y = 0; y < height; y += 2) {
			for (uint32_t x = 0; x < width; x += 2) {
				if (x < width - 2) {
					uint32_t i = Idx(x + 1, y, z);
					O[i] = X[i] = X0[i] = 0.5f * (X[Idx(x, y, z)] + X[Idx(x + 2, y, z)]);
				}
				if (y < height - 2) {
					uint32_t i = Idx(x, y + 1, z);
					O[i] = X[i] = X0[i] = 0.5f * (X[Idx(x, y, z)] + X[Idx(x, y + 2, z)]);
				}
				if (z < depth - 2) {
					uint32_t i = Idx(x, y, z + 1);
					O[i] = X[i] = X0[i] = 0.5f * (X[Idx(x, y, z)] + X[Idx(x, y, z + 2)]);
				}

				if (x < width - 2 && y < height - 2) {
					uint32_t i = Idx(x + 1, y + 1, z);
					O[i] = X[i] = X0[i] = 0.25f * (X[Idx(x, y, z)] + X[Idx(x + 2, y, z)] + X[Idx(x, y + 2, z)] + X[Idx(x + 2, y + 2, z)]);
				}
				if (y < height - 2 && z < depth - 2) {
					uint32_t i = Idx(x, y + 1, z + 1);
					O[i] = X[i] = X0[i] = 0.25f * (X[Idx(x, y, z)] + X[Idx(x, y + 2, z)] + X[Idx(x, y, z + 2)] + X[Idx(x, y + 2, z + 2)]);
				}
				if (z < depth - 2 && x < width - 2) {
					uint32_t i = Idx(x + 1, y, z + 1);
					O[i] = X[i] = X0[i] = 0.25f * (X[Idx(x, y, z)] + X[Idx(x, y, z + 2)] + X[Idx(x + 2, y, z)] + X[Idx(x + 2, y, z + 2)]);
				}

				if (x < width - 2 && y < height - 2 && z < depth - 2) {
					uint32_t i = Idx(x + 1, y + 1, z + 1);
					O[i] = X[i] = X0[i] = 0.125f *
						(X[Idx(x, y, z)] + X[Idx(x + 2, y, z)] + X[Idx(x, y + 2, z)] + X[Idx(x + 2, y + 2, z)]
						+X[Idx(x, y, z + 2)] + X[Idx(x + 2, y, z + 2)] + X[Idx(x, y + 2, z + 2)] + X[Idx(x + 2, y + 2, z + 2)]);
				}
			}
		}
	}

	// Add elements
	uint32_t ie = 0;
	for (uint32_t z = 0; z < depth - 1; z += 2) {
		for (uint32_t y = 0; y < height - 1; y += 2) {
			for (uint32_t x = 0; x < width - 1; x += 2) {
				for (uint32_t k = 0; k < 3; k++) {
					for (uint32_t j = 0; j < 3; j++) {
						for (uint32_t i = 0; i < 3; i++) {
							e[ie].i[i + j * 3 + k * 9] = (x + i) + (y + j) * width + (z + k) * width * height;
						}
					}
				}

				dvec3 dP[8] = {
					X[Idx(x, y, z)], X[Idx(x + 2, y, z)], X[Idx(x, y + 2, z)], X[Idx(x + 2, y + 2, z)],
					X[Idx(x, y, z + 2)], X[Idx(x + 2, y, z + 2)], X[Idx(x, y + 2, z + 2)], X[Idx(x + 2, y + 2, z + 2)]
				};
				vec3 P[8];
				for (uint32_t i = 0; i < 8; i++) {
					P[i] = vec3(dP[i] - dP[7]);
				}

				float m[27];
				GenerateMasses(m, density, P);
				for (uint32_t i = 0; i < 27; i++) { w[e[ie].i[i]] += m[i]; }

				GenerateElementParams(&e[ie].ep, P);

				++ie;
			}
		}
	}
	assert(ie == elemCount);

	// Flip masses into inverse masses
	for (uint32_t i = 0; i < vertCount; i++) { w[i] = 1.0f / w[i]; }
}

void H27Block::Constrain(const Settings& settings, float dt) {
	for (uint32_t i = 0; i < elemCount; i++) {
		SolveElement(dt, X, O, w, e[i].i, e[i].ep, settings);
	}
	for (uint32_t itr = 0; itr < settings.volumePasses; itr++) {
		for (uint32_t i = 0; i < elemCount; i++) {
			SolveElementVolume(dt, X, O, w, e[i].i, e[i].ep, settings);
		}
	}
}

void H27Block::Damp(float damping) {
	for (uint32_t i = 0; i < elemCount; i++) {
		PbdDamp(X, V, w, e[i].i, min(1.0f, damping / e[i].ep.surfaceArea));
	}
}

void H27Block::Render(const Settings& settings) {
	auto H27QUADLINE = [&](uint32_t i, uint32_t j, uint32_t k) {
		uint8_t f = flags[i] & flags[k];
		//if (!(f & Geo::Surface)) { return; }
		vec3 c = (f & Geo::Edge) ? color : (f & Geo::Surface) ? mix(color, vec3(1.0f), 0.2f) : mix(color, vec3(1.0f), 0.7f);
		if ((settings.flags & Settings_LockLeft) && (f & Geo::Left)) { c = vec3(0.7f, 0.0f, 0.0f); }
		if ((settings.flags & Settings_LockRight) && (f & Geo::Right)) { c = vec3(0.7f, 0.0f, 0.0f); }

		const uint32_t kCount = 7;
		vec3 P[kCount];
		for (uint32_t x = 0; x < kCount; x++) {
			float t = (float)x * (1.0f / (float)(kCount - 1)) * 2.0f - 1.0f;
			float Ni = (t * (t - 1.0f)) * 0.5f, Nj = -(t + 1.0f) * (t - 1.0f), Nk = ((t + 1.0f) * t) * 0.5f;
			P[x] = Ni * vec3(X[i]) + Nj * vec3(X[j]) + Nk * vec3(X[k]);
		}
		DLINES(P, kCount, c);
	};
	auto H27POINT = [&](uint32_t i) {
		if (!(flags[i] & Geo::Surface)) { return; }
		vec3 c = mix(color, vec3(1.0f), 0.25f);
		DPOINT(vec3(X[i]), c);
	};
	for (uint32_t z = 0; z < depth; z += 2) {
		for (uint32_t y = 0; y < height; y += 2) {
			for (uint32_t x = 0; x < width; x += 2) {
				uint32_t i = x + y * width + z * width * height;
				if (x < width - 2) { H27QUADLINE(i, i + 1, i + 2); }
				if (y < height - 2) { H27QUADLINE(i, i + width, i + 2 * width); }
				if (z < depth - 2) { H27QUADLINE(i, i + width * height, i + 2 * (width * height)); }
				if (x < width - 2 && y < height - 2) { H27POINT(i + 1 + width); }
				if (x < width - 2 && z < depth - 2) { H27POINT(i + 1 + width * height); }
				if (y < height - 2 && z < depth - 2) { H27POINT(i + width + width * height); }
			}
		}
	}
}

float H27Block::CalculateVolume() const {
	float volume = 0.0f;
	for (uint32_t i = 0; i < elemCount; i++) {
		volume += CalculateElementVolume(X, e[i].i, e[i].ep);
	}
	return volume;
}
