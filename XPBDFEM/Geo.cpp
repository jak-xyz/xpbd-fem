#include "Geo.h"
#include <assert.h>
#include <float.h>
#include <stdio.h>
#include "DebugGeo.h"
#include "Xpbd.h"

//-----------------------------------------------------------------------------
// Fem.cpp outgrowth
const uint32_t kMaxRingNodes = 256;

template<typename Dvec, typename Vec>
void SolveSimplex(uint32_t Nodes, float dt, Dvec* X, const Dvec* O, Dvec* V, const float* w, const uint32_t* is, float volume, const Settings& settings, float(&__restrict U)[2], Vec* __restrict g) {
	float I1 = U[0];
	float J = U[1];
	Vec* g0 = g;
	Vec* g1 = g + Nodes;

	uint32_t energyFunction = (settings.flags >> Settings_EnergyBit) & Settings_EnergyMask;
	uint32_t rayleighDampingType = (settings.flags >> Settings_RayleighTypeBit) & Settings_RayleighTypeMask;

	if (energyFunction == Energy_ContinuousPixar) {
		float mu = 2.0f;
		float lambda = mu / (1.0f - 2.0f * settings.poissonsRatio);
		float a = 1.0f + mu / lambda;
		float comp = settings.compliance / volume;
		float damp = settings.damping;

		// Combine the I1 and and J terms in proper proportions
		U[0] = max(0.0000001f, (mu / 2.0f) * I1 + (lambda / 2.0f) * (J - a) * (J - a));
		for (uint32_t n = 0; n < Nodes; n++) { g0[n] = (mu / 2.0f) * g0[n] + (lambda / 2.0f) * 2.0f * (J - a) * g1[n]; }

		if (V) {
			RayleighDamp(Nodes, V, w, is, U[0], g0, dt, comp, damp);
		} else {
			EnergyXpbdConstrain(Nodes, X, O, w, is, U[0], g0, dt, comp, damp, rayleighDampingType);
		}
	} else {
		float mu = 1.0f / settings.compliance;
		float lambda = (2.0f * mu * settings.poissonsRatio) / (1.0f - 2.0f * settings.poissonsRatio);
		float a = 1.0f + mu / lambda;
		float comp[2] = { 1.0f / (mu * volume), 1.0f / (lambda * volume) };
		float damp[2] = { settings.damping, settings.damping };

		if (energyFunction == Energy_ContinuousSkin) {
			const float CSkin[3] = { 0.1095f, 14.95f, 4.595f };
			const float(&C)[3] = CSkin;
			float IM;
			if constexpr (sizeof(Vec) == sizeof(vec2)) { IM = I1 - 2.0f; } else { IM = I1 - 3.0f; }
			U[0] = max(0.0001f, C[0] * IM + C[1] * IM * IM + C[2] * IM * IM * IM);
			float g0Scale = C[0] + 2.0f * C[1] * IM + 3.0f * C[2] * IM * IM;
			for (uint32_t n = 0; n < Nodes; n++) { g0[n] *= g0Scale; }
		}

		U[1] = (J - a) * (J - a);
		float g1Scale = 2.0f * (J - a);
		for (uint32_t n = 0; n < Nodes; n++) { g1[n] *= g1Scale; }

		if (V) {
			RayleighDamp(Nodes, V, w, is, U, g, dt, comp, damp);
		} else {
			EnergyXpbdConstrainSimultaneous(Nodes, X, O, w, is, U, g, dt, comp, damp, rayleighDampingType);
		}
	}
}

template<typename T, typename Q>
void SolveRing(float dt, dvec2* X, dvec2* O, dvec2* V, float* w, const Connectivity2d& con, const T* t, const Q* q, const OneRing& ring, const Settings& settings) {
	uint32_t Nodes = ring.vertCount;
	const uint32_t TVertCount = sizeof(T::i) / sizeof(T::i[0]);
	const uint32_t QVertCount = sizeof(Q::i) / sizeof(Q::i[0]);

	const uint32_t* verts = ring.Verts();
	const uint32_t* polys = ring.Polys();
	const uint8_t* sizes = ring.Sizes();
	const uint16_t* idxs = ring.Idxs();

	uint32_t energyFunction = (settings.flags >> Settings_EnergyBit) & Settings_EnergyMask;
	uint32_t rayleighDampingType = (settings.flags >> Settings_RayleighTypeBit) & Settings_RayleighTypeMask;

	if (V && settings.damping <= 0.0f) { return; }

	// Calculate the volume affected by the constraint
	float volume = 0.0f;
	for (uint32_t i = 0; i < ring.polyCount; i++) {
		volume += (sizes[i] == TVertCount) ? ((1.0f / 3.0f) * t[polys[i]].ep.volume) : ((1.0f / 4.0f) * q[polys[i] - con.triCount].ep.volume);
	}

	// Gather energies and gradients
	float U[2] = { 0.0f, 0.0f };
	vec2 g[2 * kMaxRingNodes];
	for (uint32_t n = 0; n < Nodes; n++) { g[0 * Nodes + n] = vec2(0.0f); g[1 * Nodes + n] = vec2(0.0f); }

	for (uint32_t i = 0; i < ring.polyCount; i++) {
		if (sizes[i] == TVertCount) {
			const T& ti = t[polys[i]];
			uint32_t pressurePoint = 0;
			vec2* eg[2][TVertCount];
			for (uint32_t j = 0; j < TVertCount; j++) { eg[0][j] = &g[0 * Nodes + idxs[j]]; eg[1][j] = &g[1 * Nodes + idxs[j]]; }
			CalculateElementEnergyAndGradients(dt, X, O, w, ti.i, ti.ep, pressurePoint, volume, settings, U, eg);
		} else {
			const Q& qi = q[con.triCount + polys[i]];
			uint32_t pressurePoint = 0;
			vec2* eg[2][QVertCount];
			for (uint32_t j = 0; j < QVertCount; j++) { eg[0][j] = &g[0 * Nodes + idxs[j]]; eg[1][j] = &g[1 * Nodes + idxs[j]]; }
			CalculateElementEnergyAndGradients(dt, X, O, w, qi.i, qi.ep, pressurePoint, volume, settings, U, eg);
		}
		idxs += sizes[i];
	}

	// Apply the correct constraints
	if (TVertCount == 3 && con.triCount > 0) {
		SolveSimplex(Nodes, dt, X, O, V, w, verts, volume, settings, U, (vec2*)g);
		return;
	}
	if (energyFunction == Energy_ContinuousPixar) {
		float mu = 2.0f;
		float lambda = mu / (1.0f - 2.0f * settings.poissonsRatio);
		float comp = settings.compliance / volume;
		float damp = settings.damping;

		// Combine the I1 and and J terms in proper proportions
		U[0] = max(0.000001f, (mu / 2.0f) * U[0] + (lambda / 2.0f) * U[1]);
		for (uint32_t n = 0; n < Nodes; n++) { g[0 * Nodes + n] = (mu / 2.0f) * g[0 * Nodes + n] + (lambda / 2.0f) * g[1 * Nodes + n]; }

		if (V) {
			RayleighDamp(Nodes, V, w, verts, U[0], g, dt, comp, damp);
		} else {
			EnergyXpbdConstrain(Nodes, X, O, w, verts, U[0], g, dt, comp, damp, rayleighDampingType);
		}
	} else {
		float mu = 1.0f / settings.compliance;
		float lambda = (2.0f * mu * settings.poissonsRatio) / (1.0f - 2.0f * settings.poissonsRatio);
		float comp[2] = { 1.0f / (mu * volume), 1.0f / (lambda * volume) };
		float damp[2] = { settings.damping, settings.damping };

		U[0] = max(0.0001f, U[0]);

		if (V) {
			RayleighDamp(Nodes, V, w, verts, U, g, dt, comp, damp);
		} else {
			EnergyXpbdConstrainSimultaneous(Nodes, X, O, w, verts, U, g, dt, comp, damp, rayleighDampingType);
		}
	}
}

template<typename T, typename H>
void SolveRing(float dt, dvec3* X, dvec3* O, dvec3* V, float* w, const Connectivity3d& con, const T* t, const H* h, const OneRing& ring, const Settings& settings) {
	uint32_t Nodes = ring.vertCount;
	const uint32_t TVertCount = sizeof(T::i) / sizeof(T::i[0]);
	const uint32_t HVertCount = sizeof(H::i) / sizeof(H::i[0]);

	const uint32_t* verts = ring.Verts();
	const uint32_t* polys = ring.Polys();
	const uint8_t* sizes = ring.Sizes();
	const uint16_t* idxs = ring.Idxs();

	uint32_t energyFunction = (settings.flags >> Settings_EnergyBit) & Settings_EnergyMask;
	uint32_t rayleighDampingType = (settings.flags >> Settings_RayleighTypeBit) & Settings_RayleighTypeMask;

	if (V && settings.damping <= 0.0f) { return; }

	// Calculate the volume affected by the constraint
	float volume = 0.0f;
	for (uint32_t i = 0; i < ring.polyCount; i++) {
		volume += (sizes[i] == TVertCount) ? ((1.0f / 4.0f) * t[polys[i]].ep.volume) : ((1.0f / 8.0f) * h[polys[i] - con.tetCount].ep.volume);
	}

	// Gather energies and gradients
	float U[2] = { 0.0f, 0.0f };
	vec3 g[2 * kMaxRingNodes];
	for (uint32_t n = 0; n < Nodes; n++) { g[0 * Nodes + n] = vec3(0.0f); g[1 * Nodes + n] = vec3(0.0f); }

	for (uint32_t i = 0; i < ring.polyCount; i++) {
		if (sizes[i] == TVertCount) {
			const T& ti = t[polys[i]];
			uint32_t pressurePoint = 0;
			vec3* eg[2][TVertCount];
			for (uint32_t j = 0; j < TVertCount; j++) { eg[0][j] = &g[0 * Nodes + idxs[j]]; eg[1][j] = &g[1 * Nodes + idxs[j]]; }
			CalculateElementEnergyAndGradients(dt, X, O, w, ti.i, ti.ep, pressurePoint, volume, settings, U, eg);
		} else {
			const H& hi = h[con.tetCount + polys[i]];
			uint32_t pressurePoint = 0;
			vec3* eg[2][HVertCount];
			for (uint32_t j = 0; j < HVertCount; j++) { eg[0][j] = &g[0 * Nodes + idxs[j]]; eg[1][j] = &g[1 * Nodes + idxs[j]]; }
			CalculateElementEnergyAndGradients(dt, X, O, w, hi.i, hi.ep, pressurePoint, volume, settings, U, eg);
		}
		idxs += sizes[i];
	}

	// Apply the correct constraints
	if (TVertCount == 4 && con.tetCount > 0) {
		SolveSimplex(Nodes, dt, X, O, V, w, verts, volume, settings, U, (vec3*)g);
		return;
	}
	if (energyFunction == Energy_ContinuousPixar) {
		float mu = 2.0f;
		float lambda = mu / (1.0f - 2.0f * settings.poissonsRatio);
		float comp = settings.compliance / volume;
		float damp = settings.damping;

		// Combine the I1 and and J terms in proper proportions
		U[0] = max(0.000001f, (mu / 2.0f) * U[0] + (lambda / 2.0f) * U[1]);
		for (uint32_t n = 0; n < Nodes; n++) { g[0 * Nodes + n] = (mu / 2.0f) * g[0 * Nodes + n] + (lambda / 2.0f) * g[1 * Nodes + n]; }

		if (V) {
			RayleighDamp(Nodes, V, w, verts, U[0], g, dt, comp, damp);
		} else {
			EnergyXpbdConstrain(Nodes, X, O, w, verts, U[0], g, dt, comp, damp, rayleighDampingType);
		}
	} else {
		float mu = 1.0f / settings.compliance;
		float lambda = (2.0f * mu * settings.poissonsRatio) / (1.0f - 2.0f * settings.poissonsRatio);
		float comp[2] = { 1.0f / (mu * volume), 1.0f / (lambda * volume) };
		float damp[2] = { settings.damping, settings.damping };

		U[0] = max(0.0001f, U[0]);

		if (V) {
			RayleighDamp(Nodes, V, w, verts, U, g, dt, comp, damp);
		} else {
			EnergyXpbdConstrainSimultaneous(Nodes, X, O, w, verts, U, g, dt, comp, damp, rayleighDampingType);
		}
	}
}

//-----------------------------------------------------------------------------
// Geo2d methods
void Geo2d::Substep(const Settings& settings, const Manipulator& manip, float dt) {
	// Predict
	for (uint32_t i = 0; i < vertCount; i++) {
		V[i] += dvec2(settings.gravity * dt);
		V[i] *= 1.0f - settings.timeCorrectedDrag;
		O[i] = X[i];
		X[i] += V[i] * dt;
	}

	// Apply constraints
	Constrain(settings, dt);

	// Undo any motion on locked verts
	if (settings.flags & Settings_LockLeft) {
		for (uint32_t i = 0; i < vertCount; i++) {
			if (flags[i] & Left) { X[i] = O[i]; w[i] = 0.0f; }
		}
	}
	if (settings.flags & Settings_LockRight) {
		for (uint32_t i = 0; i < vertCount; i++) {
			// Sometimes we animate the right vert to demonstrate various things
			if (flags[i] & Right) {
				X[i] = O[i] = dvec2(origin + (settings.lockedRightTransform * vec2(X0[i])));
				w[i] = 0.0f;
			}
		}
	}

	// Potentially apply the manipulation constraint
	if (manip.pickedGeo == this) {
		uint32_t i = manip.pickedPointIdx;
		vec2 target = manip.pos.xy + manip.pickDirTarget.xy * (-manip.pos.z / manip.pickDirTarget.z);
		X[i] += dvec2((target - vec2(X[i])) * (w[i] / (max(0.000001f, w[i]) + 0.1f / (dt * dt))));
	}

	// Update velocity
	for (uint32_t i = 0; i < vertCount; i++) {
		V[i] = (X[i] - O[i]) * (1.0f / dt);
	}

	// Apply damping
	uint32_t rayleighDampingType = (settings.flags >> Settings_RayleighTypeBit) & Settings_RayleighTypeMask;
	if (rayleighDampingType == Rayleigh_PostAmortized) {
		Settings amortizedSettings = settings;
		amortizedSettings.damping *= (float)AmortizationPeriod;
		amortizedSettings.areaAndTimeCorrectedPbdDamping = settings.amortizedAreaAndTimeCorrectedPbdDamping;
		Damp(amortizedSettings, dt);
	} else {
		Damp(settings, dt);
	}
}

void Geo2d::Transform(mat3 transform) {
	for (uint32_t i = 0; i < vertCount; i++) {
		O[i] = X[i] = dvec2((transform * vec3(vec2(X[i]), 1.0f)).xy);
	}
	origin = (transform * vec3(origin, 1.0f)).xy;
}

void Geo2d::Pick(vec3 rayOrigin, vec3 rayDir, vec3* outNearestPoint, uint32_t* outNearestPointIdx, float* outDistance) {
	vec2 pos = rayOrigin.xy + rayDir.xy * (-rayOrigin.z / rayDir.z);
	*outNearestPointIdx = 0;
	*outDistance = FLT_MAX;
	for (uint32_t i = 0; i < vertCount; i++) {
		if (w[i] == 0.0f) { continue; }
		float dist = distance(vec2(X[i]), pos);
		if (*outDistance > dist) {
			*outDistance = dist;
			*outNearestPoint = vec3(vec2(X[i]), 0.0f);
			*outNearestPointIdx = i;
		}
	}
}

//-----------------------------------------------------------------------------
// Geo3d methods
void Geo3d::Substep(const Settings& settings, const Manipulator& manip, float dt) {
	// Predict
	for (uint32_t i = 0; i < vertCount; i++) {
		V[i] += dvec3(settings.gravity.x * dt, settings.gravity.y * dt, 0.0f);
		V[i] *= 1.0f - settings.timeCorrectedDrag;
		O[i] = X[i];
		X[i] += V[i] * dt;
	}

	// Apply constraints
	Constrain(settings, dt);

	// Undo any motion on locked verts
	if (settings.flags & Settings_LockLeft) {
		for (uint32_t i = 0; i < vertCount; i++) {
			if (flags[i] & Left) { X[i] = O[i]; w[i] = 0.0f; }
		}
	}
	if (settings.flags & Settings_LockRight) {
		for (uint32_t i = 0; i < vertCount; i++) {
			// Sometimes we animate the right vert to demonstrate various things
			if (flags[i] & Right) {
				X[i] = O[i] = dvec3(origin + (settings.lockedRightTransform3d * vec3(X0[i])));
				w[i] = 0.0f;
			}
		}
	}

	// Potentially apply the manipulation constraint
	if (manip.pickedGeo == this) {
		uint32_t i = manip.pickedPointIdx;
		float t = dot(manip.manipPlaneNormal, manip.pick0 - manip.pos) / dot(manip.manipPlaneNormal, manip.pickDirTarget);
		vec3 target = manip.pos + t * manip.pickDirTarget; // target point is on the plane defined by the original pick location and camera normal
		X[i] += dvec3((target - vec3(X[i])) * (w[i] / (max(0.000001f, w[i]) + 1.8f / (dt * dt))));
	}

	// Update velocity
	for (uint32_t i = 0; i < vertCount; i++) {
		V[i] = (X[i] - O[i]) * (1.0f / dt);
	}

	// Apply damping
	uint32_t rayleighDampingType = (settings.flags >> Settings_RayleighTypeBit) & Settings_RayleighTypeMask;
	if (rayleighDampingType == Rayleigh_PostAmortized) {
		Settings amortizedSettings = settings;
		amortizedSettings.damping *= (float)AmortizationPeriod;
		amortizedSettings.volumeAndTimeCorrectedPbdDamping = settings.amortizedVolumeAndTimeCorrectedPbdDamping;
		Damp(amortizedSettings, dt);
	} else {
		Damp(settings, dt);
	}
}

void Geo3d::Transform(mat3 t3) {
	mat4 t4 = mat4(vec4(t3[0], 0.0f), vec4(t3[1], 0.0f), vec4(0.0f, 0.0f, 1.0f, 0.0f), vec4(t3[2].xy, 0.0f, 1.0f));
	for (uint32_t i = 0; i < vertCount; i++) {
		O[i] = X[i] = dvec3((t4 * vec4(vec3(X[i]), 1.0f)).xyz);
	}
	origin = (t4 * vec4(origin, 1.0f)).xyz;
}

void Geo3d::Pick(vec3 rayOrigin, vec3 rayDir, vec3* outNearestPoint, uint32_t* outNearestPointIdx, float* outDistance) {
	auto TestPoint = [&](vec3 p, float* dist, float* weightedDist) {
		float t = dot(rayDir, p - rayOrigin);
		if (t < 0.0f) { return; }
		*dist = distance(p, rayOrigin + t * rayDir);
		*weightedDist = t * *dist;
	};
	float weightedDist = FLT_MAX;
	for (uint32_t i = 0; i < vertCount; i++) {
		if (!(flags[i] & Pickable) || w[i] == 0.0f) { continue; }
		float iDist = FLT_MAX, iWeightedDist = FLT_MAX;
		TestPoint(vec3(X[i]), &iDist, &iWeightedDist);
		if (weightedDist >= iWeightedDist) {
			weightedDist = iWeightedDist;
			*outDistance = iDist;
			*outNearestPoint = vec3(X[i]);
			*outNearestPointIdx = i;
		}
	}
}

//-----------------------------------------------------------------------------
// GeoLinear2d methods
void GeoLinear2d::Init(Allocator* alloc, float density, const float* nodeData, uint32_t nodeDataCount, const uint32_t* idxData, uint32_t idxDataCount, bool autoResize) {
	// Build the connectivity information
	con.Init(alloc, idxData, idxDataCount);
	assert(con.vertCount == nodeDataCount / 2);

	// Allocate and initialize verts
	vertCount = con.vertCount;
	X0 = alloc->Alloc<dvec2>(vertCount);
	X = alloc->Alloc<dvec2>(vertCount);
	O = alloc->Alloc<dvec2>(vertCount);
	V = alloc->Alloc<dvec2>(vertCount);
	w = alloc->Alloc<float>(vertCount);
	flags = alloc->Alloc<uint8_t>(vertCount);

	vec2 aabbMin = vec2(1.0e24f), aabbMax(-1.0e24f);
	for (uint32_t i = 0; i < vertCount; i++) {
		O[i] = X[i] = X0[i] = dvec2(nodeData[2 * i + 0], nodeData[2 * i + 1]);
		V[i] = dvec2(0.0);
		w[i] = 0.0f; // We will accumulate un-inverted mass in here as we construct elements
		flags[i] = Geo::Pickable;
		aabbMin = min(aabbMin, vec2(X[i]));
		aabbMax = max(aabbMax, vec2(X[i]));
	}
	// Determine "left" and "right" verts for locking, and potentially resize
	vec2 aabbDims = aabbMax - aabbMin;
	double maxDim = max(aabbDims.x, aabbDims.y);
	for (uint32_t i = 0; i < vertCount; i++) {
		if (((float)X[i].x - aabbMin.x) / (aabbMax.x - aabbMin.x) < 0.02f) { flags[i] |= Geo::Left; }
		if (((float)X[i].x - aabbMin.x) / (aabbMax.x - aabbMin.x) > 0.98f) { flags[i] |= Geo::Right; }
		if (autoResize) {
			O[i] = X[i] = X0[i] = 0.085 * ((X[i] - 0.5 * dvec2(aabbMin + aabbMax)) / maxDim);
		}
	}

	// Allocate the elements
	t = alloc->Alloc<T3>(con.triCount);
	for (uint32_t i = 0; i < con.triCount; i++) {
		for (uint32_t j = 0; j < 3; j++) { t[i].i[j] = con.tris[i].v[j]; }
		InitFiniteElement(X, t[i].i, density, w, &t[i].ep);
	}
	q = alloc->Alloc<Q4>(con.quadCount);
	for (uint32_t i = 0; i < con.quadCount; i++) {
		for (uint32_t j = 0; j < 4; j++) { q[i].i[j] = con.quads[i].v[j]; }
		InitFiniteElement(X, q[i].i, density, w, &q[i].ep);
	}

	// Flip masses into inverse masses
	for (uint32_t i = 0; i < vertCount; i++) { w[i] = 1.0f / w[i]; }

	// Allocate and init the one ring connectivity info for each vert
	rings = alloc->Alloc<OneRing*>(con.vertCount);
	for (uint32_t i = 0; i < con.vertCount; i++) {
		OneRingConstructor ctor(i);
		const ConVert2d& v = con.verts[i];
		const uint32_t* triIdxs = con.triIdxs + v.trisOffset;
		const uint32_t* quadIdxs = con.quadIdxs + v.quadsOffset;
		for (uint32_t j = 0; j < v.triCount; j++) { ctor.AddPoly(triIdxs[j], t[triIdxs[j]].i, 3); }
		for (uint32_t j = 0; j < v.quadCount; j++) { ctor.AddPoly(quadIdxs[j], q[quadIdxs[j]].i, 4); }
		rings[i] = ctor.Alloc(alloc);
		assert(rings[i]->vertCount <= kMaxRingNodes);
	}

	// Generate randomized orderings for the constraints
	tOrder = alloc->Alloc<uint32_t>(con.triCount);
	qOrder = alloc->Alloc<uint32_t>(con.quadCount);
	vOrder = alloc->Alloc<uint32_t>(con.vertCount);
	for (uint32_t i = 0; i < con.triCount; i++) { tOrder[i] = i; }
	for (uint32_t i = 0; i < con.quadCount; i++) { qOrder[i] = i; }
	for (uint32_t i = 0; i < con.vertCount; i++) { vOrder[i] = i; }
	uint32_t randState = 1;
	auto Rand = [&]() { return randState = ((uint64_t)randState * 48271) % 0x7fffffff; };
	auto Swap = [&](uint32_t& a, uint32_t& b) { uint32_t c = a; a = b; b = c; };
	for (uint32_t i = 0; i < con.triCount; i++) { Swap(tOrder[i], tOrder[Rand() % con.triCount]); }
	for (uint32_t i = 0; i < con.quadCount; i++) { Swap(qOrder[i], qOrder[Rand() % con.quadCount]); }
	for (uint32_t i = 0; i < con.vertCount; i++) { Swap(vOrder[i], vOrder[Rand() % con.vertCount]); }
}

void GeoLinear2d::Constrain(const Settings& settings, float dt) {
	uint32_t energyFunction = (settings.flags >> Settings_EnergyBit) & Settings_EnergyMask;
	if (energyFunction < Energy_ContinuousPixar || energyFunction > Energy_ContinuousSkin) {
		for (uint32_t i = 0; i < con.triCount; i++) { SolveElement(dt, X, O, w, t[tOrder[i]].i, t[tOrder[i]].ep, settings); }
		for (uint32_t i = 0; i < con.quadCount; i++) { SolveElement(dt, X, O, w, q[qOrder[i]].i, q[qOrder[i]].ep, settings); }
		for (uint32_t itr = 0; itr < settings.volumePasses; itr++) {
			for (uint32_t i = 0; i < con.triCount; i++) { SolveElementVolume(dt, X, O, w, t[tOrder[i]].i, t[tOrder[i]].ep, settings); }
			for (uint32_t i = 0; i < con.quadCount; i++) { SolveElementVolume(dt, X, O, w, q[qOrder[i]].i, q[qOrder[i]].ep, settings); }
		}
	} else {
		for (uint32_t i = 0; i < con.vertCount; i++) { SolveRing(dt, X, O, 0, w, con, t, q, *rings[vOrder[i]], settings); }
	}
}

void GeoLinear2d::Damp(const Settings& settings, float dt) {
	uint32_t energyFunction = (settings.flags >> Settings_EnergyBit) & Settings_EnergyMask;
	uint32_t rayleighDampingType = (settings.flags >> Settings_RayleighTypeBit) & Settings_RayleighTypeMask;
	bool isAmortized = rayleighDampingType == Rayleigh_PostAmortized;
	auto Begin = [&](uint32_t count) { return isAmortized ? count * (settings.tickId % AmortizationPeriod) / AmortizationPeriod : 0; };
	auto End = [&](uint32_t count) { return isAmortized ? count * ((settings.tickId % AmortizationPeriod) + 1) / AmortizationPeriod : count; };
	uint32_t triBegin = Begin(con.triCount), triEnd = End(con.triCount);
	uint32_t quadBegin = Begin(con.quadCount), quadEnd = End(con.quadCount);
	uint32_t vertBegin = Begin(con.vertCount), vertEnd = End(con.vertCount);
	if (rayleighDampingType >= Rayleigh_Post) {
		if (energyFunction < Energy_ContinuousPixar || energyFunction > Energy_ContinuousSkin) {
			for (uint32_t i = triBegin; i < triEnd; i++) { DampElement(dt, X, V, w, t[tOrder[i]].i, t[tOrder[i]].ep, settings); }
			for (uint32_t i = quadBegin; i < quadEnd; i++) { DampElement(dt, X, V, w, q[qOrder[i]].i, q[qOrder[i]].ep, settings); }
		} else {
			for (uint32_t i = vertBegin; i < vertEnd; i++) { SolveRing(dt, X, 0, V, w, con, t, q, *rings[vOrder[i]], settings); }
		}
	}
	if (settings.pbdDamping > 0.0f) {
		for (uint32_t i = triBegin; i < triEnd; i++) { PbdDamp(X, V, w, t[tOrder[i]].i, min(1.0f, settings.areaAndTimeCorrectedPbdDamping / t[tOrder[i]].ep.volume)); }
		for (uint32_t i = quadBegin; i < quadEnd; i++) { PbdDamp(X, V, w, q[qOrder[i]].i, min(1.0f, settings.areaAndTimeCorrectedPbdDamping / q[qOrder[i]].ep.volume)); }
	}
}

void GeoLinear2d::Render(const Settings& settings) {
	for (uint32_t i = 0; i < con.edgeCount; i++) {
		const ConEdge2d& e = con.edges[i];
		uint8_t f = flags[e.v[0]] & flags[e.v[1]];
		vec3 c = color;
		if ((settings.flags & Settings_LockLeft) && (f & Geo::Left)) { c = vec3(0.7f, 0.0f, 0.0f); }
		if ((settings.flags & Settings_LockRight) && (f & Geo::Right)) { c = vec3(0.7f, 0.0f, 0.0f); }
		DLINE(vec3(vec2(X[e.v[0]]), 0.0f), vec3(vec2(X[e.v[1]]), 0.0f), c);
	}
}

float GeoLinear2d::CalculateVolume() const {
	float volume = 0.0f;
	for (uint32_t i = 0; i < con.triCount; i++) { volume += CalculateElementVolume(X, t[i].i, t[i].ep); }
	for (uint32_t i = 0; i < con.quadCount; i++) { volume += CalculateElementVolume(X, q[i].i, q[i].ep); }
	return volume;
}

//-----------------------------------------------------------------------------
// GeoQuadratic2d methods
void GeoQuadratic2d::Init(Allocator* alloc, float density, const float* nodeData, uint32_t nodeDataCount, const uint32_t* idxData, uint32_t idxDataCount, bool autoResize) {
	// Build the connectivity information
	con.Init(alloc, idxData, idxDataCount);
	assert(con.vertCount == nodeDataCount / 2);

	// Allocate and initialize verts
	vertCount = con.vertCount + con.edgeCount + con.quadCount;
	X0 = alloc->Alloc<dvec2>(vertCount);
	X = alloc->Alloc<dvec2>(vertCount);
	O = alloc->Alloc<dvec2>(vertCount);
	V = alloc->Alloc<dvec2>(vertCount);
	w = alloc->Alloc<float>(vertCount);
	flags = alloc->Alloc<uint8_t>(vertCount);

	// Init the "linear" vert positions
	for (uint32_t i = 0; i < con.vertCount; i++) {
		O[i] = X[i] = X0[i] = dvec2(nodeData[2 * i + 0], nodeData[2 * i + 1]);
	}
	// Position the middle verts
	uint32_t vi = con.vertCount;
	for (uint32_t i = 0; i < con.edgeCount; i++) {
		const ConEdge2d& e = con.edges[i];
		O[vi] = X[vi] = X0[vi] = 0.5 * (X[e.v[0]] + X[e.v[1]]);
		++vi;
	}
	for (uint32_t i = 0; i < con.quadCount; i++) {
		const ConQuad2d& q = con.quads[i];
		O[vi] = X[vi] = X0[vi] = 0.25 * (X[q.v[0]] + X[q.v[1]] + X[q.v[2]] + X[q.v[3]]);
		++vi;
	}
	assert(vi == vertCount);
	// Fill in other vert attributes
	vec2 aabbMin = vec2(1.0e24f), aabbMax(-1.0e24f);
	for (uint32_t i = 0; i < vertCount; i++) {
		aabbMin = min(aabbMin, vec2(X[i]));
		aabbMax = max(aabbMax, vec2(X[i]));
		V[i] = dvec2(0.0);
		w[i] = 0.0f; // We will accumulate un-inverted mass in here as we construct elements
		flags[i] = Geo::Pickable;
	}
	// Determine "left" and "right" verts for locking, and potentially resize
	vec2 aabbDims = aabbMax - aabbMin;
	double maxDim = max(aabbDims.x, aabbDims.y);
	for (uint32_t i = 0; i < vertCount; i++) {
		if (((float)X[i].x - aabbMin.x) / (aabbMax.x - aabbMin.x) < 0.02f) { flags[i] |= Geo::Left; }
		if (((float)X[i].x - aabbMin.x) / (aabbMax.x - aabbMin.x) > 0.98f) { flags[i] |= Geo::Right; }
		if (autoResize) {
			O[i] = X[i] = X0[i] = 0.085 * ((X[i] - 0.5 * dvec2(aabbMin + aabbMax)) / maxDim);
		}
	}

	// Allocate the elements
	t = alloc->Alloc<T6>(con.triCount);
	for (uint32_t i = 0; i < con.triCount; i++) {
		for (uint32_t j = 0; j < 3; j++) { t[i].i[j] = con.tris[i].v[j]; }
		t[i].i[3] = con.vertCount + con.EdgeFromVerts(t[i].i[0], t[i].i[1]);
		t[i].i[4] = con.vertCount + con.EdgeFromVerts(t[i].i[1], t[i].i[2]);
		t[i].i[5] = con.vertCount + con.EdgeFromVerts(t[i].i[2], t[i].i[0]);
		InitFiniteElement(X, t[i].i, density, w, &t[i].ep);
	}
	q = alloc->Alloc<Q9>(con.quadCount);
	for (uint32_t i = 0; i < con.quadCount; i++) {
		for (uint32_t j = 0; j < 4; j++) { q[i].i[j] = con.quads[i].v[j]; }
		q[i].i[4] = con.vertCount + con.EdgeFromVerts(q[i].i[0], q[i].i[1]);
		q[i].i[5] = con.vertCount + con.EdgeFromVerts(q[i].i[1], q[i].i[2]);
		q[i].i[6] = con.vertCount + con.EdgeFromVerts(q[i].i[2], q[i].i[3]);
		q[i].i[7] = con.vertCount + con.EdgeFromVerts(q[i].i[3], q[i].i[0]);
		q[i].i[8] = con.vertCount + con.edgeCount + i;
		InitFiniteElement(X, q[i].i, density, w, &q[i].ep);
	}

	// Flip masses into inverse masses
	for (uint32_t i = 0; i < vertCount; i++) { w[i] = 1.0f / w[i]; }

	// Allocate and init the one ring connectivity info for each vert
	rings = alloc->Alloc<OneRing*>(con.vertCount + con.edgeCount + con.quadCount);
	for (uint32_t i = 0; i < con.vertCount; i++) {
		OneRingConstructor ctor(i);
		const ConVert2d& v = con.verts[i];
		const uint32_t* triIdxs = con.triIdxs + v.trisOffset;
		const uint32_t* quadIdxs = con.quadIdxs + v.quadsOffset;
		for (uint32_t j = 0; j < v.triCount; j++) { ctor.AddPoly(triIdxs[j], t[triIdxs[j]].i, 6); }
		for (uint32_t j = 0; j < v.quadCount; j++) { ctor.AddPoly(quadIdxs[j], q[quadIdxs[j]].i, 9); }
		rings[i] = ctor.Alloc(alloc);
		assert(rings[i]->vertCount <= kMaxRingNodes);
	}

	// Generate randomized orderings for the constraints
	tOrder = alloc->Alloc<uint32_t>(con.triCount);
	qOrder = alloc->Alloc<uint32_t>(con.quadCount);
	vOrder = alloc->Alloc<uint32_t>(con.vertCount);
	for (uint32_t i = 0; i < con.triCount; i++) { tOrder[i] = i; }
	for (uint32_t i = 0; i < con.quadCount; i++) { qOrder[i] = i; }
	for (uint32_t i = 0; i < con.vertCount; i++) { vOrder[i] = i; }
	uint32_t randState = 1;
	auto Rand = [&]() { return randState = ((uint64_t)randState * 48271) % 0x7fffffff; };
	auto Swap = [&](uint32_t& a, uint32_t& b) { uint32_t c = a; a = b; b = c; };
	for (uint32_t i = 0; i < con.triCount; i++) { Swap(tOrder[i], tOrder[Rand() % con.triCount]); }
	for (uint32_t i = 0; i < con.quadCount; i++) { Swap(qOrder[i], qOrder[Rand() % con.quadCount]); }
	for (uint32_t i = 0; i < con.vertCount; i++) { Swap(vOrder[i], vOrder[Rand() % con.vertCount]); }
}

void GeoQuadratic2d::Constrain(const Settings& settings, float dt) {
	uint32_t energyFunction = (settings.flags >> Settings_EnergyBit) & Settings_EnergyMask;
	if (energyFunction < Energy_ContinuousPixar || energyFunction > Energy_ContinuousSkin) {
		for (uint32_t i = 0; i < con.triCount; i++) { SolveElement(dt, X, O, w, t[tOrder[i]].i, t[tOrder[i]].ep, settings); }
		for (uint32_t i = 0; i < con.quadCount; i++) { SolveElement(dt, X, O, w, q[qOrder[i]].i, q[qOrder[i]].ep, settings); }
		for (uint32_t itr = 0; itr < settings.volumePasses; itr++) {
			for (uint32_t i = 0; i < con.triCount; i++) { SolveElementVolume(dt, X, O, w, t[tOrder[i]].i, t[tOrder[i]].ep, settings); }
			for (uint32_t i = 0; i < con.quadCount; i++) { SolveElementVolume(dt, X, O, w, q[qOrder[i]].i, q[qOrder[i]].ep, settings); }
		}
	} else {
		for (uint32_t i = 0; i < con.vertCount; i++) { SolveRing(dt, X, O, 0, w, con, t, q, *rings[vOrder[i]], settings); }
	}
}

void GeoQuadratic2d::Damp(const Settings& settings, float dt) {
	uint32_t energyFunction = (settings.flags >> Settings_EnergyBit) & Settings_EnergyMask;
	uint32_t rayleighDampingType = (settings.flags >> Settings_RayleighTypeBit) & Settings_RayleighTypeMask;
	bool isAmortized = rayleighDampingType == Rayleigh_PostAmortized;
	auto Begin = [&](uint32_t count) { return isAmortized ? count * (settings.tickId % AmortizationPeriod) / AmortizationPeriod : 0; };
	auto End = [&](uint32_t count) { return isAmortized ? count * ((settings.tickId % AmortizationPeriod) + 1) / AmortizationPeriod : count; };
	uint32_t triBegin = Begin(con.triCount), triEnd = End(con.triCount);
	uint32_t quadBegin = Begin(con.quadCount), quadEnd = End(con.quadCount);
	uint32_t vertBegin = Begin(con.vertCount), vertEnd = End(con.vertCount);
	if (rayleighDampingType >= Rayleigh_Post) {
		if (energyFunction < Energy_ContinuousPixar || energyFunction > Energy_ContinuousSkin) {
			for (uint32_t i = triBegin; i < triEnd; i++) { DampElement(dt, X, V, w, t[tOrder[i]].i, t[tOrder[i]].ep, settings); }
			for (uint32_t i = quadBegin; i < quadEnd; i++) { DampElement(dt, X, V, w, q[qOrder[i]].i, q[qOrder[i]].ep, settings); }
		} else {
			for (uint32_t i = vertBegin; i < vertEnd; i++) { SolveRing(dt, X, 0, V, w, con, t, q, *rings[vOrder[i]], settings); }
		}
	}
	if (settings.pbdDamping > 0.0f) {
		for (uint32_t i = triBegin; i < triEnd; i++) { PbdDamp(X, V, w, t[tOrder[i]].i, min(1.0f, settings.areaAndTimeCorrectedPbdDamping / t[tOrder[i]].ep.volume)); }
		for (uint32_t i = quadBegin; i < quadEnd; i++) { PbdDamp(X, V, w, q[qOrder[i]].i, min(1.0f, settings.areaAndTimeCorrectedPbdDamping / q[qOrder[i]].ep.volume)); }
	}
}

void GeoQuadratic2d::Render(const Settings& settings) {
	auto DQLINE = [&](vec2 Xi, vec2 Xj, vec2 Xk, const vec3& color) {
		const uint32_t kCount = 7;
		vec3 P[kCount];
		for (uint32_t x = 0; x < kCount; x++) {
			float t = (float)x * (1.0f / (float)(kCount - 1)) * 2.0f - 1.0f;
			float Ni = (t * (t - 1.0f)) * 0.5f, Nj = -(t + 1.0f) * (t - 1.0f), Nk = ((t + 1.0f) * t) * 0.5f;
			P[x] = vec3(Ni * Xi + Nj * Xj + Nk * Xk, 0.0f);
		}
		DLINES(P, kCount, color);
	};

	for (uint32_t i = 0; i < con.edgeCount; i++) {
		const ConEdge2d& e = con.edges[i];
		uint8_t f = flags[e.v[0]] & flags[e.v[1]];
		vec3 c = color;
		if ((settings.flags & Settings_LockLeft) && (f & Geo::Left)) { c = vec3(0.7f, 0.0f, 0.0f); }
		if ((settings.flags & Settings_LockRight) && (f & Geo::Right)) { c = vec3(0.7f, 0.0f, 0.0f); }
		DQLINE(vec2(X[e.v[0]]), vec2(X[con.vertCount + i]), vec2(X[e.v[1]]), c);
	}

	for (uint32_t i = 0; i < con.quadCount; i++) {
		DPOINT(vec3(vec2(X[con.vertCount + con.edgeCount + i]), 0.0f), color);
	}
}

float GeoQuadratic2d::CalculateVolume() const {
	float volume = 0.0f;
	for (uint32_t i = 0; i < con.triCount; i++) { volume += CalculateElementVolume(X, t[i].i, t[i].ep); }
	for (uint32_t i = 0; i < con.quadCount; i++) { volume += CalculateElementVolume(X, q[i].i, q[i].ep); }
	return volume;
}

//-----------------------------------------------------------------------------
// GeoLinear3d methods
void GeoLinear3d::Init(Allocator* alloc, float density, const float* nodeData, uint32_t nodeDataCount, const uint32_t* idxData, uint32_t idxDataCount, bool autoResize) {
	// Build the connectivity information
	con.Init(alloc, idxData, idxDataCount);
	assert(con.vertCount == nodeDataCount / 3);

	// Allocate and initialize verts
	vertCount = con.vertCount;
	X0 = alloc->Alloc<dvec3>(vertCount);
	X = alloc->Alloc<dvec3>(vertCount);
	O = alloc->Alloc<dvec3>(vertCount);
	V = alloc->Alloc<dvec3>(vertCount);
	w = alloc->Alloc<float>(vertCount);
	flags = alloc->Alloc<uint8_t>(vertCount);

	vec3 aabbMin = vec3(1.0e24f), aabbMax(-1.0e24f);
	for (uint32_t i = 0; i < vertCount; i++) {
		O[i] = X[i] = X0[i] = dvec3(nodeData[3 * i + 0], nodeData[3 * i + 1], nodeData[3 * i + 2]);
		V[i] = dvec3(0.0);
		w[i] = 0.0f; // We will accumulate un-inverted mass in here as we construct elements
		flags[i] = Geo::Pickable;
		aabbMin = min(aabbMin, vec3(X[i]));
		aabbMax = max(aabbMax, vec3(X[i]));
	}
	// Determine "left" and "right" verts for locking, and potentially resize
	vec3 aabbDims = aabbMax - aabbMin;
	double maxDim = max(aabbDims.x, max(aabbDims.y, aabbDims.z));
	for (uint32_t i = 0; i < vertCount; i++) {
		if (((float)X[i].x - aabbMin.x) / (aabbMax.x - aabbMin.x) < 0.02f) { flags[i] |= Geo::Left; }
		if (((float)X[i].x - aabbMin.x) / (aabbMax.x - aabbMin.x) > 0.98f) { flags[i] |= Geo::Right; }
		if (autoResize) {
			O[i] = X[i] = X0[i] = 0.085 * ((X[i] - 0.5 * dvec3(aabbMin + aabbMax)) / maxDim);
		}
	}

	// Allocate the elements
	t = alloc->Alloc<T4>(con.tetCount);
	for (uint32_t i = 0; i < con.tetCount; i++) {
		for (uint32_t j = 0; j < 4; j++) { t[i].i[j] = con.tets[i].v[j]; }
		InitFiniteElement(X, t[i].i, density, w, &t[i].ep);
	}
	h = alloc->Alloc<H8>(con.hexCount);
	for (uint32_t i = 0; i < con.hexCount; i++) {
		for (uint32_t j = 0; j < 8; j++) { h[i].i[j] = con.hexs[i].v[j]; }
		InitFiniteElement(X, h[i].i, density, w, &h[i].ep);
	}

	// Flip masses into inverse masses
	for (uint32_t i = 0; i < vertCount; i++) { w[i] = 1.0f / w[i]; }

	// Allocate and init the one ring connectivity info for each vert
	rings = alloc->Alloc<OneRing*>(con.vertCount);
	for (uint32_t i = 0; i < con.vertCount; i++) {
		OneRingConstructor ctor(i);
		const ConVert3d& v = con.verts[i];
		const uint32_t* tetIdxs = con.tetIdxs + v.tetsOffset;
		const uint32_t* hexIdxs = con.hexIdxs + v.hexsOffset;
		for (uint32_t j = 0; j < v.tetCount; j++) { ctor.AddPoly(tetIdxs[j], t[tetIdxs[j]].i, 4); }
		for (uint32_t j = 0; j < v.hexCount; j++) { ctor.AddPoly(hexIdxs[j], h[hexIdxs[j]].i, 8); }
		rings[i] = ctor.Alloc(alloc);
		assert(rings[i]->vertCount <= kMaxRingNodes);
	}

	// Generate randomized orderings for the constraints
	tOrder = alloc->Alloc<uint32_t>(con.tetCount);
	hOrder = alloc->Alloc<uint32_t>(con.hexCount);
	vOrder = alloc->Alloc<uint32_t>(con.vertCount);
	for (uint32_t i = 0; i < con.tetCount; i++) { tOrder[i] = i; }
	for (uint32_t i = 0; i < con.hexCount; i++) { hOrder[i] = i; }
	for (uint32_t i = 0; i < con.vertCount; i++) { vOrder[i] = i; }
	uint32_t randState = 1;
	auto Rand = [&]() { return randState = ((uint64_t)randState * 48271) % 0x7fffffff; };
	auto Swap = [&](uint32_t& a, uint32_t& b) { uint32_t c = a; a = b; b = c; };
	for (uint32_t i = 0; i < con.tetCount; i++) { Swap(tOrder[i], tOrder[Rand() % con.tetCount]); }
	for (uint32_t i = 0; i < con.hexCount; i++) { Swap(hOrder[i], hOrder[Rand() % con.hexCount]); }
	for (uint32_t i = 0; i < con.vertCount; i++) { Swap(vOrder[i], vOrder[Rand() % (con.vertCount)]); }
}

void GeoLinear3d::Constrain(const Settings& settings, float dt) {
	uint32_t energyFunction = (settings.flags >> Settings_EnergyBit) & Settings_EnergyMask;
	if (energyFunction < Energy_ContinuousPixar || energyFunction > Energy_ContinuousSkin) {
		for (uint32_t i = 0; i < con.tetCount; i++) { SolveElement(dt, X, O, w, t[tOrder[i]].i, t[tOrder[i]].ep, settings); }
		for (uint32_t i = 0; i < con.hexCount; i++) { SolveElement(dt, X, O, w, h[hOrder[i]].i, h[hOrder[i]].ep, settings); }
		for (uint32_t itr = 0; itr < settings.volumePasses; itr++) {
			for (uint32_t i = 0; i < con.tetCount; i++) { SolveElementVolume(dt, X, O, w, t[tOrder[i]].i, t[tOrder[i]].ep, settings); }
			for (uint32_t i = 0; i < con.hexCount; i++) { SolveElementVolume(dt, X, O, w, h[hOrder[i]].i, h[hOrder[i]].ep, settings); }
		}
	} else {
		for (uint32_t i = 0; i < con.vertCount; i++) {
			SolveRing(dt, X, O, 0, w, con, t, h, *rings[vOrder[i]], settings);
		}
	}
}

void GeoLinear3d::Damp(const Settings& settings, float dt) {
	uint32_t energyFunction = (settings.flags >> Settings_EnergyBit) & Settings_EnergyMask;
	uint32_t rayleighDampingType = (settings.flags >> Settings_RayleighTypeBit) & Settings_RayleighTypeMask;
	bool isAmortized = rayleighDampingType == Rayleigh_PostAmortized;
	auto Begin = [&](uint32_t count) { return isAmortized ? count * (settings.tickId % AmortizationPeriod) / AmortizationPeriod : 0; };
	auto End = [&](uint32_t count) { return isAmortized ? count * ((settings.tickId % AmortizationPeriod) + 1) / AmortizationPeriod : count; };
	uint32_t tetBegin = Begin(con.tetCount), tetEnd = End(con.tetCount);
	uint32_t hexBegin = Begin(con.hexCount), hexEnd = End(con.hexCount);
	uint32_t vertBegin = Begin(con.vertCount), vertEnd = End(con.vertCount);
	if (rayleighDampingType >= Rayleigh_Post) {
		if (energyFunction < Energy_ContinuousPixar || energyFunction > Energy_ContinuousSkin) {
			for (uint32_t i = tetBegin; i < tetEnd; i++) { DampElement(dt, X, V, w, t[tOrder[i]].i, t[tOrder[i]].ep, settings); }
			for (uint32_t i = hexBegin; i < hexEnd; i++) { DampElement(dt, X, V, w, h[hOrder[i]].i, h[hOrder[i]].ep, settings); }
		} else {
			for (uint32_t i = vertBegin; i < vertEnd; i++) { SolveRing(dt, X, 0, V, w, con, t, h, *rings[vOrder[i]], settings); }
		}
	}
	if (settings.pbdDamping > 0.0f) {
		for (uint32_t i = tetBegin; i < tetEnd; i++) { PbdDamp(X, V, w, t[tOrder[i]].i, min(1.0f, settings.volumeAndTimeCorrectedPbdDamping / t[tOrder[i]].ep.surfaceArea)); }
		for (uint32_t i = hexBegin; i < hexEnd; i++) { PbdDamp(X, V, w, h[hOrder[i]].i, min(1.0f, settings.volumeAndTimeCorrectedPbdDamping / h[hOrder[i]].ep.surfaceArea)); }
	}
}

void GeoLinear3d::Render(const Settings& settings) {
	bool tooBig = con.edgeCount > 10000;
	for (uint32_t i = 0; i < con.edgeCount; i++) {
		const ConEdge3d& e = con.edges[i];
		bool isSurface = con.surfaceEdges[i];
		if (!isSurface && tooBig) { continue; }
		uint8_t f = flags[e.v[0]] & flags[e.v[1]];
		vec3 c = isSurface ? color : mix(color, vec3(1.0f), 0.7f);
		if ((settings.flags & Settings_LockLeft) && (f & Geo::Left)) { c = vec3(0.7f, 0.0f, 0.0f); }
		if ((settings.flags & Settings_LockRight) && (f & Geo::Right)) { c = vec3(0.7f, 0.0f, 0.0f); }
		DLINE(vec3(X[e.v[0]]), vec3(X[e.v[1]]), c);
	}
}

float GeoLinear3d::CalculateVolume() const {
	float volume = 0.0f;
	for (uint32_t i = 0; i < con.tetCount; i++) { volume += CalculateElementVolume(X, t[i].i, t[i].ep); }
	for (uint32_t i = 0; i < con.hexCount; i++) { volume += CalculateElementVolume(X, h[i].i, h[i].ep); }
	return volume;
}

//-----------------------------------------------------------------------------
// GeoQuadratic3d methods
void GeoQuadratic3d::Init(Allocator* alloc, float density, const float* nodeData, uint32_t nodeDataCount, const uint32_t* idxData, uint32_t idxDataCount, bool autoResize) {
	// Build the connectivity information
	con.Init(alloc, idxData, idxDataCount);
	assert(con.vertCount == nodeDataCount / 3);

	// Allocate and initialize verts
	vertCount = con.vertCount + con.edgeCount + con.quadCount + con.hexCount;
	X0 = alloc->Alloc<dvec3>(vertCount);
	X = alloc->Alloc<dvec3>(vertCount);
	O = alloc->Alloc<dvec3>(vertCount);
	V = alloc->Alloc<dvec3>(vertCount);
	w = alloc->Alloc<float>(vertCount);
	flags = alloc->Alloc<uint8_t>(vertCount);

	// Init the "linear" vert positions
	for (uint32_t i = 0; i < con.vertCount; i++) {
		O[i] = X[i] = X0[i] = dvec3(nodeData[3 * i + 0], nodeData[3 * i + 1], nodeData[3 * i + 2]);
	}
	// Position the middle verts
	uint32_t vi = con.vertCount;
	for (uint32_t i = 0; i < con.edgeCount; i++) {
		const ConEdge3d& e = con.edges[i];
		O[vi] = X[vi] = X0[vi] = 0.5 * (X[e.v[0]] + X[e.v[1]]);
		++vi;
	}
	for (uint32_t i = 0; i < con.quadCount; i++) {
		const ConQuad3d& q = con.quads[i];
		O[vi] = X[vi] = X0[vi] = 0.25 * (X[q.v[0]] + X[q.v[1]] + X[q.v[2]] + X[q.v[3]]);
		++vi;
	}
	for (uint32_t i = 0; i < con.hexCount; i++) {
		const ConHex3d& h = con.hexs[i];
		O[vi] = X[vi] = X0[vi] = 0.125 * (X[h.v[0]] + X[h.v[1]] + X[h.v[2]] + X[h.v[3]] + X[h.v[4]] + X[h.v[5]] + X[h.v[6]] + X[h.v[7]]);
		++vi;
	}
	assert(vi == vertCount);
	// Fill in other vert attributes
	vec3 aabbMin = vec3(1.0e24f), aabbMax(-1.0e24f);
	for (uint32_t i = 0; i < vertCount; i++) {
		aabbMin = min(aabbMin, vec3(X[i]));
		aabbMax = max(aabbMax, vec3(X[i]));
		V[i] = dvec3(0.0);
		w[i] = 0.0f; // We will accumulate un-inverted mass in here as we construct elements
		flags[i] = Geo::Pickable;
	}
	// Determine "left" and "right" verts for locking, and potentially resize
	vec3 aabbDims = aabbMax - aabbMin;
	double maxDim = max(aabbDims.x, max(aabbDims.y, aabbDims.z));
	for (uint32_t i = 0; i < vertCount; i++) {
		if (((float)X[i].x - aabbMin.x) / (aabbMax.x - aabbMin.x) < 0.02f) { flags[i] |= Geo::Left; }
		if (((float)X[i].x - aabbMin.x) / (aabbMax.x - aabbMin.x) > 0.98f) { flags[i] |= Geo::Right; }
		if (autoResize) {
			O[i] = X[i] = X0[i] = 0.085 * ((X[i] - 0.5 * dvec3(aabbMin + aabbMax)) / maxDim);
		}
	}

	// Allocate the elements
	t = alloc->Alloc<T10>(con.tetCount);
	for (uint32_t i = 0; i < con.tetCount; i++) {
		for (uint32_t j = 0; j < 4; j++) { t[i].i[j] = con.tets[i].v[j]; }
		t[i].i[4] = con.vertCount + con.EdgeFromVerts(t[i].i[0], t[i].i[1]);
		t[i].i[5] = con.vertCount + con.EdgeFromVerts(t[i].i[0], t[i].i[2]);
		t[i].i[6] = con.vertCount + con.EdgeFromVerts(t[i].i[0], t[i].i[3]);
		t[i].i[7] = con.vertCount + con.EdgeFromVerts(t[i].i[1], t[i].i[2]);
		t[i].i[8] = con.vertCount + con.EdgeFromVerts(t[i].i[1], t[i].i[3]);
		t[i].i[9] = con.vertCount + con.EdgeFromVerts(t[i].i[2], t[i].i[3]);
		InitFiniteElement(X, t[i].i, density, w, &t[i].ep);
	}
	h = alloc->Alloc<H27>(con.hexCount);
	for (uint32_t i = 0; i < con.hexCount; i++) {
		for (uint32_t j = 0; j < 8; j++) { h[i].i[j] = con.hexs[i].v[j]; }
		h[i].i[0] = con.hexs[i].v[0];
		h[i].i[1] = con.vertCount + con.EdgeFromVerts(con.hexs[i].v[0], con.hexs[i].v[1]);
		h[i].i[2] = con.hexs[i].v[1];
		h[i].i[3] = con.vertCount + con.EdgeFromVerts(con.hexs[i].v[0], con.hexs[i].v[2]);
		h[i].i[4] = con.vertCount + con.edgeCount + con.QuadFromVerts(con.hexs[i].v[0], con.hexs[i].v[1], con.hexs[i].v[3], con.hexs[i].v[2]);
		h[i].i[5] = con.vertCount + con.EdgeFromVerts(con.hexs[i].v[1], con.hexs[i].v[3]);
		h[i].i[6] = con.hexs[i].v[2];
		h[i].i[7] = con.vertCount + con.EdgeFromVerts(con.hexs[i].v[2], con.hexs[i].v[3]);
		h[i].i[8] = con.hexs[i].v[3];

		h[i].i[9] = con.vertCount + con.EdgeFromVerts(con.hexs[i].v[0], con.hexs[i].v[4]);
		h[i].i[10] = con.vertCount + con.edgeCount + con.QuadFromVerts(con.hexs[i].v[4], con.hexs[i].v[5], con.hexs[i].v[1], con.hexs[i].v[0]);
		h[i].i[11] = con.vertCount + con.EdgeFromVerts(con.hexs[i].v[1], con.hexs[i].v[5]);
		h[i].i[12] = con.vertCount + con.edgeCount + con.QuadFromVerts(con.hexs[i].v[4], con.hexs[i].v[0], con.hexs[i].v[2], con.hexs[i].v[6]);
		h[i].i[13] = con.vertCount + con.edgeCount + con.quadCount + i;
		h[i].i[14] = con.vertCount + con.edgeCount + con.QuadFromVerts(con.hexs[i].v[1], con.hexs[i].v[5], con.hexs[i].v[7], con.hexs[i].v[3]);
		h[i].i[15] = con.vertCount + con.EdgeFromVerts(con.hexs[i].v[2], con.hexs[i].v[6]);
		h[i].i[16] = con.vertCount + con.edgeCount + con.QuadFromVerts(con.hexs[i].v[2], con.hexs[i].v[3], con.hexs[i].v[7], con.hexs[i].v[6]);
		h[i].i[17] = con.vertCount + con.EdgeFromVerts(con.hexs[i].v[3], con.hexs[i].v[7]);

		h[i].i[18] = con.hexs[i].v[4];
		h[i].i[19] = con.vertCount + con.EdgeFromVerts(con.hexs[i].v[4], con.hexs[i].v[5]);
		h[i].i[20] = con.hexs[i].v[5];
		h[i].i[21] = con.vertCount + con.EdgeFromVerts(con.hexs[i].v[4], con.hexs[i].v[6]);
		h[i].i[22] = con.vertCount + con.edgeCount + con.QuadFromVerts(con.hexs[i].v[5], con.hexs[i].v[4], con.hexs[i].v[6], con.hexs[i].v[7]);
		h[i].i[23] = con.vertCount + con.EdgeFromVerts(con.hexs[i].v[5], con.hexs[i].v[7]);
		h[i].i[24] = con.hexs[i].v[6];
		h[i].i[25] = con.vertCount + con.EdgeFromVerts(con.hexs[i].v[6], con.hexs[i].v[7]);
		h[i].i[26] = con.hexs[i].v[7];

		InitFiniteElement(X, h[i].i, density, w, &h[i].ep);
	}

	// Flip masses into inverse masses
	for (uint32_t i = 0; i < vertCount; i++) { w[i] = 1.0f / w[i]; }

	// Allocate and init the one ring connectivity info for each vert
	rings = alloc->Alloc<OneRing*>(con.vertCount + con.edgeCount + con.quadCount);
	for (uint32_t i = 0; i < con.vertCount; i++) {
		OneRingConstructor ctor(i);
		const ConVert3d& v = con.verts[i];
		const uint32_t* tetIdxs = con.tetIdxs + v.tetsOffset;
		const uint32_t* hexIdxs = con.hexIdxs + v.hexsOffset;
		for (uint32_t j = 0; j < v.tetCount; j++) { ctor.AddPoly(tetIdxs[j], t[tetIdxs[j]].i, 10); }
		for (uint32_t j = 0; j < v.hexCount; j++) { ctor.AddPoly(hexIdxs[j], h[hexIdxs[j]].i, 27); }
		rings[i] = ctor.Alloc(alloc);
		assert(rings[i]->vertCount <= kMaxRingNodes);
	}

	// Generate randomized orderings for the constraints
	tOrder = alloc->Alloc<uint32_t>(con.tetCount);
	hOrder = alloc->Alloc<uint32_t>(con.hexCount);
	vOrder = alloc->Alloc<uint32_t>(con.vertCount);
	for (uint32_t i = 0; i < con.tetCount; i++) { tOrder[i] = i; }
	for (uint32_t i = 0; i < con.hexCount; i++) { hOrder[i] = i; }
	for (uint32_t i = 0; i < con.vertCount; i++) { vOrder[i] = i; }
	uint32_t randState = 1;
	auto Rand = [&]() { return randState = ((uint64_t)randState * 48271) % 0x7fffffff; };
	auto Swap = [&](uint32_t& a, uint32_t& b) { uint32_t c = a; a = b; b = c; };
	for (uint32_t i = 0; i < con.tetCount; i++) { Swap(tOrder[i], tOrder[Rand() % con.tetCount]); }
	for (uint32_t i = 0; i < con.hexCount; i++) { Swap(hOrder[i], hOrder[Rand() % con.hexCount]); }
	for (uint32_t i = 0; i < con.vertCount; i++) { Swap(vOrder[i], vOrder[Rand() % (con.vertCount)]); }
}

void GeoQuadratic3d::Constrain(const Settings& settings, float dt) {
	uint32_t energyFunction = (settings.flags >> Settings_EnergyBit) & Settings_EnergyMask;
	if (energyFunction < Energy_ContinuousPixar || energyFunction > Energy_ContinuousSkin) {
		for (uint32_t i = 0; i < con.tetCount; i++) { SolveElement(dt, X, O, w, t[tOrder[i]].i, t[tOrder[i]].ep, settings); }
		for (uint32_t i = 0; i < con.hexCount; i++) { SolveElement(dt, X, O, w, h[hOrder[i]].i, h[hOrder[i]].ep, settings); }
		for (uint32_t itr = 0; itr < settings.volumePasses; itr++) {
			for (uint32_t i = 0; i < con.tetCount; i++) { SolveElementVolume(dt, X, O, w, t[tOrder[i]].i, t[tOrder[i]].ep, settings); }
			for (uint32_t i = 0; i < con.hexCount; i++) { SolveElementVolume(dt, X, O, w, h[hOrder[i]].i, h[hOrder[i]].ep, settings); }
		}
	} else {
		for (uint32_t i = 0; i < con.vertCount; i++) {
			SolveRing(dt, X, O, 0, w, con, t, h, *rings[vOrder[i]], settings);
		}
	}
}

void GeoQuadratic3d::Damp(const Settings& settings, float dt) {
	uint32_t energyFunction = (settings.flags >> Settings_EnergyBit) & Settings_EnergyMask;
	uint32_t rayleighDampingType = (settings.flags >> Settings_RayleighTypeBit) & Settings_RayleighTypeMask;
	bool isAmortized = rayleighDampingType == Rayleigh_PostAmortized;
	auto Begin = [&](uint32_t count) { return isAmortized ? count * (settings.tickId % AmortizationPeriod) / AmortizationPeriod : 0; };
	auto End = [&](uint32_t count) { return isAmortized ? count * ((settings.tickId % AmortizationPeriod) + 1) / AmortizationPeriod : count; };
	uint32_t tetBegin = Begin(con.tetCount), tetEnd = End(con.tetCount);
	uint32_t hexBegin = Begin(con.hexCount), hexEnd = End(con.hexCount);
	uint32_t vertBegin = Begin(con.vertCount), vertEnd = End(con.vertCount);
	if (rayleighDampingType >= Rayleigh_Post) {
		if (energyFunction < Energy_ContinuousPixar || energyFunction > Energy_ContinuousSkin) {
			for (uint32_t i = tetBegin; i < tetEnd; i++) { DampElement(dt, X, V, w, t[tOrder[i]].i, t[tOrder[i]].ep, settings); }
			for (uint32_t i = hexBegin; i < hexEnd; i++) { DampElement(dt, X, V, w, h[hOrder[i]].i, h[hOrder[i]].ep, settings); }
		} else {
			for (uint32_t i = vertBegin; i < vertEnd; i++) { SolveRing(dt, X, 0, V, w, con, t, h, *rings[vOrder[i]], settings); }
		}
	}
	if (settings.pbdDamping > 0.0f) {
		for (uint32_t i = tetBegin; i < tetEnd; i++) { PbdDamp(X, V, w, t[tOrder[i]].i, min(1.0f, settings.volumeAndTimeCorrectedPbdDamping / t[tOrder[i]].ep.surfaceArea)); }
		for (uint32_t i = hexBegin; i < hexEnd; i++) { PbdDamp(X, V, w, h[hOrder[i]].i, min(1.0f, settings.volumeAndTimeCorrectedPbdDamping / h[hOrder[i]].ep.surfaceArea)); }
	}
}

void GeoQuadratic3d::Render(const Settings& settings) {
	auto DQLINE = [&](vec3 Xi, vec3 Xj, vec3 Xk, const vec3& color) {
		const uint32_t kCount = 7;
		vec3 P[kCount];
		for (uint32_t x = 0; x < kCount; x++) {
			float t = (float)x * (1.0f / (float)(kCount - 1)) * 2.0f - 1.0f;
			float Ni = (t * (t - 1.0f)) * 0.5f, Nj = -(t + 1.0f) * (t - 1.0f), Nk = ((t + 1.0f) * t) * 0.5f;
			P[x] = Ni * Xi + Nj * Xj + Nk * Xk;
		}
		DLINES(P, kCount, color);
	};

	bool tooBig = con.edgeCount > 2000;
	for (uint32_t i = 0; i < con.edgeCount; i++) {
		const ConEdge3d& e = con.edges[i];
		bool isSurface = con.surfaceEdges[i];
		if (!isSurface && tooBig) { continue; }
		uint8_t f = flags[e.v[0]] & flags[e.v[1]];
		vec3 c = isSurface ? color : mix(color, vec3(1.0f), 0.7f);
		if ((settings.flags & Settings_LockLeft) && (f & Geo::Left)) { c = vec3(0.7f, 0.0f, 0.0f); }
		if ((settings.flags & Settings_LockRight) && (f & Geo::Right)) { c = vec3(0.7f, 0.0f, 0.0f); }
		DQLINE(vec3(X[e.v[0]]), vec3(X[con.vertCount + i]), vec3(X[e.v[1]]), c);
	}

	for (uint32_t i = 0; i < con.quadCount; i++) {
		bool isSurface = con.quads[i].p[1] == CON_NULL_IDX;
		if (!isSurface && tooBig) { continue; }
		uint8_t f = flags[con.vertCount + con.edgeCount + i];
		vec3 c = isSurface ? color : mix(color, vec3(1.0f), 0.7f);
		if ((settings.flags & Settings_LockLeft) && (f & Geo::Left)) { c = vec3(0.7f, 0.0f, 0.0f); }
		if ((settings.flags & Settings_LockRight) && (f & Geo::Right)) { c = vec3(0.7f, 0.0f, 0.0f); }
		DPOINT(vec3(X[con.vertCount + con.edgeCount + i]), c);
	}
	for (uint32_t i = 0; i < con.hexCount; i++) {
		DPOINT(vec3(X[con.vertCount + con.edgeCount + con.quadCount + i]), mix(color, vec3(1.0f), 0.7f));
	}
}

float GeoQuadratic3d::CalculateVolume() const {
	float volume = 0.0f;
	for (uint32_t i = 0; i < con.tetCount; i++) { volume += CalculateElementVolume(X, t[i].i, t[i].ep); }
	for (uint32_t i = 0; i < con.hexCount; i++) { volume += CalculateElementVolume(X, h[i].i, h[i].ep); }
	return volume;
}
