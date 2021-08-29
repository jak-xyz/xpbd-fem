#include "Geo.h"
#include <float.h>

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
		vec2 target = manip.pos.xy + manip.pickDirTarget.xy * (-manip.pos.z / manip.pickDirTarget.z);
		X[manip.pickedPointIdx] += dvec2((target - vec2(X[manip.pickedPointIdx])) * (1.0f / (1.0f + 0.000001f / (dt * dt))));
	}

	// Update velocity
	for (uint32_t i = 0; i < vertCount; i++) {
		V[i] = (X[i] - O[i]) * (1.0f / dt);
	}

	// Apply damping
	if (settings.pbdDamping > 0.0f) {
		Damp(settings.areaAndTimeCorrectedPbdDamping);
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
	*outDistance = distance(vec2(X[0]), pos);
	for (uint32_t i = 1; i < vertCount; i++) {
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
		float t = dot(manip.manipPlaneNormal, manip.pick0 - manip.pos) / dot(manip.manipPlaneNormal, manip.pickDirTarget);
		vec3 target = manip.pos + t * manip.pickDirTarget; // target point is on the plane defined by the original pick location and camera normal
		X[manip.pickedPointIdx] += dvec3((target - vec3(X[manip.pickedPointIdx])) * (1.0f / (1.0f + 0.000001f / (dt * dt))));
	}

	// Update velocity
	for (uint32_t i = 0; i < vertCount; i++) {
		V[i] = (X[i] - O[i]) * (1.0f / dt);
	}

	// Apply damping
	if (settings.pbdDamping > 0.0f) {
		Damp(settings.volumeAndTimeCorrectedPbdDamping);
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
		if (!(flags[i] & Pickable)) { continue; }
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
