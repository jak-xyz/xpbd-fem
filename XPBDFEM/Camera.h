//------------------------------------------------------------------------------
// A unified place for systems to affect the camera.
//------------------------------------------------------------------------------
#pragma once

#include "vectormath.h"
#include "Device.h"

class Camera {
public:
	const vec3& Pos() const { return pos; }
	vec3 Dir() const { return -z; }
	const vec3& Right() const { return x; }
	const vec3& Up() const { return y; }

	void Update(const vec3& pos, const vec3& dir, const vec3& targetUp, float verticalFovDegrees) {
		this->pos = pos;
		z = -normalize(dir);
		x = normalize(cross(targetUp, z));
		y = normalize(cross(z, x));
		fov = verticalFovDegrees;
	}

	mat4 GetViewFromWorld() const {
		// Affine inverse:
		mat3 u33 = transpose(mat3(x, y, z));
		return mat4(vec4(u33[0], 0.0f), vec4(u33[1], 0.0f), vec4(u33[2], 0.0f), vec4(-(u33 * pos), 1.0f));
	}

	mat4 GetProjFromView() const {
		return projection(0.05f, fov, GetDevice().viewportDims.x, GetDevice().viewportDims.y);
	}

	vec3 NdcToDir(vec2 ndc) const {
		float t = tan(0.5f * fov * (3.1415926f / 180.0f));
		float aspectRatio = GetDevice().viewportDims.x / GetDevice().viewportDims.y;
		return safeNormalize(-z + t * (aspectRatio * ndc.x * x + ndc.y * y));
	}

private:
	vec3 pos = vec3(0.0f);
	vec3 x = vec3(1.0f, 0.0f, 0.0f);
	vec3 y = vec3(0.0f, 1.0f, 0.0f);
	vec3 z = vec3(0.0f, 0.0f, -1.0f);
	float fov = 60.0f;
};
