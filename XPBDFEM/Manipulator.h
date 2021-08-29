//-----------------------------------------------------------------------------
// Allows the user to manipulate individual verts of the physics objects
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include "vectormath.h"

struct Manipulator {
	vec3 pos, manipPlaneNormal, pick0, pickDir, pickDirOld, pickDirTarget;
	struct Geo* pickedGeo = nullptr;
	uint32_t pickedPointIdx = 0;
};
