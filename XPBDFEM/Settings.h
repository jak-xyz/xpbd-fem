//-----------------------------------------------------------------------------
// Global program tweak values
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include "vectormath.h"

const uint32_t Settings_ElementTypeBit = 0;
const uint32_t Settings_ElementTypeMask = (1 << 4) - 1;
const uint32_t Settings_ConstructBlockFromSettings = (1 << 4);
const uint32_t Settings_ForceResetBlock = (1 << 5);
const uint32_t Settings_EnergyBit = 6;
const uint32_t Settings_EnergyMask = (1 << 4) - 1;
const uint32_t Settings_ShapeBit = 10;
const uint32_t Settings_ShapeMask = (1 << 4) - 1;
const uint32_t Settings_ConstraintOrderBit = 14;
const uint32_t Settings_ConstraintOrderMask = (1 << 2) - 1;
const uint32_t Settings_Rotate90Degrees = (1 << 17);
const uint32_t Settings_LockLeft = (1 << 18);
const uint32_t Settings_LockRight = (1 << 19);
const uint32_t Settings_RotateLock = (1 << 20);

const uint32_t Element_Null = 0;
const uint32_t Element_T3 = 1;
const uint32_t Element_Q4 = 2;
const uint32_t Element_Q9 = 3;
const uint32_t Element_H8 = 4;
const uint32_t Element_H27 = 5;

const uint32_t Energy_Pixar = 0;
const uint32_t Energy_PixarReduced = 1;
const uint32_t Energy_PixarSel = 2;
const uint32_t Energy_Mixed = 3;
const uint32_t Energy_MixedSerial = 4;
const uint32_t Energy_MixedSub = 5;
const uint32_t Energy_YeohRubber = 6;
const uint32_t Energy_YeohSkin = 7;

const uint32_t Shape_Single = 0;
const uint32_t Shape_Line = 1;
const uint32_t Shape_BeamL = 2;
const uint32_t Shape_BeamM = 3;
const uint32_t Shape_BeamH = 4;
const uint32_t Shape_BeamL1x2 = 5;
const uint32_t Shape_BeamL2x1 = 6;
const uint32_t Shape_BeamL4x1 = 7;
const uint32_t Shape_BeamL8x1 = 8;
const uint32_t Shape_BoxL = 9;
const uint32_t Shape_BoxM = 10;
const uint32_t Shape_BoxH = 11;

struct Settings {
	float timeScale = 1.0f;
	float substepsPerSecond = 1800.0f;
	float overRelaxation = 1.0f;
	uint32_t volumePasses = 0;
	vec2 gravity = vec2(0.0f, -9.81f);
	float compliance = 0.01f;
	float damping = 0.0f;
	float pbdDamping = 0.0f;
	float drag = 0.0f;
	float poissonsRatio = 0.45f;
	float wonkiness = 0.0f;
	float leftRightSeparation = 1.0f;
	uint32_t flags = 0;

	// Automatically updated vars
	float areaAndTimeCorrectedPbdDamping;
	float volumeAndTimeCorrectedPbdDamping;
	float timeCorrectedDrag;
	mat2 lockedRightTransform;
	mat3 lockedRightTransform3d;
	uint32_t tickId = 0;
};
