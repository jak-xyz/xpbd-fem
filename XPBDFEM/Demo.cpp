#define _USE_MATH_DEFINES

#include "Demo.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "DebugGeo.h"
#include "Input.h"

#include "T3.h"
#include "Q4.h"
#include "Q9.h"
#include "H8.h"
#include "H27.h"
extern "C" double performance_now();

//-----------------------------------------------------------------------------
const float kSpacing = (20.0f / 100.0f) / (float)(31);

//-----------------------------------------------------------------------------
void Sim::Initialize() {
	settings.gravity = vec2(0.0f, -0.5f);
	settings.drag = 0.005f;
	settings.flags =
		(Energy_MixedSub << Settings_EnergyBit) |
		Settings_LockLeft;

	Reset();
}

void Sim::Update(float dt, float medianFrameTime, Manipulator* manip) {
	// Calculate how many substeps to take
	float sdt = 1.0f / settings.substepsPerSecond; // Substep dt
	dtResidual += settings.timeScale * dt;
	uint32_t substeps = (uint32_t)(dtResidual / sdt);
	uint32_t maxSubsteps = (uint32_t)ceil((float)medianFrameTime / sdt);
	//@TODO: If the processor can't keep up, just slow down time
	if (substeps > maxSubsteps) {
		substeps = maxSubsteps;
		dtResidual = 0.0f;
	} else {
		dtResidual -= sdt * (float)substeps;
	}

	// Update frame-dependant settings
	// The user set damping assuming 1000 steps a second, so we have to adjust based on actual step count
	float timeCorrectedPbdDamping = 1.0f - pow(1.0f - settings.pbdDamping, 1000.0f * sdt);
	// Furthermore, the damping is affected by element size, so we pre-scale based on average element size,
	// and then divide out the actual element size when applying it
	settings.areaAndTimeCorrectedPbdDamping = timeCorrectedPbdDamping * kSpacing * kSpacing;
	settings.volumeAndTimeCorrectedPbdDamping = timeCorrectedPbdDamping * 6.0f * kSpacing * kSpacing;
	settings.timeCorrectedDrag = 1.0f - pow(1.0f - settings.drag, 1000.0f * sdt);

	// Substep
	double simStart = performance_now();
	for (uint32_t substep = 0; substep < substeps; substep++) {
		// If we're automatically animating the right sides of blocks, update the transform
		if (settings.flags & Settings_LockRight) {
			mat2 scale = mat2(vec2(settings.leftRightSeparation * 2.0f - 1.0f, 0.0f), vec2(0.0f, 1.0f));
			if (settings.flags & Settings_RotateLock) {
				rightRotationTheta += sdt;
				if (rightRotationTheta > (float)(2.0 * M_PI)) { rightRotationTheta -= (float)(2.0 * M_PI); }
			} else {
				rightRotationTheta = (settings.flags & Settings_Rotate90Degrees) ? (float)(1.5f * M_PI) : 0.0f;
			}
			settings.lockedRightTransform = rotation(rightRotationTheta) * scale;
			settings.lockedRightTransform3d = mat3(vec3(settings.lockedRightTransform[0], 0.0f), vec3(settings.lockedRightTransform[1], 0.0f), vec3(0.0f, 0.0f, 1.0f));
		}
		settings.tickId = tickId;

		// Move the pick manipulator with substep resolution to avoid janky movement when interacting with geo
		manip->pickDirTarget = mix(manip->pickDirOld, manip->pickDir, (float)(substep + 1) / (float)substeps);

		for (uint32_t i = 0; i < geoCount; i++) {
			geos[i]->Substep(settings, *manip, sdt);
		}
		++tickId;
	}
	double simEnd = performance_now();
	updateTimeAccum += simEnd - simStart;
	updateCount++;
	if (lastReportTimestamp == 0.0) { lastReportTimestamp = simStart; }
	if (lastReportTimestamp + 500.0 < simStart) {
		savedAverageUpdate = (float)updateTimeAccum / (float)updateCount;
		savedFrameTime = 1000.0f * medianFrameTime;
		updateTimeAccum = 0.0;
		updateCount = 0;
		lastReportTimestamp = simStart;
	}
}

void Sim::Render() {
	for (uint32_t i = 0; i < geoCount; i++) {
		geos[i]->Render(settings);

		float volumePercent = 100.0f * geos[i]->CalculateVolume() / geos[i]->volume0;
		DTEXT(vec3(debugTextOffset - vec2(0.0, 0.01f * (float)i), 0.0f), vec2(0.005f), mat3(1.0f), geos[i]->color, "volume: %4.1f%%", volumePercent);
	}
	if (geoCount) {
		DTEXT(vec3(-0.1125f, debugTextOffset.y, 0.0f), vec2(0.005f), mat3(1.0f), vec3(0.0f), "%.2f / %.2f (%.1f%%)", savedAverageUpdate, savedFrameTime, 100.0f * savedAverageUpdate / savedFrameTime);
		//DTEXT(vec3(-0.1125f, 0.1075f, 0.0f), vec2(0.005f), mat3(1.0f), vec3(0.0f), "%.2f / %.2f (%.1f%%)", savedAverageUpdate, savedFrameTime, 100.0f * savedAverageUpdate / savedFrameTime);
	}
}

void Sim::AddBlock(uint32_t type, uint32_t width, uint32_t height, float xScale, float yScale) {
	assert(geoCount < kMaxGeos);

	vec2 elemScale = 0.7f * vec2(kSpacing) * vec2(xScale, yScale);
	float mass = 0.0025f;

	switch (type) {
	case Element_Null: {} break;
	case Element_T3: {
		T3Block* t3block = alloc.New<T3Block>();
		t3block->InitBlock(settings, &alloc, width, height, elemScale, mass);
		geos[geoCount++] = t3block;
		break; }
	case Element_Q4: {
		Q4Block* q4block = alloc.New<Q4Block>();
		q4block->InitBlock(settings, &alloc, width, height, elemScale, mass);
		geos[geoCount++] = q4block;
		break; }
	case Element_Q9: {
		Q9Block* q9block = alloc.New<Q9Block>();
		q9block->InitBlock(settings, &alloc, width, height, elemScale, mass);
		geos[geoCount++] = q9block;
		break; }
	case Element_H8: {
		H8Block* h8block = alloc.New<H8Block>();
		h8block->InitBlock(settings, &alloc, width, height, height, elemScale.xyy, mass);
		geos[geoCount++] = h8block;
		break; }
	case Element_H27: {
		H27Block* h27block = alloc.New<H27Block>();
		h27block->InitBlock(settings, &alloc, width, height, height, elemScale.xyy, mass);
		geos[geoCount++] = h27block;
		break; }
	default:
		assert(!"Invalid type given");
	}
}

void Sim::FinishAddingBlocks() {
	for (uint32_t i = 0; i < geoCount; i++) {
		float yOff = (float)i * -0.03f + (float)(geoCount - 1) * 0.015f;
		bool shouldRotate = settings.flags & Settings_Rotate90Degrees;
		mat2 rot = shouldRotate ? mat2(vec2(0.0f, -1.0f), vec2(1.0f, 0.0f)) : mat2(1.0f);
		geos[i]->Transform(transform(rot, geoOffset + rot * vec2(0.0f, yOff)));
		geos[i]->color = vec3(0.0f);
		geos[i]->volume0 = geos[i]->CalculateVolume();
	}

	uint32_t memUsed = (uint32_t)(alloc.head - alloc.mem);
	printf("Mem used: %u / %u (%u%%)\n", memUsed, (uint32_t)alloc.capacity, 100 * memUsed / (uint32_t)alloc.capacity);
}

void Sim::Reset() {
	alloc.Initialize(allocMem, sizeof(allocMem));
	geoCount = 0;
	dtResidual = 0.0f;
	tickId = 0;
	rightRotationTheta = 0.0;
}

void Sim::SetGeoOffset(vec2 offset) {
	if (geoOffset.x == offset.x && geoOffset.y == offset.y) { return; }
	vec2 delta = offset - geoOffset;
	for (uint32_t i = 0; i < geoCount; i++) {
		geos[i]->Transform(transform(mat2(1.0f), delta));
	}
	geoOffset = offset;
}

//-----------------------------------------------------------------------------
void Demo::Initialize() {
	debugGeo.Initialize();
	sims[0].Initialize();
	sims[1].Initialize();
	sims[0].debugTextOffset = vec2(0.07125f, 0.1075f);
	sims[1].debugTextOffset = vec2(0.07125f, -0.1125f);
}

void Demo::Update(float dt, float medianFrameTime) {
	debugGeo.BeginFrame();

	// Update the camera
	vec3 camTarget = vec3(0.0f, 0.0f, 0.0f);
	static bool s_initCam = true;
	if (s_initCam) {
		vec3 camPos = vec3(0.0f, 0.0f, 0.2f);
		camera.Update(camPos, camTarget - camPos, vec3(0.0f, 1.0f, 0.0f), 60.0f);
		s_initCam = false;
	}
	if (input.mouse.button[1]) {
		vec2 move = (input.mouse.pos - input.mouse.posOld) * vec2(2.0f, -2.0f) / device.viewportDims;
		vec3 camPos = camera.Pos();
		camPos += 1.5f * (camera.Right() * -move.x + camera.Up() * move.y);
		camPos = camTarget + 0.2f * normalize(camPos - camTarget);
		camera.Update(camPos, camTarget - camPos, vec3(0.0f, 1.0f, 0.0f), 60.0f);
	}

	// Update the picker/manipulator
	if (input.mouse.button[0]) {
		manip.pos = camera.Pos();
		manip.manipPlaneNormal = camera.Dir();
		vec2 ndc = input.mouse.pos * vec2(2.0f, -2.0f) / GetDevice().viewportDims - vec2(1.0f, -1.0f);
		manip.pickDirOld = manip.pickDir;
		manip.pickDir = camera.NdcToDir(ndc);
		if (!input.mouse.buttonOld[0]) {
			float distance = 1.0e20f;
			manip.pickedGeo = nullptr;
			for (uint32_t simIdx = 0; simIdx < 2; simIdx++) {
				for (uint32_t i = 0; i < sims[simIdx].geoCount; i++) {
					vec3 iNearestPoint;
					uint32_t iNearestPointIdx = 0;
					float iDistance = 1.0e20f;
					sims[simIdx].geos[i]->Pick(manip.pos, manip.pickDir, &iNearestPoint, &iNearestPointIdx, &iDistance);
					if (distance > iDistance) {
						distance = iDistance;
						manip.pickedGeo = sims[simIdx].geos[i];
						manip.pick0 = iNearestPoint;
						manip.pickedPointIdx = iNearestPointIdx;
					}
				}
			}
			manip.pickDirOld = manip.pickDir;
		}
	} else {
		manip.pickedGeo = nullptr;
	}

	sims[0].Update(dt, medianFrameTime, &manip);
	sims[1].Update(dt, medianFrameTime, &manip);

	input.Flush();
}

void Demo::Render() {
	sims[0].Render();
	sims[1].Render();
}

void Demo::UpdateSettings(uint32_t simIdx, const Settings& settings) {
	// See if the settings object is trying to set the current block
	if (settings.flags & Settings_ConstructBlockFromSettings) {
		uint32_t curType = (sims[simIdx].settings.flags >> Settings_ElementTypeBit) & Settings_ElementTypeMask;
		uint32_t newType = (settings.flags >> Settings_ElementTypeBit) & Settings_ElementTypeMask;
		uint32_t curShape = (sims[simIdx].settings.flags >> Settings_ShapeBit) & Settings_ShapeMask;
		uint32_t newShape = (settings.flags >> Settings_ShapeBit) & Settings_ShapeMask;
		bool curRotate = sims[simIdx].settings.flags & Settings_Rotate90Degrees;
		bool newRotate = settings.flags & Settings_Rotate90Degrees;
		uint32_t curLock = sims[simIdx].settings.flags & (Settings_LockLeft | Settings_LockRight | Settings_RotateLock);
		uint32_t newLock = settings.flags & (Settings_LockLeft | Settings_LockRight | Settings_RotateLock);
		bool shapeChanged = curType != newType || curShape != newShape || curRotate != newRotate || curLock != newLock || sims[simIdx].settings.wonkiness != settings.wonkiness;
		if (shapeChanged || (settings.flags & Settings_ForceResetBlock)) {
			manip.pickedGeo = nullptr;
			sims[simIdx].Reset();

			bool isQuadratic = newType == Element_Q9 || newType == Element_H27;
			bool is3D = newType == Element_H8 || newType == Element_H27;
			uint32_t width, height;
			float dimX, dimY;
			switch (newShape) {
			case Shape_Single:
				width = height = 1;
				dimX = dimY = 2.0f;
				break;
			case Shape_Line:
				width = 24; height = 1;
				dimX = dimY = 1.0f;
				break;
			case Shape_BeamL:
				width = 16; height = 4;
				dimX = dimY = 1.5f;
				break;
			case Shape_BeamM:
				width = 32; height = 8;
				dimX = dimY = 0.75f;
				break;
			case Shape_BeamH:
				width = 48; height = 12;
				dimX = dimY = 0.5f;
				break;
			case Shape_BeamL1x2:
				width = 16; height = 8;
				dimX = 1.5f; dimY = 0.75f;
				break;
			case Shape_BeamL2x1:
				width = 32; height = 4;
				dimX = 0.75f; dimY = 1.5f;
				break;
			case Shape_BeamL4x1:
				width = 64; height = 4;
				dimX = 0.375f; dimY = 1.5f;
				break;
			case Shape_BeamL8x1:
				width = 128; height = 4;
				dimX = 0.1875f; dimY = 1.5f;
				break;
			case Shape_BoxL:
				width = height = 16;
				dimX = dimY = 1.5f;
				break;
			case Shape_BoxM:
				width = height = 24;
				dimX = dimY = 1.0f;
				break;
			case Shape_BoxH: [[fallthrough]];
			default:
				width = height = 48;
				dimX = dimY = 0.5f;
				break;
			}
			if (is3D) {
				width = width == 1 ? 1 : (width / 2);
				height = height == 1 ? 1 : (height / 2);
				dimX *= 2.0f; dimY *= 2.0f;
			}
			if (isQuadratic && width < 2) { width = 2; }
			if (isQuadratic && height < 2) { height = 2; }
			sims[simIdx].settings = settings;
			sims[simIdx].AddBlock(newType, width + 1, height + 1, dimX, dimY);
			sims[simIdx].FinishAddingBlocks();
			bool bothHaveGeos = sims[0].geoCount && sims[1].geoCount;
			bool rot0 = sims[0].settings.flags & Settings_Rotate90Degrees;
			bool rot1 = sims[1].settings.flags & Settings_Rotate90Degrees;
			vec2 off0, off1;
			if (!bothHaveGeos) { off0 = vec2(0.0f, 0.0f); off1 = vec2(0.0f, 0.0f); }
			else if (rot0 && rot1) { off0 = vec2(-0.0325f, 0.0f); off1 = vec2(0.0325f, 0.0f); }
			else if (!rot0 && rot1) { off0 = vec2(-0.02f, 0.0f); off1 = vec2(0.06f, 0.0f); }
			else if (rot0 && !rot1) { off0 = vec2(-0.06f, 0.0f); off1 = vec2(0.02f, 0.0f); }
			else /*(!rot0 && !rot1)*/ { off0 = vec2(0.0f, 0.03f); off1 = vec2(0.0f, -0.03f); }
			sims[0].SetGeoOffset(off0);
			sims[1].SetGeoOffset(off1);
		}
	}

	sims[simIdx].settings = settings;
}

void Demo::AddBlock(uint32_t type, uint32_t width, uint32_t height, float xScale, float yScale) {
	sims[0].AddBlock(type, width, height, xScale, yScale);
}

void Demo::FinishAddingBlocks() {
	sims[0].FinishAddingBlocks();
	sims[1].FinishAddingBlocks();
}

void Demo::Reset() {
	manip.pickedGeo = nullptr;
	sims[0].Reset();
	sims[1].Reset();
}
