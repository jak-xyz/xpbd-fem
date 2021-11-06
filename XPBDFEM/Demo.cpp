#define _USE_MATH_DEFINES

#include "Demo.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "Connectivity.h"
#include "DebugGeo.h"
#include "Input.h"
#include "MeshGen.h"

extern "C" double performance_now();

const float kSpacing = (20.0f / 100.0f) / (float)(31);

//-----------------------------------------------------------------------------
// Sim methods
void Sim::Initialize() {
	//settings.timeScale = 0.005f;
	settings.substepsPerSecond = 20000.0f;
	settings.gravity = vec2(0.0f, -0.602f);
	settings.compliance = 3.2f;
	settings.pbdDamping = 0.0f;
	settings.damping = 0.00165f;
	settings.poissonsRatio = 0.5f;
	settings.wonkiness = 0.0f;
	settings.flags =
		(Energy_YeohSkinFast << Settings_EnergyBit) |
		Settings_Rotate90Degrees |
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
	// Now we need similar values for the amortized damping, which only gets run every AmortizationPeriod supsteps
	float amortizedTimeCorrectedPbdDamping = 1.0f - pow(1.0f - settings.pbdDamping, 1000.0f * (float)AmortizationPeriod * sdt);
	settings.amortizedAreaAndTimeCorrectedPbdDamping = amortizedTimeCorrectedPbdDamping * kSpacing * kSpacing;
	settings.amortizedVolumeAndTimeCorrectedPbdDamping = amortizedTimeCorrectedPbdDamping * 6.0f * kSpacing * kSpacing;
	// Similar idea for time corrected drag, too
	settings.timeCorrectedDrag = 1.0f - pow(1.0f - settings.drag, 1000.0f * sdt);

	// Substep
	double simStart = performance_now();
	for (uint32_t substep = 0; substep < substeps; substep++) {
		// If we're automatically animating the right sides of blocks, update the transform
		if (settings.flags & Settings_LockRight) {
			float leftRightSeparation = mix(leftRightSeparationOld, settings.leftRightSeparation, (float)substep / (float)substeps);
			mat2 scale = mat2(vec2(leftRightSeparation * 2.0f - 1.0f, 0.0f), vec2(0.0f, 1.0f));
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
	leftRightSeparationOld = settings.leftRightSeparation;
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
		DTEXT(vec3(debugTextOffset - vec2(0.0f, 0.01f * (float)i), 0.0f), vec2(0.005f), mat3(1.0f), geos[i]->color, "volume: %4.1f%%", volumePercent);
		vec3 reducedAttentionColor = mix(geos[i]->color, vec3(1.0f), 0.45f);
		DTEXT(vec3(debugTextOffset - vec2(0.0485f, 0.01f * (float)i), 0.0f), vec2(0.005f), mat3(1.0f), reducedAttentionColor, "elements: %d", geos[i]->ElementCount());
		DTEXT(vec3(debugTextOffset - vec2(0.0875f, 0.01f * (float)i), 0.0f), vec2(0.005f), mat3(1.0f), reducedAttentionColor, "nodes: %d", geos[i]->VertCount());
	}
	if (geoCount) {
		DTEXT(vec3(-0.1125f, debugTextOffset.y, 0.0f), vec2(0.005f), mat3(1.0f), vec3(0.0f), "%.2f / %.2f (%.1f%%)", savedAverageUpdate, savedFrameTime, 100.0f * savedAverageUpdate / savedFrameTime);
	}
}

void Sim::AddBlock(uint32_t type, const float* nodeData, uint32_t nodeDataCount, const uint32_t* idxData, uint32_t idxDataCount, bool autoResize) {
	assert(geoCount < kMaxGeos);

	float density = autoResize ? 2.0f : 1.0f; //@HACK: Make the armadillo heavy!

	switch (type) {
	case Element_Null: {} break;
	case Element_T3: [[fallthrough]];
	case Element_Q4: {
		GeoLinear2d* linear2d = alloc.New<GeoLinear2d>();
		linear2d->Init(&alloc, density, nodeData, nodeDataCount, idxData, idxDataCount, autoResize);
		geos[geoCount++] = linear2d;
		break; }
	case Element_T6: [[fallthrough]];
	case Element_Q9: {
		GeoQuadratic2d* quadratic2d = alloc.New<GeoQuadratic2d>();
		quadratic2d->Init(&alloc, density, nodeData, nodeDataCount, idxData, idxDataCount, autoResize);
		geos[geoCount++] = quadratic2d;
		break; }
	case Element_T4: [[fallthrough]];
	case Element_H8: {
		GeoLinear3d* linear3d = alloc.New<GeoLinear3d>();
		linear3d->Init(&alloc, density, nodeData, nodeDataCount, idxData, idxDataCount, autoResize);
		geos[geoCount++] = linear3d;
		break; }
	case Element_T10: [[fallthrough]];
	case Element_H27: {
		GeoQuadratic3d* quadratic3d = alloc.New<GeoQuadratic3d>();
		quadratic3d->Init(&alloc, density, nodeData, nodeDataCount, idxData, idxDataCount, autoResize);
		geos[geoCount++] = quadratic3d;
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
	leftRightSeparationOld = 1.0f;
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
// Demo methods
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
		vec3 camPos = vec3(0.0f, 0.0f, 31.0f * kSpacing);
		camera.Update(camPos, camTarget - camPos, vec3(0.0f, 1.0f, 0.0f), 60.0f);
		s_initCam = false;
	}
	if (input.mouse.button[1]) {
		vec2 move = (input.mouse.pos - input.mouse.posOld) * vec2(2.0f, -2.0f) / device.viewportDims;
		move *= 7.5f * distance(camTarget, camera.Pos());
		vec3 camPos = camera.Pos();
		camPos += camera.Right() * -move.x + camera.Up() * move.y;
		camPos = camTarget + (31.0f * kSpacing) * normalize(camPos - camTarget);
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
		uint32_t curPattern = (sims[simIdx].settings.flags >> Settings_ElementPatternBit) & Settings_ElementPatternMask;
		uint32_t newPattern = (settings.flags >> Settings_ElementPatternBit) & Settings_ElementPatternMask;
		bool curRotate = sims[simIdx].settings.flags & Settings_Rotate90Degrees;
		bool newRotate = settings.flags & Settings_Rotate90Degrees;
		bool curCentered = sims[simIdx].settings.flags & Settings_Centered;
		bool newCentered = settings.flags & Settings_Centered;
		uint32_t curLock = sims[simIdx].settings.flags & (Settings_LockLeft | Settings_LockRight | Settings_RotateLock);
		uint32_t newLock = settings.flags & (Settings_LockLeft | Settings_LockRight | Settings_RotateLock);
		bool shapeChanged = curType != newType || curShape != newShape || curPattern != newPattern || curRotate != newRotate || curCentered != newCentered || curLock != newLock || sims[simIdx].settings.wonkiness != settings.wonkiness;
		if (shapeChanged || (settings.flags & Settings_ForceResetBlock)) {
			manip.pickedGeo = nullptr;
			sims[simIdx].Reset();

			const float* nodeData = nullptr;
			uint32_t nodeDataCount = 0;
			const uint32_t* idxData = nullptr;
			uint32_t idxDataCount = 0;

			if (newShape == Shape_Armadillo) {
				GenerateArmadillo(newType, &nodeData, &nodeDataCount, &idxData, &idxDataCount);
			} else {
				uint32_t width, height;
				float dimX, dimY;
				switch (newShape) {
				case Shape_Single: width = 1; height = 1; dimX = dimY = 2.0f; break;
				case Shape_Line:   width = 12; height = 1; dimX = dimY = 2.0f; break;
				case Shape_BeamL:  width = 16; height = 4; dimX = dimY = 1.5f; break;
				case Shape_BeamM:  width = 32; height = 8; dimX = dimY = 0.75f; break;
				case Shape_BeamH:  width = 48; height = 12; dimX = dimY = 0.5f; break;
				case Shape_BeamL1x2: width = 16; height = 8; dimX = 1.5f; dimY = 0.75f; break;
				case Shape_BeamL2x1: width = 32; height = 4; dimX = 0.75f; dimY = 1.5f; break;
				case Shape_BeamL4x1: width = 64; height = 4; dimX = 0.375f; dimY = 1.5f; break;
				case Shape_BeamL8x1: width = 128; height = 4; dimX = 0.1875f; dimY = 1.5f; break;
				case Shape_BoxL: width = height = 16; dimX = dimY = 1.5f; break;
				case Shape_BoxM: width = height = 24; dimX = dimY = 1.0f; break;
				case Shape_BoxH: [[fallthrough]];
				default:         width = height = 48; dimX = dimY = 0.5f; break;
				}
				//width = height = 1; //@HACK:
				bool isQuadratic = newType == Element_T6 || newType == Element_Q9 || newType == Element_T10 || newType == Element_H27;
				bool is3D = newType == Element_T4 || newType == Element_T10 || newType == Element_H8 || newType == Element_H27;
				if (isQuadratic && newShape > Shape_Line) {
					width /= 2; height /= 2;
					dimX *= 2.0f; dimY *= 2.0f;
				}
				if (is3D) {
					width /= 2; height /= 2;
					dimX *= 2.0f; dimY *= 2.0f;
				}
				if (width < 1) { width = 1; }
				if (height < 1) { height = 1; }

				GenerateBlock(newType, &sims[simIdx].alloc, width, height, (0.7f * kSpacing) * vec2(dimX, dimY), newPattern, settings.wonkiness, &nodeData, &nodeDataCount, &idxData, &idxDataCount);
			}

			sims[simIdx].settings = settings;
			sims[simIdx].AddBlock(newType, nodeData, nodeDataCount, idxData, idxDataCount, newShape == Shape_Armadillo);
			sims[simIdx].FinishAddingBlocks();
			bool bothHaveGeos = sims[0].geoCount && sims[1].geoCount;
			bool rot0 = sims[0].settings.flags & Settings_Rotate90Degrees;
			bool rot1 = sims[1].settings.flags & Settings_Rotate90Degrees;
			vec2 off0, off1;
			if (!bothHaveGeos) { off0 = vec2(0.0f, 0.0f); off1 = vec2(0.0f, 0.0f); }
			else if (rot0 && rot1) { off0 = vec2(-5.05f * kSpacing, 0.0f); off1 = vec2(5.05f * kSpacing, 0.0f); }
			else if (!rot0 && rot1) { off0 = vec2(-3.1f * kSpacing, 0.0f); off1 = vec2(9.3f * kSpacing, 0.0f); }
			else if (rot0 && !rot1) { off0 = vec2(-9.3f * kSpacing, 0.0f); off1 = vec2(3.1f * kSpacing, 0.0f); }
			else /*(!rot0 && !rot1)*/ { off0 = vec2(0.0f, 4.65f * kSpacing); off1 = vec2(0.0f, -4.65f * kSpacing); }
			bool overlap0 = sims[0].settings.flags & Settings_Centered;
			bool overlap1 = sims[1].settings.flags & Settings_Centered;
			if (overlap0) { off0 = vec2(0.0f, 0.0f); }
			if (overlap1) { off1 = vec2(0.0f, 0.0f); }
			sims[0].SetGeoOffset(off0);
			sims[1].SetGeoOffset(off1);
		}
	}

	sims[simIdx].settings = settings;
}

void Demo::AddBlock(uint32_t type, uint32_t width, uint32_t height, float xScale, float yScale) {
	uint32_t pattern = (sims[0].settings.flags >> Settings_ElementPatternBit) & Settings_ElementPatternMask;
	const float* nodeData = nullptr;
	uint32_t nodeDataCount = 0;
	const uint32_t* idxData = nullptr;
	uint32_t idxDataCount = 0;
	GenerateBlock(type, &sims[0].alloc, width, height, (0.7f * kSpacing) * vec2(xScale, yScale), pattern, 0.0f, &nodeData, &nodeDataCount, &idxData, &idxDataCount);
	sims[0].AddBlock(type, nodeData, nodeDataCount, idxData, idxDataCount, false);
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
