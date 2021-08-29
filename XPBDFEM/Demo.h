//-----------------------------------------------------------------------------
// The main app orchestrator. Container for all the demo data and handler for
// all the demo logic.
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include "Allocator.h"
#include "Camera.h"
#include "DebugGeo.h"
#include "Device.h"
#include "Geo.h"
#include "Input.h"
#include "Manipulator.h"
#include "Settings.h"

//-----------------------------------------------------------------------------
struct Sim {
	void Initialize();
	void Update(float dt, float medianFrameTime, Manipulator* manip);
	void Render();

	void AddBlock(uint32_t type, uint32_t width, uint32_t height, float xScale, float yScale);
	void FinishAddingBlocks();
	void Reset();

	void SetGeoOffset(vec2 offset);

	Settings settings;

	Allocator alloc;
	uint8_t allocMem[20 * 1024 * 1024];

	static const uint32_t kMaxGeos = 32;
	Geo* geos[kMaxGeos];
	uint32_t geoCount = 0;

	vec2 geoOffset = vec2(0.0f);

	float dtResidual = 0.0f;
	uint32_t tickId = 0;
	float rightRotationTheta = 0.0;

	double updateTimeAccum = 0.0;
	uint32_t updateCount = 0;
	double lastReportTimestamp = 0.0;
	float savedAverageUpdate = 0.0f;
	float savedFrameTime = 1.0f;

	vec2 debugTextOffset = vec2(0.0f);
};

//-----------------------------------------------------------------------------
struct Demo {
	void Initialize();
	void Update(float dt, float medianFrameTime);
	void Render();

	void UpdateSettings(uint32_t simIdx, const Settings& settings);

	void AddBlock(uint32_t type, uint32_t width, uint32_t height, float xScale, float yScale);
	void FinishAddingBlocks();
	void Reset();

	Device device;
	Input input;
	Camera camera;
	DebugGeo debugGeo;
	Manipulator manip;

	Sim sims[2];
};
