#include "wasm-compat.h"
#include "debugout.h"

#include "Demo.h"

extern uint8_t __heap_base;

extern "C" __attribute__((visibility ("default"))) void Init();
extern "C" __attribute__((visibility ("default"))) void Update(double timestamp, double medianFrameTime);
extern "C" __attribute__((visibility ("default"))) void* GetPtrArray();
extern "C" __attribute__((visibility ("default"))) void Resize(uint32_t width, uint32_t height, float ppi);
extern "C" __attribute__((visibility ("default"))) void OnTouchEvent(uint32_t eventType, uint32_t device, double timeStamp, uint32_t osId, double clientX, double clientY);
extern "C" __attribute__((visibility ("default"))) void OnFocus();
extern "C" __attribute__((visibility ("default"))) void AddBlock(uint32_t type, uint32_t width, uint32_t height, float xScale, float yScale);
extern "C" __attribute__((visibility ("default"))) void FinishAddingBlocks();
extern "C" __attribute__((visibility ("default"))) void ResetDemo();
extern "C" __attribute__((visibility ("default"))) void UpdateSettings(uint32_t simIdx, float timeScale, float substepsPerSecond, float overRelaxation, uint32_t volumePasses, float gravity, float compliance, float damping, float pbdDamping, float drag, float poissonsRatio, float wonkiness, float leftRightSeparation, uint32_t flags);

struct State {
	double timestampOld = 0.0;
	Demo demo;
	mat4 viewFromWorld;
	mat4 projFromView;
	void* uniformPtr;
	void* pointPtr;
	void* vertPtr;
	void* idxPtr;
};
State* state = nullptr;

void Init() {
	// Global ctors aren't run, so we need to make sure we allocate anything with a ctor manually
	state = new (&__heap_base) State();
	state->demo.Initialize();
	state->uniformPtr = &state->viewFromWorld;
	state->pointPtr = &state->demo.debugGeo.pointCount;
	state->vertPtr = &state->demo.debugGeo.vertCount;
	state->idxPtr = &state->demo.debugGeo.idxCount;
}

void Update(double timestamp, double medianFrameTime) {
	float dt = 0.001f * (float)(timestamp - state->timestampOld);
	double preUpdate = performance_now();
	state->demo.Update(dt, 0.001f * (float)medianFrameTime);
	double postUpdate = performance_now();
	state->demo.Render();
	state->viewFromWorld = state->demo.camera.GetViewFromWorld();
	state->projFromView = state->demo.camera.GetProjFromView();

	state->timestampOld = timestamp;
}

void* GetPtrArray() {
	return &state->uniformPtr;
}

void Resize(uint32_t width, uint32_t height, float ppi) {
	state->demo.device.viewportDims = vec2(width, height);
	state->demo.device.ppi = ppi;
}

void OnTouchEvent(uint32_t eventType, uint32_t device, double timeStamp, uint32_t osId, double clientX, double clientY) {
	state->demo.input.OnTouchEvent(eventType, device, timeStamp, osId, clientX, clientY);
}

void OnFocus() {
	state->demo.input.OnFocus();
}

void AddBlock(uint32_t type, uint32_t width, uint32_t height, float xScale, float yScale) {
	state->demo.AddBlock(type, width, height, xScale, yScale);
}

void FinishAddingBlocks() {
	state->demo.FinishAddingBlocks();
}

void ResetDemo() {
	state->demo.Reset();
}

void UpdateSettings(uint32_t simIdx, float timeScale, float substepsPerSecond, float overRelaxation, uint32_t volumePasses, float gravity, float compliance, float damping, float pbdDamping, float drag, float poissonsRatio, float wonkiness, float leftRightSeparation, uint32_t flags) {
	state->demo.UpdateSettings(simIdx, Settings {
		timeScale,
		substepsPerSecond,
		overRelaxation,
		volumePasses,
		gravity * vec2(0.0f, -9.81f),
		compliance,
		damping,
		pbdDamping,
		drag,
		poissonsRatio,
		wonkiness,
		leftRightSeparation,
		flags,
	});
}
