//------------------------------------------------------------------------------
// Retains system input messages into an easy-to-use format.
//------------------------------------------------------------------------------
#pragma once

#include "debugout.h"
#include "vectormath.h"

const uint32_t TouchEvent_Move = 0;
const uint32_t TouchEvent_Start = 1;
const uint32_t TouchEvent_End = 2;
const uint32_t TouchEvent_Cancel = 3;

const uint32_t TouchDevice_Touch = 0;
const uint32_t TouchDevice_Mouse = 1;

struct Touch {
	uint32_t id = 0;
	uint32_t device = 0;
	vec2 pos = vec2(0.0f);
	vec2 posOld = vec2(0.0f);
	vec2 posStart = vec2(0.0f);
	float maxDistFromStart = 0.0f;
	double startTime = 0.0;
	double lastUpdateTime = 0.0;
	bool isNew = true;
	bool isDead = false;
	bool isCanceled = false;

	vec2 GetPosAtTime(double time) const {
		// Case 1: Given time is after our history
		if (time >= historyTime[0]) {
			return history[0];
		}
		// Case 2: Given time is within our history
		for (uint32_t i = 1; i < historyCount; i++) {
			if (time < historyTime[i - 1] && time >= historyTime[i]) {
				float t = (float)((time - historyTime[i]) / (historyTime[i - 1] - historyTime[i]));
				return mix(history[i], history[i - 1], t);
			}
		}
		// Case 3: Given time is before our history
		return history[historyCount - 1];
	}
	vec2 GetRecentSpeed(double lastXSeconds = 0.033) const {
		lastXSeconds = max(lastXSeconds, 0.001);
		// The OS can take a little bit of time to report a dead touch, so we
		// ignore the last few milliseconds for those
		double start = historyTime[0] - (isDead ? 0.016 : 0.0);
		return (GetPosAtTime(start) - GetPosAtTime(start - lastXSeconds)) / (float)lastXSeconds;
	}

	// For internal use:
	uint32_t osId = 0;
	vec2 history[16];
	double historyTime[16];
	uint32_t historyCount = 0;
};

struct Mouse {
	static const uint32_t kMaxButtons = 3;
	vec2 pos = vec2(0.0f);
	vec2 posOld = vec2(0.0f);
	bool button[kMaxButtons] = {};
	bool buttonOld[kMaxButtons] = {};
};

class Input {
public:
	const Touch* GetTouchById(uint32_t id) {
		for (uint32_t i = 0; i < touchCount; i++) {
			if (touches[i].id == id) {
				return &touches[i];
			}
		}
		return nullptr;
	}

	void CaptureTouch(const Touch* t) {
		for (uint32_t i = 0; i < touchCount; i++) {
			if (&touches[i] != t) { continue; }

			// Move the touch from the regular touch list to the captured list
			capturedTouches[capturedTouchCount] = *t;
			++capturedTouchCount;
			--touchCount;
			for (; i < touchCount; i++) { touches[i] = touches[i + 1]; }
			return;
		}
		debugout("Given touch not active");
	}

	void OnTouchEvent(uint32_t type, uint32_t device, double timestamp, uint32_t osId, float x, float y) {
		// DVAL(type, device, timestamp, osId, x, y);
		// See if we're already tracking this touch
		Touch* t = nullptr;
		for (uint32_t i = 0; i < touchCount + capturedTouchCount; i++) {
			Touch* touches = i < touchCount ? this->touches : capturedTouches;
			uint32_t j = i < touchCount ? i : (i - touchCount);
			if (touches[j].device == device && touches[j].osId == osId && !touches[j].isDead) {
				t = &touches[j];
				break;
			}
		}

		// Handle starting and stopping touches
		if (type == TouchEvent_Start && !t) {
			t = &touches[touchCount];
			++touchCount;

			*t = Touch();
			++touchIdCounter;
			t->id = touchIdCounter;
			t->device = device;
			t->osId = osId;
			t->posOld = vec2(x, y);
			t->posStart = vec2(x, y);
			t->startTime = timestamp;
			t->isNew = true;
			t->isDead = false;
			t->historyCount = 0;

			if (t->device == TouchDevice_Mouse && t->osId < Mouse::kMaxButtons) {
				bool anyPressed = false;
				for (uint32_t i = 0; i < Mouse::kMaxButtons; i++) { anyPressed = anyPressed || mouse.button[i]; }
				if (!anyPressed) { mouse.posOld = t->posOld; }
				mouse.button[t->osId] = true;
			}
		}
		if ((type == TouchEvent_End || type == TouchEvent_Cancel) && t) {
			t->isDead = true;
			t->isCanceled = type == TouchEvent_Cancel;

			if (t->device == TouchDevice_Mouse && t->osId < Mouse::kMaxButtons) {
				mouse.button[t->osId] = false;
			}
		}

		// Do general state update
		if (t) {
			float maxDistFromStartOld = t->maxDistFromStart;
			t->pos = vec2(x, y);
			t->maxDistFromStart = max(distance(t->pos, t->posStart), t->maxDistFromStart);
			t->lastUpdateTime = timestamp;

			// Hack to eliminate the hitch when the touch first moves out of the dead zone
			bool wasAtStart = maxDistFromStartOld < 0.001f;
			float age = (float)(timestamp - t->startTime);
			if (wasAtStart && t->maxDistFromStart >= 0.001f && t->maxDistFromStart < 200.0f * age) {
				t->posOld = t->pos;
			}

			if (!t->historyCount || timestamp - t->historyTime[0] > 0.004) {
				const uint32_t maxHistory = sizeof(t->history) / sizeof(t->history[0]);
				for (uint32_t i = maxHistory - 1; i > 0; i--) {
					t->history[i] = t->history[i - 1];
					t->historyTime[i] = t->historyTime[i - 1];
				}
				if (t->historyCount < maxHistory) { ++t->historyCount; }

				t->history[0] = vec2(x, y);
				t->historyTime[0] = timestamp;
			}

			if (t->device == TouchDevice_Mouse) {
				mouse.pos = t->pos;
			}
		}

		// Emulate the mouse with touches, the number of fingers on the screen
		// corresponding with the currently pressed mouse button
		Touch* liveTouches[Mouse::kMaxButtons] = {};
		uint32_t liveTouchCount = 0;
		bool anyRealMouse = false;
		for (uint32_t i = 0; i < touchCount; i++) {
			Touch* mt = &touches[i];
			if (mt->device == TouchDevice_Mouse) { anyRealMouse = true; break; }
			if (mt->isDead) { continue; }
			if (liveTouchCount < Mouse::kMaxButtons) { liveTouches[liveTouchCount++] = mt; }
		}
		if (!anyRealMouse) {
			if (liveTouchCount > 0) {
				mouse.pos = vec2(0.0f);
				mouse.posOld = vec2(0.0f);
				for (uint32_t i = 0; i < liveTouchCount; i++) {
					mouse.pos += liveTouches[i]->pos;
					mouse.posOld += liveTouches[i]->posOld;
				}
				mouse.pos /= (float)liveTouchCount;
				mouse.posOld /= (float)liveTouchCount;
			}
			for (uint32_t i = 0; i < Mouse::kMaxButtons; i++) { mouse.button[i] = (i + 1) == liveTouchCount; }
		}
	}

	void OnFocus() {
		touchCount = 0;
		for (uint32_t i = 0; i < capturedTouchCount; i++) {
			capturedTouches[i].isDead = true;
			capturedTouches[i].isCanceled = true;
		}
		capturedTouchCount = 0;
		for (uint32_t i = 0; i < Mouse::kMaxButtons; i++) { mouse.button[i] = false; }
	}

	void Flush() {
		// Update the isNew flags and clear out dead touches
		for (uint32_t touchType = 0; touchType < 2; touchType++) {
			Touch* touches = touchType == 0 ? this->touches : capturedTouches;
			uint32_t* count = touchType == 0 ? &touchCount : &capturedTouchCount;
			for (int32_t i = (int32_t)*count - 1; i >= 0; i--) {
				touches[i].posOld = touches[i].pos;
				touches[i].isNew = false;
				if (touches[i].isDead) {
					--* count;
					for (uint32_t j = i; j < *count; j++) { touches[j] = touches[j + 1]; }
				}
			}
		}
		// Update mouse params
		mouse.posOld = mouse.pos;
		for (uint32_t i = 0; i < Mouse::kMaxButtons; i++) { mouse.buttonOld[i] = mouse.button[i]; }
	}

	static const uint32_t kMaxTouches = 20;

	Touch touches[kMaxTouches];
	uint32_t touchCount = 0;
	Touch capturedTouches[kMaxTouches];
	uint32_t capturedTouchCount = 0;
	uint32_t touchIdCounter = 0;

	Mouse mouse;
};
