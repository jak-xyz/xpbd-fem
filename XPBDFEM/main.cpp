//-----------------------------------------------------------------------------
// All the code for hosting the Demo in Windows.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Header files
//-----------------------------------------------------------------------------
#include <stdint.h>
#include <stdio.h>

#define NOMINMAX 1
#define WIN32_LEAN_AND_MEAN 1
#include "glad/glad.h"
#include "GLFW/glfw3.h"

#include "Demo.h"
#include "GlRenderer.h"

//-----------------------------------------------------------------------------
// Globals
//-----------------------------------------------------------------------------
Demo demo;
GlRenderer glRenderer;
extern "C" double performance_now() { return glfwGetTime(); }

//-----------------------------------------------------------------------------
// GLFW callbacks
//-----------------------------------------------------------------------------
void GlfwErrorCallback(int error, const char* description) {
	fprintf(stderr, "GLFW Error: %s\n", description);
}
void GlfwKeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, GLFW_TRUE);
	}
	if (key == GLFW_KEY_1 && action == GLFW_PRESS) {
		demo.Initialize();
	}
}
void GlfwMouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
	double timestamp = glfwGetTime();
	if (action != GLFW_PRESS && action != GLFW_RELEASE) { return; }
	uint32_t type = action == GLFW_PRESS ? TouchEvent_Start : TouchEvent_End;
	double mouseX, mouseY;
	glfwGetCursorPos(window, &mouseX, &mouseY);
	demo.input.OnTouchEvent(type, TouchDevice_Mouse, timestamp, (uint32_t)button, (float)mouseX, (float)mouseY);
}
void GlfwCursorPositionCallback(GLFWwindow* window, double xpos, double ypos) {
	double timestamp = glfwGetTime();
	for (uint32_t i = 0; i < 5; i++) {
		if (!glfwGetMouseButton(window, (int)i) == GLFW_PRESS) { continue; }
		demo.input.OnTouchEvent(TouchEvent_Move, TouchDevice_Mouse, timestamp, i, (float)xpos, (float)ypos);
	}
}
void GlfwFocusCallback(GLFWwindow* window, int focused) {
	if (focused) { demo.input.OnFocus(); }
}

//-----------------------------------------------------------------------------
// Program entry
//-----------------------------------------------------------------------------
int main() {
	//-------------------------------------------------------------------------
	// GLFW init
	//-------------------------------------------------------------------------
	if (!glfwInit()) {
		fprintf(stderr, "glfwInit failed.\n");
	}
	glfwSetErrorCallback(GlfwErrorCallback);

	GLFWmonitor* monitor = glfwGetPrimaryMonitor();
	const GLFWvidmode* mode = glfwGetVideoMode(monitor);
	glfwWindowHint(GLFW_RED_BITS, mode->redBits);
	glfwWindowHint(GLFW_GREEN_BITS, mode->greenBits);
	glfwWindowHint(GLFW_BLUE_BITS, mode->blueBits);
	glfwWindowHint(GLFW_REFRESH_RATE, mode->refreshRate);
	glfwWindowHint(GLFW_DEPTH_BITS, 24);
	glfwWindowHint(GLFW_STENCIL_BITS, 8);
	glfwWindowHint(GLFW_SAMPLES, 8);
//	GLFWwindow* window = glfwCreateWindow(mode->width, mode->height, "XPBD FEM", monitor, NULL);
	GLFWwindow* window = glfwCreateWindow(704, 704, "XPBD FEM", NULL, NULL);
	if (!window) {
		fprintf(stderr, "GLFWwindow creation failed.\n");
	}
	glfwSetKeyCallback(window, GlfwKeyCallback);
	glfwSetMouseButtonCallback(window, GlfwMouseButtonCallback);
	glfwSetCursorPosCallback(window, GlfwCursorPositionCallback);
	glfwSetWindowFocusCallback(window, GlfwFocusCallback);

	glfwMakeContextCurrent(window);
	gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);

	printf("GL Version: %d.%d \n", GLVersion.major, GLVersion.minor);

	int widthMM, heightMM;
	glfwGetMonitorPhysicalSize(monitor, &widthMM, &heightMM);
	printf("Monitor dimensions: %d, %d\n", widthMM, heightMM);

	glfwSwapInterval(1);

	//-------------------------------------------------------------------------
	// Sim object init
	//-------------------------------------------------------------------------
	demo.Initialize();
	glRenderer.Initialize(32 * 1024, 128 * 1024, 256 * 1024);

	//
	demo.Reset();
	//demo.AddBlock(Element_T3, 25, 7, 1.0f, 1.0f);
	//demo.AddBlock(Element_Q4, 25, 7, 1.0f, 1.0f);
	//demo.AddBlock(Element_Q5, 25, 7, 1.0f, 1.0f);
	//demo.AddBlock(Element_Q9, 25, 7, 1.0f, 1.0f);
	demo.AddBlock(Element_H8, 3, 3, 1.0f, 1.0f);
	//demo.AddBlock(Element_H27, 3, 3, 1.0f, 1.0f);
	demo.FinishAddingBlocks();

	//-------------------------------------------------------------------------
	// Polling update
	//-------------------------------------------------------------------------
	static double timeOld = glfwGetTime() - 1.0 / 60.0;
	while (!glfwWindowShouldClose(window)) {
		double time = glfwGetTime();
		float dt = (float)(time - timeOld);
		timeOld = time;

		//---------------------------------------------------------------------
		// Update sim and graphics
		//---------------------------------------------------------------------
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		glViewport(0, 0, width, height);
		demo.device.viewportDims = vec2((float)width, (float)height);
		//device.ppi = GLFWmonitor * glfwGetPrimaryMonitor(void);

		// Tick the demo
		double simStart = glfwGetTime();
		demo.Update(dt, 1.0f / 60.0f);
		double simEnd = glfwGetTime();
		static double s_accum = 0.0;
		static int s_count = 0;
		s_accum += simEnd - simStart;
		s_count += 1;
		if (s_accum > 0.25) {
			printf("sim time: %g ms\n", 1000.0 * s_accum / (double)s_count);
			s_accum = 0.0;
			s_count = 0;
		}
		
		// Render
		glClearColor(0.95f, 0.95f, 0.95f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_DEPTH_TEST);

		demo.Render();
		glRenderer.Render(
			(GlSimpleVertex*)demo.debugGeo.points, demo.debugGeo.pointCount,
			(GlSimpleVertex*)demo.debugGeo.verts, demo.debugGeo.vertCount,
			demo.debugGeo.idxs, demo.debugGeo.idxCount,
			demo.camera.GetViewFromWorld(), demo.camera.GetProjFromView()
		);

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	//-------------------------------------------------------------------------
	// Cleanup
	//-------------------------------------------------------------------------
	glRenderer.Terminate();
	
	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}
