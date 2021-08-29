#include "GlRenderer.h"

#define NOMINMAX 1
#define WIN32_LEAN_AND_MEAN 1
#include <Windows.h>

#include "glad/glad.h"
#include "GLFW/glfw3.h"

#include <assert.h>
#include <stdio.h>

//------------------------------------------------------------------------------
// Shader helpers
#define VERT_POS 0
#define VERT_COLOR 1
#define UNIFORM_ENV 0

const char* kShaderPreamble =
"#version 440 core \n"
"#define VERT_POS 0 \n"
"#define VERT_COLOR 1 \n"
"#define UNIFORM_ENV 0 \n"
;

bool CompileShader(GLuint* shader, const GLenum type, const GLchar* source, const char* name) {
	*shader = glCreateShader(type);
	const GLchar* sources[2] = { kShaderPreamble, source };
	int32_t lengths[2] = { (int32_t)strlen(kShaderPreamble), (int32_t)strlen(source) };
	glShaderSource(*shader, 2, sources, lengths);
	glCompileShader(*shader);

	GLint status;
	glGetShaderiv(*shader, GL_COMPILE_STATUS, &status);
	if (status == 0) {
		GLint logLength;
		glGetShaderiv(*shader, GL_INFO_LOG_LENGTH, &logLength);
		if (logLength > 0) {
			GLchar* log = (GLchar*)malloc(logLength);
			glGetShaderInfoLog(*shader, logLength, &logLength, log);
			fprintf(stderr, "%s shader '%s' compile log:\n%s", type == GL_FRAGMENT_SHADER ? "Fragment" : "Vertex", name, log);
			free(log);
		}

		glDeleteShader(*shader);
		return false;
	}

	return true;
}

bool LinkShaderProgram(GLuint prog) {
	GLint status;

	glLinkProgram(prog);

	glGetProgramiv(prog, GL_LINK_STATUS, &status);
	if (status == 0) {
		fprintf(stderr, "Program link failed\n");

		GLint logLength;
		glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &logLength);
		if (logLength > 0) {
			GLchar* log = (GLchar*)malloc(logLength);
			glGetProgramInfoLog(prog, logLength, &logLength, log);
			fprintf(stderr, "Program link log:\n%s", log);
			free(log);
		}

		return false;
	}

	return true;
}

//------------------------------------------------------------------------------
void GlRenderer::Initialize(uint32_t maxPointVerts, uint32_t maxLineVerts, uint32_t maxLineIdxs) {
	this->maxPointVerts = maxPointVerts;
	this->maxLineVerts = maxLineVerts;
	this->maxLineIdxs = maxLineIdxs;

	glEnable(GL_PROGRAM_POINT_SIZE);

	// Allocate the buffers
	glGenBuffers(kBufferCount, pvbos);
	glGenBuffers(kBufferCount, lvbos);
	glGenBuffers(kBufferCount, libos);
	glGenBuffers(kBufferCount, ubos);
	for (uint32_t i = 0; i < kBufferCount; i++) {
		glBindBuffer(GL_ARRAY_BUFFER, pvbos[i]);
		glBufferData(GL_ARRAY_BUFFER, maxPointVerts * 16, 0, GL_DYNAMIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, lvbos[i]);
		glBufferData(GL_ARRAY_BUFFER, maxLineVerts * 16, 0, GL_DYNAMIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, libos[i]);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, maxLineIdxs * 4, 0, GL_DYNAMIC_DRAW);
		glBindBuffer(GL_UNIFORM_BUFFER, ubos[i]);
		glBufferData(GL_UNIFORM_BUFFER, 2 * 16 * 4, 0, GL_DYNAMIC_DRAW);
	}
	bufferIdx = 2;

	// Create the shader
	const char* vertSource =
		"layout(location = VERT_POS) in vec3 pos; \n"
		"layout(location = VERT_COLOR) in vec4 color; \n"
		" \n"
		"out vec4 v_color; \n"
		" \n"
		"layout(binding = UNIFORM_ENV) uniform EnvBlock { \n"
		"	mat4 viewFromWorld; \n"
		"	mat4 projFromView; \n"
		"} u_env; \n"
		" \n"
		"void main() { \n"
		"	v_color = color; \n"
		"	gl_Position = u_env.projFromView * (u_env.viewFromWorld * vec4(pos, 1.0)); \n"
		"	gl_PointSize = 4.0;"
		"} \n"
		"";
	const char* fragSource =
		"in vec4 v_color; \n"
		" \n"
		"out vec4 o_color; \n"
		" \n"
		"void main() { \n"
		"	o_color = v_color; \n"
		"} \n"
		"";

	GLuint vertShader = 0;
	GLuint fragShader = 0;
	bool vertSuccess = CompileShader(&vertShader, GL_VERTEX_SHADER, (char*)vertSource, "GL Line Renderer");
	bool fragSuccess = CompileShader(&fragShader, GL_FRAGMENT_SHADER, (char*)fragSource, "GL Line Renderer");
	assert(vertSuccess && fragSuccess);

	program = glCreateProgram();
	glAttachShader(program, vertShader);
	glAttachShader(program, fragShader);
	bool linkSuccess = LinkShaderProgram(program);
	assert(linkSuccess);

	glDeleteShader(vertShader);
	glDeleteShader(fragShader);
}

void GlRenderer::Terminate() {
	glDeleteBuffers(kBufferCount, pvbos);
	glDeleteBuffers(kBufferCount, lvbos);
	glDeleteBuffers(kBufferCount, libos);
	glDeleteBuffers(kBufferCount, ubos);
	glDeleteProgram(program);
}

void GlRenderer::Render(
	GlSimpleVertex* pointVerts, uint32_t pointVertCount,
	GlSimpleVertex* lineVerts, uint32_t lineVertCount,
	uint32_t* lineIdxs, uint32_t lineIdxCount,
	const mat4& viewFromWorld, const mat4& projFromView
) {
	assert(pointVertCount <= maxPointVerts);
	assert(lineVertCount <= maxLineVerts);
	assert(lineIdxCount <= maxLineIdxs);

	// Advance the buffer index and upload new data
	bufferIdx = bufferIdx + 1 < kBufferCount ? bufferIdx + 1 : 0;
	glBindBuffer(GL_ARRAY_BUFFER, pvbos[bufferIdx]);
	glBufferSubData(GL_ARRAY_BUFFER, 0, pointVertCount * 16, pointVerts);
	glBindBuffer(GL_ARRAY_BUFFER, lvbos[bufferIdx]);
	glBufferSubData(GL_ARRAY_BUFFER, 0, lineVertCount * 16, lineVerts);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, libos[bufferIdx]);
	glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, lineIdxCount * 4, lineIdxs);
	struct UniformBuffer { mat4 vfw, pfv; } ub = { viewFromWorld, projFromView };
	glBindBuffer(GL_UNIFORM_BUFFER, ubos[bufferIdx]);
	glBufferSubData(GL_UNIFORM_BUFFER, 0, 2 * 16 * 4, &ub);

	// Set up the shader state
	glEnableVertexAttribArray(VERT_POS);
	glEnableVertexAttribArray(VERT_COLOR);
	glBindBufferBase(GL_UNIFORM_BUFFER, UNIFORM_ENV, ubos[bufferIdx]);
	glUseProgram(program);

	// Render points
	glBindBuffer(GL_ARRAY_BUFFER, pvbos[bufferIdx]);
	glVertexAttribPointer(VERT_POS, 3, GL_FLOAT, GL_FALSE, 16, (void*)0);
	glVertexAttribPointer(VERT_COLOR, 4, GL_UNSIGNED_BYTE, GL_TRUE, 16, (void*)12);
	glDrawArrays(GL_POINTS, 0, pointVertCount);

	// Render lines
	glBindBuffer(GL_ARRAY_BUFFER, lvbos[bufferIdx]);
	glVertexAttribPointer(VERT_POS, 3, GL_FLOAT, GL_FALSE, 16, (void*)0);
	glVertexAttribPointer(VERT_COLOR, 4, GL_UNSIGNED_BYTE, GL_TRUE, 16, (void*)12);
	glDrawElements(GL_LINES, lineIdxCount, GL_UNSIGNED_INT, (void*)0);

	// Unbind everything
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glBindBuffer(GL_UNIFORM_BUFFER, 0);
	glDisableVertexAttribArray(VERT_POS);
	glDisableVertexAttribArray(VERT_COLOR);
	glUseProgram(0);
}
