// Kicks off a web worker thread to run a wasm-based physics simulation, and
// then handles input and rendering for that simulation.

//-----------------------------------------------------------------------------
// Wasm worker
const MSG_RPC = 0;
const MSG_FRAME = 1;
const MSG_SIM_RESULT = 2;

var simWorker = new Worker("sim-worker.js");
var simResults = [null, null];
simWorker.onmessage = function(e) {
	let msg = e.data;
	// Just store off the sim result until we can use it in our anim update
	simResults[!simResults[0] ? 0 : 1] = msg.simResult;
}
function rpc(name, ...args) {
	simWorker.postMessage({ type: MSG_RPC, name, args });
}

//-----------------------------------------------------------------------------
// WebGL interface
var gl = null;
var gfx = null;

class Shader {
	static compileShader(gl, shaderSource, shaderType) {
		let shader = gl.createShader(shaderType);
		gl.shaderSource(shader, shaderSource);
		gl.compileShader(shader);
		let success = gl.getShaderParameter(shader, gl.COMPILE_STATUS);
		if (!success) {
			throw 'Could not compile shader:'+gl.getShaderInfoLog(shader);
		}
		return shader;
	}

	static createProgram(gl, vertexShader, fragmentShader) {
		let program = gl.createProgram();
		gl.attachShader(program, vertexShader);
		gl.attachShader(program, fragmentShader);
		gl.linkProgram(program);
		let success = gl.getProgramParameter(program, gl.LINK_STATUS);
		if (!success) {
			throw ('Program filed to link:'+gl.getProgramInfoLog(program));
		}
		return program;
	}
}

class Gfx {
	constructor(parentElement) {
		this.canvas = document.createElement('canvas');
		parentElement.appendChild(this.canvas);
		gl = this.canvas.getContext('webgl', {
			alpha: false,
			depth: true,
			desynchronized: true,
			preserveDrawingBuffer: true,
			// antialias: Boolean that indicates whether or not to perform anti-aliasing.
		});
		this.OES_element_index_uint = gl.getExtension('OES_element_index_uint');
		if (!this.OES_element_index_uint) {
			console.error('OES_element_index_uint extension is required, but not supported by this platform');
		}

		this.viewportWidth = 100.0;
		this.viewportHeight = 100.0;
		this.ppi = 96.0;

		const vertSrc = `
		#version 100
		uniform mat4 viewFromWorld;
		uniform mat4 projFromView;
		attribute vec3 a_pos;
		attribute vec3 a_color;
		varying vec3 v_color;
		void main()
		{
			v_color = a_color;
			gl_Position = projFromView * viewFromWorld * vec4(a_pos, 1.0);
			gl_PointSize = 3.0;
		}
		`;
		const fragSrc = `
			#version 100
			varying highp vec3 v_color;
			void main() {
				gl_FragData[0] = vec4(v_color, 1.0);
			}
		`;
		this.vertShader = Shader.compileShader(gl, vertSrc, gl.VERTEX_SHADER);
		this.fragShader = Shader.compileShader(gl, fragSrc, gl.FRAGMENT_SHADER);
		this.program = Shader.createProgram(gl, this.vertShader, this.fragShader);
		this.attributes = {
			pos: gl.getAttribLocation(this.program, 'a_pos'),
			color: gl.getAttribLocation(this.program, 'a_color'),
		};
		this.uniforms = {
			viewFromWorld: gl.getUniformLocation(this.program, 'viewFromWorld'),
			projFromView: gl.getUniformLocation(this.program, 'projFromView'),
		};
		this.uniformData = {
			viewFromWorld: new Float32Array(16),
			projFromView: new Float32Array(16),
		};

		this.pointBufSize = 16 * 32 * 1024;
		this.vertBufSize = 16 * 128 * 1024;
		this.idxBufSize = 2 * 256 * 1024;
		this.pointBufs = [];
		this.vertBufs = [];
		this.idxBufs = [];
		for (let i = 0; i < 3; i++) {
			this.pointBufs[i] = gl.createBuffer();
			gl.bindBuffer(gl.ARRAY_BUFFER, this.pointBufs[i]);
			gl.bufferData(gl.ARRAY_BUFFER, this.pointBufSize, gl.DYNAMIC_DRAW);
			this.vertBufs[i] = gl.createBuffer();
			gl.bindBuffer(gl.ARRAY_BUFFER, this.vertBufs[i]);
			gl.bufferData(gl.ARRAY_BUFFER, this.vertBufSize, gl.DYNAMIC_DRAW);
			this.idxBufs[i] = gl.createBuffer();
			gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.idxBufs[i]);
			gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, this.idxBufSize, gl.DYNAMIC_DRAW);
		}
		this.bufferIdx = 2; // Keeps track of which of our triple buffered buffer we're actually using

		this.pointCount = 0;
		this.vertCount = 0;
		this.idxCount = 0;

		this.resize();
	}

	resize() {
		let canvas = this.canvas;
		// let w = canvas.parentElement.clientWidth;
		let w = canvas.clientWidth;
		let h = w;//canvas.parentElement.clientHeight;
		this.viewportWidth = Math.floor(w * window.devicePixelRatio);
		this.viewportHeight = Math.floor(h * window.devicePixelRatio);
		this.ppi = 96.0 * window.devicePixelRatio;
		canvas.width = this.viewportWidth;
		canvas.height = this.viewportHeight;
	}

	ingestSimResult(simResult) {
		let rU32 = new Uint32Array(simResult);
		let head = 0;

		this.uniformData.viewFromWorld.set(new Float32Array(simResult, head, 16));
		head += 4 * 16;
		this.uniformData.projFromView.set(new Float32Array(simResult, head, 16));
		head += 4 * 16;

		this.bufferIdx = this.bufferIdx < 2 ? this.bufferIdx + 1 : 0;

		this.pointCount = rU32[head / 4];
		let pointData = new DataView(simResult, head + 4, 16 * this.pointCount);
		head += 4 + 16 * this.pointCount;
		gl.bindBuffer(gl.ARRAY_BUFFER, this.pointBufs[this.bufferIdx]);
		gl.bufferSubData(gl.ARRAY_BUFFER, 0, pointData);

		this.vertCount = rU32[head / 4];
		let vertData = new DataView(simResult, head + 4, 16 * this.vertCount);
		head += 4 + 16 * this.vertCount;
		gl.bindBuffer(gl.ARRAY_BUFFER, this.vertBufs[this.bufferIdx]);
		gl.bufferSubData(gl.ARRAY_BUFFER, 0, vertData);

		this.idxCount = rU32[head / 4];
		let idxData = new DataView(simResult, head + 4, 4 * this.idxCount);
		head += 4 + 4 * this.idxCount;
		gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.idxBufs[this.bufferIdx]);
		gl.bufferSubData(gl.ELEMENT_ARRAY_BUFFER, 0, idxData);
	}

	render() {
		gl.viewport(0, 0, this.viewportWidth, this.viewportHeight);
		gl.clearColor(0.95, 0.95, 0.95, 1.0);
		gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

		gl.enable(gl.DEPTH_TEST);

		gl.useProgram(this.program);
		gl.enableVertexAttribArray(this.attributes.pos);
		gl.enableVertexAttribArray(this.attributes.color);
		gl.uniformMatrix4fv(this.uniforms.viewFromWorld, false, this.uniformData.viewFromWorld);
		gl.uniformMatrix4fv(this.uniforms.projFromView, false, this.uniformData.projFromView);

		gl.bindBuffer(gl.ARRAY_BUFFER, this.pointBufs[this.bufferIdx]);
		gl.vertexAttribPointer(this.attributes.pos, 3, gl.FLOAT, false, 16, 0);
		gl.vertexAttribPointer(this.attributes.color, 4, gl.UNSIGNED_BYTE, true, 16, 12);
		gl.drawArrays(gl.POINTS, 0, this.pointCount);

		gl.bindBuffer(gl.ARRAY_BUFFER, this.vertBufs[this.bufferIdx]);
		gl.vertexAttribPointer(this.attributes.pos, 3, gl.FLOAT, false, 16, 0);
		gl.vertexAttribPointer(this.attributes.color, 4, gl.UNSIGNED_BYTE, true, 16, 12);
		gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.idxBufs[this.bufferIdx]);
		gl.drawElements(gl.LINES, this.idxCount, gl.UNSIGNED_INT, 0);

		gl.disableVertexAttribArray(this.attributes.color);
		gl.disableVertexAttribArray(this.attributes.pos);
	}
}

//-----------------------------------------------------------------------------
// Update loop
function update(timestamp) {
	simWorker.postMessage({ type: MSG_FRAME, timestamp });

	// Ingest the sim result
	if (simResults[0]) {
		gfx.ingestSimResult(simResults[0]);

		// Return the large result buffer to the thread for reuse
		simWorker.postMessage({ type: MSG_SIM_RESULT, simResult: simResults[0] }, [simResults[0]]);
		simResults[0] = simResults[1];
		simResults[1] = null;
	}

	gfx.render();

	window.requestAnimationFrame(update);
}

//-----------------------------------------------------------------------------
// Initialization
document.addEventListener('DOMContentLoaded', function(event) {
	let content = document.getElementById('content');
	gfx = new Gfx(content);

	let resize = () => {
		gfx.resize();
		rpc('Resize', gfx.viewportWidth, gfx.viewportHeight, gfx.ppi);
	};
	resize();
	window.addEventListener('resize', resize, false);

	//-----------------------------
	// Pass along input to the native code
	const TouchEvent_Move = 0;
	const TouchEvent_Start = 1;
	const TouchEvent_End = 2;
	const TouchEvent_Cancel = 3;

	const TouchDevice_Touch = 0;
	const TouchDevice_Mouse = 1;

	let onTouchShared = (event, type) => {
		if (type != TouchEvent_Cancel) {
			event.preventDefault();
		}

		let time = 0.001 * (event.timeStamp || performance.now());
		let canvasRect = gfx.canvas.getBoundingClientRect();
		for (let i = 0; i < event.changedTouches.length; i++) {
			let t = event.changedTouches[i];
			let id = t.identifier + 5; // Allow room for virtual mouse IDs
			let x = (t.clientX - canvasRect.x) * (gfx.viewportWidth / gfx.canvas.clientWidth);   // Convert from device pixels
			let y = (t.clientY - canvasRect.y) * (gfx.viewportHeight / gfx.canvas.clientHeight); // to viewport pixels
			rpc('OnTouchEvent', type, TouchDevice_Touch, time, id, x, y);
		}
	};
	let onTouchMove = (event) => { onTouchShared(event, TouchEvent_Move); };
	let onTouchStart = (event) => { onTouchShared(event, TouchEvent_Start); };
	let onTouchEnd = (event) => { onTouchShared(event, TouchEvent_End); };
	let onTouchCancel = (event) => { onTouchShared(event, TouchEvent_Cancel); };
	gfx.canvas.addEventListener("touchmove", onTouchMove, false);
	gfx.canvas.addEventListener("touchstart", onTouchStart, false);
	gfx.canvas.addEventListener("touchend", onTouchEnd, false);
	gfx.canvas.addEventListener("touchcancel", onTouchCancel, false);
	
	// Simulate touches with the mouse
	let onMouseShared = (event, type) => {
		let time = 0.001 * (event.timeStamp || performance.now());
		let canvasRect = gfx.canvas.getBoundingClientRect();
		let x = (event.clientX - canvasRect.x) * (gfx.viewportWidth / gfx.canvas.clientWidth);
		let y = (event.clientY - canvasRect.y) * (gfx.viewportHeight / gfx.canvas.clientHeight);

		if (type == TouchEvent_Move) {
			for (let i = 0; i < 5; i++) {
				if (event.buttons & (1 << i)) {
					rpc('OnTouchEvent', type, TouchDevice_Mouse, time, i, x, y);
				}
			}
		} else {
			let map = [0, 2, 1, 3, 4];
			rpc('OnTouchEvent', type, TouchDevice_Mouse, time, map[event.button], x, y);
		}

		return false;
	};
	let onMouseMove = (event) => { onMouseShared(event, TouchEvent_Move); };
	let onMouseDown = (event) => {
		onMouseShared(event, TouchEvent_Start);
		window.getSelection().removeAllRanges();
	};
	let onMouseUp = (event) => { onMouseShared(event, TouchEvent_End); };
	window.addEventListener("mousemove", onMouseMove, false);
	gfx.canvas.addEventListener("mousedown", onMouseDown, false);
	window.addEventListener("mouseup", onMouseUp, false);
	gfx.canvas.addEventListener("contextmenu", (event) => { event.preventDefault(); }, false);

	let onWindowFocus = () => {
		rpc('OnFocus');
		console.log('OnFocus');
	};
	window.addEventListener('focus', onWindowFocus, false);

	//-----------------------------
	// Kick off an update loop
	window.requestAnimationFrame(update);
});
