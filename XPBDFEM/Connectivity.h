//-----------------------------------------------------------------------------
// Connectivity information for area and volume meshes.
// Connectivity2d and Connectivity3d can ingest element idx arrays, and then
// make it relatively efficient to iterate over edges, find surface elements,
// etc. In practice, the interface is quite unstructured, so watch out!
//-----------------------------------------------------------------------------
#pragma once

#include "Allocator.h"
#include "vectormath.h"

//-----------------------------------------------------------------------------
const uint32_t CON_TRI = 3;
const uint32_t CON_QUAD = 4;

const uint32_t CON_TET = 4;
const uint32_t CON_HEX = 8;

const uint32_t CON_NULL_IDX = 0xffffffff;

//-----------------------------------------------------------------------------
struct ConVert2d {
	uint32_t edgesOffset, edgeCount;
	uint32_t trisOffset, triCount;
	uint32_t quadsOffset, quadCount;
};
struct ConEdge2d { uint32_t v[2], p[2]; };
struct ConTri2d { uint32_t v[3]; };
struct ConQuad2d { uint32_t v[4]; };
struct alignas(uint32_t) ConOneRing2d {
	uint8_t vertCount, triCount, quadCount, idxCount;
	const uint32_t* Verts() const { return (const uint32_t*)(((uint32_t*)this) + 1); }
	const uint32_t* Tris() const { return (const uint32_t*)(((uint32_t*)this) + 1 + vertCount); }
	const uint32_t* Quads() const { return (const uint32_t*)(((uint32_t*)this) + 1 + vertCount + triCount); }
	const uint8_t* Idxs() const { return (const uint8_t*)(((uint32_t*)this) + 1 + vertCount + triCount + quadCount); }
};

struct Connectivity2d {
	void Init(struct Allocator* alloc, const uint32_t* idxData, uint32_t idxDataCount);

	uint32_t EdgeFromVerts(uint32_t v0Idx, uint32_t v1Idx) const;
	ConOneRing2d* AllocOneRing(struct Allocator* alloc, uint32_t vIdx) const;

	ConVert2d* verts;
	ConEdge2d* edges;
	ConTri2d* tris;
	ConQuad2d* quads;
	uint32_t* triIdxs;
	uint32_t* quadIdxs;

	uint32_t vertCount = 0;
	uint32_t edgeCount = 0;
	uint32_t triCount = 0;
	uint32_t quadCount = 0;
};

//-----------------------------------------------------------------------------
struct ConVert3d {
	uint32_t edgesOffset, edgeCount;
	uint32_t trisOffset, triCount;
	uint32_t quadsOffset, quadCount;
	uint32_t tetsOffset, tetCount;
	uint32_t hexsOffset, hexCount;
};
struct ConEdge3d { uint32_t v[2]; };
struct ConTri3d { uint32_t v[3], p[2]; };
struct ConQuad3d { uint32_t v[4], p[2]; };
struct ConTet3d { uint32_t v[4]; };
struct ConHex3d { uint32_t v[8]; };
struct alignas(uint32_t) ConOneRing3d {
	uint16_t vertCount, tetCount, hexCount, idxCount;
	const uint32_t* Verts() const { return (const uint32_t*)(((uint32_t*)this) + 2); }
	const uint32_t* Tets() const { return (const uint32_t*)(((uint32_t*)this) + 2 + vertCount); }
	const uint32_t* Hexs() const { return (const uint32_t*)(((uint32_t*)this) + 2 + vertCount + tetCount); }
	const uint8_t* Idxs() const { return (const uint8_t*)(((uint32_t*)this) + 2 + vertCount + tetCount + hexCount); }
};

struct Connectivity3d {
	void Init(struct Allocator* alloc, const uint32_t* idxData, uint32_t idxDataCount);

	uint32_t EdgeFromVerts(uint32_t v0Idx, uint32_t v1Idx) const;
	uint32_t QuadFromVerts(uint32_t v0Idx, uint32_t v1Idx, uint32_t v2Idx, uint32_t v3Idx) const;
	void EdgesFromTri(uint32_t triIdx, uint32_t(&outEdgeIdxs)[3]) const;
	void EdgesFromQuad(uint32_t quadIdx, uint32_t(&outEdgeIdxs)[4]) const;
	ConOneRing3d* AllocOneRing(struct Allocator* alloc, uint32_t vIdx) const;

	ConVert3d* verts;
	ConEdge3d* edges;
	ConTri3d* tris;
	ConQuad3d* quads;
	ConTet3d* tets;
	ConHex3d* hexs;
	uint32_t* tetIdxs;
	uint32_t* hexIdxs;
	uint8_t* surfaceEdges;

	uint32_t vertCount = 0;
	uint32_t edgeCount = 0;
	uint32_t triCount = 0;
	uint32_t quadCount = 0;
	uint32_t tetCount = 0;
	uint32_t hexCount = 0;
};

//-----------------------------------------------------------------------------
struct alignas(uint32_t) OneRing {
	uint16_t vertCount, polyCount, idxCount, padding;
	const uint32_t* Verts() const { return (const uint32_t*)(((uint32_t*)this) + 2); }
	const uint32_t* Polys() const { return (const uint32_t*)(((uint32_t*)this) + 2 + vertCount); }
	const uint8_t* Sizes() const { return (const uint8_t*)(((uint32_t*)this) + 2 + vertCount + polyCount); }
	const uint16_t* Idxs() const { return (const uint16_t*)(((uint32_t*)this) + 2 + vertCount + polyCount + (polyCount + 3) / 4); }
};

struct OneRingConstructor {
	OneRingConstructor(uint32_t vIdx0) {
		verts[0] = vIdx0;
		vertCount = 1;
	}

	void AddPoly(uint32_t polyIdx, uint32_t* polyVerts, uint32_t polyVertCount) {
		polys[polyCount] = polyIdx;
		sizes[polyCount] = polyVertCount;
		++polyCount;
		for (uint32_t i = 0; i < polyVertCount; i++) {
			uint32_t vertIdx = vertCount;
			for (uint32_t j = 0; j < vertCount; j++) {
				if (polyVerts[i] == verts[j]) { vertIdx = j; break; }
			}
			if (vertIdx == vertCount) {
				verts[vertCount++] = polyVerts[i];
			}
			idxs[idxCount++] = vertIdx;
		}
	}

	OneRing* Alloc(Allocator* alloc) {
		// Alloc and init the one ring struct, and then allocate arrays off the end of it
		OneRing* r = alloc->New<OneRing>();
		r->vertCount = vertCount;
		r->polyCount = polyCount;
		r->idxCount = idxCount;
		r->padding = 0;
		uint32_t* rverts = alloc->Alloc<uint32_t>(r->vertCount);
		uint32_t* rpolys = alloc->Alloc<uint32_t>(r->polyCount);
		uint8_t* rsizes = alloc->Alloc<uint8_t>(((polyCount + 3) / 4) * 4);
		uint16_t* ridxs = alloc->Alloc<uint16_t>(r->idxCount);
		for (uint32_t i = 0; i < r->vertCount; i++) { rverts[i] = verts[i]; }
		for (uint32_t i = 0; i < r->polyCount; i++) { rpolys[i] = polys[i]; rsizes[i] = sizes[i]; }
		for (uint32_t i = 0; i < r->idxCount; i++) { ridxs[i] = idxs[i]; }
		return r;
	}

	uint32_t verts[768];
	uint32_t vertCount = 0;
	uint32_t polys[128];
	uint32_t polyCount = 0;
	uint8_t sizes[128];
	uint16_t idxs[6144];
	uint32_t idxCount = 0;
};

//-----------------------------------------------------------------------------
void GmshPrintTriBorder(struct Allocator* alloc, const float* nodeData, uint32_t nodeDataCount, const uint32_t* idxData, uint32_t idxDataCount);
