#include "Connectivity.h"

#include <assert.h>
#include <stdio.h>
#include "Allocator.h"

//-----------------------------------------------------------------------------
// Connectivity2d methods
void Connectivity2d::Init(struct Allocator* alloc, const uint32_t* idxData, uint32_t idxDataCount) {
	// Count the number of each type of polygon, so we can allocate space for them
	triCount = 0;
	quadCount = 0;
	// Also figure out the max vert idx, so we can allocate verts
	uint32_t maxVert = 0;
	for (uint32_t i = 0; i < idxDataCount; i++) {
		auto CheckMaxVert = [&](uint32_t count) {
			for (uint32_t j = 0; j < count; j++) {
				uint32_t idx = idxData[i + j + 1];
				if (maxVert < idx) { maxVert = idx; }
			}
		};
		switch (idxData[i]) {
		case CON_TRI: ++triCount; CheckMaxVert(3); i += 3; break;
		case CON_QUAD: ++quadCount; CheckMaxVert(4); i += 4; break;
		default:
			printf("Init(): Unknown polygon type '%d' at index %d.", idxData[i], i);
			assert(false);
		}
	}
	tris = alloc->Alloc<ConTri2d>(triCount);
	quads = alloc->Alloc<ConQuad2d>(quadCount);
	triIdxs = alloc->Alloc<uint32_t>(3 * triCount);
	quadIdxs = alloc->Alloc<uint32_t>(4 * quadCount);
	vertCount = 1 + maxVert;
	verts = alloc->Alloc<ConVert2d>(vertCount);

	// Count how many edges/tris/quads each vert is connected to
	for (uint32_t i = 0; i < vertCount; i++) {
		verts[i].edgesOffset = 0;
		verts[i].edgeCount = 0;
		verts[i].trisOffset = 0;
		verts[i].triCount = 0;
		verts[i].quadsOffset = 0;
		verts[i].quadCount = 0;
	}
	for (uint32_t i = 0; i < idxDataCount; i++) {
		auto ReserveEdge = [&](uint32_t v0, uint32_t v1) { verts[v0 < v1 ? v0 : v1].edgeCount++; };
		const uint32_t* vs = &idxData[i + 1];
		switch (idxData[i]) {
		case CON_TRI:
			ReserveEdge(vs[0], vs[1]);
			ReserveEdge(vs[1], vs[2]);
			ReserveEdge(vs[2], vs[0]);

			for (uint32_t j = 0; j < 3; j++) {
				verts[vs[j]].triCount++;
			}

			i += 3;
			break;
		case CON_QUAD:
			ReserveEdge(vs[0], vs[1]);
			ReserveEdge(vs[1], vs[2]);
			ReserveEdge(vs[2], vs[3]);
			ReserveEdge(vs[3], vs[0]);

			for (uint32_t j = 0; j < 4; j++) {
				verts[vs[j]].quadCount++;
			}

			i += 4;
			break;
		default: break; // Already checked this earlier ...
		}
	}

	// Now that we know many objects each vert points to, we can allot each vert a portion of the object arrays
	uint32_t edgeHead = 0;
	uint32_t triHead = 0;
	uint32_t quadHead = 0;
	for (uint32_t i = 0; i < vertCount; i++) {
		verts[i].edgesOffset = edgeHead;
		edgeHead += verts[i].edgeCount;
		verts[i].edgeCount = 0;
		verts[i].trisOffset = triHead;
		triHead += verts[i].triCount;
		verts[i].triCount = 0;
		verts[i].quadsOffset = quadHead;
		quadHead += verts[i].quadCount;
		verts[i].quadCount = 0;
	}

	// Allocate non-unique edges, tris and quads and link them to verts (we'll uniquify them later)
	// These non-unique arrays are temporary, so we allocate them in a way we can easily roll back
	Allocator tempAlloc = *alloc;
	tempAlloc.Alloc<ConEdge2d>(3 * triCount + 4 * quadCount); // Leave room for our unique edges. We'll fill this in later
	uint32_t* edgeIdxs = tempAlloc.Alloc<uint32_t>(3 * triCount + 4 * quadCount);
	ConEdge2d* nuEdges = tempAlloc.Alloc<ConEdge2d>(3 * triCount + 4 * quadCount);
	edgeHead = 0;
	triHead = 0;
	quadHead = 0;
	for (uint32_t i = 0; i < idxDataCount; i++) {
		auto AddEdge = [&](uint32_t v0, uint32_t v1, uint32_t parentPolygon) {
			nuEdges[edgeHead].v[0] = v0 < v1 ? v0 : v1;
			nuEdges[edgeHead].v[1] = v0 < v1 ? v1 : v0;
			nuEdges[edgeHead].p[0] = parentPolygon;
			nuEdges[edgeHead].p[1] = CON_NULL_IDX;
			ConVert2d& v = verts[v0 < v1 ? v0 : v1];
			edgeIdxs[v.edgesOffset + v.edgeCount++] = edgeHead++;
		};
		const uint32_t* vs = &idxData[i + 1];
		switch (idxData[i]) {
		case CON_TRI:
			AddEdge(vs[0], vs[1], triHead);
			AddEdge(vs[1], vs[2], triHead);
			AddEdge(vs[2], vs[0], triHead);

			for (uint32_t j = 0; j < 3; j++) {
				triIdxs[verts[vs[j]].trisOffset + verts[vs[j]].triCount++] = triHead;
				tris[triHead].v[j] = vs[j];
			}

			++triHead;
			i += 3;
			break;
		case CON_QUAD:
			AddEdge(vs[0], vs[1], triCount + quadHead);
			AddEdge(vs[1], vs[2], triCount + quadHead);
			AddEdge(vs[2], vs[3], triCount + quadHead);
			AddEdge(vs[3], vs[0], triCount + quadHead);

			for (uint32_t j = 0; j < 4; j++) {
				quadIdxs[verts[vs[j]].quadsOffset + verts[vs[j]].quadCount++] = quadHead;
				quads[quadHead].v[j] = vs[j];
			}

			++quadHead;
			i += 4;
			break;
		default: break; // Already checked this earlier ...
		}
	}

	// Uniquify the edges
	// First go through each vert and uniquify the pointers into the non-unique edge array
	edgeCount = 0;
	for (uint32_t i = 0; i < vertCount; i++) {
		ConVert2d& v = verts[i];
		for (uint32_t j = 0; j < v.edgeCount; j++) {
			ConEdge2d& ej = nuEdges[edgeIdxs[v.edgesOffset + j]];
			for (uint32_t k = j + 1; k < v.edgeCount; k++) {
				ConEdge2d& ek = nuEdges[edgeIdxs[v.edgesOffset + k]];
				// Overwrite the pointers to all duplicate edges
				if (ej.v[1] == ek.v[1]) {
					assert(ej.p[1] == CON_NULL_IDX); // There shouldn't be more than 1 duplicate edge
					ej.p[1] = ek.p[0];
					edgeIdxs[v.edgesOffset + k] = edgeIdxs[v.edgesOffset + v.edgeCount - 1];
					--v.edgeCount;
					--k;
				}
			}
		}
		edgeCount += v.edgeCount;
	}
	// Now compact the unique edges into the final edge array
	edges = alloc->Alloc<ConEdge2d>(edgeCount); // This is a permanent buffer, so use the real allocator
	edgeHead = 0;
	for (uint32_t i = 0; i < vertCount; i++) {
		ConVert2d& v = verts[i];
		for (uint32_t j = 0; j < v.edgeCount; j++) {
			edges[edgeHead + j] = nuEdges[edgeIdxs[v.edgesOffset + j]];
		}
		v.edgesOffset = edgeHead;
		edgeHead += v.edgeCount;
	}
}

uint32_t Connectivity2d::EdgeFromVerts(uint32_t v0Idx, uint32_t v1Idx) const {
	uint32_t i0 = v0Idx < v1Idx ? v0Idx : v1Idx;
	uint32_t i1 = v0Idx > v1Idx ? v0Idx : v1Idx;
	const ConVert2d& v0 = verts[i0];
	for (uint32_t j = 0; j < v0.edgeCount; j++) {
		const ConEdge2d& e = edges[v0.edgesOffset + j];
		if (e.v[1] == i1) { return v0.edgesOffset + j; }
	}
	assert(false);
	return CON_NULL_IDX;
}

ConOneRing2d* Connectivity2d::AllocOneRing(struct Allocator* alloc, uint32_t vIdx) const {
	// Find the unique verts in the surrounding tris and quads and generate local indexing
	const ConVert2d& v = verts[vIdx];
	uint32_t ringVerts[64];
	uint32_t ringVertCount = 0;
	ringVerts[ringVertCount++] = vIdx;
	uint8_t idxs[256];
	uint32_t idxCount = 0;
	for (uint32_t i = 0; i < v.triCount; i++) {
		const ConTri2d& t = tris[triIdxs[v.trisOffset + i]];
		for (uint32_t j = 0; j < 3; j++) {
			uint32_t triVertIdx = ringVertCount;
			for (uint32_t k = 0; k < ringVertCount; k++) {
				if (t.v[j] == ringVerts[k]) { triVertIdx = k; break; }
			}
			if (triVertIdx == ringVertCount) {
				ringVerts[ringVertCount++] = t.v[j];
			}
			idxs[idxCount++] = triVertIdx;
		}
	}
	for (uint32_t i = 0; i < v.quadCount; i++) {
		const ConQuad2d& q = quads[quadIdxs[v.quadsOffset + i]];
		for (uint32_t j = 0; j < 4; j++) {
			uint32_t quadVertIdx = ringVertCount;
			for (uint32_t k = 0; k < ringVertCount; k++) {
				if (q.v[j] == ringVerts[k]) { quadVertIdx = k; break; }
			}
			if (quadVertIdx == ringVertCount) {
				ringVerts[ringVertCount++] = q.v[j];
			}
			idxs[idxCount++] = quadVertIdx;
		}
	}

	// Alloc and init the one ring struct, and then allocate arrays off the end of it
	ConOneRing2d* r = alloc->New<ConOneRing2d>();
	r->vertCount = ringVertCount;
	r->triCount = v.triCount;
	r->quadCount = v.quadCount;
	r->idxCount = idxCount;
	uint32_t* rverts = alloc->Alloc<uint32_t>(r->vertCount);
	uint32_t* rtris = alloc->Alloc<uint32_t>(r->triCount);
	uint32_t* rquads = alloc->Alloc<uint32_t>(r->quadCount);
	uint8_t* ridxs = alloc->Alloc<uint8_t>(r->idxCount);
	for (uint32_t i = 0; i < r->vertCount; i++) { rverts[i] = ringVerts[i]; }
	for (uint32_t i = 0; i < r->triCount; i++) { rtris[i] = triIdxs[v.trisOffset + i]; }
	for (uint32_t i = 0; i < r->quadCount; i++) { rquads[i] = quadIdxs[v.quadsOffset + i]; }
	for (uint32_t i = 0; i < r->idxCount; i++) { ridxs[i] = idxs[i]; }

	return r;
}

//-----------------------------------------------------------------------------
// Connectivity3d methods
void Connectivity3d::Init(struct Allocator* alloc, const uint32_t* idxData, uint32_t idxDataCount) {
	// Count the number of each type of polyhedron, so we can allocate space for them
	tetCount = 0;
	hexCount = 0;
	// Also figure out the max vert idx, so we can allocate those
	uint32_t maxVert = 0;
	for (uint32_t i = 0; i < idxDataCount; i++) {
		auto CheckMaxVert = [&](uint32_t count) {
			for (uint32_t j = 0; j < count; j++) {
				uint32_t idx = idxData[i + j + 1];
				if (maxVert < idx) { maxVert = idx; }
			}
		};
		switch (idxData[i]) {
		case CON_TET: ++tetCount; CheckMaxVert(4); i += 4; break;
		case CON_HEX: ++hexCount; CheckMaxVert(8); i += 8; break;
		default:
			printf("Init(): Unknown polyhedron type '%d' at index %d.", idxData[i], i);
			assert(false);
		}
	}
	tets = alloc->Alloc<ConTet3d>(tetCount);
	hexs = alloc->Alloc<ConHex3d>(hexCount);
	tetIdxs = alloc->Alloc<uint32_t>(4 * tetCount);
	hexIdxs = alloc->Alloc<uint32_t>(8 * hexCount);
	vertCount = 1 + maxVert;
	verts = alloc->Alloc<ConVert3d>(vertCount);

	// Count how many edges/tris/quads/tets/hexes each vert is connected to
	for (uint32_t i = 0; i < vertCount; i++) {
		verts[i].edgesOffset = 0;
		verts[i].edgeCount = 0;
		verts[i].trisOffset = 0;
		verts[i].triCount = 0;
		verts[i].quadsOffset = 0;
		verts[i].quadCount = 0;
		verts[i].tetsOffset = 0;
		verts[i].tetCount = 0;
		verts[i].hexsOffset = 0;
		verts[i].hexCount = 0;
	}
	for (uint32_t i = 0; i < idxDataCount; i++) {
		auto ReserveEdge = [&](uint32_t v0, uint32_t v1) { verts[v0 < v1 ? v0 : v1].edgeCount++; };
		auto ReserveTri = [&](uint32_t v0, uint32_t v1, uint32_t v2) { verts[v0 < v1 && v0 < v2 ? v0 : v1 < v2 ? v1 : v2].triCount++; };
		auto ReserveQuad = [&](uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3) {
			uint32_t vminA = v0 < v1 ? v0 : v1;
			uint32_t vminB = v2 < v3 ? v2 : v3;
			uint32_t vmin = vminA < vminB ? vminA : vminB;
			verts[vmin].quadCount++;
		};
		const uint32_t* vs = &idxData[i + 1];
		switch (idxData[i]) {
		case CON_TET:
			ReserveEdge(vs[0], vs[1]);
			ReserveEdge(vs[0], vs[2]);
			ReserveEdge(vs[0], vs[3]);
			ReserveEdge(vs[1], vs[2]);
			ReserveEdge(vs[1], vs[3]);
			ReserveEdge(vs[2], vs[3]);

			ReserveTri(vs[0], vs[1], vs[2]);
			ReserveTri(vs[0], vs[1], vs[3]);
			ReserveTri(vs[0], vs[2], vs[3]);
			ReserveTri(vs[1], vs[2], vs[3]);

			for (uint32_t j = 0; j < 4; j++) {
				verts[vs[j]].tetCount++;
			}

			i += 4;
			break;
		case CON_HEX:
			ReserveEdge(vs[0], vs[1]);
			ReserveEdge(vs[0], vs[2]);
			ReserveEdge(vs[0], vs[4]);
			ReserveEdge(vs[1], vs[3]);
			ReserveEdge(vs[1], vs[5]);
			ReserveEdge(vs[2], vs[3]);
			ReserveEdge(vs[2], vs[6]);
			ReserveEdge(vs[3], vs[7]);
			ReserveEdge(vs[4], vs[5]);
			ReserveEdge(vs[4], vs[6]);
			ReserveEdge(vs[5], vs[7]);
			ReserveEdge(vs[6], vs[7]);

			ReserveQuad(vs[0], vs[1], vs[3], vs[2]);
			ReserveQuad(vs[4], vs[5], vs[7], vs[6]);
			ReserveQuad(vs[0], vs[1], vs[5], vs[4]);
			ReserveQuad(vs[2], vs[3], vs[7], vs[6]);
			ReserveQuad(vs[0], vs[2], vs[6], vs[4]);
			ReserveQuad(vs[1], vs[3], vs[7], vs[5]);

			for (uint32_t j = 0; j < 8; j++) {
				verts[vs[j]].hexCount++;
			}

			i += 8;
			break;
		default: break; // Already checked this earlier ...
		}
	}

	// Now that we know many objects each vert points to, we can allot each vert a portion of the object arrays
	uint32_t edgeHead = 0;
	uint32_t triHead = 0;
	uint32_t quadHead = 0;
	uint32_t tetHead = 0;
	uint32_t hexHead = 0;
	for (uint32_t i = 0; i < vertCount; i++) {
		verts[i].edgesOffset = edgeHead;
		edgeHead += verts[i].edgeCount;
		verts[i].edgeCount = 0;
		verts[i].trisOffset = triHead;
		triHead += verts[i].triCount;
		verts[i].triCount = 0;
		verts[i].quadsOffset = quadHead;
		quadHead += verts[i].quadCount;
		verts[i].quadCount = 0;
		verts[i].tetsOffset = tetHead;
		tetHead += verts[i].tetCount;
		verts[i].tetCount = 0;
		verts[i].hexsOffset = hexHead;
		hexHead += verts[i].hexCount;
		verts[i].hexCount = 0;
	}

	// Allocate non-unique edges, tris and quads and link them to verts (we'll uniquify them later)
	// These non-unique arrays are temporary, so we allocate them in a way we can easily roll back
	Allocator tempAlloc = *alloc;
	tempAlloc.Alloc<ConEdge3d>(6 * tetCount + 12 * hexCount); // Leave room for our unique edges, tris and quads. We'll fill this in later
	tempAlloc.Alloc<ConTri3d>(4 * tetCount);
	tempAlloc.Alloc<ConQuad3d>(6 * hexCount);
	uint32_t* edgeIdxs = tempAlloc.Alloc<uint32_t>(6 * tetCount + 12 * hexCount);
	uint32_t* triIdxs = tempAlloc.Alloc<uint32_t>(4 * tetCount);
	uint32_t* quadIdxs = tempAlloc.Alloc<uint32_t>(6 * hexCount);
	ConEdge3d* nuEdges = tempAlloc.Alloc<ConEdge3d>(6 * tetCount + 12 * hexCount);
	ConTri3d* nuTris = tempAlloc.Alloc<ConTri3d>(4 * tetCount);
	ConQuad3d* nuQuads = tempAlloc.Alloc<ConQuad3d>(6 * hexCount);
	edgeHead = 0;
	triHead = 0;
	quadHead = 0;
	tetHead = 0;
	hexHead = 0;
	for (uint32_t i = 0; i < idxDataCount; i++) {
		auto AddEdge = [&](uint32_t v0, uint32_t v1) {
			nuEdges[edgeHead].v[0] = v0 < v1 ? v0 : v1;
			nuEdges[edgeHead].v[1] = v0 < v1 ? v1 : v0;
			ConVert3d& v = verts[v0 < v1 ? v0 : v1];
			edgeIdxs[v.edgesOffset + v.edgeCount++] = edgeHead++;
		};
		auto AddTri = [&](uint32_t v0, uint32_t v1, uint32_t v2) {
			uint32_t vl = v0 < v1&& v0 < v2 ? v0 : v1 < v2 ? v1 : v2;
			uint32_t vh = v0 > v1 && v0 > v2 ? v0 : v1 > v2 ? v1 : v2;
			uint32_t vm = v0 != vl && v0 != vh ? v0 : v1 != vl && v1 != vh ? v1 : v2;
			nuTris[triHead].v[0] = vl;
			nuTris[triHead].v[1] = vm;
			nuTris[triHead].v[2] = vh;
			nuTris[triHead].p[0] = tetHead;
			nuTris[triHead].p[1] = CON_NULL_IDX;
			ConVert3d& v = verts[vl];
			triIdxs[v.trisOffset + v.triCount++] = triHead++;
		};
		auto AddQuad = [&](uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3) {
			uint32_t vs[4] = { v0, v1, v2, v3 };
			for (uint32_t j = 0; j < 4; j++) {
				if (vs[j] > vs[(j + 1) % 4] || vs[j] > vs[(j + 2) % 4] || vs[j] > vs[(j + 3) % 4]) { continue; }
				nuQuads[quadHead].v[0] = vs[j];
				nuQuads[quadHead].v[1] = vs[(j + 1) % 4];
				nuQuads[quadHead].v[2] = vs[(j + 2) % 4];
				nuQuads[quadHead].v[3] = vs[(j + 3) % 4];
				nuQuads[quadHead].p[0] = hexHead;
				nuQuads[quadHead].p[1] = CON_NULL_IDX;
				ConVert3d& v = verts[vs[j]];
				quadIdxs[v.quadsOffset + v.quadCount++] = quadHead++;
				break;
			}
		};
		const uint32_t* vs = &idxData[i + 1];
		switch (idxData[i]) {
		case CON_TET:
			AddEdge(vs[0], vs[1]);
			AddEdge(vs[0], vs[2]);
			AddEdge(vs[0], vs[3]);
			AddEdge(vs[1], vs[2]);
			AddEdge(vs[1], vs[3]);
			AddEdge(vs[2], vs[3]);

			AddTri(vs[0], vs[1], vs[2]);
			AddTri(vs[0], vs[1], vs[3]);
			AddTri(vs[0], vs[2], vs[3]);
			AddTri(vs[1], vs[2], vs[3]);

			for (uint32_t j = 0; j < 4; j++) {
				tetIdxs[verts[vs[j]].tetsOffset + verts[vs[j]].tetCount++] = tetHead;
				tets[tetHead].v[j] = vs[j];
			}

			++tetHead;
			i += 4;
			break;
		case CON_HEX:
			AddEdge(vs[0], vs[1]);
			AddEdge(vs[0], vs[2]);
			AddEdge(vs[0], vs[4]);
			AddEdge(vs[1], vs[3]);
			AddEdge(vs[1], vs[5]);
			AddEdge(vs[2], vs[3]);
			AddEdge(vs[2], vs[6]);
			AddEdge(vs[3], vs[7]);
			AddEdge(vs[4], vs[5]);
			AddEdge(vs[4], vs[6]);
			AddEdge(vs[5], vs[7]);
			AddEdge(vs[6], vs[7]);

			AddQuad(vs[0], vs[1], vs[3], vs[2]);
			AddQuad(vs[4], vs[5], vs[7], vs[6]);
			AddQuad(vs[0], vs[1], vs[5], vs[4]);
			AddQuad(vs[2], vs[3], vs[7], vs[6]);
			AddQuad(vs[0], vs[2], vs[6], vs[4]);
			AddQuad(vs[1], vs[3], vs[7], vs[5]);

			for (uint32_t j = 0; j < 8; j++) {
				hexIdxs[verts[vs[j]].hexsOffset + verts[vs[j]].hexCount++] = hexHead;
				hexs[hexHead].v[j] = vs[j];
			}

			++hexHead;
			i += 8;
			break;
		default: break; // Already checked this earlier ...
		}
	}

	// Uniquify the edges, tris and quads
	// First go through each vert and uniquify the pointers into the non-unique element arrays
	edgeCount = 0;
	triCount = 0;
	quadCount = 0;
	for (uint32_t i = 0; i < vertCount; i++) {
		ConVert3d& v = verts[i];
		// Edges
		for (uint32_t j = 0; j < v.edgeCount; j++) {
			ConEdge3d& ej = nuEdges[edgeIdxs[v.edgesOffset + j]];
			for (uint32_t k = j + 1; k < v.edgeCount; k++) {
				ConEdge3d& ek = nuEdges[edgeIdxs[v.edgesOffset + k]];
				// Overwrite the pointers to all duplicate edges
				if (ej.v[1] == ek.v[1]) {
					edgeIdxs[v.edgesOffset + k] = edgeIdxs[v.edgesOffset + v.edgeCount - 1];
					--v.edgeCount;
					--k;
				}
			}
		}
		edgeCount += v.edgeCount;
		// Tris
		for (uint32_t j = 0; j < v.triCount; j++) {
			ConTri3d& tj = nuTris[triIdxs[v.trisOffset + j]];
			for (uint32_t k = j + 1; k < v.triCount; k++) {
				ConTri3d& tk = nuTris[triIdxs[v.trisOffset + k]];
				// Overwrite the pointers to all duplicate tris
				if (tj.v[1] == tk.v[1] && tj.v[2] == tk.v[2]) {
					assert(tj.p[1] == CON_NULL_IDX); // There shouldn't be more than 1 duplicate tri
					tj.p[1] = tk.p[0];
					triIdxs[v.trisOffset + k] = triIdxs[v.trisOffset + v.triCount - 1];
					--v.triCount;
					--k;
				}
			}
		}
		triCount += v.triCount;
		// Quads
		for (uint32_t j = 0; j < v.quadCount; j++) {
			ConQuad3d& qj = nuQuads[quadIdxs[v.quadsOffset + j]];
			for (uint32_t k = j + 1; k < v.quadCount; k++) {
				ConQuad3d& qk = nuQuads[quadIdxs[v.quadsOffset + k]];
				// Overwrite the pointers to all duplicate quads
				if ((qj.v[1] == qk.v[1] && qj.v[2] == qk.v[2] && qj.v[3] == qk.v[3]) ||
					(qj.v[1] == qk.v[3] && qj.v[2] == qk.v[2] && qj.v[3] == qk.v[1])) { // <-- Check for backwards winding
					assert(qj.p[1] == CON_NULL_IDX); // There shouldn't be more than 1 duplicate quad
					qj.p[1] = qk.p[0];
					quadIdxs[v.quadsOffset + k] = quadIdxs[v.quadsOffset + v.quadCount - 1];
					--v.quadCount;
					--k;
				}
			}
		}
		quadCount += v.quadCount;
	}

	// Now compact the unique elements into the final element arrays
	edges = alloc->Alloc<ConEdge3d>(edgeCount); // These are permanent buffers, so use the real allocator
	tris = alloc->Alloc<ConTri3d>(triCount);
	quads = alloc->Alloc<ConQuad3d>(quadCount);
	edgeHead = 0;
	triHead = 0;
	quadHead = 0;
	for (uint32_t i = 0; i < vertCount; i++) {
		ConVert3d& v = verts[i];
		// Edges
		for (uint32_t j = 0; j < v.edgeCount; j++) {
			edges[edgeHead + j] = nuEdges[edgeIdxs[v.edgesOffset + j]];
		}
		v.edgesOffset = edgeHead;
		edgeHead += v.edgeCount;
		// Tris
		for (uint32_t j = 0; j < v.triCount; j++) {
			tris[triHead + j] = nuTris[triIdxs[v.trisOffset + j]];
		}
		v.trisOffset = triHead;
		triHead += v.triCount;
		// Quads
		for (uint32_t j = 0; j < v.quadCount; j++) {
			quads[quadHead + j] = nuQuads[quadIdxs[v.quadsOffset + j]];
		}
		v.quadsOffset = quadHead;
		quadHead += v.quadCount;
	}

	// At this point the vert offsets point directly into the edge, tri and quad arrays

	// Mark edges on the surface
	surfaceEdges = alloc->Alloc<uint8_t>(edgeCount);
	for (uint32_t i = 0; i < edgeCount; i++) { surfaceEdges[i] = 0; }
	for (uint32_t i = 0; i < triCount; i++) {
		if (tris[i].p[1] != CON_NULL_IDX) { continue; }
		uint32_t edgeIdxs[3];
		EdgesFromTri(i, edgeIdxs);
		surfaceEdges[edgeIdxs[0]] = 1;
		surfaceEdges[edgeIdxs[1]] = 1;
		surfaceEdges[edgeIdxs[2]] = 1;
	}
	for (uint32_t i = 0; i < quadCount; i++) {
		if (quads[i].p[1] != CON_NULL_IDX) { continue; }
		uint32_t edgeIdxs[4];
		EdgesFromQuad(i, edgeIdxs);
		surfaceEdges[edgeIdxs[0]] = 1;
		surfaceEdges[edgeIdxs[1]] = 1;
		surfaceEdges[edgeIdxs[2]] = 1;
		surfaceEdges[edgeIdxs[3]] = 1;
	}
}

uint32_t Connectivity3d::EdgeFromVerts(uint32_t v0Idx, uint32_t v1Idx) const {
	uint32_t i0 = v0Idx < v1Idx ? v0Idx : v1Idx;
	uint32_t i1 = v0Idx > v1Idx ? v0Idx : v1Idx;
	const ConVert3d& v0 = verts[i0];
	for (uint32_t j = 0; j < v0.edgeCount; j++) {
		const ConEdge3d& e = edges[v0.edgesOffset + j];
		if (e.v[1] == i1) { return v0.edgesOffset + j; }
	}
	assert(false);
	return CON_NULL_IDX;
}

uint32_t Connectivity3d::QuadFromVerts(uint32_t v0Idx, uint32_t v1Idx, uint32_t v2Idx, uint32_t v3Idx) const {
	uint32_t vs[4] = { v0Idx, v1Idx, v2Idx, v3Idx };
	for (uint32_t j = 0; j < 4; j++) {
		if (vs[j] > vs[(j + 1) % 4] || vs[j] > vs[(j + 2) % 4] || vs[j] > vs[(j + 3) % 4]) { continue; }
		const ConVert3d& v0 = verts[vs[j]];
		for (uint32_t k = 0; k < v0.quadCount; k++) {
			const ConQuad3d& q = quads[v0.quadsOffset + k];
			if ((q.v[1] == vs[(j + 1) % 4] && q.v[2] == vs[(j + 2) % 4] && q.v[3] == vs[(j + 3) % 4]) ||
				(q.v[1] == vs[(j + 3) % 4] && q.v[2] == vs[(j + 2) % 4] && q.v[3] == vs[(j + 1) % 4])) { return v0.quadsOffset + k; }
		}
		break;
	}
	assert(false);
	return CON_NULL_IDX;
}

void Connectivity3d::EdgesFromTri(uint32_t triIdx, uint32_t(&outEdgeIdxs)[3]) const {
	uint32_t found = 0;
	const ConTri3d& t = tris[triIdx];
	const ConVert3d& v0 = verts[t.v[0]];
	for (uint32_t j = 0; j < v0.edgeCount; j++) {
		const ConEdge3d& e = edges[v0.edgesOffset + j];
		if (e.v[1] == t.v[1] || e.v[1] == t.v[2]) { outEdgeIdxs[found++] = v0.edgesOffset + j; }
	}
	const ConVert3d& v1 = verts[t.v[1]];
	for (uint32_t j = 0; j < v1.edgeCount; j++) {
		const ConEdge3d& e = edges[v1.edgesOffset + j];
		if (e.v[1] == t.v[2]) { outEdgeIdxs[found++] = v1.edgesOffset + j; break; }
	}
	assert(found == 3);
}

void Connectivity3d::EdgesFromQuad(uint32_t quadIdx, uint32_t(&outEdgeIdxs)[4]) const {
	uint32_t found = 0;
	const ConQuad3d& q = quads[quadIdx];
	const ConVert3d& v0 = verts[q.v[0]];
	for (uint32_t j = 0; j < v0.edgeCount; j++) {
		const ConEdge3d& e = edges[v0.edgesOffset + j];
		if (e.v[1] == q.v[1] || e.v[1] == q.v[3]) { outEdgeIdxs[found++] = v0.edgesOffset + j; }
	}
	// Now the problem is the verts are out of order
	{
		uint32_t vi1 = q.v[1] < q.v[2] ? q.v[1] : q.v[2];
		uint32_t vi2 = q.v[1] > q.v[2] ? q.v[1] : q.v[2];
		const ConVert3d& v1 = verts[vi1];
		for (uint32_t j = 0; j < v1.edgeCount; j++) {
			const ConEdge3d& e = edges[v1.edgesOffset + j];
			if (e.v[1] == vi2) { outEdgeIdxs[found++] = v1.edgesOffset + j; }
		}
	}
	{
		uint32_t vi1 = q.v[2] < q.v[3] ? q.v[2] : q.v[3];
		uint32_t vi2 = q.v[2] > q.v[3] ? q.v[2] : q.v[3];
		const ConVert3d& v1 = verts[vi1];
		for (uint32_t j = 0; j < v1.edgeCount; j++) {
			const ConEdge3d& e = edges[v1.edgesOffset + j];
			if (e.v[1] == vi2) { outEdgeIdxs[found++] = v1.edgesOffset + j; }
		}
	}

	assert(found == 4);
}

ConOneRing3d* Connectivity3d::AllocOneRing(struct Allocator* alloc, uint32_t vIdx) const {
	// Find the unique verts in the surrounding tets and hexs and generate local indexing
	const ConVert3d& v = verts[vIdx];
	uint32_t ringVerts[256];
	uint32_t ringVertCount = 0;
	ringVerts[ringVertCount++] = vIdx;
	uint8_t idxs[2048];
	uint32_t idxCount = 0;
	for (uint32_t i = 0; i < v.tetCount; i++) {
		const ConTet3d& t = tets[tetIdxs[v.tetsOffset + i]];
		for (uint32_t j = 0; j < 4; j++) {
			uint32_t tetVertIdx = ringVertCount;
			for (uint32_t k = 0; k < ringVertCount; k++) {
				if (t.v[j] == ringVerts[k]) { tetVertIdx = k; break; }
			}
			if (tetVertIdx == ringVertCount) {
				ringVerts[ringVertCount++] = t.v[j];
			}
			idxs[idxCount++] = tetVertIdx;
		}
	}
	for (uint32_t i = 0; i < v.hexCount; i++) {
		const ConHex3d& h = hexs[hexIdxs[v.hexsOffset + i]];
		for (uint32_t j = 0; j < 8; j++) {
			uint32_t hexVertIdx = ringVertCount;
			for (uint32_t k = 0; k < ringVertCount; k++) {
				if (h.v[j] == ringVerts[k]) { hexVertIdx = k; break; }
			}
			if (hexVertIdx == ringVertCount) {
				ringVerts[ringVertCount++] = h.v[j];
			}
			idxs[idxCount++] = hexVertIdx;
		}
	}

	// Alloc and init the one ring struct, and then allocate arrays off the end of it
	ConOneRing3d* r = alloc->New<ConOneRing3d>();
	r->vertCount = ringVertCount;
	r->tetCount = v.tetCount;
	r->hexCount = v.hexCount;
	r->idxCount = idxCount;
	uint32_t* rverts = alloc->Alloc<uint32_t>(r->vertCount);
	uint32_t* rtets = alloc->Alloc<uint32_t>(r->tetCount);
	uint32_t* rhexs = alloc->Alloc<uint32_t>(r->hexCount);
	uint8_t* ridxs = alloc->Alloc<uint8_t>(r->idxCount);
	for (uint32_t i = 0; i < r->vertCount; i++) { rverts[i] = ringVerts[i]; }
	for (uint32_t i = 0; i < r->tetCount; i++) { rtets[i] = tetIdxs[v.tetsOffset + i]; }
	for (uint32_t i = 0; i < r->hexCount; i++) { rhexs[i] = hexIdxs[v.hexsOffset + i]; }
	for (uint32_t i = 0; i < r->idxCount; i++) { ridxs[i] = idxs[i]; }

	return r;
}

//-----------------------------------------------------------------------------
// Helper functions
void GmshPrintTriBorder(struct Allocator* alloc, const float* nodeData, uint32_t nodeDataCount, const uint32_t* idxData, uint32_t idxDataCount) {
	Connectivity2d con;
	con.Init(alloc, idxData, idxDataCount);

	// Trace the border of the mesh (assuming no holes)
	uint32_t borderVerts[1024];
	uint32_t borderVertCount = 0;

	// Grab a list of all edges on the border
	uint32_t borderEdges[1024] = { 0 };
	uint32_t borderEdgeCount = 0;
	for (uint32_t i = 0; i < con.edgeCount; i++) {
		if (con.edges[i].p[1] == CON_NULL_IDX) { borderEdges[borderEdgeCount++] = i; }
	}

	// Sort them so one connects to the next
	borderVerts[borderVertCount++] = con.edges[borderEdges[0]].v[0];
	borderVerts[borderVertCount++] = con.edges[borderEdges[0]].v[1];
	while (borderVertCount < borderEdgeCount) {
		uint32_t lastVert = borderVerts[borderVertCount - 1];
		uint32_t lastLastVert = borderVerts[borderVertCount - 2];
		for (uint32_t i = 1; i < borderEdgeCount; i++) {
			const auto& e = con.edges[borderEdges[i]];
			if (e.v[0] == lastVert && e.v[1] != lastLastVert) { borderVerts[borderVertCount++] = e.v[1]; break; }
			if (e.v[1] == lastVert && e.v[0] != lastLastVert) { borderVerts[borderVertCount++] = e.v[0]; break; }
		}
	}

	// Print the gmsh commands for building this border polygon
	for (uint32_t i = 0; i < borderVertCount; i++) {
		float x = nodeData[2 * borderVerts[i] + 0];
		float y = nodeData[2 * borderVerts[i] + 1];
		printf("Point(%d) = {%f, %f, 0, 1.0};\n", (i + 1), x, y);
	}
	for (uint32_t i = 0; i < borderVertCount; i++) {
		printf("Line(%d) = {%d, %d};\n", (i + 1), (i + 1), ((i + 1) % borderVertCount) + 1);
	}
	printf("Curve Loop(1) = {");
	for (uint32_t i = 0; i < borderVertCount; i++) {
		printf("%d%s", (i + 1), i < borderVertCount - 1 ? ", " : "");
	}
	printf("};\n");
	printf("Plane Surface(1) = {1};\n");
}
