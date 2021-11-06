//-----------------------------------------------------------------------------
// Functions for generating node and index data that can be fed to Geo to
// create various shapes.
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include "Allocator.h"
#include "Settings.h"
#include "vectormath.h"

//-----------------------------------------------------------------------------
void GenerateArmadillo(uint32_t elementType, const float** outNodeData, uint32_t* outNodeDataCount, const uint32_t** outIdxData, uint32_t* outIdxDataCount);
void GenerateBlock(uint32_t elementType, Allocator* alloc, uint32_t width, uint32_t height, vec2 scale, uint32_t elementPattern, float wonkiness, const float** outNodeData, uint32_t* outNodeDataCount, const uint32_t** outIdxData, uint32_t* outIdxDataCount);
