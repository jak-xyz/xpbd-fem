#include "Allocator.h"
#include <assert.h>

//-----------------------------------------------------------------------------
void* Allocator::Alloc(size_t align, size_t size) {
	uint8_t* p = (uint8_t*)((((size_t)head + align - 1) / align) * align);
	head = p + size;
	assert((size_t)(head - mem) <= capacity);
	return p;
}

void Allocator::Initialize(void* buffer, size_t size) {
	mem = (uint8_t*)buffer;
	head = mem;
	capacity = size;
}
