//-----------------------------------------------------------------------------
// Push allocator
//-----------------------------------------------------------------------------
#pragma once

#include <stdint.h>

//-----------------------------------------------------------------------------
struct Allocator {
	void* Alloc(size_t align, size_t size);

	template <typename T>
	T* Alloc(size_t count) { return (T*)Alloc(alignof(T), sizeof(T) * count); }

	template <typename T>
	T* New() { return new (Alloc<T>(1)) T(); }

	void Initialize(void* buffer, size_t size);
	
	uint8_t* mem;
	uint8_t* head;
	size_t capacity;
};
