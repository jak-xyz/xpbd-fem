//------------------------------------------------------------------------------
// debugout function, a minimal WASM printf alternative, as well as the
// DVAL function, used to quickly output a list of values.
//------------------------------------------------------------------------------
#pragma once

#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define DEBUG_OUT 1

inline void debugout(const char* format, ...) {
	char buffer[1024];
	va_list args;
	va_start(args, format);
	vsnprintf(buffer, sizeof(buffer), format, args);
	puts(buffer);
	va_end(args);
}

namespace DVAL_ADL {
	struct Count {
		Count(size_t count_) : count(count_) {}
		size_t count;
		operator size_t() const { return count; }
	};
}

namespace DVAL_ADL {
	inline int Stringify(char* buffer, size_t count, bool b) { return snprintf(buffer, count, "%s", b ? "true" : "false"); }
	inline int Stringify(char* buffer, size_t count, char c) { return snprintf(buffer, count, "%c", c); }
	inline int Stringify(char* buffer, size_t count, int16_t i) { return snprintf(buffer, count, "%i", i); }
	inline int Stringify(char* buffer, size_t count, int32_t i) { return snprintf(buffer, count, "%i", i); }
	inline int Stringify(char* buffer, size_t count, uint8_t u) { return snprintf(buffer, count, "%u", u); }
	inline int Stringify(char* buffer, size_t count, uint16_t u) { return snprintf(buffer, count, "%u", u); }
	inline int Stringify(char* buffer, size_t count, uint32_t u) { return snprintf(buffer, count, "%u", u); }
	inline int Stringify(char* buffer, size_t count, float f) { return snprintf(buffer, count, "%f", f); }
	inline int Stringify(char* buffer, size_t count, double f) { return snprintf(buffer, count, "%f", f); }
	inline int Stringify(char* buffer, size_t count, const char* s) { return snprintf(buffer, count, "%s", s); }
}

template <typename T>
inline size_t DVALFmtKernel(char* buffer, size_t count, T t) {
	int used = Stringify(buffer, DVAL_ADL::Count(count), t);
	return (size_t)(used > 0 ? used : 0);
}
template <typename T, typename... Targs>
inline size_t DVALFmtKernel(char* buffer, size_t count, T t, Targs... args) {
	size_t used = DVALFmtKernel(buffer, count, t);
	used = used < count ? used : count;
	if (used < count) { buffer[used++] = '\t'; }
	return used + DVALFmtKernel(buffer + used, count - used, args...);
}

template <typename... Targs>
inline void DVAL(Targs ... args) {
	char buffer[1024];
	size_t used = DVALFmtKernel(buffer, sizeof(buffer), args...);
	if (used >= sizeof(buffer)) { used = sizeof(buffer) - 1; }
	buffer[used++] = '\0';
	puts(buffer);
}

#define DVAL_ONCE(...) \
	static bool once##__LINE__ = false; \
	if (!once##__LINE__) { \
		once##__LINE__ = true; \
		DVAL(__VA_ARGS__); \
	} while(0)
