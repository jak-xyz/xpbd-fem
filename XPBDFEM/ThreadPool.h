//-----------------------------------------------------------------------------
// Quick and dirty thread pool to begin testing multithreading physics.
// Using pImpl idiom, to avoid contaminating the rest of the project with
// giant std headers.
//-----------------------------------------------------------------------------
#pragma once

//-----------------------------------------------------------------------------
#include <stdint.h>

//-----------------------------------------------------------------------------
class ThreadPool {
public:
	void Initialize();
	void Terminate();

	template <typename Fn, class... Args>
	void Kick(Fn& l, Args... args);

	uint32_t GetThreadCount() const;

private:
	void KickTask(struct Task* task);

	alignas(16) uint8_t data[1024];
	struct ThreadPoolImpl* impl = 0;
};

struct Task {
	virtual void Do(uint32_t slot) = 0;
};

template <typename Fn, class... Args>
void ThreadPool::Kick(Fn& l, Args... args) {
	if (!impl) { return; }

	auto ll = [&](uint32_t slot) { l(slot, args...); };
	struct CurrentTask : Task {
		decltype(ll)& ll;
		CurrentTask(decltype(ll)& ll_) : ll(ll_) {}
		virtual void Do(uint32_t slot) {
			ll(slot);
		}
	} currentTask(ll);
	KickTask(&currentTask);
}
