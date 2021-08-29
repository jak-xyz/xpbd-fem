#include "ThreadPool.h"
#include <condition_variable>
#include <mutex>
#include <new>
#include <thread>

//-----------------------------------------------------------------------------

struct ThreadPoolImpl {
	ThreadPoolImpl() {
		threadCount = std::thread::hardware_concurrency();
		threads = new std::thread[threadCount];
		kickedFlags = new bool[threadCount];

		// Launch all the threads, but they'll immediately just start waiting to be kicked
		for (uint32_t i = 0; i < threadCount; i++) {
			kickedFlags[i] = false;
			threads[i] = std::thread(ThreadMainStub, this, i);
		}
	}

	~ThreadPoolImpl() {
		// Wake up all the tasks threads so they can die
		std::unique_lock<std::mutex> kickLock(kickMutex);
		die = true;
		for (uint32_t i = 0; i < threadCount; i++) { kickedFlags[i] = true; }
		kickLock.unlock();
		kickCv.notify_all();

		// Now just wait for the threads to end, so we can delete them
		for (uint32_t i = 0; i < threadCount; i++) {
			threads[i].join();
		}
		delete[] kickedFlags;
		delete[] threads;
	}

	void Kick(Task* kickedTask) {
		std::unique_lock<std::mutex> doneLock(doneMutex);
		doneCount = 0;

		std::unique_lock<std::mutex> kickLock(kickMutex);
		task = kickedTask;
		for (uint32_t i = 0; i < threadCount; i++) { kickedFlags[i] = true; }
		kickLock.unlock();
		kickCv.notify_all();

		// Block until all the workers finish
		doneCv.wait(doneLock, [&] { return doneCount == threadCount; });
		doneLock.unlock();
	}

	static uint32_t ThreadMainStub(ThreadPoolImpl* self, uint32_t slot) { return self->ThreadMain(slot); }
	uint32_t ThreadMain(uint32_t slot) {
		while (true) {
			// Mostly, we just wait around to be kicked
			std::unique_lock<std::mutex> kickLock(kickMutex);
			kickCv.wait(kickLock, [&] { return kickedFlags[slot]; });
			kickedFlags[slot] = false;
			kickLock.unlock();

			// See if we were just woken up to exit
			if (die) { return 0; }

			// Otherwise, run our task to completion
			task->Do(slot);

			// If we're the last task, wake up the kicking thread
			std::unique_lock<std::mutex> doneLock(doneMutex);
			++doneCount;
			doneLock.unlock();
			if (doneCount == threadCount) {
				doneCv.notify_one();
			}
		}
	}

	std::thread* threads = 0;
	uint32_t threadCount = 0;

	std::mutex kickMutex;
	std::condition_variable kickCv;
	bool* kickedFlags = 0;

	Task* task = 0;

	std::mutex doneMutex;
	std::condition_variable doneCv;
	uint32_t doneCount = 0;

	bool die = false;
};

//-----------------------------------------------------------------------------
void ThreadPool::Initialize() {
	static_assert(sizeof(data) >= sizeof(ThreadPoolImpl));
	impl = new (data) ThreadPoolImpl();
}

void ThreadPool::Terminate() {
	if (impl) { impl->~ThreadPoolImpl(); }
}

uint32_t ThreadPool::GetThreadCount() const {
	return impl ? impl->threadCount : 0;
}

void ThreadPool::KickTask(struct Task* task) {
	if (impl) { impl->Kick(task); }
}
