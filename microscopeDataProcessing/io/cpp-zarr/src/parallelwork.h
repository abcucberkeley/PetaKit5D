// How the parallel reader and writer share out their chunks and set up their
// threads (internal to the library, not installed)
#ifndef CPPZARR_PARALLELWORK_H
#define CPPZARR_PARALLELWORK_H
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>
#ifdef __linux__
#include <sys/syscall.h>
#include <unistd.h>
#endif

// Contiguous blocks of work items with work stealing. Each thread starts with its
// own contiguous block of the items, so neighboring chunks stay with one thread and
// its writes to the array stay together (and on the thread's NUMA node); a thread
// that runs out takes the back half of the largest block left. Fixed blocks left up
// to a third of the threads' time idle when chunks took different times to decode or
// compress, and handing chunks out one at a time from a shared counter split
// neighboring chunks between threads, which made C-order reads up to 38% slower.
class WorkRanges {
    struct alignas(64) Range {
        std::mutex m;
        std::atomic<int64_t> begin{0}, end{0};
    };
    std::vector<Range> ranges;
public:
    WorkRanges(const int threads, const int64_t items) : ranges(threads > 0 ? threads : 1) {
        const int64_t t = (int64_t)ranges.size();
        for(int64_t i = 0; i < t; i++){
            ranges[i].begin.store(items*i/t, std::memory_order_relaxed);
            ranges[i].end.store(items*(i+1)/t, std::memory_order_relaxed);
        }
    }
    // The item thread t takes next from its own block, or -1 (another thread may still
    // take it first)
    int64_t peek(const int t) const {
        const int64_t b = ranges[t].begin.load(std::memory_order_relaxed);
        return b < ranges[t].end.load(std::memory_order_relaxed) ? b : -1;
    }
    // The next item for thread t, or -1 once none is left anywhere
    int64_t next(const int t){
        {
            std::lock_guard<std::mutex> g(ranges[t].m);
            const int64_t b = ranges[t].begin.load(std::memory_order_relaxed);
            if(b < ranges[t].end.load(std::memory_order_relaxed)){
                ranges[t].begin.store(b+1, std::memory_order_relaxed);
                return b;
            }
        }
        while(true){
            int victim = -1;
            int64_t most = 0;
            for(int v = 0; v < (int)ranges.size(); v++){
                const int64_t left = ranges[v].end.load(std::memory_order_relaxed)-ranges[v].begin.load(std::memory_order_relaxed);
                if(left > most){
                    most = left;
                    victim = v;
                }
            }
            if(victim < 0) return -1;
            int64_t first, last;
            {
                std::lock_guard<std::mutex> g(ranges[victim].m);
                const int64_t b = ranges[victim].begin.load(std::memory_order_relaxed);
                last = ranges[victim].end.load(std::memory_order_relaxed);
                if(last <= b) continue;
                first = last-(last-b+1)/2;
                ranges[victim].end.store(first, std::memory_order_relaxed);
            }
            {
                std::lock_guard<std::mutex> g(ranges[t].m);
                ranges[t].begin.store(first+1, std::memory_order_relaxed);
                ranges[t].end.store(last, std::memory_order_relaxed);
            }
            return first;
        }
    }
};

// A pool of the buffers and coders a thread uses (T, with a shrink() that lets go of
// anything too big to keep), kept between calls. Allocating them in every call had
// glibc grow and trim each thread's heap in every call, and small calls spent most of
// their time faulting that memory back in and in mprotect and munmap calls, which
// stall every thread of the process. A thread gets back the ones it gave back last
// (place: its thread number) unless another call took them: handed out in any order,
// they were in other cores' caches and on the other socket's memory, which made small
// writes up to 8% slower.
template <typename T>
class ScratchPool {
    std::mutex m;
    std::vector<std::unique_ptr<T>> byPlace, idle;
public:
    // place < 0: any
    std::unique_ptr<T> take(const int place){
        {
            std::lock_guard<std::mutex> g(m);
            if(place >= 0 && place < (int)byPlace.size() && byPlace[place]) return std::move(byPlace[place]);
            if(!idle.empty()){
                std::unique_ptr<T> s = std::move(idle.back());
                idle.pop_back();
                return s;
            }
            for(auto &b : byPlace) if(b) return std::move(b);
        }
        return std::unique_ptr<T>(new T());
    }
    void give(const int place, std::unique_ptr<T> s){
        if(!s) return;
        s->shrink();
        std::lock_guard<std::mutex> g(m);
        if(place >= 0){
            if(place >= (int)byPlace.size()) byPlace.resize(place+1);
            if(!byPlace[place]){
                byPlace[place] = std::move(s);
                return;
            }
        }
        idle.push_back(std::move(s));
    }
};

// Wait until done() is true (and return true), or until stop() is (and return false).
// The thread yields its CPU at first, then sleeps a little longer each time up to 0.16
// ms: a long wait leaves the CPU to others (the thread being waited for among them, when
// the system has paused it), while a short one ends soon after.
template <typename Done, typename Stop>
bool waitFor(const Done &done, const Stop &stop){
    for(int i = 0; !done(); i++){
        if(stop()) return false;
        if(i < 32) std::this_thread::yield();
        else std::this_thread::sleep_for(std::chrono::microseconds(10 << std::min(i-32, 4)));
    }
    return true;
}

// For the length of a call, a thread with the default memory policy gets an explicit
// local one: its memory still comes from its own NUMA node, but Linux's automatic NUMA
// balancing skips threads with an explicit policy. Otherwise long-running threads (as
// in MATLAB) have the kernel mark the process's pages, and writes took up to millions
// of hinting faults and page migrations a call (8-30% of their time). A thread with a
// policy of its own (numactl, say) keeps it.
class LocalMemoryPolicy {
#if defined(__linux__) && defined(SYS_get_mempolicy) && defined(SYS_set_mempolicy)
    bool set = false;
public:
    LocalMemoryPolicy(){
        int mode = -1;
        if(syscall(SYS_get_mempolicy, &mode, nullptr, 0UL, nullptr, 0UL) == 0 && mode == 0 /* MPOL_DEFAULT */){
            set = syscall(SYS_set_mempolicy, 4 /* MPOL_LOCAL */, nullptr, 0UL) == 0;
        }
    }
    ~LocalMemoryPolicy(){
        if(set) syscall(SYS_set_mempolicy, 0 /* MPOL_DEFAULT */, nullptr, 0UL);
    }
#else
public:
    LocalMemoryPolicy(){}
#endif
    LocalMemoryPolicy(const LocalMemoryPolicy&) = delete;
    LocalMemoryPolicy& operator=(const LocalMemoryPolicy&) = delete;
};

#endif
