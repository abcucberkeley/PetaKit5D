// Internal to cpp/fft: this module's own implementations behind SharedApi (backend.cpp, driver.cpp). Everything
// else calls the functions of backend.h, which use the attached state when there is one.
#pragma once

#include <cstdint>
#include <mutex>
#include <string>

#include "backend.h"

namespace pk5d::fft::own {

const Backend& backend(int id);
bool exact_loaded();
void load_exact(const std::string& libf, const std::string& libd);
void unload_exact();
int default_backend();
void set_default_backend(int id);
std::mutex& planner_lock();
void* cached_plan(int be, int prec, int kind, std::int64_t n, int howmany);
void clear_plans(int be);

}  // namespace pk5d::fft::own
