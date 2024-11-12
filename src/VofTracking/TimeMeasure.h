#pragma once

#ifndef MEASURE_TIME

#define TIME_MEASURE_RESET(C)
#define TIME_MEASURE_START(N)
#define TIME_MEASURE_END(N)
#define TIME_MEASURE_SYNC_START(N)
#define TIME_MEASURE_SYNC_END(N)
#define TIME_MEASURE_BARRIER()
#define TIME_MEASURE_JSON(J)

#else

#define TIME_MEASURE_RESET(C) TimeMeasure::getInstance().reset(C)
#define TIME_MEASURE_START(N) TimeMeasure::getInstance().start(N)
#define TIME_MEASURE_END(N) TimeMeasure::getInstance().end(N)
#define TIME_MEASURE_SYNC_START(N) TimeMeasure::getInstance().syncStart(N)
#define TIME_MEASURE_SYNC_END(N) TimeMeasure::getInstance().syncEnd(N)
#define TIME_MEASURE_BARRIER() TimeMeasure::getInstance().barrier()
#define TIME_MEASURE_JSON(J) TimeMeasure::getInstance().toJson(J)

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <nlohmann/json.hpp>
#include <vtkMPIController.h>

class TimeMeasure {
public:
    using json = nlohmann::json;

    enum class Name {
        All = 0,
        Seeding,
        Advection,
        Advection_Init,
        Advection_Advection,
        Advection_OOB,
        Advection_Corrector,
        Advection_PhaseChange,
        Advection_Exchange,
        Components,
        ParticleLabeling,
        SeedLabeling,
        Boundary,
        Output,
    };

    static inline TimeMeasure& getInstance() {
        if (instance_ == nullptr) {
            instance_ = std::unique_ptr<TimeMeasure>(new TimeMeasure());
        }
        return *instance_;
    }

    inline void reset(vtkMPIController* mpiController = nullptr) {
        mpiController_ = mpiController;
        // Store a high precision time stamp for measuring nanoseconds relative to this and a system clock to print
        // actual date and time of measurement start. The actual date and time is meant as metadata of the measurement
        // and should not be combined with the relative nanoseconds. Therefore, offset between them is ignored here.
        barrier();
        resetTime_ = clock_type::now();
        resetTimeSystem_ = std::chrono::system_clock::now();
        start_.clear();
        end_.clear();
    }

    inline void barrier() {
        if (mpiController_ != nullptr && mpiController_->GetCommunicator() != nullptr) {
            mpiController_->Barrier();
        }
    }

    inline void start(Name name) {
        start_[name].push_back(clock_type::now());
    }

    inline void end(Name name) {
        end_[name].push_back(clock_type::now());
    }

    inline void syncStart(Name name) {
        barrier();
        start(name);
    }

    inline void syncEnd(Name name) {
        barrier();
        end(name);
    }

    void toJson(json& j) {
        std::stringstream s;
        auto time = std::chrono::system_clock::to_time_t(resetTimeSystem_);
        s << std::put_time(std::localtime(&time), "%FT%T%z");
        j["ResetTime"] = s.str();

        json timings;
        timings["All"] = getTimingsJson(Name::All);
        timings["Seeding"] = getTimingsJson(Name::Seeding);
        timings["Advection"] = getTimingsJson(Name::Advection);
        timings["AdvectionInit"] = getTimingsJson(Name::Advection_Init);
        timings["AdvectionAdvection"] = getTimingsJson(Name::Advection_Advection);
        timings["AdvectionOOB"] = getTimingsJson(Name::Advection_OOB);
        timings["AdvectionCorrector"] = getTimingsJson(Name::Advection_Corrector);
        timings["AdvectionPhaseChange"] = getTimingsJson(Name::Advection_PhaseChange);
        timings["AdvectionExchange"] = getTimingsJson(Name::Advection_Exchange);
        timings["Components"] = getTimingsJson(Name::Components);
        timings["ParticleLabeling"] = getTimingsJson(Name::ParticleLabeling);
        timings["SeedLabeling"] = getTimingsJson(Name::SeedLabeling);
        timings["Boundary"] = getTimingsJson(Name::Boundary);
        timings["Output"] = getTimingsJson(Name::Output);
        j["Timings"] = timings;
    }

private:
    typedef std::chrono::high_resolution_clock clock_type;
    typedef std::unordered_map<Name, std::vector<clock_type::time_point>> map_type;

    TimeMeasure() : mpiController_(nullptr) {
        reset();
    }

    inline std::vector<std::uint64_t> nanoSecondsSinceReset(const std::vector<clock_type::time_point>& time) const {
        std::vector<uint64_t> result(time.size());
        for (std::size_t i = 0; i < time.size(); i++) {
            result[i] = std::chrono::duration_cast<std::chrono::nanoseconds>(time[i] - resetTime_).count();
        }
        return result;
    }

    inline json getTimingsJson(Name name) const {
        auto start_it = start_.find(name);
        auto end_it = end_.find(name);
        if (start_it == start_.end() || end_it == end_.end()) {
            return {};
        }
        auto start = start_it->second;
        auto end = end_it->second;
        if (start.size() != end.size()) {
            throw std::runtime_error("Invalid time measure!");
        }
        return {
            {"start", nanoSecondsSinceReset(start)},
            {"end", nanoSecondsSinceReset(end)},
        };
    }

    static inline std::unique_ptr<TimeMeasure> instance_ = nullptr;
    vtkMPIController* mpiController_;
    clock_type::time_point resetTime_;
    std::chrono::system_clock::time_point resetTimeSystem_;
    map_type start_;
    map_type end_;
};

#endif
