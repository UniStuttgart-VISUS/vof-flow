#pragma once

namespace VofFlow {
    class ArrayNames {
    public:
        static constexpr const char* POINTS = "Points"; // Name as used by vtkPoints internally.
        static constexpr const char* PHASE_IDX = "Phase";
        static constexpr const char* SEED_IDX = "SeedIdx";
        static constexpr const char* LABELS = "Labels";
        static constexpr const char* ID = "Id";
        static constexpr const char* PROCESS_ID = "ProcessId";
        static constexpr const char* CHANGED_PHASE_STEP = "ChangedPhaseStep";
        static constexpr const char* UNCERTAINTY = "Uncertainty";
        static constexpr const char* TARGET_PHASE = "TargetPhase";
        static constexpr const char* SAME_END_PHASE_RATIO = "SameEndPhaseRatio";
        static constexpr const char* STAYED_IN_PHASE_RATIO = "StayedInPhaseRatio";
    };

    class MPITags {
    public:
        static constexpr int NEIGHBOR_SEEDS = 101;
        static constexpr int COMPONENT_EXCHANGE_SIZE = 110;
        static constexpr int COMPONENT_EXCHANGE_DATA = 111;
        static constexpr int COMPONENT_EXCHANGE_NUM_LABELS = 112;
        static constexpr int COMPONENT_EXCHANGE_PAIR_LIST = 113;
        static constexpr int COMPONENT_EXCHANGE_MAPPING = 114;
        static constexpr int PARTICLE_EXCHANGE = 120;
        static constexpr int PARTICLE_TO_SEED_ID = 150;
        static constexpr int PARTICLE_TO_SEED_POS = 151;
        static constexpr int PARTICLE_TO_SEED_LABEL = 152;
        static constexpr int PARTICLE_TO_SEED_CHANGED_PHASE_STEP = 153;
        static constexpr int PARTICLE_TO_SEED_UNCERTAINTY = 154;
        static constexpr int PARTICLE_TO_SEED_TARGET_PHASE = 155;
        static constexpr int PARTICLE_TO_SEED_OUT_OF_BOUNDS_ID = 156;
        static constexpr int PARTICLE_TO_SEED_OUT_OF_BOUNDS_CHANGED_PHASE_STEP = 157;
    };

    class ErrorLabels {
    public:
        static constexpr int EMPTY_COMPONENT = -1;
        static constexpr int BAD_COMPONENT_MAP = -2;
        static constexpr int BAD_PARTICLE = -3;
        static constexpr int PARTICLE_OUT_OF_BOUNDS = -4;
        static constexpr int BAD_SEED = -5;

        // Always keep last and below all other values, is used for boundary grid init.
        static constexpr int MAX_ERROR_LABEL = -6;
    };
} // namespace VofFlow
