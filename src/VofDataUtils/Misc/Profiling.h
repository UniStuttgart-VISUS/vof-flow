#pragma once

#ifdef VOFFLOW_USE_TRACY
#include <tracy/Tracy.hpp>
#else

// Use empty macro definitions from Tracy.hpp
#define ZoneScoped
#define ZoneName(x, y)

#endif
