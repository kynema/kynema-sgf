#ifndef GODUNOV_H
#define GODUNOV_H

#include <cstdint>

namespace godunov {

enum class scheme : std::uint8_t {
    PLM,
    PPM,
    PPM_NOLIM,
    BDS,
    WENO_JS,
    WENOZ,
    MINMOD,
    UPWIND
};

} // namespace godunov

#endif /* Godunov_H */
