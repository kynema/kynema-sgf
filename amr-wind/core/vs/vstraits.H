#ifndef VSTRAITS_H
#define VSTRAITS_H

#include <limits>

namespace amr_wind::vs {

template <typename T>
struct DTraits
{};

template <>
struct DTraits<int>
{
    static constexpr int zero() { return 0; }
    static constexpr int one() { return 1; }
    static constexpr int max() { return std::numeric_limits<int>::max(); }
    static constexpr int min() { return std::numeric_limits<int>::min(); }
};

template <>
struct DTraits<double>
{
    static constexpr double zero() { return 0.0; }
    static constexpr double one() { return 1.0; }
    static constexpr double max() { return std::numeric_limits<double>::max(); }
    static constexpr double min() { return std::numeric_limits<double>::min(); }
    static constexpr double eps()
    {
        return std::numeric_limits<double>::epsilon();
    }
};

template <>
struct DTraits<float>
{
    static constexpr float zero() { return 0.0F; }
    static constexpr float one() { return 1.0F; }
    static constexpr float max() { return std::numeric_limits<float>::max(); }
    static constexpr float min() { return std::numeric_limits<float>::min(); }
    static constexpr float eps()
    {
        return std::numeric_limits<float>::epsilon();
    }
};

} // namespace amr_wind::vs

#endif /* VSTRAITS_H */
