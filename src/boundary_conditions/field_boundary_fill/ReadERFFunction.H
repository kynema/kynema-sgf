#ifndef READERFFUNCTION_H
#define READERFFUNCTION_H
#include <functional>

class MultiBlockContainer;
namespace kynema_sgf {
class InletData;
}

using ReadERFFunction = std::function<void(
    const amrex::Real time,
    amrex::Vector<amrex::Real>&,
    kynema_sgf::InletData&,
    const amrex::Vector<kynema_sgf::Field*>&,
    MultiBlockContainer*)>;

#endif /* READERFFUNCTION_H */
