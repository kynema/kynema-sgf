#ifndef SGFTIOGAIFACE_H
#define SGFTIOGAIFACE_H

#include <memory>
#include <string>
#include <vector>

namespace kynema_sgf {
class CFDSim;
}

namespace TIOGA {
class tioga;
struct AMRMeshInfo;
} // namespace TIOGA

namespace kynema_sgf {

class SGFTiogaIface
{
public:
    // Modular coupling TODO: let the SGF-owned module interface own this
    // optional bridge and its registered field lists. Export these hooks for
    // driver-coordinated connectivity/interpolation; do not perform a second
    // exchange inside advance_states(). Keep JSON inputs/outputs independent
    // of Tioga and usable in standalone mode. The bridge borrows CFDSim and
    // Tioga, so their lifetimes must outlast it.
    SGFTiogaIface(kynema_sgf::CFDSim&, TIOGA::tioga& tg);
    ~SGFTiogaIface();

    void pre_overset_conn_work();

    void post_overset_conn_work();

    void register_mesh();

    void register_solution(
        const std::vector<std::string>& cell_vars,
        const std::vector<std::string>& node_vars);

    void update_solution();

private:
    kynema_sgf::CFDSim& m_sim;
    TIOGA::tioga& m_tg;

    std::unique_ptr<TIOGA::AMRMeshInfo> m_info;
};

} // namespace kynema_sgf

#endif /* SGFTIOGAIFACE_H */
