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
