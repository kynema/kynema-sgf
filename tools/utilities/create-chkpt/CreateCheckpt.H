#ifndef CREATECHECKPT_H
#define CREATECHECKPT_H

#include "src/incflo.H"

namespace kynema_sgf {
namespace tools {

class CreateCheckpt : public incflo
{
public:
    CreateCheckpt();

    virtual ~CreateCheckpt() = default;

    void interpolate();

    void read_fvcom_data();

    void run_utility();

private:

};

} // namespace tools
} // namespace kynema_sgf

#endif /* CREATECHECKPT_H */
