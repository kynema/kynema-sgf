#include "src/equation_systems/source_terms/DampingLayerSource.H"

namespace kynema_sgf::pde {

template class DampingLayerSource<MomentumSource>;
template class DampingLayerSource<TemperatureSource>;
template class DampingLayerSource<DensitySource>;
template class DampingLayerSource<TKESource>;
template class DampingLayerSource<SDRSource>;

// Currently, cannot identify field for generic SourceTerm, so there is no
// specialization for it. SourceTerm is associated with passive_scalar,
// levelset, and vof. To apply a damping layer to these fields, they would
// need a specifically named source term instead or this template would need
// more information to discern the field name.

} // namespace kynema_sgf::pde