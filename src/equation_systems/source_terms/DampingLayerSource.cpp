#include "src/equation_systems/source_terms/DampingLayerSource.H"

namespace kynema_sgf::pde {

template class DampingLayerSource<MomentumSource>;
template class DampingLayerSource<TemperatureSource>;
template class DampingLayerSource<DensitySource>;
template class DampingLayerSource<TKESource>;
template class DampingLayerSource<SDRSource>;
template class DampingLayerSource<SourceTerm>;

} // namespace kynema_sgf::pde