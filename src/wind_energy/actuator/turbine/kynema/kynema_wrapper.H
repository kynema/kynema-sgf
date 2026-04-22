#ifndef KYNEMA_WRAPPER_H
#define KYNEMA_WRAPPER_H

#ifdef AMR_WIND_USE_KYNEMA
#include <array>
#include <numbers>
#include <ranges>
#include <yaml-cpp/yaml.h>

#include "elements/beams/hollow_circle_properties.hpp"
#include "interfaces/components/aerodynamics.hpp"
#include "interfaces/components/inflow.hpp"
#include "interfaces/turbine/turbine_interface.hpp"
#include "interfaces/turbine/turbine_interface_builder.hpp"
#else
#endif

#endif /* KYNEMA_WRAPPER_H */
