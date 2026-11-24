#ifndef KYNEMA_FMB_WRAPPER_H
#define KYNEMA_FMB_WRAPPER_H

#ifdef KYNEMA_SGF_USE_KYNEMA_FMB
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

#endif /* KYNEMA_FMB_WRAPPER_H */
