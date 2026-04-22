.. _multiphase:

Multiphase flow modeling
------------------------

Kynema-SGF employs the volume-of-fluid method for simulating two-phase (water-air) flows.
More specifically, the volume fraction field is advected explicitly using a
directional split geometric approach, and the advection of momentum is
discretized in a mass-consistent manner. Overall, this approach conserves mass
and momentum while remaining stable at high density ratios (typically 1000).
Viscosities can be specified for each fluid independently, but surface tension
is not modeled by Kynema-SGF currently. For further detail, see
`Kuhn, Deskos, Sprague (Computers & Fluids 2023)
<https://doi.org/10.1016/j.compfluid.2022.105770>`_.
