.. _alm_fllc:

Filtered lifting line correction (FLLC)
=======================================

To accurately represent wind turbine using the actuator line model (ALM) with OpenFAST, the user must select ALM parameters that are appropriate for the turbine of interest and for the near-turbine grid resolution.

The standard ALM has uniformly spaced points and thus requires a fine value of epsilon and a large number of points along the blade to properly model it. In recent years, a correction, called filtered lifting line correction (FLLC), has been developed :cite:p:`martinez-tossas_meneveau_2019`. Such a correction also includes the use of points that are not uniformly spaced, thus lifting the requirement of small values of epsilon. It is recommended that you run the ALM with the FLLC on. 

To enable the ALM with OpenFAST coupling, ``Actuator`` should be added to ``incflo.physics``, ``ActuatorForcing`` should be added to ``ICNS.source_terms``, and ``Actuator.type = TurbineFastLine`` should be set. Next, the user should enable FLLC, choose its type, and set the option for a non-uniform point distribution:

.. code-block:: none

    Actuator.TurbineFastLine.fllc = 1
    Actuator.TurbineFastLine.fllc_nonuniform = 1
    Actuator.TurbineFastLine.fllc_type = variable_chord

When selecting the ``Actuator.TurbineFastLine.fllc_nonuniform`` as ``1``, a new distribution of points is calculated internally, based on the ``Actuator.TurbineFastLine.num_points_blade`` and ``Actuator.TurbineFastLine.fllc_epsilon_dr_ratio`` entries. The ``Actuator.TurbineFastLine.fllc_epsilon_dr_ratio`` value should be at least 1 with a maximum value of 3. It is recommended to use the value 3. This parameter controls the distribution of points, and its value should be selected based on the desired accuracy of the correction, according to Table 1 given in :cite:`MartinezFLLCimpl`.

.. code-block:: none

    Actuator.TurbineFastLine.fllc_epsilon_dr_ratio = 3

.. note::

    The entries ``Actuator.TurbineFastLine.num_points_blade`` and ``Actuator.TurbineFastLine.num_points_tower`` should match the entries ``NumBlNds`` from the AeroDyn blade file  and ``NumTwrNds`` from the AeroDyn input file.

Two values of epsilon should be set for the FLLC: the regular epsilon and an epsilon for the chord. The regular ``epsilon`` should be set to twice the grid resolution. If you have a very fine grid near the turbine, e.g. 0.625 m, then you can set ``epsilon`` to be 3 to 4 times the local grid resolution. The smaller the epsilon value, the more sensitive to the resolution the ALM will be. The ``epsilon_chord`` parameter represents the ideal epsilon value if FLLC were not to be used. However, we set both epsilon parameters for FLLC. One of the key details of the correction is that it uses the ``epsilon_chord`` on the non-uniform points, and this should be set to a target value regardless of the underlying grid resolution. The target values again depend on the desired accuracy of the correction and guidelines are also given on Table 1 of :cite:`MartinezFLLCimpl`. The correction itself uses one value of each input: the first value of ``epsilon``, the first value of ``epsilon_min`` when that input is given (the force is spread with the larger of the two, so the correction takes the larger of the two as the width that the grid resolves) and the second value of ``epsilon_chord`` (see :input_param:`Actuator.TurbineFastLine.fllc`). Give each of them the same value in all three directions.

.. dropdown:: More details on implementation and warnings

    When the force points and the velocity points differ in number or in position, the FLLC locates every point by its distance along the blade from the first point of its own set. It interpolates the sampled velocity linearly to the force points to form the relative velocity, computes the correction at the force points, and interpolates only the correction back to the velocity points, where it is added to the sampled velocity.
    When the correction starts, a warning is printed for each blade whose force points do not resolve the optimal epsilon: at some force point, ``epsilon_chord`` times the chord is smaller than the larger of the distances to the neighboring force points. The message gives the number of such points and the lowest ratio of the two. The lift is known only at the force points, so the correction then depends on their number: use more points or a larger ``epsilon_chord``. With ``Actuator.TurbineFastLine.fllc_nonuniform`` set to ``0`` the correction is wrong at these points, so keep the default there.
    A second warning is printed when the velocity points are not the force points and do not resolve the optimal epsilon either: at some velocity point, ``epsilon_chord`` times the chord is smaller than the larger of the distances to the neighboring velocity points. The correction is computed at the force points and handed over at the velocity points, so where it varies over the length of the optimal epsilon, coarser velocity points show only a few values of it. This is the case at the tip of a blade whose chord, and with it the lift, does not go to zero there: within a few optimal epsilons of the tip the correction typically turns strongly negative and then positive at the tip point itself. The velocity handed over at the last points then depends on where these points are and on the number of force points. Use more blade nodes near the tip in the external solver, or a larger ``epsilon_chord``.


.. code-block:: none

    Actuator.TurbineFastLine.epsilon       = 10 10 10        # for a near-turbine grid of 5 m
    Actuator.TurbineFastLine.epsilon_chord = 0.25 0.25 0.25  # for any grid resolution

.. tip::

    If your near-turbine grid is not isotropic, then consider your "grid resolution" to be the largest of :math:`\Delta x`,  :math:`\Delta y`, and  :math:`\Delta z`. 

Lastly, in certain cases, a numerical instability can arise from the application of the correction on the first few time steps. When a simulation begins, there is often an initial spike in some of the ALM quantities that affect the correction. Thus, in these scenarios, the user can choose a delay, in seconds, to start the correction. In cases where instabilities are observed, a 5-second delay has been shown to be sufficient:

.. code-block:: none

    Actuator.TurbineFastLine.fllc_start_time = 5

More details on other parameters are available in the description of each entry on the :ref:`input file reference page <inputs_actuator>`.