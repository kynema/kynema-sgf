.. _inputs_channel_builder:

Section: ChannelBuilder
~~~~~~~~~~~~~~~~~~~~~~~

This section is for setting parameters for the channel builder,
which is a physics module for initializing channel flows with
different geometries. It can be used in single-phase and multiphase
contexts. It can be used for closed channels,
enabling internal flows, or open channels, enabling free surface flows.
This setup is comprised of defining channel segments.

The channel builder is exclusively for initializing the flow and
does not perform any operation during the progression of the simulation.
The channel builder capability allows for initializing a variety of velocity
profiles, which can be awkward to replicate as standard boundary conditions
for a simulation. As a result, the BoundaryPlane field boundary type is
often best suited for these cases by leveraging the
`BoundaryPlane.output_and_use_initial_plane`
capability.
The channel builder is activated by including `ChannelBuilder` in the list
of physics under :input_param:`incflo.physics`. Finally, the channel builder
sets up channels by defining the terrain blank variable. This variable
must be used in the momentum forcing in order to actually affect the flow.
Currently, this is done by including the `DragForcing` source term in
:input_param:`ICNS.source_terms`.

.. input_param:: ChannelBuilder.segment_labels

   **type:** List of strings, required

   After listing the segment labels, the details of each
   channel segment are defined through additional input arguments,
   which describe the geometry and velocity profile of each segment.
   The example segment label used for the input arguments below is `s1`.

.. input_param:: ChannelBuilder.s1.segment_start

   **type:** Vector<Real>, mandatory

   Spatial location of the start of the segment.

.. input_param:: ChannelBuilder.s1.segment_end

   **type:** Vector<Real>, mandatory

   Spatial location of the end of the segment. The velocity profile
   (specified through other input arguments) is oriented along the
   line from the segment start to the segment end. The cross-sectional
   geometry (i.e., trapezoidal or elliptical) is perpendicular to this line.

.. input_param:: ChannelBuilder.s1.type

   **type:** String, optional, default = "ellipse"

   The type of each channel segment. Currently implemented types include
   "trapezoid" and "ellipse". The geometry parameters of a given segment
   depend on the type of that segment.

.. input_param:: ChannelBuilder.s1.horizontal_axis_start

   **type:** Real, optional

   For an elliptical segment, the length of the horizontal axis of the
   ellipse at the start of the segment.

.. input_param:: ChannelBuilder.s1.horizontal_axis_end

   **type:** Real, optional

   For an elliptical segment, the length of the horizontal axis of the
   ellipse at the end of the segment.

.. input_param:: ChannelBuilder.s1.vertical_axis_start

   **type:** Real, optional

   For an elliptical segment, the length of the vertical axis of the
   ellipse at the start of the segment.

.. input_param:: ChannelBuilder.s1.vertical_axis_end

   **type:** Real, optional

   For an elliptical segment, the length of the vertical axis of the
   ellipse at the end of the segment.

.. input_param:: ChannelBuilder.s1.horizontal_axis

   **type:** Real, optional

   For an elliptical segment, if the cross-section does not change
   along the segment, this parameter can be used to specify the length of
   the horizontal axis of the ellipse along the entire segment instead of
   specifying both the start and end values.

.. input_param:: ChannelBuilder.s1.vertical_axis

   **type:** Real, optional

   For an elliptical segment, if the cross-section does not change
   along the segment, this parameter can be used to specify the length of
   the vertical axis of the ellipse along the entire segment instead of
   specifying both the start and end values.

.. input_param:: ChannelBuilder.s1.diameter_start

   **type:** Real, optional

   For an elliptical segment, if the vertical and horizontal axes are the same,
   this parameter can be used to specify the diameter of the ellipse at the segment start.

.. input_param:: ChannelBuilder.s1.diameter_end

   **type:** Real, optional

   For an elliptical segment, if the vertical and horizontal axes are the same,
   this parameter can be used to specify the diameter of the ellipse at the segment end.

.. input_param:: ChannelBuilder.s1.diameter

   **type:** Real, optional

   For an elliptical segment, if the vertical and horizontal axes are the same
   and they do not change along the length of the segment,
   this parameter can be used to specify the diameter of the ellipse. Although
   all of these elliptical parameters are individually optional, at least one
   set of parameters must be provided for an elliptical segment or the code will fail.

.. input_param:: ChannelBuilder.s1.top_width_start

   **type:** Real, optional

   For a trapezoidal segment, the top width of the trapezoid
   at the start of the segment.

.. input_param:: ChannelBuilder.s1.top_width_end

   **type:** Real, optional

   For a trapezoidal segment, the top width of the trapezoid
   at the end of the segment.

.. input_param:: ChannelBuilder.s1.bottom_width_start

   **type:** Real, optional

   For a trapezoidal segment, the bottom width of the trapezoid
   at the start of the segment.

.. input_param:: ChannelBuilder.s1.bottom_width_end

   **type:** Real, optional

   For a trapezoidal segment, the bottom width of the trapezoid
   at the end of the segment.

.. input_param:: ChannelBuilder.s1.height_start

   **type:** Real, optional

   For a trapezoidal segment, the height of the trapezoid
   at the start of the segment.

.. input_param:: ChannelBuilder.s1.height_end

   **type:** Real, optional

   For a trapezoidal segment, the height of the trapezoid
   at the end of the segment.

.. input_param:: ChannelBuilder.s1.top_width

   **type:** Real, optional

   For a trapezoidal segment, if the cross-section does not change
   along the segment, this parameter can be used to specify the top width
   of the trapezoid along the entire segment instead of specifying both
   the start and end values.

.. input_param:: ChannelBuilder.s1.bottom_width

   **type:** Real, optional

   For a trapezoidal segment, if the cross-section does not change
   along the segment, this parameter can be used to specify the bottom width
   of the trapezoid along the entire segment.

.. input_param:: ChannelBuilder.s1.height

   **type:** Real, optional

   For a trapezoidal segment, if the cross-section does not change
   along the segment, this parameter can be used to specify the height
   of the trapezoid along the entire segment.

.. input_param:: ChannelBuilder.velocity_profile

   **type:** String, optional, default = "uniform"

   The initial velocity profile along the channel segment. Currently implemented
   profiles include "uniform", "linear", and "parabolic". The velocity
   profile is oriented along the line from the segment start to the segment end.
   For elliptical segments, the velocity profile is based on the distance (in any direction) from the
   center of the ellipse, but for trapezoidal segments, only the vertical distance is
   used to determine the velocity profile. For multiphase cases, the velocity profile
   is only applied to the liquid phase, setting the velocity above the
   :input_param:`ChannelBuilder.water_level` to zero.

.. input_param:: ChannelBuilder.s1.flow_speed

   **type:** Real, mandatory

   The flow speed associated with the velocity profile for the segment.
   This is the maximum flow speed of the velocity profile regardless of the type.

.. input_param:: ChannelBuilder.water_level

   **type:** Real, optional, default = 0.

   The water level determines the location of the free surface for
   multiphase channel flows. This setup parameter applies to the entire domain.
   This parameter is not used in single-phase cases.

.. input_param:: ChannelBuilder.land_level

   **type:** Real, mandatory when :input_param:`ChannelBuilder.water_level` is specified

   The land level determines where the channel opens up for free surface flows.
   Above this vertical location, no terrain blank is applied, creating an open channel.
   This parameter is not used in single-phase cases.