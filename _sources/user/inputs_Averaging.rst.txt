.. _inputs_averaging:
   
Section: Averaging
~~~~~~~~~~~~~~~~~~

This section controls data-averaging actions supported within
AMR-wind. 
Averaging is included as one of the identifiers of 
post-processing by
``incflo.post_processing = averaging``.
Time-averaging can be done to compute the means of any variable and 
is denoted by ``ReAveraging``. The code will add the prefix ``_mean`` 
to any variable that is time-averaged.
Reynolds stresses can also be computed using the ``ReynoldsStress`` averaging 
type.
Note that to compute the Reynolds stresses, the mean of the velocity
field is required.

.. input_param:: averaging.averaging_window

   **type:** Real, optional, default = huge number (no window)
   
   Specify the averaging window over which the time-averaging is done.

.. input_param:: averaging.averaging_start_time

   **type:** Real, optional, default = 0

   Specify the time to start computing averages.

.. input_param:: averaging.averaging_stop_time

   **type:** Real, optional, default = simulation end time

   Specify the time to stop time-averaging.

.. input_param:: averaging.averaging_time_interval

   **type:** Real, optional, default = -1. (off)

   Specify an interval of simulation time at which to do averaging.
   This allows for the averaging to take place only at instances
   that correspond to the averaging interval, which can be useful
   for phase averaging in contexts like periodic domains, regular
   waves, or fixed rotor speed simulations. When used with an
   averaging start time, the interval starts relative to the
   averaging start time.

Example::

   incflo.post_processing = averaging

   averaging.type = TimeAveraging
   averaging.labels = means  stress

   averaging.averaging_window = 10.0
   averaging.averaging_start_time = 0.0

   averaging.means.fields = velocity
   averaging.means.averaging_type = ReAveraging

   averaging.stress.fields = velocity
   averaging.stress.averaging_type = ReynoldsStress

  

