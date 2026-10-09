KLAxell precursor/inflow homogeneity check
========================================

Status
------

These are experiment inputs, not measured results. The checkout was based on
``main`` revision ``733cea204e07738fdbe37c43ef239fb248c4e126``. Configuration
of the existing solver build failed because the AMReX and AMReX-Hydro
submodule directories were empty; no simulations have been executed.
Profile agreement and the cause of the reported heterogeneity remain unverified.

Start with the high inversion case
----------------------------------

The inputs under ``docs/sphinx/walkthrough/rans_inflow`` use a
2048 x 128 x 1024 m domain, 64 x 4 x 64 cells, dx = dy = 32 m and dz = 16 m.
This is a quasi-2D calculation using the existing 3D solver, with periodic y,
no imposed perturbations, flat terrain and no turbines.

The high inversion starts at 768 m and ends at 832 m. The assumed jump is
8 K over 64 m; potential temperature is constant below and above it, so the
background lapse rate is zero. Medium and low inputs are prepared at
512--576 m and 256--320 m, respectively, but should not be run until the high
case is understood. All cases have zero surface heat flux and z0 = 0.1 m.
The seed profile has five columns, as required by the implementation:
height, u, v, an unused fourth value, and TKE. It is only an initial guess.

Both runs use ``BoussinesqBuoyancy GeostrophicForcing``; neither uses
``CoriolisForcing``, ``ABLForcing`` or ``BodyForce``. Latitude is 45 degrees N
and Earth's default rotation period, 86164.091 s, is left untouched.
The resulting f is approximately 1.0313e-4 per second.

**Coordinate convention:** the 10 m/s geostrophic vector is (0, -10, 0).
``GeostrophicForcing`` applies f*(-Vg, Ug, 0), so this choice drives positive x
with acceleration approximately 0.0010313 m/s squared. A vector (10, 0, 0)
would instead drive positive y; rotate the domain/inflow accordingly if that
is the production convention. The initial streamwise wind is 10 m/s, but,
without Coriolis, geostrophic forcing does not constrain the final wind speed
to 10 m/s. Its equilibrium is determined by the imposed pressure gradient
and stress balance.

Execution
---------

Use an existing Kynema-SGF executable. Native boundary output and ASCII
sampling avoid a NetCDF dependency. All paths below are absolute; ``FILE``
and auxiliary paths inside the decks are resolved from the run directory.
Work outside the repository so outputs are not committed.

First copy the shared inputs and the high-case deck into a clean spinup
directory::

   ROOT=/home/runner/work/kynema-sgf/kynema-sgf
   INPUT="$ROOT/docs/sphinx/walkthrough/rans_inflow"
   RUN=/tmp/kynema-rans-high
   EXE=/absolute/path/to/kynema_sgf
   mkdir -p "$RUN/spinup" "$RUN/record" "$RUN/inflow"
   cp "$INPUT/common.inp" "$INPUT/precursor.inp" "$INPUT/seed.info" "$INPUT/probes.txt" "$RUN/spinup/"
   cp "$INPUT/high/precursor.inp" "$RUN/spinup/case.inp"
   cd "$RUN/spinup"
   "$EXE" "$RUN/spinup/case.inp" > "$RUN/spinup/run.log" 2>&1

Require the convergence monitor's complete window and hold to pass, and
check the full temperature and velocity profiles for continuing drift.
Reaching the 60000 s backstop alone is not convergence. Extend the run if
needed. Check that spinup is horizontally uniform at every height.

Identify the final checkpoint and its physical time T from the log/checkpoint.
Set ``CHK`` to its absolute path and ``END`` to T + 4096 s, rounded consistently
with the 2 s timestep. This recording interval is twenty nominal domain
transits at 10 m/s; slower near-wall flow may need more time. Record a periodic
continuation from that checkpoint::

   cp "$RUN/spinup/"*.inp "$INPUT/seed.info" "$INPUT/probes.txt" "$RUN/record/"
   cd "$RUN/record"
   "$EXE" "$RUN/record/case.inp" io.restart_file="$CHK" time.stop_time="$END" \
       ABL.bndry_io_mode=0 convergence.stop_on_convergence=false \
       > "$RUN/record/run.log" 2>&1

Then restart the inflow calculation from the **same spinup checkpoint**, not
the end of the recording run::

   cp "$INPUT/common.inp" "$INPUT/inflow.inp" "$INPUT/seed.info" "$INPUT/probes.txt" "$RUN/inflow/"
   cp "$INPUT/high/inflow.inp" "$RUN/inflow/case.inp"
   cd "$RUN/inflow"
   "$EXE" "$RUN/inflow/case.inp" io.restart_file="$CHK" time.stop_time="$END" \
       > "$RUN/inflow/run.log" 2>&1

Confirm boundary times span the entire restart interval before interpreting
results. Velocity, temperature **and TKE** are recorded/replayed. The
turbulent viscosity is computed by KLAxell rather than imposed as inflow data.
The grid, wall treatment, forcing, sponge setting and initial state match.
``ABL.inflow_outflow_mode=false`` intentionally retains the same wall-statistics
path in both runs; changing it to true requires consistent precursor-derived
``wf_velocity``, ``wf_vmag`` and ``wf_theta``, not their defaults.

Measure spatial heterogeneity, not just a domain mean
----------------------------------------------------

The ASCII ``profiles`` sampler covers every x-z cell centre at y = 16 m,
every 60 s and at the final step. Pair files using their ``*_info.txt`` physical
times, and match points by coordinates, not row order. Compare all 64 vertical
levels at x = 16, 272, 528, 1008, 1520 and 2032 m with the simultaneous
periodic continuation. Also compare later averaging windows after flushing
the slowest relevant near-wall transit.

For each sampled quantity q, report:

* profile mismatch: q_inflow(x,z,t) - q_periodic(z,t);
* spatial range at each height: max_x(q) - min_x(q);
* maximum absolute mismatch, RMS mismatch and downstream evolution;
* time drift separately from spatial drift.

Check u, v, w, temperature, TKE, turbulent viscosity and length scale,
not just horizontal speed. Suggested diagnostic tolerances are 0.01 m/s for
velocity components, 0.01 K for temperature, and max(0.001 m squared/s squared,
1 percent of precursor TKE) for TKE. These are proposed acceptance criteria,
not validated accuracy bounds. Report results above/below/inside the inversion
separately and retain full x-z fields to detect an outlet-localized disturbance.

The coupled sponge and length-scale switch
-----------------------------------------

``KLAxell.cpp`` assigns the local ``lengthscale_switch`` from
``ABL.meso_sponge_start``. It is **not a separate input parameter**. For zero
surface heat flux, below that height it overrides the stability-dependent
length scale with the neutral shear length and sets Rt to zero.
``KransAxell.cpp`` uses the same input to relax TKE above that height toward
the five-column reference profile, with a coefficient of 1/dt.

Consequently, moving ``meso_sponge_start`` to 1024 m disables the vertical TKE
sponge, but extends the neutral length-scale override throughout the domain,
including the inversion. Omitting the parameter has the same qualitative
effect here because its default is 2000 m. Setting it to zero instead enables
the TKE sponge almost everywhere; that is not a sponge-off test.
The common input explicitly uses 1024 m for a matched, sponge-free baseline;
it does **not** claim to disable the override. No temperature sponge source
is selected in these inputs.

For the high case, after the baseline, compare a separate matched pair using
``ABL.meso_sponge_start=752``. This moves the switch below the inversion but
also activates the TKE sponge, so use a converged reference TKE profile rather
than the seed. This is a coupled sensitivity experiment, not an isolated
length-scale test. If it changes the result, independently configurable
sponge enablement and length-scale-switch height would be needed to isolate
the two effects. Do not change only the inflow run.

Source-based hypotheses, not established causes
-----------------------------------------------

* Mismatched momentum forcing or missing TKE inflow changes the downstream
  stress balance even when inlet velocity looks correct.
* Changing wall statistics or frozen wall reference values can change shear.
* Changing the coupled sponge/switch or its reference profile changes TKE,
  length scale and viscosity. An inversion shifted by mixing can amplify
  downstream differences.
* A restart pressure adjustment or pressure-outflow boundary can create a
  localized disturbance; locate it in the x-z fields before blaming turbulence.
* An unconverged precursor or inadequate flush time can resemble persistent
  heterogeneity. The monitor tests speed/TKE, not temperature or direction.

At exactly zero stratification the unstable/neutral branch gives the neutral
shear length (Rt = 0); the epsilon floor in the unused buoyancy length is not
by itself evidence of neutral length-scale suppression. Inversion **height**
and inversion **strength** are distinct; these cases vary height only.
Do not attribute the colleague's discrepancy to a model defect without
measured matched-run evidence.
