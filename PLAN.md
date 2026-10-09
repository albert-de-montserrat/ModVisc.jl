# Plan: rewrite as `ModVisc.jl`

## Progress

- [x] Audit the 2012 paper appendix and original viscous iteration.
- [x] Establish concrete geometry, mixture, solver-option, and result types.
- [x] Implement allocation-free scalar viscous kernels and fixed-point solve.
- [x] Add the equation 5 spheroid fit and profile API with collection dispatch.
- [x] Add initial validation, numerical, convergence, inference, and allocation tests.
- [x] Add mathematical notes and Figure 2/3 reproduction scripts.
- [x] Validate and refine the generated figures against the published plots.
- [x] Implement the elastic unrelaxed/relaxed model and seismic velocities.
- [x] Add a Documenter-based manual and API reference.
- [x] Add continuous integration for tests and documentation builds.
- [x] Audit equations A1-A15 against the paper and MATLAB reference.
- [x] Add limiting-case and near-spherical numerical validation.
- [ ] Add documentation deployment once the repository URL and hosting branch are known.
- [ ] Assess GPU execution with an actual backend after the CPU reference is fixed.

### Decisions recorded during implementation

- Fixed tuples, rather than `StaticArrays.jl`, represent geometry mixtures. They
  preserve compile-time size and concrete element types without a core dependency.
- The scalar solver is separate from curve allocation and plotting.
- Both scalar solvers preserve `Float32` and `Float64`, are inference stable,
  and allocate zero bytes after compilation for concrete model types.
- `Tube(Inf)` is supported through its analytical circular limit.
- Solver results carry iteration count and convergence state. Near and above
  disaggregation the fixed-point iteration approaches zero asymptotically, so
  plotting uses the original code's scale-aware `1e-3` cutoff.
- Figure 2's external experimental and grain-contiguity overlays are deferred
  because their numerical data are not present in the supplied folder.
- The PDF equations are authoritative where the MATLAB implementation differs.
  The audit corrected the isolated-pore auxiliary fraction and the exact
  incompressible tube limit; details are in `docs/src/validation.md`.

The paper’s core is a self-consistent poroelastic model for melt inclusions. Its geometry terms cover films, tubes, and oblate spheroids; the elastic formulation separates isolated and connected inclusions and calculates unrelaxed and relaxed moduli. The viscous formulation applies the correspondence principle, assumes an incompressible matrix and zero pore-fluid viscosity, and calculates relaxed bulk and shear viscosity. The MATLAB code also calculates Poisson ratios and seismic velocities.

The rewrite should follow the published equations and physical concepts, using the MATLAB implementation as a cross-check rather than a design template. The first reproduction targets are the paper’s viscosity figures: **Figure 2** compares geometries and mixtures, and **Figure 3** compares spheroid model results with the published fitting equations. Figure 1 is a MATLAB GUI example, so it is not a numerical reproduction target.

## Package structure

- `src/ModVisc.jl`: public API and exports.
- `src/geometries.jl`: concrete geometry types for films, tubes, and spheroids.
- `src/parameters.jl`: validated material properties, geometry fractions, and solver settings.
- `src/kernels/`: scalar formulas for the paper’s geometry terms (`θ`, `A`, and auxiliary quantities).
- `src/solver.jl`: self-consistent iteration for one melt fraction, with convergence status and iteration count.
- `src/elastic.jl`: unrelaxed and relaxed moduli, Poisson ratios, density, and seismic velocities.
- `src/viscous.jl`: relaxed effective viscosities and viscous Poisson ratio.
- `src/curves.jl`: profile APIs to calculate a model over a melt-fraction grid.
- `docs/`: model assumptions, equations, API examples, and figure reproduction instructions.
- `test/`: equation, solver, and published-result checks.
- `examples/`: scripts that generate the paper figures and export their data.

## Design for type stability and GPU readiness

- Represent geometry choices with concrete types and multiple dispatch, rather than integer flags or dynamically typed collections.
- Put each model’s fixed parameters in immutable structs. Use tuples or statically sized containers for the three geometry contributions, whose size is known at compile time.
- Keep the scalar solver kernels separate from grid generation, plotting, I/O, and result storage. The scalar path should be usable independently inside CPU or GPU kernels.
- Avoid heap allocation in the per-porosity iteration. Preallocate result arrays at the outer level; use a compact immutable result for each scalar solve.
- Keep backend-specific concerns out of the mathematical core. Start with CPU-compatible generic scalar code, then assess GPU execution against the chosen solver control flow and transcendental functions.
- Consider `StaticArrays.jl` for genuinely small fixed-size values, while keeping large or variable-length curve data in ordinary arrays.

## Implementation stages

### 1. Equation and reference audit

- Map the paper’s equations and appendix terms to their model roles and source locations.
- Resolve notation and special cases, especially isolated versus connected fractions, spheroid limits, zero melt-fluid modulus, and the incompressible matrix limit.
- Record the paper’s defaults and figure-specific parameter sets, including the Figure 2 geometry mixtures and Figure 3 spheroid aspect ratios.
- Establish which outputs and conventions the first release supports.

### 2. Julia package foundation

- Replace the scaffold with the package API, model types, input validation, and documentation setup.
- Define explicit result types for elastic and viscous calculations, including solver convergence information.
- Decide a consistent policy for non-convergence and singular or invalid parameter combinations.

### 3. Scalar mathematical kernels

- Implement geometry functions and auxiliary quantities from the appendix.
- Add focused tests against hand-calculated cases, limiting cases, and independent evaluations of the published equations.
- Pay particular attention to spheroid aspect ratio 1, thin films, large tube shape parameter, and divisions that become ill-conditioned near physical limits.

### 4. Viscous solver

- Implement the relaxed bulk and shear iteration for a single melt fraction.
- Test pure-geometry cases and mixtures, convergence behavior, and the expected drop to zero at finite melt fraction for geometries that disaggregate.
- Compare against the MATLAB implementation for selected parameter sets while treating the paper as the authoritative specification.

### 5. Elastic solver

- Add unrelaxed and relaxed bulk/shear moduli, including the isolated-pore auxiliary solve.
- Derive Poisson ratios and P- and S-wave velocities from the resulting moduli and density.
- Test against the paper’s cited agreement with Schmeling (1985), as well as elastic limiting cases.

### 6. Curve API and first figure reproductions

- Add log-spaced melt-fraction grids and batch curve calculation.
- Generate Figure 2 viscosity curves and Figure 3 spheroid curves plus the published fit (equation 5).
- Save plotted data and figure-generation parameters so discrepancies can be traced to inputs or numerical behavior.

### 7. Stability, performance, and documentation

- Check inference and allocations for the scalar solver and curve path; keep any plotting or export allocations outside the numerical kernel.
- Add API documentation, assumptions and limitations, equation references, and a guide for regenerating Figures 2 and 3.
- Once CPU results and numerical behavior are established, assess GPU compatibility with representative kernels and supported precisions.

## Initial acceptance criteria

- The API expresses geometry, material parameters, and isolated/connected fractions without runtime type instability.
- Scalar solver results report whether they converged and do not silently hide invalid values.
- Figure 2 and Figure 3 can be regenerated from documented Julia scripts, with figure data and input parameters available for review.
- Tests cover equation helpers, solver behavior, elastic/viscous outputs, and the published fit.

## Initial technical risk

Some appendix formulas are dense and special-case-heavy. Their transcription and limiting behavior should be resolved before optimizing or promising allocation-free GPU execution.
