# Mathematical validation

The implementation was audited against the appendix of Schmeling, Kruse &
Richard (2012), with the PDF equations treated as the mathematical authority
and the supplied MATLAB program used as a second transcription reference.

## Equation checks

- Equations A1-A7 were mapped to the three elastic fixed-point solves and the
  final relaxed bulk-modulus correction.
- Equations A8-A15 were checked term by term for films, tubes, and spheroids.
- The viscous correspondence was checked by setting pore bulk and shear
  viscosity to zero and the matrix bulk viscosity to infinity.
- Equation 5 was checked against its spherical low-porosity limit,
  `eta_b / eta_s0 -> 4/(3phi)`.
- Density and P- and S-wave speed expressions were checked against the original
  GUI calculation.

## Corrections relative to the MATLAB program

Two MATLAB expressions do not agree with the published equations. The Julia
implementation follows the paper:

1. The auxiliary isolated-pore fraction is
   `phi_mi / (1 - phi_mi)`, as stated immediately before equation A6. The
   MATLAB loop divides by `1 - phi_mc`.
2. The exact `nu = 1/2` limit of equation A9 is
   `2(1-nu) * ratio / mu`. The MATLAB special branch has twice this value and
   also contains dimensionally inconsistent parentheses. The Julia branch is
   the continuous limit of the published equation.

The direct spheroid auxiliary expressions lose precision as the aspect ratio
approaches one. `ModVisc.jl` evaluates sixth-order series in
`x = 1 - alpha^2` near the spherical limit. This preserves Float32 and Float64
accuracy and joins continuously to the analytical spherical formulas.

## Numerical checks

- Scalar Float32 and Float64 return types are inferred and allocate zero bytes
  after compilation for concrete models.
- The spherical low-porosity bulk-viscosity asymptote is recovered.
- The Figure 2 disaggregation locations are recovered near melt fractions 0.20
  for tapered tubes and 0.50 for spheres.
- Spheroid numerical curves track the published equation 5 fits and reproduce
  Figure 3.
- Purely isolated elastic pores give identical relaxed and unrelaxed results.
- Fixed-point residuals, circular-tube limits, incompressible limits, and
  near-spherical continuity are covered by tests.
- Default stopping tolerances reproduce the original values: `1e-10` for the
  viscous solver and `1e-11` for the elastic solver.

## Remaining external validation

The supplied folder does not contain the numerical tables from Schmeling
(1985), nor the experimental and external-model data plotted in Figure 2.
Those data therefore cannot be used as independent numerical fixtures here.
The elastic solver has been checked against the printed equations and the
supplied MATLAB algorithm, but a published tabulated benchmark would provide a
stronger independent validation.
