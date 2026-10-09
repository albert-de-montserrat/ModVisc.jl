# Mathematical model

## Scope

The implementation covers the relaxed viscous and full elastic formulations in
the appendix of Schmeling, Kruse & Richard (2012). In the viscous case the matrix has intrinsic shear
viscosity `eta_s0`, infinite intrinsic bulk viscosity, and the pores have zero
bulk and shear viscosities. Effective values are normalized by `eta_s0` in the
paper figures.

The elastic model distinguishes isolated and hydraulically connected pores and
returns unrelaxed and relaxed moduli, Poisson ratios, and P- and S-wave speeds.

## Equation map

| Paper expression | Package implementation | Role |
|---|---|---|
| A7 | `poisson_ratio` | Viscous Poisson ratio |
| A8-A10 | `_bulk_compliance` | Geometry-dependent bulk compliance |
| A11 | `_spheroid_terms`, `resistance_ratio` | Spheroid auxiliaries |
| A12-A15 | `_shear_compliance` | Geometry-dependent shear compliance |
| A4 and A6 correspondence | `viscosity` | Coupled self-consistent fixed point |
| A1-A3 | `elasticity` | Unrelaxed bulk/shear moduli and Poisson ratio |
| A4-A7 | `elasticity` | Auxiliary and relaxed elastic properties |
| Equation 5 | `spheroid_fit` | Published spheroid fitting law |

The source paper uses `phi` both for total melt fraction and for an auxiliary
spheroid function in its appendix. The implementation uses `porosity` for the
total melt fraction and keeps the auxiliary internal, avoiding that ambiguity.

## Geometry conventions

- `Film(alpha)` uses the thin oblate-film expressions. The paper recommends
  `alpha << 1`.
- `Tube(kappa)` uses the Mavko cross-section parameter. `kappa = 0` is tapered;
  `Tube(Inf)` evaluates the circular limit without an `Inf/Inf` operation.
- `Spheroid(alpha)` accepts `0 < alpha <= 1`; `alpha = 1` uses the analytical
  spherical special case.
- `MeltFraction(geometry, fraction)` is a fraction of total pore volume. All
  fractions in a `ViscousModel` must sum to one.

## Solver behavior

`viscosity` starts both effective viscosities at the intrinsic matrix shear
viscosity and performs a fixed-point iteration for the coupled effective
viscosities. Convergence is tested against

```text
absolute_tolerance + relative_tolerance * matrix_shear
```

for both effective viscosities. The default relative tolerance is `1e-10`.

The elastic API defaults to a relative tolerance of `1e-11`.

Above the disaggregation porosity, the physical fixed point is zero. The
iteration generally approaches a small positive number, so callers should use
the `converged` field and a scale-aware plotting cutoff. The figure
reproductions use `1e-3 * eta_s0`.

The elastic path performs three coupled solves.
Moduli must use a consistent unit system. Wave speeds have the corresponding
derived units; using pascals and kilograms per cubic metre produces metres per
second.

## Figure reproduction limits

The included Julia example reproduces the self-consistent model curves in
Figures 2 and 3 and the equation 5 fits. Figure 2 also contains experimental
and external-model curves from Kohlstedt et al. (2000), Hirth & Kohlstedt
(1995), and Takei & Holtzman (2009). Their underlying numerical data are not
distributed with this package, so those overlays are absent.
