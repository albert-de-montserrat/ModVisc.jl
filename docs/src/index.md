# ModVisc.jl

`ModVisc.jl` calculates effective elastic moduli, seismic velocities, and
viscosities for partially molten materials using the self-consistent inclusion
model of Schmeling, Kruse & Richard (2012).

## Viscous example

```@example viscous
using ModVisc

model = ViscousModel(1.0, MeltFraction(Tube(0.0), 1.0))
viscosity(model, 0.05)
```

`viscosity_profile` also dispatches on a scalar porosity and returns a
`ViscosityResult`:

```@example viscous
result = viscosity_profile(model, 0.05)
result.bulk
```

For several porosities, use the explicit profile API or collection dispatch:

```@example viscous
porosities = 10.0 .^ range(-4, -1; length=20)
profile = viscosity_profile(model, porosities)
profile.bulk == viscosity(model, porosities).bulk
```

Use `viscosity_profile!` to write into preallocated arrays:

```@example viscous
bulk, shear, poisson = similar(porosities), similar(porosities), similar(porosities)
viscosity_profile!(bulk, shear, poisson, model, porosities)
bulk == profile.bulk
```

## Elastic example

```@example elastic
using ModVisc

model = ElasticModel(
    66e9,
    40e9,
    20e9,
    PoreFraction(Spheroid(0.1), 0.0, 1.0);
    matrix_density=3300.0,
    density_contrast=300.0,
)
elasticity(model, 0.01)
```

```@example elastic
profile = elasticity_profile(model, [0.001, 0.01, 0.1])
profile.relaxed_bulk == elasticity(model, [0.001, 0.01, 0.1]).relaxed_bulk
```

`elasticity_profile!` and `spheroid_fit!` likewise fill preallocated arrays;
see their docstrings for the output-array order.
