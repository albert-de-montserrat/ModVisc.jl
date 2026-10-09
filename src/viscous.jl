"""
    viscosity(model, porosity[, options]) -> ViscosityResult

Solve the relaxed self-consistent equations for effective bulk and shear
viscosity at one melt fraction. This is the viscous correspondence of equations
A4 and A6 in Schmeling, Kruse & Richard (2012), with an incompressible matrix
and zero pore-fluid viscosities.
"""
function viscosity(model::ViscousModel, porosity::Real,
                   options::SolverOptions=SolverOptions())
    eta0, phi = promote(model.matrix_shear, float(porosity))
    rtol = typeof(eta0)(options.relative_tolerance)
    atol = typeof(eta0)(options.absolute_tolerance)
    isfinite(phi) && zero(phi) < phi <= one(phi) ||
        throw(ArgumentError("porosity must lie in (0, 1]"))

    bulk = eta0
    shear = eta0
    tolerance = muladd(rtol, eta0, atol)

    for iteration in 1:options.max_iterations
        nu = poisson_ratio(bulk, shear)
        R = resistance_ratio(bulk, shear)
        bulk_compliance, shear_compliance = _component_compliances(
            model.components, bulk, shear, phi, nu, R)
        next_bulk = inv(bulk_compliance)
        next_shear = inv(inv(eta0) + shear_compliance)

        if !(isfinite(next_bulk) && isfinite(next_shear)) ||
           next_bulk <= zero(next_bulk) || next_shear <= zero(next_shear)
            z = zero(next_bulk)
            return ViscosityResult(z, z, typeof(z)(0.5), iteration, false)
        end

        error = max(abs(next_bulk - bulk), abs(next_shear - shear))
        bulk, shear = next_bulk, next_shear
        if error <= tolerance
            return ViscosityResult(bulk, shear, poisson_ratio(bulk, shear),
                                   iteration, true)
        end
    end

    return ViscosityResult(bulk, shear, poisson_ratio(bulk, shear),
                           options.max_iterations, false)
end
