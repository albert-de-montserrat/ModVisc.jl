@inline function _elastic_sums(::Tuple{}, bulk, shear, phi, model, mode)
    return (zero(bulk), zero(shear))
end

@inline function _elastic_sums(components::Tuple, bulk, shear, phi, model, mode)
    component = first(components)
    tail_bulk, tail_shear = _elastic_sums(
        Base.tail(components), bulk, shear, phi, model, mode)
    nu = poisson_ratio(bulk, shear)
    R = resistance_ratio(bulk, shear)
    theta = _bulk_compliance(component.geometry, bulk, shear, nu, R)
    Afluid = _elastic_shear_compliance(component.geometry, bulk, shear, nu, R,
                                       model.fluid_bulk, model.matrix_bulk, theta)

    isolated = component.isolated * phi
    connected = component.connected * phi
    if mode === Val(:unrelaxed)
        delta = inv(model.fluid_bulk) - inv(model.matrix_bulk)
        bulk_term = delta * (isolated + connected) /
                    (one(bulk) + inv(theta) *
                     (inv(model.fluid_bulk) - inv(bulk)))
        shear_term = Afluid * (isolated + connected)
    elseif mode === Val(:bar)
        delta = inv(model.fluid_bulk) - inv(model.matrix_bulk)
        bulk_term = delta * isolated /
                    (one(bulk) + inv(theta) *
                     (inv(model.fluid_bulk) - inv(bulk))) + theta * connected
        Avoid = _elastic_shear_compliance(component.geometry, bulk, shear, nu, R,
                                           zero(model.fluid_bulk), model.matrix_bulk, theta)
        shear_term = Afluid * isolated + Avoid * connected
    else
        # Paper appendix: phi'_mi = phi_mi / (1 - phi_mi).
        modified = isolated / (one(phi) - isolated)
        delta = inv(model.fluid_bulk) - inv(model.matrix_bulk)
        bulk_term = delta * modified /
                    (one(bulk) + inv(theta) *
                     (inv(model.fluid_bulk) - inv(bulk)))
        shear_term = Afluid * modified
    end
    return bulk_term + tail_bulk, shear_term + tail_shear
end

function _elastic_fixed_point(model, phi, options, mode)
    bulk, shear = model.matrix_bulk, model.matrix_shear
    tolerance_bulk = muladd(options.relative_tolerance, model.matrix_bulk,
                            options.absolute_tolerance)
    tolerance_shear = muladd(options.relative_tolerance, model.matrix_shear,
                             options.absolute_tolerance)
    for iteration in 1:options.max_iterations
        bulk_sum, shear_sum = _elastic_sums(model.components, bulk, shear,
                                            phi, model, mode)
        next_bulk = inv(inv(model.matrix_bulk) + bulk_sum)
        next_shear = inv(inv(model.matrix_shear) + shear_sum)
        converged = abs(next_bulk - bulk) <= tolerance_bulk &&
                    abs(next_shear - shear) <= tolerance_shear
        bulk, shear = next_bulk, next_shear
        converged && return bulk, shear, iteration, true
    end
    return bulk, shear, options.max_iterations, false
end

"""Solve unrelaxed and relaxed effective elastic properties at one porosity."""
function elasticity(model::ElasticModel, porosity::Real,
                    options::SolverOptions=SolverOptions(relative_tolerance=1e-11))
    phi = typeof(model.matrix_bulk)(porosity)
    isfinite(phi) && zero(phi) < phi < one(phi) ||
        throw(ArgumentError("porosity must lie in (0, 1)"))
    model.fluid_bulk > zero(model.fluid_bulk) ||
        throw(ArgumentError("elastic calculations require a positive fluid bulk modulus"))

    Ku, muu, iu, cu = _elastic_fixed_point(model, phi, options, Val(:unrelaxed))
    Kbar, mur, ib, cb = _elastic_fixed_point(model, phi, options, Val(:bar))
    Ki, _, ii, ci = _elastic_fixed_point(model, phi, options, Val(:isolated))

    connected = phi * sum(c -> c.connected, model.components)
    if iszero(connected)
        Kr = Kbar
    else
        F = model.fluid_bulk * (Ki - Kbar) /
            (connected * (Ki - model.fluid_bulk))
        Kr = Ki * (Kbar + F) / (Ki + F)
    end

    nuu = poisson_ratio(Ku, muu)
    nur = poisson_ratio(Kr, mur)
    rho = model.matrix_density - model.density_contrast * phi
    rho > zero(rho) || throw(DomainError(rho, "mixture density must be positive"))
    vpu = sqrt((Ku + 4muu / 3) / rho)
    vsu = sqrt(muu / rho)
    vpr = sqrt((Kr + 4mur / 3) / rho)
    vsr = sqrt(mur / rho)
    return ElasticResult(Ku, muu, Kr, mur, nuu, nur, vpu, vsu, vpr, vsr,
                         (iu, ib, ii), cu && cb && ci)
end
