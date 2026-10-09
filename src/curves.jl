"""
    viscosity_profile!(bulk, shear, poisson, model, porosities[, options])
    viscosity!(bulk, shear, poisson, model, porosities[, options])

Fill preallocated result arrays for the supplied porosity values.
"""
function viscosity_profile!(bulk, shear, poisson, model::ViscousModel,
                            porosities::AbstractArray,
                            options::SolverOptions=SolverOptions())
    axes(bulk) == axes(shear) == axes(poisson) == axes(porosities) ||
        throw(DimensionMismatch("output arrays and porosities must have matching axes"))
    for i in eachindex(bulk, shear, poisson, porosities)
        result = viscosity(model, porosities[i], options)
        bulk[i] = result.bulk
        shear[i] = result.shear
        poisson[i] = result.poisson
    end
    return bulk, shear, poisson
end

viscosity!(bulk, shear, poisson, model::ViscousModel, porosities::AbstractArray,
           options::SolverOptions=SolverOptions()) =
    viscosity_profile!(bulk, shear, poisson, model, porosities, options)

"""
    viscosity_profile(model, porosity_or_porosities[, options])

Calculate effective viscosities at one porosity or over an array of porosities.
Scalar input returns a `ViscosityResult`; array input returns arrays in a named
tuple.
"""
@inline viscosity_profile(model::ViscousModel, porosity::Real,
                          options::SolverOptions=SolverOptions()) =
    viscosity(model, porosity, options)

function viscosity_profile(model::ViscousModel, porosities::AbstractArray,
                           options::SolverOptions=SolverOptions())
    T = promote_type(typeof(model.matrix_shear), eltype(porosities))
    bulk = similar(porosities, T)
    shear = similar(porosities, T)
    poisson = similar(porosities, T)
    viscosity_profile!(bulk, shear, poisson, model, porosities, options)
    return (; porosity=copy(porosities), bulk, shear, poisson)
end

viscosity(model::ViscousModel, porosities::AbstractArray,
          options::SolverOptions=SolverOptions()) =
    viscosity_profile(model, porosities, options)

"""
    elasticity_profile!(unrelaxed_bulk, unrelaxed_shear, relaxed_bulk,
                        relaxed_shear, unrelaxed_p_speed, unrelaxed_s_speed,
                        relaxed_p_speed, relaxed_s_speed, converged,
                        model, porosities[, options])

Fill preallocated arrays with the elastic properties for each porosity.
"""
function elasticity_profile!(unrelaxed_bulk, unrelaxed_shear, relaxed_bulk,
                             relaxed_shear, unrelaxed_p_speed, unrelaxed_s_speed,
                             relaxed_p_speed, relaxed_s_speed, converged,
                             model::ElasticModel, porosities::AbstractArray,
                             options::SolverOptions=SolverOptions(relative_tolerance=1e-11))
    axes(unrelaxed_bulk) == axes(unrelaxed_shear) == axes(relaxed_bulk) ==
        axes(relaxed_shear) == axes(unrelaxed_p_speed) == axes(unrelaxed_s_speed) ==
        axes(relaxed_p_speed) == axes(relaxed_s_speed) == axes(converged) ==
        axes(porosities) ||
        throw(DimensionMismatch("output arrays and porosities must have matching axes"))
    for i in eachindex(porosities)
        result = elasticity(model, porosities[i], options)
        unrelaxed_bulk[i] = result.unrelaxed_bulk
        unrelaxed_shear[i] = result.unrelaxed_shear
        relaxed_bulk[i] = result.relaxed_bulk
        relaxed_shear[i] = result.relaxed_shear
        unrelaxed_p_speed[i] = result.unrelaxed_p_speed
        unrelaxed_s_speed[i] = result.unrelaxed_s_speed
        relaxed_p_speed[i] = result.relaxed_p_speed
        relaxed_s_speed[i] = result.relaxed_s_speed
        converged[i] = result.converged
    end
    return (unrelaxed_bulk, unrelaxed_shear, relaxed_bulk, relaxed_shear,
            unrelaxed_p_speed, unrelaxed_s_speed, relaxed_p_speed,
            relaxed_s_speed, converged)
end

"""
    elasticity_profile!(unrelaxed_bulk, unrelaxed_shear, relaxed_bulk,
                        relaxed_shear, unrelaxed_p_speed, unrelaxed_s_speed,
                        relaxed_p_speed, relaxed_s_speed, converged,
                        model, porosities[, options])
    elasticity!(...)

Fill preallocated arrays with the elastic properties for each porosity.
"""
elasticity!(unrelaxed_bulk, unrelaxed_shear, relaxed_bulk, relaxed_shear,
            unrelaxed_p_speed, unrelaxed_s_speed, relaxed_p_speed, relaxed_s_speed,
            converged, model::ElasticModel, porosities::AbstractArray,
            options::SolverOptions=SolverOptions(relative_tolerance=1e-11)) =
    elasticity_profile!(unrelaxed_bulk, unrelaxed_shear, relaxed_bulk, relaxed_shear,
                        unrelaxed_p_speed, unrelaxed_s_speed, relaxed_p_speed,
                        relaxed_s_speed, converged, model, porosities, options)

"""
    elasticity_profile(model, porosities[, options])

Calculate unrelaxed and relaxed elastic properties over porosity values.
"""
function elasticity_profile(model::ElasticModel, porosities::AbstractArray,
                            options::SolverOptions=SolverOptions(relative_tolerance=1e-11))
    T = promote_type(typeof(model.matrix_bulk), eltype(porosities))
    unrelaxed_bulk = similar(porosities, T)
    unrelaxed_shear = similar(porosities, T)
    relaxed_bulk = similar(porosities, T)
    relaxed_shear = similar(porosities, T)
    unrelaxed_p_speed = similar(porosities, T)
    unrelaxed_s_speed = similar(porosities, T)
    relaxed_p_speed = similar(porosities, T)
    relaxed_s_speed = similar(porosities, T)
    converged = similar(porosities, Bool)
    elasticity_profile!(unrelaxed_bulk, unrelaxed_shear, relaxed_bulk, relaxed_shear,
                        unrelaxed_p_speed, unrelaxed_s_speed, relaxed_p_speed,
                        relaxed_s_speed, converged, model, porosities, options)
    return (; porosity=copy(porosities), unrelaxed_bulk, unrelaxed_shear,
            relaxed_bulk, relaxed_shear, unrelaxed_p_speed, unrelaxed_s_speed,
            relaxed_p_speed, relaxed_s_speed, converged)
end

elasticity(model::ElasticModel, porosities::AbstractArray,
           options::SolverOptions=SolverOptions(relative_tolerance=1e-11)) =
    elasticity_profile(model, porosities, options)

@inline _curve_argument(value, i) = value isa AbstractArray ? value[i] : value

"""
    spheroid_fit!(bulk, shear, porosities, aspect_ratio, matrix_shear=1)

Fill preallocated bulk and shear viscosity arrays using `spheroid_fit`.
`aspect_ratio` and `matrix_shear` may be scalars or arrays with axes matching
`porosities`.
"""
function spheroid_fit!(bulk, shear, porosities::AbstractArray, aspect_ratio,
                       matrix_shear=1)
    axes(bulk) == axes(shear) == axes(porosities) ||
        throw(DimensionMismatch("output arrays and porosities must have matching axes"))
    aspect_ratio isa AbstractArray && axes(aspect_ratio) != axes(porosities) &&
        throw(DimensionMismatch("aspect ratio and porosities must have matching axes"))
    matrix_shear isa AbstractArray && axes(matrix_shear) != axes(porosities) &&
        throw(DimensionMismatch("matrix shear and porosities must have matching axes"))
    for i in eachindex(bulk, shear, porosities)
        result = spheroid_fit(porosities[i], _curve_argument(aspect_ratio, i),
                              _curve_argument(matrix_shear, i))
        bulk[i] = result.bulk
        shear[i] = result.shear
    end
    return bulk, shear
end

"""
    spheroid_fit(porosity, aspect_ratio, matrix_shear=1)

Equation 5 of Schmeling, Kruse & Richard (2012). Returns the fitted `(bulk,
shear)` viscosities for spheroidal inclusions. Both are zero at and above the
critical porosity.
"""
function spheroid_fit(porosity::Real, aspect_ratio::Real, matrix_shear::Real=1)
    phi, alpha, eta0 = promote(float(porosity), float(aspect_ratio), float(matrix_shear))
    zero(alpha) < alpha <= one(alpha) ||
        throw(ArgumentError("aspect ratio must lie in (0, 1]"))
    zero(phi) <= phi <= one(phi) || throw(ArgumentError("porosity must lie in [0, 1]"))
    eta0 > zero(eta0) || throw(ArgumentError("matrix shear viscosity must be positive"))

    T = typeof(alpha)
    a1, a2 = T(0.97), T(0.8)
    b1, b2 = T(2.2455), T(3.45)
    k2, k3, c3 = T(1.25), T(1.29), T(2.4)
    k1 = a1 * (a2 + alpha * (one(alpha) - a2))
    critical = b1 * alpha / (one(alpha) + b2 * alpha^k3)
    phi >= critical && return (bulk=zero(eta0), shear=zero(eta0))
    c2 = (T(4) / T(3)) * alpha * critical^(-k2) *
         (c3 * (one(alpha) - alpha) + alpha)
    shear = eta0 * (one(phi) - phi / critical)^k1
    bulk = eta0 * c2 * (critical - phi)^k2 / phi
    return (; bulk, shear)
end
