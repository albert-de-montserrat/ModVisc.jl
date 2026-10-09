"""
    ViscousModel(matrix_shear, components...)

Self-consistent viscous model. `matrix_shear` is the intrinsic matrix shear
viscosity and the `MeltFraction` component fractions must sum to one.
"""
struct ViscousModel{T<:AbstractFloat,C<:Tuple}
    matrix_shear::T
    components::C
end

function ViscousModel(matrix_shear::Real, components::MeltFraction...)
    η0 = float(matrix_shear)
    isfinite(η0) && η0 > zero(η0) ||
        throw(ArgumentError("matrix shear viscosity must be finite and positive"))
    isempty(components) && throw(ArgumentError("at least one geometry is required"))

    fractions = map(c -> c.fraction, components)
    total = sum(fractions)
    T = promote_type(typeof(η0), map(typeof, fractions)...)
    tolerance = 8 * eps(T) * max(one(T), T(length(components)))
    isapprox(total, one(total); rtol=tolerance, atol=tolerance) ||
        throw(ArgumentError("geometry fractions must sum to one (received $total)"))

    promoted = map(c -> MeltFraction(_convert_geometry(c.geometry, T), T(c.fraction)), components)
    return ViscousModel{T,typeof(promoted)}(T(η0), promoted)
end

"""Controls fixed-point convergence for the self-consistent solver."""
struct SolverOptions{T<:AbstractFloat}
    relative_tolerance::T
    absolute_tolerance::T
    max_iterations::Int
    function SolverOptions{T}(relative_tolerance::T, absolute_tolerance::T,
                              max_iterations::Int) where {T<:AbstractFloat}
        relative_tolerance >= zero(T) ||
            throw(ArgumentError("relative tolerance must be nonnegative"))
        absolute_tolerance >= zero(T) ||
            throw(ArgumentError("absolute tolerance must be nonnegative"))
        max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
        new{T}(relative_tolerance, absolute_tolerance, max_iterations)
    end
end

function SolverOptions(; relative_tolerance::Real=1e-10,
                         absolute_tolerance::Real=zero(relative_tolerance),
                         max_iterations::Integer=1000)
    rtol, atol = promote(float(relative_tolerance), float(absolute_tolerance))
    return SolverOptions{typeof(rtol)}(rtol, atol, Int(max_iterations))
end

"""Effective viscosities and convergence metadata for one porosity."""
struct ViscosityResult{T<:AbstractFloat}
    bulk::T
    shear::T
    poisson::T
    iterations::Int
    converged::Bool
end

"""
    ElasticModel(matrix_bulk, matrix_shear, fluid_bulk, components...;
                 matrix_density=1, density_contrast=0)

Material and pore geometry for the elastic model. Moduli must use consistent
units. Densities are used only for wave speeds.
"""
struct ElasticModel{T<:AbstractFloat,C<:Tuple}
    matrix_bulk::T
    matrix_shear::T
    fluid_bulk::T
    matrix_density::T
    density_contrast::T
    components::C
end


function ElasticModel(matrix_bulk::Real, matrix_shear::Real, fluid_bulk::Real,
                      components::PoreFraction...;
                      matrix_density::Real=one(float(matrix_bulk)),
                      density_contrast::Real=zero(float(matrix_bulk)))
    isempty(components) && throw(ArgumentError("at least one geometry is required"))
    values = promote(float(matrix_bulk), float(matrix_shear), float(fluid_bulk),
                     float(matrix_density), float(density_contrast))
    K0, mu0, Kf, rho0, delta_rho = values
    all(isfinite, values) || throw(ArgumentError("material values must be finite"))
    K0 > 0 && mu0 > 0 && Kf >= 0 || throw(ArgumentError("moduli must be positive, with fluid bulk modulus nonnegative"))
    rho0 > 0 || throw(ArgumentError("matrix density must be positive"))

    T = typeof(K0)
    total = sum(c -> c.isolated + c.connected, components)
    tolerance = 8eps(T) * max(one(T), T(length(components)))
    isapprox(total, one(total); rtol=tolerance, atol=tolerance) ||
        throw(ArgumentError("pore fractions must sum to one (received $total)"))
    promoted = map(c -> PoreFraction(_convert_geometry(c.geometry, T),
                                     T(c.isolated), T(c.connected)), components)
    return ElasticModel{T,typeof(promoted)}(K0, mu0, Kf, rho0, delta_rho, promoted)
end

"""Effective elastic properties at one porosity."""
struct ElasticResult{T<:AbstractFloat}
    unrelaxed_bulk::T
    unrelaxed_shear::T
    relaxed_bulk::T
    relaxed_shear::T
    unrelaxed_poisson::T
    relaxed_poisson::T
    unrelaxed_p_speed::T
    unrelaxed_s_speed::T
    relaxed_p_speed::T
    relaxed_s_speed::T
    iterations::NTuple{3,Int}
    converged::Bool
end
