"""Supertype for all supported pore inclusion geometries."""
abstract type AbstractGeometry end

"""Film-like oblate inclusions with aspect ratio `aspect_ratio ≪ 1`."""
struct Film{T<:AbstractFloat} <: AbstractGeometry
    aspect_ratio::T
    function Film{T}(aspect_ratio::T) where {T<:AbstractFloat}
        isfinite(aspect_ratio) && zero(T) < aspect_ratio <= one(T) ||
            throw(ArgumentError("film aspect ratio must lie in (0, 1]"))
        new{T}(aspect_ratio)
    end
end

Film(aspect_ratio::T) where {T<:AbstractFloat} = Film{T}(aspect_ratio)
Film(aspect_ratio::Real) = Film(float(aspect_ratio))

"""
Tubular inclusions with Mavko's cross-section parameter `shape`.

`shape = 0` describes tapered tubes; increasing values approach a circular
cross-section.
"""
struct Tube{T<:AbstractFloat} <: AbstractGeometry
    shape::T
    function Tube{T}(shape::T) where {T<:AbstractFloat}
        !isnan(shape) && shape >= zero(T) ||
            throw(ArgumentError("tube shape parameter must be nonnegative"))
        new{T}(shape)
    end
end

Tube(shape::T) where {T<:AbstractFloat} = Tube{T}(shape)
Tube(shape::Real) = Tube(float(shape))

"""Oblate spheroidal inclusions with aspect ratio in `(0, 1]`."""
struct Spheroid{T<:AbstractFloat} <: AbstractGeometry
    aspect_ratio::T
    function Spheroid{T}(aspect_ratio::T) where {T<:AbstractFloat}
        isfinite(aspect_ratio) && zero(T) < aspect_ratio <= one(T) ||
            throw(ArgumentError("spheroid aspect ratio must lie in (0, 1]"))
        new{T}(aspect_ratio)
    end
end

Spheroid(aspect_ratio::T) where {T<:AbstractFloat} = Spheroid{T}(aspect_ratio)
Spheroid(aspect_ratio::Real) = Spheroid(float(aspect_ratio))

@inline _convert_geometry(geometry::Film, ::Type{T}) where {T} = Film(T(geometry.aspect_ratio))
@inline _convert_geometry(geometry::Tube, ::Type{T}) where {T} = Tube(T(geometry.shape))
@inline _convert_geometry(geometry::Spheroid, ::Type{T}) where {T} = Spheroid(T(geometry.aspect_ratio))

"""A geometry and its fraction of the total pore volume."""
struct MeltFraction{G<:AbstractGeometry,T<:AbstractFloat}
    geometry::G
    fraction::T
    function MeltFraction(geometry::G, fraction::T) where {G<:AbstractGeometry,T<:AbstractFloat}
        isfinite(fraction) && zero(T) <= fraction <= one(T) ||
            throw(ArgumentError("geometry fraction must lie in [0, 1]"))
        new{G,T}(geometry, fraction)
    end
end

MeltFraction(geometry::AbstractGeometry, fraction::Real) =
    MeltFraction(geometry, float(fraction))

"""A geometry with isolated and hydraulically connected pore fractions."""
struct PoreFraction{G<:AbstractGeometry,T<:AbstractFloat}
    geometry::G
    isolated::T
    connected::T
    function PoreFraction(geometry::G, isolated::T, connected::T) where
            {G<:AbstractGeometry,T<:AbstractFloat}
        isfinite(isolated) && isfinite(connected) &&
            isolated >= zero(T) && connected >= zero(T) &&
            isolated + connected <= one(T) ||
            throw(ArgumentError("isolated and connected fractions must be nonnegative and sum to at most one"))
        new{G,T}(geometry, isolated, connected)
    end
end

function PoreFraction(geometry::AbstractGeometry, isolated::Real, connected::Real)
    i, c = promote(float(isolated), float(connected))
    return PoreFraction(geometry, i, c)
end
