"""Return Poisson's ratio from isotropic bulk and shear moduli or viscosities."""
@inline function poisson_ratio(bulk::T, shear::T) where {T<:AbstractFloat}
    isinf(bulk) && return T(0.5)
    return (3bulk - 2shear) / (6bulk + 2shear)
end

@inline resistance_ratio(bulk, shear) = 3shear / (3bulk + 4shear)

@inline function _spheroid_terms(alpha::T) where {T}
    alpha == one(T) && return (zero(T), zero(T))
    alpha2 = alpha * alpha
    x = one(T) - alpha2
    if x < T(0.1)
        q = T(2) / T(3) + x * (-T(2) / T(15) + x * (-T(8) / T(105) +
            x * (-T(16) / T(315) + x * (-T(128) / T(3465) +
            x * (-T(256) / T(9009) - x * T(1024) / T(45045))))))
        g = -T(2) / T(5) + x * (T(6) / T(35) + x * (T(8) / T(105) +
            x * (T(16) / T(385) + x * (T(128) / T(5005) +
            x * (T(256) / T(15015) + x * T(1024) / T(85085))))))
        return q, g
    end
    q = alpha * x^(-T(1.5)) * (acos(alpha) - alpha * sqrt(x))
    g = alpha2 / x * (3q - 2)
    return q, g
end

@inline function _tube_ratios(kappa)
    isinf(kappa) && return (one(kappa), one(kappa))
    h = (2 + kappa)^2
    return ((h + 2) / (h - 2), h / (h - 2))
end

@inline function _bulk_compliance(geometry::Film, bulk, shear, nu, R)
    alpha = geometry.aspect_ratio
    if nu == typeof(nu)(0.5)
        return 2(one(nu) - nu) / (typeof(nu)(pi) * alpha * shear)
    end
    return 4 / (3 * typeof(nu)(pi) * bulk) *
           (one(nu) - nu * nu) / (one(nu) - 2nu) / alpha
end

@inline function _bulk_compliance(geometry::Tube, bulk, shear, nu, R)
    kappa = geometry.shape
    cross_section, _ = _tube_ratios(kappa)
    if nu == typeof(nu)(0.5)
        return 2(one(nu) - nu) * cross_section / shear +
               (one(nu) - 2nu) / (3bulk)
    end
    return 2 / (3bulk) *
           (2(one(nu) - nu * nu) / (one(nu) - 2nu) * cross_section +
            (one(nu) - 2nu) / 2)
end

@inline function _bulk_compliance(geometry::Spheroid, bulk, shear, nu, R)
    q, g = _spheroid_terms(geometry.aspect_ratio)
    if iszero(q) && iszero(g)
        nu == typeof(nu)(0.5) &&
            return 9(one(nu) - nu) / (4shear * (one(nu) + nu))
        return 3 / (2bulk) * (one(nu) - nu) / (one(nu) - 2nu)
    end
    four_thirds = oftype(bulk, 4) / oftype(bulk, 3)
    d1 = one(bulk) - 3(g + q) / 2 + R * (3g / 2 + 5q / 2 - four_thirds)
    d2 = one(bulk) - (one(bulk) + 3(g + q) / 2 - R * (3g / 2 + 5q / 2)) +
         (3 - 4R) / 2 * (g + q - R * (g - q + 2q * q))
    return d1 / (bulk * d2)
end

@inline function _shear_compliance(geometry::Film, bulk, shear, nu, R)
    alpha = geometry.aspect_ratio
    return 8 / (15 * typeof(nu)(pi) * shear) * (one(nu) - nu) / (2 - nu) *
           ((2 - nu) + 3) / alpha
end

@inline function _shear_compliance(geometry::Tube, bulk, shear, nu, R)
    kappa = geometry.shape
    cross_section, circular = _tube_ratios(kappa)
    two_fifteenths = oftype(nu, 2) / oftype(nu, 15)
    eight_fifths = oftype(nu, 8) / oftype(nu, 5)
    three_halves = oftype(nu, 3) / oftype(nu, 2)
    return (two_fifteenths * (one(nu) + nu) +
            two_fifteenths * (one(nu) - nu) * cross_section +
            eight_fifths * (three_halves - nu) * circular) / shear
end

@inline function _shear_compliance(geometry::Spheroid, bulk, shear, nu, R)
    return _spheroid_shear_compliance(geometry, bulk, shear, R, zero(bulk))
end

@inline function _spheroid_shear_compliance(geometry::Spheroid, bulk, shear, R, B)
    alpha = geometry.aspect_ratio
    alpha == one(alpha) && return 5 / shear * (4shear + 3bulk) / (8shear + 9bulk)

    q, g = _spheroid_terms(alpha)
    four_thirds = oftype(alpha, 4) / oftype(alpha, 3)
    common = one(alpha) - (3q + g - R * (g - q)) / 4
    t1 = 2 / (one(alpha) - (-(one(alpha) + alpha^2) * g / alpha^2 +
                         R * (2 - q + (one(alpha) + alpha^2) * g / alpha^2)) / 2)
    t2 = inv(common)
    t3 = common
    t4 = B * q * (3 - 4R) + g - R * (g + q - four_thirds)
    t5 = 2 * (R * (g + q) - g + B * (one(B) - q) * (3 - 4R))
    t6 = one(alpha) - (9q + 3g - R * (5q + 3g)) / 8 + B * q * (3 - 4R) / 2
    t7 = 2 * (B * (one(B) - q) * (3 - 4R) - one(alpha) + 3q / 2 +
              g / 2 - R * (5q + g - 4) / 2)
    t8 = (B * q * (3 - 4R) + g - R * (g - q)) / 2
    t10 = inv(R * (3g / 2 + 5q / 2) - 3(g + q) / 2 +
              B * (3 - 4R) - (3B - one(B)) * (3 - 4R) *
              (g + q - R * (g - q + 2q^2)) / 2)
    return (t1 + t2 + (t3 * t4 + t5 * t6 - t7 * t8) * t2 * t10) /
           (5shear)
end

@inline _fluid_B(::Film, fluid_bulk, bulk, nu, shear) = zero(bulk)
@inline _fluid_B(geometry::Tube, fluid_bulk, bulk, nu, shear) =
    _shear_compliance(geometry, bulk, shear, nu, resistance_ratio(bulk, shear))
@inline _fluid_B(::Spheroid, fluid_bulk, bulk, nu, shear) = fluid_bulk / (3bulk)

@inline function _elastic_shear_compliance(geometry::Film, bulk, shear, nu, R,
                                           fluid_bulk, matrix_bulk, theta)
    D = iszero(fluid_bulk) ? one(bulk) :
        (inv(fluid_bulk) - inv(matrix_bulk)) / (theta + inv(fluid_bulk))
    alpha = geometry.aspect_ratio
    return 8 / (15 * typeof(nu)(pi) * shear) * (one(nu) - nu) / (2 - nu) *
           ((2 - nu) * D + 3) / alpha
end

@inline function _elastic_shear_compliance(geometry::Tube, bulk, shear, nu, R,
                                           fluid_bulk, matrix_bulk, theta)
    B = _fluid_B(geometry, fluid_bulk, bulk, nu, shear)
    iszero(fluid_bulk) && return B
    kappa = geometry.shape
    h = (2 + kappa)^2
    invh = isinf(h) ? zero(h) : inv(h)
    plus = one(h) + 2invh
    minus = one(h) - 2invh
    s1 = 2 * (one(nu) - nu) * plus / minus - one(nu) + 2nu
    s2_over_h = 2 * (one(nu) - nu) * plus - (one(nu) - 2nu) * minus
    denominator_over_h = -2 * (one(nu) - nu) * plus + shear * minus *
        (inv(bulk) - inv(fluid_bulk) -
         (one(nu) - 2nu)^2 / (2shear * (one(nu) + nu)))
    return B + s1 * s2_over_h / denominator_over_h / (15shear)
end

@inline function _elastic_shear_compliance(geometry::Spheroid, bulk, shear, nu, R,
                                           fluid_bulk, matrix_bulk, theta)
    B = fluid_bulk / (3bulk)
    return _spheroid_shear_compliance(geometry, bulk, shear, R, B)
end

@inline _component_compliances(::Tuple{}, bulk, shear, porosity, nu, R) =
    (zero(bulk), zero(shear))

@inline function _component_compliances(components::Tuple, bulk, shear,
                                        porosity, nu, R)
    component = first(components)
    tail_bulk, tail_shear = _component_compliances(
        Base.tail(components), bulk, shear, porosity, nu, R)
    weight = component.fraction * porosity
    return (
        muladd(weight, _bulk_compliance(component.geometry, bulk, shear, nu, R), tail_bulk),
        muladd(weight, _shear_compliance(component.geometry, bulk, shear, nu, R), tail_shear),
    )
end
