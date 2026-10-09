using Test
using ModVisc

@testset "geometry and model validation" begin
    @test_throws ArgumentError Film(0.0)
    @test_throws ArgumentError Film(1.01)
    @test_throws ArgumentError Tube(-1.0)
    @test_throws ArgumentError Spheroid(1.01)
    @test_throws ArgumentError ViscousModel(
        1.0,
        MeltFraction(Film(0.01), 0.4),
        MeltFraction(Spheroid(1.0), 0.5),
    )

    model = ViscousModel(
        1.0,
        MeltFraction(Tube(0.0), 0.5),
        MeltFraction(Spheroid(1.0), 0.5),
    )
    @test isbitstype(typeof(model))
end

@testset "elastic self-consistent solver" begin
    connected = ElasticModel(
        0.66,
        0.4,
        0.2,
        PoreFraction(Spheroid(0.1), 0.0, 1.0);
        matrix_density=3300.0,
        density_contrast=300.0,
    )
    result = @inferred elasticity(connected, 0.01)
    @test result.converged
    @test result.unrelaxed_bulk ≈ 0.6482651401218136 rtol=1e-10
    @test result.unrelaxed_shear ≈ 0.38505945312590206 rtol=1e-10
    @test result.relaxed_bulk ≈ 0.6482199586600882 rtol=1e-10
    @test result.relaxed_shear ≈ 0.38139936055754264 rtol=1e-10
    @test result.unrelaxed_p_speed > result.unrelaxed_s_speed > 0
    @test result.relaxed_p_speed > result.relaxed_s_speed > 0

    Kbar, mubar, _, _ = ModVisc._elastic_fixed_point(
        connected, 0.01, SolverOptions(relative_tolerance=1e-11), Val(:bar)
    )
    for (mode, bulk, shear) in ((Val(:unrelaxed), result.unrelaxed_bulk,
                                 result.unrelaxed_shear),
                                (Val(:bar), Kbar, mubar))
        bulk_sum, shear_sum = ModVisc._elastic_sums(
            connected.components, bulk, shear, 0.01, connected, mode
        )
        @test bulk ≈ inv(inv(connected.matrix_bulk) + bulk_sum) rtol=2e-10
        @test shear ≈ inv(inv(connected.matrix_shear) + shear_sum) rtol=2e-10
    end

    isolated = ElasticModel(0.66, 0.4, 0.2,
                            PoreFraction(Spheroid(0.1), 1.0, 0.0))
    isolated_result = elasticity(isolated, 0.01)
    @test isolated_result.relaxed_bulk ≈ isolated_result.unrelaxed_bulk
    @test isolated_result.relaxed_shear ≈ isolated_result.unrelaxed_shear

    circular = ElasticModel(0.66, 0.4, 0.2,
                            PoreFraction(Tube(Inf), 0.0, 1.0))
    @test elasticity(circular, 0.01).converged

    mixed = ElasticModel(
        0.66,
        0.4,
        0.2,
        PoreFraction(Film(0.01), 0.1, 0.1),
        PoreFraction(Tube(0.0), 0.1, 0.2),
        PoreFraction(Spheroid(0.1), 0.2, 0.3),
    )
    @test elasticity(mixed, 0.01).converged

    # Appendix A defines the auxiliary isolated fraction as phi_mi/(1-phi_mi),
    # independently of the connected fraction.
    auxiliary_model = ElasticModel(
        0.66, 0.4, 0.2, PoreFraction(Spheroid(1.0), 0.25, 0.75)
    )
    phi = 0.1
    K, mu = auxiliary_model.matrix_bulk, auxiliary_model.matrix_shear
    bulk_sum, shear_sum = ModVisc._elastic_sums(
        auxiliary_model.components, K, mu, phi, auxiliary_model, Val(:isolated)
    )
    component = first(auxiliary_model.components)
    nu = poisson_ratio(K, mu)
    R = ModVisc.resistance_ratio(K, mu)
    theta = ModVisc._bulk_compliance(component.geometry, K, mu, nu, R)
    A = ModVisc._elastic_shear_compliance(
        component.geometry, K, mu, nu, R,
        auxiliary_model.fluid_bulk, auxiliary_model.matrix_bulk, theta,
    )
    modified = component.isolated * phi / (1 - component.isolated * phi)
    delta = inv(auxiliary_model.fluid_bulk) - inv(auxiliary_model.matrix_bulk)
    @test bulk_sum ≈ delta * modified /
          (1 + inv(theta) * (inv(auxiliary_model.fluid_bulk) - inv(K)))
    @test shear_sum ≈ A * modified

    model32 = ElasticModel(0.66f0, 0.4f0, 0.2f0,
                           PoreFraction(Spheroid(0.1), 0f0, 1f0))
    options32 = SolverOptions(relative_tolerance=1f-6)
    result32 = @inferred elasticity(model32, 0.01f0, options32)
    @test result32 isa ElasticResult{Float32}
    elasticity(model32, 0.01f0, options32)
    @test @allocated(elasticity(model32, 0.01f0, options32)) == 0

    profile = elasticity_profile(connected, [0.001, 0.01, 0.1])
    @test all(profile.converged)
    @test all(diff(profile.unrelaxed_bulk) .< 0)
    @test all(profile.unrelaxed_p_speed .> profile.unrelaxed_s_speed)
    @test elasticity(connected, [0.001, 0.01]).relaxed_bulk ==
          elasticity_profile(connected, [0.001, 0.01]).relaxed_bulk
end


@testset "geometry kernel limits" begin
    # Continuous incompressible limit of equation A9.
    tube = Tube(0.0)
    theta_limit = ModVisc._bulk_compliance(tube, Inf, 1.0, 0.5, 0.0)
    @test theta_limit ≈ 3.0

    # Near-spherical auxiliaries approach q=2/3 and g=-2/5 without cancellation.
    q, g = ModVisc._spheroid_terms(1 - 1e-12)
    @test q ≈ 2 / 3 rtol=1e-11
    @test g ≈ -2 / 5 rtol=1e-11

    q32, g32 = ModVisc._spheroid_terms(0.99f0)
    qbig, gbig = ModVisc._spheroid_terms(BigFloat(0.99f0))
    @test q32 ≈ Float32(qbig) rtol=2eps(Float32)
    @test g32 ≈ Float32(gbig) rtol=2eps(Float32)

    sphere = Spheroid(1.0)
    near_sphere = Spheroid(1 - 1e-10)
    K, mu = 2.0, 1.0
    nu = poisson_ratio(K, mu)
    R = ModVisc.resistance_ratio(K, mu)
    @test ModVisc._bulk_compliance(near_sphere, K, mu, nu, R) ≈
          ModVisc._bulk_compliance(sphere, K, mu, nu, R) rtol=1e-8
    @test ModVisc._shear_compliance(near_sphere, K, mu, nu, R) ≈
          ModVisc._shear_compliance(sphere, K, mu, nu, R) rtol=1e-8
end

@testset "viscous self-consistent solver" begin
    sphere = ViscousModel(1.0, MeltFraction(Spheroid(1.0), 1.0))
    result = @inferred viscosity(sphere, 0.01)
    @test result.converged
    @test result.bulk ≈ 129.79264214046762 rtol=1e-12
    @test result.shear ≈ 0.9832775919732443 rtol=1e-12
    @test result.poisson ≈ poisson_ratio(result.bulk, result.shear)
    nu = poisson_ratio(result.bulk, result.shear)
    R = ModVisc.resistance_ratio(result.bulk, result.shear)
    cb, cs = ModVisc._component_compliances(
        sphere.components, result.bulk, result.shear, 0.01, nu, R
    )
    @test result.bulk ≈ inv(cb) rtol=2e-10
    @test result.shear ≈ inv(1 + cs) rtol=2e-10

    # The spherical low-porosity bulk asymptote is 4(1-phi)/(3phi).
    @test result.bulk ≈ 4 * (1 - 0.01) / (3 * 0.01) rtol=0.02

    tapered = ViscousModel(1.0, MeltFraction(Tube(0.0), 1.0))
    circular = ViscousModel(1.0, MeltFraction(Tube(Inf), 1.0))
    @test viscosity(tapered, 0.01).converged
    @test viscosity(circular, 0.01).converged

    # Figure 2: tapered tubes disaggregate near 0.20; spheres near 0.50.
    @test viscosity(tapered, 0.19).shear > 0
    @test viscosity(tapered, 0.21).shear < 1e-8
    @test viscosity(sphere, 0.49).shear > 0
    @test viscosity(sphere, 0.51).shear < 1e-8

    spheroid = ViscousModel(1.0, MeltFraction(Spheroid(0.1), 1.0))
    critical = 2.2455 * 0.1 / (1 + 3.45 * 0.1^1.29)
    @test viscosity(spheroid, 0.95critical).shear > 0
    @test viscosity(spheroid, 1.05critical).shear < 1e-8

    viscosity(sphere, 0.01) # compile before measuring
    @test @allocated(viscosity(sphere, 0.01)) == 0

    sphere32 = ViscousModel(1.0f0, MeltFraction(Spheroid(0.1f0), 1.0f0))
    options32 = SolverOptions(relative_tolerance=1.0f-6)
    result32 = @inferred viscosity(sphere32, 0.01f0, options32)
    @test result32 isa ViscosityResult{Float32}
    @test @inferred(viscosity(sphere32, 0.01f0)) isa ViscosityResult{Float32}
    viscosity(sphere32, 0.01f0, options32)
    @test @allocated(viscosity(sphere32, 0.01f0, options32)) == 0
end

@testset "profiles and published spheroid fit" begin
    sphere = ViscousModel(1.0, MeltFraction(Spheroid(1.0), 1.0))
    scalar_profile = @inferred viscosity_profile(sphere, 0.01)
    @test scalar_profile isa ViscosityResult{Float64}
    @test scalar_profile == viscosity(sphere, 0.01)
    @test @allocated(viscosity_profile(sphere, 0.01)) == 0

    porosity = 10.0 .^ range(-4, -1; length=8)
    profile = viscosity_profile(sphere, porosity)
    @test size(profile.bulk) == size(porosity)
    @test all(profile.bulk .> profile.shear)
    @test all(diff(profile.bulk) .< 0)
    @test all(diff(profile.shear) .< 0)
    @test viscosity(sphere, porosity).bulk == profile.bulk

    fit = spheroid_fit(1e-6, 1.0)
    @test fit.bulk ≈ 4 / (3e-6) rtol=1e-4
    @test fit.shear ≈ 1 rtol=1e-4
    @test spheroid_fit(0.51, 1.0) == (bulk=0.0, shear=0.0)
    @test @inferred(spheroid_fit(0.01f0, 0.1f0, 1.0f0)) isa
          @NamedTuple{bulk::Float32, shear::Float32}
end
