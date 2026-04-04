@testset "Test1" begin
    jd0, _ = datevec2jdate([2020.0, 2, 23, 0, 0, 0])
    # Might need to add a few seconds there to keep the day on the 23rd and not go back to the 22nd when converting to UT1
    lat = 55.0
    lon = 45.0

    # results are cm^-3
    alt = 100.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 6.563e-10, atol = 1.0e-11)
    # @test isapprox(density, 6.56e-10, atol = 1e-9)

    alt = 150.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.648e-12, atol = 1.0e-13)

    alt = 200.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.679e-13, atol = 1.0e-14)

    alt = 250.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 2.874e-14, atol = 1.0e-15)

    alt = 300.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 6.391e-15, atol = 3.0e-16)

    alt = 350.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.652e-15, atol = 2.0e-16)

    alt = 400.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 4.723e-16, atol = 4.0e-17)

    alt = 450.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.504e-16, atol = 2.0e-17)

    alt = 500.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 5.603e-17, atol = 5.0e-18)

    alt = 550.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 2.578e-17, atol = 2.0e-18)

    alt = 600.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.462e-17, atol = 1.0e-18)

    alt = 650.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 9.617e-18, atol = 3.0e-19)

    alt = 700.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 6.878e-18, atol = 3.0e-19)

    alt = 750.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 5.134e-18, atol = 3.0e-19)

    alt = 800.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.921e-18, atol = 1.0e-19)

    alt = 850.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.038e-18, atol = 1.0e-19)

    alt = 900.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 2.382e-18, atol = 1.0e-19)

    alt = 950.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.886e-18, atol = 1.0e-19)

    alt = 1000.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.509e-18, atol = 1.0e-19)
end
