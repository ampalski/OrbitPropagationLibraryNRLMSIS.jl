@testset "Alt Sweep" begin
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

@testset "Northern Longitude Sweep" begin
    jd0, _ = datevec2jdate([2020.0, 3, 21, 8, 0, 0])
    # Might need to add a few seconds there to keep the day on the 23rd and not go back to the 22nd when converting to UT1
    lat = 55.0
    alt = 250.0

    # results are cm^-3
    lon = 0.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.72e-14, atol = 1.0e-14)
    # @test isapprox(density, 6.56e-10, atol = 1e-9)

    lon = 20.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 4.026e-14, atol = 1.0e-14)

    lon = 40.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 4.306e-14, atol = 1.0e-14)

    lon = 60.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 4.538e-14, atol = 1.0e-14)

    lon = 80.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 4.652e-14, atol = 1.0e-14)

    lon = 100.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 4.554e-14, atol = 1.0e-14)

    lon = 120.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 4.264e-14, atol = 1.0e-14)

    lon = 140.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.968e-14, atol = 1.0e-14)

    lon = 160.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.809e-14, atol = 1.0e-14)

    lon = 180.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.739e-14, atol = 1.0e-14)

    lon = 200.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.659e-14, atol = 1.0e-14)

    lon = 220.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.562e-14, atol = 1.0e-14)

    lon = 240.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.489e-14, atol = 1.0e-14)

    lon = 260.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.432e-14, atol = 1.0e-14)

    lon = 280.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.354e-14, atol = 1.0e-14)

    lon = 300.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.271e-14, atol = 1.0e-14)

    lon = 320.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.275e-14, atol = 1.0e-14)

    lon = 340.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.439e-14, atol = 1.0e-14)
end

@testset "Southern Longitude Sweep" begin
    jd0, _ = datevec2jdate([2020.0, 9, 21, 16, 0, 0])
    # Might need to add a few seconds there to keep the day on the 23rd and not go back to the 22nd when converting to UT1
    lat = -25.0
    alt = 450.0

    # results are cm^-3
    lon = 0.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 4.818e-16, atol = 1.0e-16)
    # @test isapprox(density, 6.56e-10, atol = 1e-9)

    lon = 20.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 4.187e-16, atol = 1.0e-16)

    lon = 40.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.305e-16, atol = 1.0e-16)

    lon = 60.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 2.456e-16, atol = 1.0e-16)

    lon = 80.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.835e-16, atol = 1.0e-16)

    lon = 100.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.514e-16, atol = 1.0e-16)

    lon = 120.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.41e-16, atol = 1.0e-16)

    lon = 140.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.374e-16, atol = 1.0e-16)

    lon = 160.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.309e-16, atol = 1.0e-16)

    lon = 180.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.238e-16, atol = 1.0e-16)

    lon = 200.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.255e-16, atol = 1.0e-16)

    lon = 220.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.438e-16, atol = 1.0e-16)

    lon = 240.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 1.817e-16, atol = 1.0e-16)

    lon = 260.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 2.357e-16, atol = 1.0e-16)

    lon = 280.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.001e-16, atol = 1.0e-16)

    lon = 300.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 3.715e-16, atol = 1.0e-16)

    lon = 320.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 4.424e-16, atol = 1.0e-16)

    lon = 340.0
    density = getDensity(jd0, lat, lon, alt) / 1000
    @test isapprox(density, 4.893e-16, atol = 1.0e-16)
end
