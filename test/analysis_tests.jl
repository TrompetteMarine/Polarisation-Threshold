@testset "Analysis" begin
    params = ModelParams()
    mock_data = Dict(:G => rand(20), :variance_u => abs.(randn(20)))
    κ_est = Analysis.kappa_star_from_abm(mock_data, params; window=5)
    @test isfinite(κ_est)

    orbit = (times=collect(0:0.1:1.0), states=rand(1, 11))
    stability = Analysis.stability_fraction(orbit, params)
    @test 0.0 <= stability <= 1.0

    comp = Analysis.compare_thresholds(mock_data, orbit, params)
    @test haskey(comp, :abm)
    @test haskey(comp, :deterministic)
    @test haskey(comp, :difference)
end
