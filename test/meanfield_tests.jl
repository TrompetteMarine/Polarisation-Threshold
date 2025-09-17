@testset "MeanField" begin
    params = ModelParams()
    system = MeanField.mean_field_system(params)
    @test isa(system, ContinuousDynamicalSystem)

    orbit = MeanField.orbit_data(system; steps=20, discard=5)
    @test haskey(orbit, :times)
    @test haskey(orbit, :states)
    @test length(orbit.times) == size(orbit.states, 2)

    factory = κ -> MeanField.mean_field_system(update(params; κ=κ))
    data = MeanField.bifurcation_data(factory, 0.1:0.1:0.3; steps=50, discard=10)
    @test length(data.parameters) == length(data.observable)
end
