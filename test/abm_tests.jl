@testset "ABM" begin
    params = ModelParams(N=10, T=1.0, dt=0.1, seed=42)
    model = ABM.initialize_model(params)
    @test model.params == params
    @test !isempty(model.agents)

    result = ABM.run_abm(params; steps=5)
    @test haskey(result.data, :G)
    @test length(result.data[:G]) == 5
    @test haskey(result.data, :variance_u)
end
