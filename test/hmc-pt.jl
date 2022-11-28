using DrWatson
@quickactivate

using main
using Distributions, BenchmarkTools, MCMCChains, Plots


begin
    Σ = [1.0 0.5; 0.5 1.0]
    R = [0.0 -1.0; 1.0 0.0]
    model = Model(
        ξ=MixtureModel([MvNormal([3.0, 3.0], Σ), MvNormal([-3.0, -3.0], R * Σ * R')])
    )
end

begin
    num = 4
    met = HMC(ϵ=0.5, L=50)
    tmp = InitializeState(rand(model.ξ), met, model)
    pt_tmp = fill(tmp, num)
    temps = 5.0 .^ (-1 .* [0:(num-1)...])
    pt_tmp = map((x, y) -> OnePath(x, met, model, τ=y)[1], pt_tmp, temps)
    paths = path2pts.(pt_tmp)

    cls = palette(:tab10)[1:num]

    a = @animate for i in eachindex(paths[1])
        plt = contourplot(model, l=10.0)
        for j in eachindex(paths)
            plt = plot(plt, paths[j][1:i], c=cls[j], lw=3)
            plt = scatter(plt, paths[j][i], c=cls[j], ms=5)
        end
        plt
    end
    gif(a, fps=10)
end