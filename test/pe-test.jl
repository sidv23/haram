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
    s1, a1 = mcmc(HaRAM(ϵ=0.5, L=20, γ=0.1), model; n=1e4, n_burn=1e4)
    X1 = s1[a1, :]
end

begin
    s2, w2, a2 = mcmc(PEHMC(ϵ=0.65, L=20, N=3), model; n=1e3, n_burn=1e2)
    X2, w2 = [s2[a2, :]...] |> a2m, [w2[a2, :]...]
end

plot(
    scatter(contourplot(model), X1 |> m2t, c=:violet, ma=0.1),
    scatter(contourplot(model), X2 |> m2t, ms=[w2...] * 10, c=:green, ma=0.2),
    size=(800, 400)
)

plot(
    plot(X1[:, 1], bins=50, normalize=true),
    plot(X2[:, 1], weights=w2, bins=50, normalize=true),
    layout=(1, 2),
    size=(800, 400)
)

plot(
    begin
        histogram(X1[:, 1], bins=50, normalize=true)
        plot!(x -> 0.5 * (pdf(Normal(-3), x) + pdf(Normal(3), x)), lw=3)
    end,
    begin
        histogram(X2[:, 1], weights=w2, bins=50, normalize=true)
        plot!(x -> 0.5 * (pdf(Normal(-3), x) + pdf(Normal(3), x)), lw=3)
    end,
    size=(800, 400)
)

W2_minibatch(X1, Matrix(rand(model.ξ, size(X1, 1))'), iters=10, k=100, N=32)
W2_minibatch(X2, Matrix(rand(model.ξ, size(X2, 1))'), iters=10, k=100, N=32)