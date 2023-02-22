using DrWatson
@quickactivate


begin
  using main
  using BlockDiagonals, LinearAlgebra, Plots, ProgressMeter
  using Distributions, MCMCChains, Random
end

begin
  function generate_model(; m=5.0, s=0.5, d=10)
    R = prod(
    # [BlockDiagonal([diagm(ones(j)), [0.0 -1.0; 1.0 0.0], diagm(ones(d - j - 2))]) for j in 0:d-2]
      [
      BlockDiagonal(
        [diagm(ones(2j)),
        [0.0 -1.0; 1.0 0.0],
        diagm(ones(d - 2j - 2))]
      ) for j in 0:round(Int, d / 2 - 1)
    ]
    )

    Σ₁ = [s^abs(i - j) for i in 1:d, j in 1:d]
    Σ₂ = R * Σ₁ * R'
    μ = [-m .* ones(d), m .* ones(d)]
    Σ = [Σ₁, Σ₂]

    ξ = MixtureModel(
      [MvNormal(x, y) for (x, y) in zip(μ, Σ)],
      [0.5, 0.5]
    )
    return Model(ξ=ξ)
  end

  function acfplots(chains, names, lags=0:2:20; kwargs...)
    plt = plot(0, 0)
    for (x, n) in zip(chains, names)
      plt = plot(plt, mean(1 .* (autocor(x, lags=[lags...])[:, :]), dims=1)', label=n; kwargs...)
    end
    return plt
  end

  function scatterplots(xs, names; baseplt=plot(0, 0, label=nothing), l=400, kwargs...)
    L = length(xs)
    ds = sample(1:size(xs[1], 2), 2, replace=false)
    plts = [scatter(baseplt, x[:, ds[1]], x[:, ds[2]], label=names[i]) for (x, i) in zip(xs, eachindex(xs))]
    plot(plts..., axes=false, ticks=false; kwargs...)
  end

  function traceplots(xs, names, args...; baseplt=plot(0, 0, label=nothing), l=400, kwargs...)
    L = length(xs)
    ds = sample(1:size(xs[1], 2), 1, replace=false)
    plts = [plot(x[:, ds[1]], label=names[i]; args...) for (x, i) in zip(xs, eachindex(xs))]
    plot(plts...; kwargs...)
  end

  function scatterplot(x; baseplt=plot(0, 0, label=""), kwargs...)
    plt = plot(baseplt, x |> m2t, c=:black, lw=0.1, la=0.25, label="")
    plt = scatter(plt, x |> m2t, c=:orange; kwargs...)
    return plt
  end

  function mean_ess(chains)
    return [ess_rhat(chn)[:, 2] |> mean for chn in chains]
  end

  function w2_minibatch(xs, model; eps=0.1, iters=100, k=100, N=64)
    results = zeros(length(xs))
    for (i, x) in zip(eachindex(xs), xs)
      z = Matrix(rand(model.ξ, size(x, 1))')
      results[i] = W2_minibatch(x, z, eps=eps, iters=iters, k=k, N=N)
    end
    return results
  end
end

d = 2
model = generate_model(d=d);
state = (; p=randn(d), q=randn(d), m=1.0)

@profview for _ in 1:1000 OneStep(state, method, model) end





@code_warntype mcmc(
  DualAverage(λ=10, δ=0.65),
  HMC(),
  model; n=5e3, n_burn=1e3
)