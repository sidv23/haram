using DrWatson, Revise;
@quickactivate "haram";

begin
  using main
  using DynamicPPL, Bijectors
  using Zygote
  using LinearAlgebra, Plots, Pipe
  using CSV, DataFrames
  using Turing
  using Libtask: TArray
  using Libtask
  using ForwardDiff
  using ProgressMeter

  using Distributions
  using FillArrays
  using StatsPlots

  using LinearAlgebra
  using Random
end

Random.seed!(2022)

begin
  N = 1000
  K = 15
  r = 1.5
  σ = 0.2

  function f(k, μ=0.0)
    θ = range(2π / k, 2π, round(Int, k))
    return [μ .+ (cos(t), sin(t)) for t in θ]
  end

  μ = [f(K / 3, -r); f(2K / 3, r)]
  v = sample(1:K, N, replace=true)
  w = v .<= K / 3

  x = [rand(MvNormal([μ[v[i]]...], σ / sqrt(v[i]))) for i in eachindex(v)]
  x = x |> a2m

  scatter(x |> m2t, c=1 .+ w)
end



@model function gaussian_mixture_model(x)
  D, N = size(x)

  K ~ Uniform(1, N)
  μ ~ MvNormal(Zeros(K), I)
  w ~ Dirichlet(K, 1.0)


  distribution_assignments = Categorical(w)
  distribution_clusters = [MvNormal(Fill(μₖ, D), I) for μₖ in μ]

  k = Vector{Int}(undef, N)
  for i in 1:N
    k[i] ~ distribution_assignments
    x[:, i] ~ distribution_clusters[k[i]]
  end

  return k
end

model = gaussian_mixture_model(x);
sampler = Gibbs(PG(100, :k), Turing.HMC(0.05, 10, :μ, :w))

nsamples = 100
nchains = 3
chains = sample(model, sampler, MCMCThreads(), nsamples, nchains);
plot(chains)