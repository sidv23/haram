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
end

v2n(V::Vector{T}, S::NTuple) where {T<:Real} = NamedTuple{S}(V)
s2v(s::Symbol, v::Any) = @eval (($s) = ($v))
s2v(s::Symbol, v::D) where {D<:Distribution} = @eval (($s) ~ ($v))

Dists = (;
    Δ=Uniform(-1178.939, 1179.939),
    β=Uniform(-60, 60),
    μ=Uniform(-39, 30),
    σ2=InverseGamma(1.0, 2e-7),
    τ=InverseGamma(1.0, 1.0)
);
K = keys(Dists)
B = map(x -> bijector(x), Dists)
Binv = map(x -> inverse(x), B)

df = CSV.File(datadir("q0957usno.csv")) |> DataFrame

@model function π_prior()
    Δ ~ Dists.Δ |> transformed
    β ~ Dists.β |> transformed
    μ ~ Dists.μ |> transformed
    σ2 ~ Dists.σ2 |> transformed
    τ ~ Dists.τ |> transformed
end

function log_likelihood(Binv, θ::Vector{T}, df) where {T<:Real}
    Δ = θ[1] |> Binv[1]
    β = θ[2] |> Binv[2]
    μ = θ[3] |> Binv[3]
    σ2 = θ[4] |> Binv[4]
    τ = θ[5] |> Binv[5]

    n, m = size(df)
    var0 = τ * σ2 / 2

    T2n = [df.time; df.time .- Δ]
    ord = sortperm(T2n)

    tΔ = T2n[ord]
    ind = [fill(0, n); fill(1, n)][ord]
    lc = [df.lca; df.lcb .- β][ord]
    se = [df.sea; df.seb][ord]

    x = lc .- μ

    a = exp.(-1 .* (tΔ[2:end] .- tΔ[1:end-1]) ./ τ)
    B = fill(T(+1), 2n)
    m = fill(T(-1), 2n)
    M = fill(T(-1), 2n)
end




foo(x) = ForwardDiff.gradient(x_ -> logprior(π_prior(), v2n(x_, K)), x)

tmp1 = [df.time; df.time .- 1000]
tmp2 = sortperm(tmp1)
tmp1[tmp2]'

[df.lca; df.lcb]'
[df.lca; df.lcb][tmp2]'


mod = astro(df) | (; X=0.0)
logjoint(mod, (; Δ=1.0, μ=1.0, β=1.0, σ2=1.0, τ=1.0))


tmp = (; abc=Uniform(0, 1))
map((k, v) -> s2v(k, v), keys(tmp), values(tmp))




