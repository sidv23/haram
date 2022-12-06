### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 72323508-742e-11ed-14e5-9f08721c77e0
import Pkg; Pkg.activate("/storage/work/s/suv87/julia/haram/")

# ╔═╡ 85bd3a37-cbe3-4bd1-9302-c88ef0399370
begin
	using PlutoUI
	using main
	using DrWatson
	using Random
    using DynamicPPL, Bijectors
    using Zygote
    using LinearAlgebra, Plots, Pipe
    using CSV, DataFrames
    using Turing
    using Libtask: TArray
    using Libtask
    using ForwardDiff
	using ProgressMeter
	using Optim
	using ProgressLogging
end

# ╔═╡ ad0b9f2d-a9a9-48aa-bc18-c383ccb30845
Pkg.add("PlutoUI")

# ╔═╡ 8860958a-b685-4899-97c6-5b8d64c43c86
begin
	Dists = (;
	    Δ=Uniform(-1278.939, 1279.939),
	    β=Uniform(-60, 60),
	    μ=Uniform(-39, 30),
	    τ=InverseGamma(1.0, 1.0),
	    σ2=InverseGamma(1.0, 2e-7)
	);
	K = keys(Dists)
	B = map(x -> bijector(x), Dists)
	Binv = map(x -> inverse(x), B)
end

# ╔═╡ dd09b8f9-bbfa-49b0-8b2c-1c92aca940b0
begin
	fmap(x, B) = map((f, x) -> f(x), B, x)
	v2n(V::Vector{T}, S::NTuple) where {T<:Real} = Zygote.@ignore NamedTuple{S}(V)
	v2n(V::Vector{T}, S::Vector{Symbol}) where {T<:Real} = NamedTuple{(; S)}(V)
	df = CSV.File(datadir("q0957usno.csv")) |> DataFrame
end

# ╔═╡ a9e274b4-337b-4522-815f-f563a35c88d5
# @model function π_delay(df::DataFrame, ::Type{T}=Float64) where {T}
#     Δ ~ Dists.Δ |> transformed
#     β ~ Dists.β |> transformed
#     μ ~ Dists.μ |> transformed
#     τ ~ Dists.τ |> transformed
#     σ2 ~ Dists.σ2 |> transformed
#     Δ0, β0, μ0, τ0, σ20 = fmap([Δ, β, μ, τ, σ2], Binv)

#     n, m = size(df)
#     var0 = τ0 * σ20 / 2

#     T2n = [df.time; df.time .- Δ0]
#     ord = sortperm(T2n)
#     ord_inv = invperm(ord)

#     se = [df.sea; df.seb][ord]
#     x = [df.lca; df.lcb][ord]

#     tΔ = T2n[ord]
#     a = [0; exp.(-1 .* (tΔ[2:end] .- tΔ[1:end-1]) ./ τ0)]

#     x_hat = fill(T(0.0), 2n)
#     Ω = fill(T(var0), 2n)

#     for i in 2:2n
#         Ω[i] = var0 * (1 - a[i]^2) +
#                (Ω[i-1] * a[i]^2 * (1 - (Ω[i-1] / Ω[i-1] + se[i]^2)))

#         x_hat[i] = (a[i] * x_hat[i-1]) +
#                    ((Ω[i-1] / Ω[i-1] + se[i]^2) * a[i] * (x[i-1] - x_hat[i-1]))
#     end

#     df.lca ~ MvNormal(x_hat[ord_inv][1:n] .+ μ, Ω[ord_inv][1:n])
#     df.lcb ~ MvNormal(x_hat[ord_inv][n+1:2n] .+ β0, Ω[ord_inv][n+1:2n])
# end

# ╔═╡ 5f3e2382-48ca-4807-924a-0e8a19e4f998
@model function π_delay(df::DataFrame, ::Type{T}=Float64) where {T}
    Δ ~ Dists.Δ |> transformed
    β ~ Dists.β |> transformed
    μ ~ Dists.μ |> transformed
    τ ~ Dists.τ |> transformed
    σ2 ~ Dists.σ2 |> transformed
    Δ0, β0, μ0, τ0, σ20 = fmap([Δ, β, μ, τ, σ2], Binv)

    n, m = size(df)
    ω = τ0 * σ20 / 2

    T2n = [df.time; df.time .- Δ0]
    ord = sortperm(T2n)
    ord_inv = invperm(ord)

    se = [df.sea; df.seb][ord]
    x = [df.lca; df.lcb][ord]

    tΔ = T2n[ord]
    a = [0; exp.(-diff(tΔ) ./ τ0)]

    Ω = fill(T(se[1]^2 / (se[1]^2 + ω)), 2n)
    s = fill(T((se[1]^2 + ω) / se[1]), 2n)

    for i in 2:2n
        Ω[i] = (se[i]^2) / (
			(se[i]^2) + (a[i-1]^2 * se[i-1]^2 * (1-Ω[i])) +(ω * (1 - a[i]^2))
		)
		s[i] = (se[i]^2) / (Ω[i]^2)
    end

	x = [0; x[1:end-1]]
    df.lca ~ MvNormal((a .* x)[ord_inv][1:n] .+ μ, s[ord_inv][1:n])
    df.lcb ~ MvNormal((a .* x)[ord_inv][n+1:2n] .+ μ .+ β0, s[ord_inv][n+1:2n])
end

# ╔═╡ 38077261-3763-4f54-93a9-2c9bd8de4a45
logjoint(π_delay(df), v2n(randn(5), K))

# ╔═╡ b978fed2-1787-4150-ada8-cbd8bb1a7f8a
K_small = K[2:end]

# ╔═╡ db2f93d4-62b2-441c-a109-bd73b7d3f2ab
begin
	cond1(y) = π_delay(df) | (; v2n([y...], K_small)...)
	cond2(x) = π_delay(df) | (; Δ = [x...])
	ℓ1(x, y) = loglikelihood(cond1(y), (;Δ=x))
	ℓ2(x, y) = loglikelihood(cond2(x), (; v2n([y...], K_small)...))
end

# ╔═╡ d1610159-1589-47c7-8edd-a5bdadd9ad06
begin
c = sample(
	π_delay(df), 
	Gibbs(
		Turing.NUTS(2, 0.5, :Δ),
		Turing.NUTS(2, 0.5, :β),
		Turing.NUTS(2, 0.5, :μ),
		Turing.NUTS(2, 0.5, :τ),
		Turing.NUTS(2, 0.5, :σ2),
		# PG(10, :μ, :σ2, :τ),
	),
	1000
);
z = map((f, x) -> f(c[x]), Binv, keys(c))
plot(histogram.(z)..., layout=(5, 1), size=(1000, 1000))
end

# ╔═╡ 4cdb2cd0-eb0f-4c6c-a695-f49b80246eeb
θ = map((f, x) -> f(c[x]), Binv, K)

# ╔═╡ 074ca61f-566e-4e32-9dd4-eba00e538e6c
begin
	plot(df.time, df.lca)
	plot!(df.time .+ mean(θ[1]), df.lcb)
end

# ╔═╡ Cell order:
# ╠═72323508-742e-11ed-14e5-9f08721c77e0
# ╠═ad0b9f2d-a9a9-48aa-bc18-c383ccb30845
# ╠═85bd3a37-cbe3-4bd1-9302-c88ef0399370
# ╠═8860958a-b685-4899-97c6-5b8d64c43c86
# ╠═dd09b8f9-bbfa-49b0-8b2c-1c92aca940b0
# ╠═a9e274b4-337b-4522-815f-f563a35c88d5
# ╠═5f3e2382-48ca-4807-924a-0e8a19e4f998
# ╠═38077261-3763-4f54-93a9-2c9bd8de4a45
# ╠═b978fed2-1787-4150-ada8-cbd8bb1a7f8a
# ╠═db2f93d4-62b2-441c-a109-bd73b7d3f2ab
# ╠═ba39f8dd-c3bd-468a-9395-04c00138d825
# ╠═81659fd0-fcc6-465d-9dbd-e6f1d0c9e0a3
# ╠═d1610159-1589-47c7-8edd-a5bdadd9ad06
# ╠═4cdb2cd0-eb0f-4c6c-a695-f49b80246eeb
# ╠═074ca61f-566e-4e32-9dd4-eba00e538e6c
