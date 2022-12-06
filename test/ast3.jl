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
end

Dists = (;
    Δ=Uniform(-1178.939, 1179.939),
    β=Uniform(-60, 60),
    μ=Uniform(-39, 30),
    τ=InverseGamma(1.0, 1.0),
    σ2=InverseGamma(1.0, 2e-7)
);
K = keys(Dists)
B = map(x -> bijector(x), Dists)
Binv = map(x -> inverse(x), B)

fmap(x, B) = map((f, x) -> f(x), B, x)
v2n(V::Vector{T}, S::NTuple) where {T<:Real} = Zygote.@ignore NamedTuple{S}(V)
v2n(V::Vector{T}, S::Vector{Symbol}) where {T<:Real} = NamedTuple{(; S)}(V)

df = CSV.File(datadir("q0957usno.csv")) |> DataFrame

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

function make_model(Dists::NamedTuple, df::DataFrame)
    K = keys(Dists)
    B = map(x -> bijector(x), Dists)
    Binv = map(x -> inverse(x), B)

    ℓ(x) = logjoint(π_delay(df), v2n(x, K))
    U(x) = isfinite(ℓ(x)) ? -ℓ(x) : 1e100
    dU(x) = ForwardDiff.gradient(U, x)
    f(x) = max(exp(-U(x)), 1e-100)
    g(x) = ForwardDiff.gradient(f, x)
    return main.Model(ξ=Dists, d=length(Dists), f=f, g=g, U=U, dU=dU)
end

model = make_model(Dists, df)

model.dU(randn(5))

begin
    mass = @pipe [3.0, 2.0, 1.5, -2.0, -7.0] .|> 10^(_ * +1.0)
    S, q0 = main.DualAveraging(
        main.DualAverage(λ=1e-8, δ=0.65),
        # main.HMC(),
        main.HaRAM(),
        model, n_burn=1000,
        m = mass
    )
    fmap(q0.q, Binv)
end

begin
    s, a = mcmc(
        # PT(τ=[1.0, 2.5, 5.0, 7.5, 10.0]),
        # PT(τ=[100.0, 100.0]),
        # main.HMC(ϵ=0.9, L=10),
        main.HaRAM(ϵ=5e-4, L=12, γ=0.15),
        # S,
        model; n=10000, n_burn=100,
        # m=[1e-1, 1e0, 1e0, 2e7, 2e5],
        init=fmap([420; 0.1; rand(3)], B),
        # m=q0.m
    )
    params = map((f, x) -> f(x), Binv, eachcol(s[a, :]))
    # params = map((x) -> x, eachcol(s[a, :]))
    plot([histogram(p, bins=200) for p in params]..., layout=(length(params), 1), size=(5, length(params)) .* 200)
end


c = sample(log_likelihood(df), NUTS(), 5000)

map(k -> Binv[k](c[k]), K)
plot(histogram.(map(k -> Binv[k](c[k]), K))..., layout=(length(K), 1), size=(7, length(K)) .* 150)


begin
    plot(df.time, df.lca, marker=:o, label="")
    plot!(df.time .+ 686, df.lcb, marker=:o, label="")
end
