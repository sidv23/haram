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
    σ2=InverseGamma(1.0, 2e-7),
);
K = keys(Dists)
B = map(x -> bijector(x), Dists)
Binv = map(x -> inverse(x), B)

fmap(x, B) = map((f, x) -> f(x), B, x)
v2n(V::Vector{T}, S::NTuple) where {T<:Real} = Zygote.@ignore NamedTuple{S}(V)
v2n(V::Vector{T}, S::Vector{Symbol}) where {T<:Real} = NamedTuple{(; S)}(V)

df = CSV.File(datadir("q0957usno.csv")) |> DataFrame

@model function π_posterior(Dists::NamedTuple, df::DataFrame, ::Type{T}=Float64) where {T}
    Δ ~ Dists.Δ |> transformed
    β ~ Dists.β |> transformed
    μ ~ Dists.μ |> transformed
    τ ~ Dists.τ |> transformed
    σ2 ~ Dists.σ2 |> transformed
    Δ0, β0, μ0, τ0, σ20 = fmap([Δ, β, μ, τ, σ2], Binv)

    σ20 = 1e-7
    Δ0, β0, μ0, τ0 = Δ, β, μ, τ

    n, m = size(df)
    var0 = τ0 * σ20 / 2

    T2n = [df.time; df.time .- Δ0]
    ord = sortperm(T2n)
    ord_inv = invperm(ord)

    tΔ = T2n[ord]
    a = [0; exp.(-1 .* (tΔ[2:end] .- tΔ[1:end-1]) ./ τ0)]

    X = fill(T(0.0), 2n)
    X[1] = rand(Normal(0, var0))

    for j in 2:n
        X[j] = rand(Normal(μ + a[j] * (X[j-1] - μ0), var0 * (1 - a[j]^2)))
    end

    df.lca ~ MvNormal(X[ord_inv][1:n], df.sea)
    df.lcb ~ MvNormal(X[ord_inv][n+1:2n] .+ β0, df.sea)
end

function make_model(Dists::NamedTuple, df::DataFrame)
    K = keys(Dists)
    B = map(x -> bijector(x), Dists)
    Binv = map(x -> inverse(x), B)

    U(x) = -max(logjoint(π_posterior(Dists, df), ), 1/eps())
    dU(x) = ForwardDiff.gradient(U, x)
    f(x) = max(exp(-U(x)), floatmin(Float64))
    g(x) = ForwardDiff.gradient(f, x)
    return main.Model(ξ=Dists, d=length(Dists), f=f, g=g, U=U, dU=dU)
end

model = make_model(Dists, df)

S, q0 = main.DualAveraging(
    main.DualAverage(λ=1, δ=0.75),
    # main.HMC(),
    main.HaRAM(),
    model, n_burn=100,
    ϵ0=1e-5,
    # m=[1e7, 1e7, 1e7, 1e7, 1e7]
    m=fill(1e10, 5)
);

state = InitializeState(randn(5), S, model, m=1.0)
newstate, _ = OneStep(state, main.HMC(ϵ=1, L=10), model)
fmap(newstate.q, Binv)

begin
    s, a = mcmc(
        # PT(τ=[1.0, 2.5, 5.0, 7.5, 10.0]),
        # main.HMC(ϵ=0.9, L=10),
        main.HaRAM(ϵ=1.5, L=12, γ=0.75),
        model; n=10000, n_burn=100,
        m=[1e-1, 1e0, 1e0, 2e7, 2e5]
    )
    params = map((f, x) -> f(x), Binv, eachcol(s[a, :]))
    # params = map((x) -> x, eachcol(s[a, :]))
    plot([histogram(p, bins=200) for p in params]..., layout=(length(params), 1), size=(5, length(params)) .* 200)
end


begin
    # del = -700
    del = -495.5
    T1 = map(x -> mod(round.(Int, x), range(extrema(round.(Int, df.time))...)), df.time)
    T2 = map(x -> mod(round.(Int, x), range(extrema(round.(Int, df.time))...)), df.time .- del)
    ord = sortperm(T2)

    plot(T1, df.lca, ribbon=df.sea, m=:o)
    plot!(T2[ord], df.lcb[ord], ribbon=df.sea, m=:o)
    vline!([minimum(df.time), 4825, 4940, 5185, 5360, 5550, 5590])
end



foo(x) = mean([model.U([x; values(rand(log_likelihood(df)))[2:end]...]) for _ in 1:1000])

X = 1000:2:1200
y = @showprogress [foo(x) for x in X]
plot(X, y)