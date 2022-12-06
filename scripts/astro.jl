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

s2v(s::Symbol, v::Any) = @eval (($s) = ($v))
s2v(s::Symbol, v::D) where {D<:Distribution} = @eval (($s) ~ ($v))

Dists = (;
    # Δ=Uniform(-1178.939, 1179.939),
    # β=Uniform(-60, 60),
    # μ=Uniform(-39, 30),
    Δ=Truncated(Normal(0, 1e5), -1178.939, 1179.939),
    β=Truncated(Normal(0, 1e5), -60, 60),
    μ=Truncated(Normal(0, 1e5), -39, 30)
    # σ2=InverseGamma(1.0, 2e-7),
    # τ=InverseGamma(1.0, 1.0)
);
K = keys(Dists)
B = map(x -> bijector(x), Dists)
Binv = map(x -> inverse(x), B)

B_map(x) = map((f, x) -> f(x), B, x)
Binv_map(x) = map((f, x) -> f(x), Binv, x)

v2n(V::Vector{T}, S::NTuple) where {T<:Real} = NamedTuple{S}(V)

df = CSV.File(datadir("q0957usno.csv")) |> DataFrame

@model function π_prior(Dists)
    Δ ~ Dists.Δ
    β ~ Dists.β
    μ ~ Dists.μ
    # σ2 ~ Dists.σ2 |> transformed
    # τ ~ Dists.τ |> transformed
end

function log_likelihood(Binv, θ::Vector{T}, df) where {T<:Real}
    Δ = θ[1]
    β = θ[2]
    μ = θ[3]
    # σ2 = θ[4] |> Binv[4]
    # τ = θ[5] |> Binv[5]

    σ2 = 2e-7
    τ = 0.5

    n, m = size(df)
    var0 = τ * σ2 / 2

    T2n = [df.time; df.time .- Δ]
    ord = sortperm(T2n)

    tΔ = T2n[ord]
    ind = [fill(0, n); fill(1, n)][ord]
    lc = [df.lca; df.lcb .- β][ord]
    se = [df.sea; df.seb][ord]

    x = lc .- μ

    a = [0; exp.(-1 .* (tΔ[2:end] .- tΔ[1:end-1]) ./ τ)]
    B = fill(T(+1), 2n)
    m = fill(T(-1), 2n)

    B[1] = se[1]^2 / (se[1]^2 + var0)
    m[1] = (1 - B[1]) * x[1]

    for i in 2:length(T2n)
        B[i] = se[i]^2 / (
            se[i]^2 +
            a[i-1]^2 * (1 - B[i-1]) * se[i-1]^2 +
            var0 * (1 - a[i-1]^2)
        )
        m[i] = (1 - B[i]) * x[i] + B[i] * a[i-1] * m[i-1]
    end

    D = MvNormal(a .* m, se ./ B)
    return logpdf(D, x)
end

function make_model(Dists, log_likelihood)
    K = keys(Dists)
    B = map(x -> bijector(x), Dists)
    Binv = map(x -> inverse(x), B)

    U(x) = min(-logprior(π_prior(Dists), v2n(x, K)) - log_likelihood(Binv, x, df), floatmax())
    # dU(x) = Zygote.gradient(x_ -> Zygote.forwarddiff(U, x_), x)[1]
    dU(x) = ForwardDiff.gradient(U, x)
    f(x) = max(exp(-U(x)), floatmin(Float64))
    # g(x) = Zygote.gradient(x_ -> Zygote.forwarddiff(f, x_), x)
    g(x) = ForwardDiff.gradient(f, x)
    return main.Model(ξ=Dists, d=length(Dists), f=f, g=g, U=U, dU=dU)
end


model = make_model(Dists, log_likelihood)

y = [1000.0, 10.0, 10.0];
model.dU(y);


S, q0 = main.DualAveraging(
    main.DualAverage(λ=1e-2, δ=0.75),
    main.HaRAM(),
    # main.HMC(),
    model, n_burn=1000
)

begin
    s, a = mcmc(
        PT(τ=[1.0, 2.5, 5.0, 7.5, 10.0]),
        # DualAverage(λ=0.0005, δ=0.65),
        # DualAverage(λ=5, δ=0.5),
        # S,
        # main.HMC(ϵ=0.001, L=5),
        HaRAM(ϵ=5e-4, L=100, γ=50.25),
        # main.HaRAM(ϵ=0.0005, L=20, γ=0.85),
        # main.HaRAM(ϵ=0.0005, L=20, γ=0.85),
        model; n=1000, n_burn=100,
        init = [1100, 5.0, 5.0],
        # init= map((f, x) -> f(x), B, [420, 0.6, -18, 0.01^2, 200])
        # init= map((f, x) -> f(x), Binv, [420, 0.6, -18, 0.01^2, 200])
    )
    # params = map((f, x) -> f(x), Binv, eachcol(s[a, :]))
    params = map((x) -> x, eachcol(s[a, :]))
    plot([histogram(p, bins=200) for p in params]..., layout=(length(params), 1), size=(5, length(params)) .* 200)
end


foo(x, a) = model.U(map((f, x) -> f(x), Binv, [x; a]))
begin
    a =
        plot(x -> foo(x, a), -2000, 2000)
end
