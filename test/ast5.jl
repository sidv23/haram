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
end

begin
    Dists = (;
        Δ=Uniform(-1178.939, 1179.939),
        β=Uniform(-60, 60),
        μ=Uniform(-39, 30),
        τ=InverseGamma(1.0, 1.0),
        σ2=InverseGamma(1.0, 2e-7)
    )
    K = keys(Dists)
    B = map(x -> bijector(x), Dists)
    Binv = map(x -> inverse(x), B)

    fmap(x, B) = map((f, x) -> f(x), B, x)
    v2n(V::Vector{T}, S::NTuple) where {T<:Real} = Zygote.@ignore NamedTuple{S}(V)
    v2n(V::Vector{T}, S::Vector{Symbol}) where {T<:Real} = NamedTuple{(; S)}(V)

    df = CSV.File(datadir("q0957usno.csv")) |> DataFrame
end

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
            (se[i]^2) + (a[i-1]^2 * se[i-1]^2 * (1 - Ω[i])) + (ω * (1 - a[i]^2))
        )
        s[i] = (se[i]^2) / (Ω[i]^2)
    end

    x = [0; x[1:end-1]]
    df.lca ~ MvNormal((a.*x)[ord_inv][1:n] .+ μ0, s[ord_inv][1:n])
    df.lcb ~ MvNormal((a.*x)[ord_inv][n+1:2n] .+ μ0 .+ β0, s[ord_inv][n+1:2n])
end


function p1(df::DataFrame, z::Vector{T}) where {T<:Float64}
    K = keys(Dists)
    K2 = K[1:end.!=1]
    ℓ(x) = logjoint(π_delay(df) | v2n([z...], K2), (; Δ=x))
    U(x) = isfinite(ℓ(x)) ? -ℓ(x) : 1e100
    dU(x) = ForwardDiff.derivative(U, x)
    f(x) = max(exp(-U(x)), 1e-100)
    g(x) = ForwardDiff.derivative(f, x)
    return main.Model(ξ=Dists, d=1, f=f, g=g, U=U, dU=dU)
end

function p2(df::DataFrame, z::Vector{T}) where {T<:Float64}
    K = keys(Dists)
    K2 = K[1:end.!=2]
    ℓ(x) = logjoint(π_delay(df) | v2n([z...], K2), (; β=x))
    U(x) = isfinite(ℓ(x)) ? -ℓ(x) : 1e100
    dU(x) = ForwardDiff.derivative(U, x)
    f(x) = max(exp(-U(x)), 1e-100)
    g(x) = ForwardDiff.derivative(f, x)
    return main.Model(ξ=Dists, d=1, f=f, g=g, U=U, dU=dU)
end


function p3(df::DataFrame, z::Vector{T}) where {T<:Float64}
    K = keys(Dists)
    K1 = K[1:end.∈Ref(1:2)]
    K2 = K[1:end.∉Ref(1:2)]
    ℓ(x) = logjoint(π_delay(df) | v2n(z, K1), v2n([x...], K2))
    U(x) = isfinite(ℓ(x)) ? -ℓ(x) : 1e100
    dU(x) = ForwardDiff.gradient(U, x)
    f(x) = max(exp(-U(x)), 1e-100)
    g(x) = ForwardDiff.gradient(f, x)
    return main.Model(ξ=Dists, d=length(Dists) - 2, f=f, g=g, U=U, dU=dU)
end

z1
S1
m1.dU(rand())


function gibbs_step3(df, S, z)
    S1, S2, S3 = S
    z1, z2, z3 = z

    m1 = p1(df, [z2.q; z3.q])
    z1 = RefreshState(z1.q, z1, S1, m1)
    z1 = (; z1..., p=only(z1.p), m=1e5)
    n1, a1 = OneStep(z1, S1, m1)
    b1 = rand() < a1
    z1 = b1 ? n1 : z1

    m2 = p2(df, [z1.q; z3.q])
    z2 = RefreshState(z2.q, z2, S2, m2)
    z2 = (; z2..., p=only(z2.p), m=1e5)
    n2, a2 = OneStep(z2, S2, m2)
    b2 = rand() < a2
    z2 = b2 ? n2 : z2

    m3 = p3(df, [z1.q; z2.q])
    z3 = RefreshState(z3.q, z3, S3, m3)
    n3, a3 = OneStep(z3, S3, m3)
    b3 = rand() < a3
    z3 = b3 ? n3 : z3

    # println(fmap([z1.q; z2.q...], Binv), (a1, a2))

    # return (z1, z2, z3, b1 || b2 || b3)
    return (z1, z2, z3, b1 && b2 && b3)
end



function gibbs_sample3(df, S; n=1000, n_burn=1, init=B[1](420.0))
    S1, S2, S3 = S
    q1, q2, q3 = init, 0.0, randn(3)

    z1 = InitializeState(q1, S1, p1(df, [q2; q3]))
    z1 = (; z1..., q=only(z1.q), p=only(z1.p), m=1.0)

    z2 = InitializeState(q2, S2, p2(df, [q1; q3]))
    z2 = (; z2..., q=only(z2.q), p=only(z2.p), m=1.0)

    z3 = InitializeState(q3, S3, p3(df, [q1; q2]))

    N = Int(n + n_burn)
    samples = repeat([z1.q; z2.q; z3.q...]', N + 1)
    accepts = fill(false, N + 1)

    p = Progress(N)
    generate_showvalues(x) = () -> [("Gibbs", x)]

    @showprogress for i ∈ 1:N
        z1, z2, z3, α_MH = gibbs_step3(df, S, (z1, z2, z3))
        samples[i+1, :] = [z1.q; z2.q; z3.q...]
        accepts[i+1] = α_MH

        next!(p; showvalues=generate_showvalues(mean(accepts[1:i+1])))
    end

    samples = samples[Int(n_burn):end, :]
    accepts = accepts[Int(n_burn):end]
    @info "Acceptance Ratio" ratio = round(mean(accepts); digits=4)
    return (samples, accepts)
end

# S1 = main.HaRAM(ϵ=0.01, L=4, γ=250.0)
S1 = main.HaRAM(ϵ=0.01, L=5, γ=100.0)
S2 = main.HaRAM(ϵ=2e-2, L=5, γ=0.5)
S2 = main.HMC(ϵ=2e-2, L=5)
S = (S1, S2, S3)

s, a = gibbs_sample3(
    df,
    (
        main.HaRAM(ϵ=0.01, L=10, γ=50.1),
        main.HaRAM(ϵ=2e-2, L=5, γ=0.2),
        main.HaRAM(ϵ=2e-5, L=3, γ=0.1)
    ),
    n=2000
)
θ = fmap(eachcol(s[a, :]), Binv)
# θ = fmap(eachcol(s[:, :]), Binv)
plot(histogram(θ[1], bins=100), histogram(θ[2], bins=100), layout=(2, 1), size=(1000, 1000))

plt = plot(df.time, df.lca, marker=:o, label="x(t)", lw=2)
for i in sample(eachindex(θ[1]), 200)
    plt = plot(plt, df.time .- θ[1][i], df.lcb, label="", lw=2, la=0.05, c=:grey)
end
plt = plot(plt, df.time .- median(θ[1]), df.lcb, marker=:o, label="y(t-Δ)", c=:firebrick1, lw=2)
vline!([minimum(df.time), 4825, 4940, 5185, 5360, 5550, 5590], label="", legend=:bottomright, ls=:dash)
plt