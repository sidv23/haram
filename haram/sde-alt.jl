### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 72323508-742e-11ed-14e5-9f08721c77e0
using Pkg; Pkg.activate("/storage/work/s/suv87/julia/haram")

# ╔═╡ 254a3c4e-3006-42c9-8d4c-29bb7f008b40
begin
	using DrWatson
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
	using ProgressLogging
	using Pipe
	
	ProgressMeter.ijulia_behavior(:clear);
end

# ╔═╡ 221d1eab-6657-4d21-bd0e-6d7f3a92998a
begin
Dists = (;
    Δ=Uniform(-1178.939, 1179.939),
    β=Uniform(-0.60, 0.60),
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
end

# ╔═╡ 0949a1e1-0775-43b7-b9d2-81feabae506f
begin
	df = CSV.File(datadir("q0957usno.csv")) |> DataFrame
	first(df, 5)
end

# ╔═╡ 9a9346b0-2d45-4073-99c9-4e1a86f9e190
Turing.@model function π_delay(df::DataFrame, ::Type{T}=Float64) where {T}
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

# ╔═╡ 0161373c-ec02-4d78-8206-ec0a2f524229
begin
	function p1(df::DataFrame, z::Vector{T}) where {T<:Float64}
	    K = keys(Dists)
	    K1 = K[1:end.∈Ref(1:2)]
	    K2 = K[1:end.∉Ref(1:2)]
	    ℓ(x) = logjoint(π_delay(df) | v2n([z...], K2), v2n([x...], K1))
	    U(x) = isfinite(ℓ(x)) ? -ℓ(x) : 1e100
	    dU(x) = ForwardDiff.gradient(U, x)
	    f(x) = max(exp(-U(x)), 1e-100)
	    g(x) = ForwardDiff.gradient(f, x)
	    return main.Model(ξ=Dists, d=2, f=f, g=g, U=U, dU=dU)
	end
	
	function p2(df::DataFrame, z::Vector{T}) where {T<:Float64}
	    K = keys(Dists)
	    K1 = K[1:end.∈Ref(1:2)]
	    K2 = K[1:end.∉Ref(1:2)]
	    ℓ(x) = logjoint(π_delay(df) | v2n([z...], K1), v2n([x...], K2))
	    U(x) = isfinite(ℓ(x)) ? -ℓ(x) : 1e100
	    dU(x) = ForwardDiff.gradient(U, x)
	    f(x) = max(exp(-U(x)), 1e-100)
	    g(x) = ForwardDiff.gradient(f, x)
	    return main.Model(ξ=Dists, d=length(Dists) - 2, f=f, g=g, U=U, dU=dU)
	end
end

# ╔═╡ 4420364d-a996-44cb-8b5f-4069adff7247
begin
	function gibbs_step2(df, S, z)
	    S1, S2 = S
	    z1, z2 = z
	
	    m1 = p1(df, z2.q)
	    z1 = RefreshState(z1.q, z1, S1, m1)
	    n1, a1 = OneStep(z1, S1, m1)
	    b1 = rand() < a1
	    z1 = b1 ? n1 : z1
	
	    m2 = p2(df, z1.q)
	    z2 = RefreshState(z2.q, z2, S2, m2)
	    n2, a2 = OneStep(z2, S2, m2)
	    b2 = rand() < a2
	    z2 = b1 && b2 ? n2 : z2
	
	    return (z1, z2, [b1, b2])
	end
	
	
	
	function gibbs_sample2(df, S; n=1000, n_burn=1, init=[B[1](420.0); B[2](0.0)], m=ones(5))
	    S1, S2 = S
	    q1, q2 = init, randn(3)
	
	    z1 = InitializeState(q1, S1, p1(df, q2))
	    z1 = (; z1..., m=m[1:2])
	
	    z2 = InitializeState(q2, S2, p2(df, q1))
	    z2 = (; z2..., m=m[3:end])
	
	    N = Int(n + n_burn)
	    samples = repeat([z1.q...; z2.q...]', N + 1)
	    accepts = fill(false, N + 1)
	
	    p = Progress(N)
	    generate_showvalues(x) = () -> [("Gibbs", x)]

		B = zeros(N+1, 2)
	
	    @progress for i ∈ 1:N
	        Z1, Z2, x = gibbs_step2(df, S, (z1, z2))
	        samples[i+1, :] = [z1.q...; z2.q...]
			B[i+1, :] = x
	        
			if prod(x)
				accepts[i+1] = true
				z1, z2 = Z1, Z2
			end

			varinfo = (i=i+1, x=Binv[1](Z1.q[1]), s=std(Binv[1](samples[accepts, 1])), B=[mean(B[1:i+1, :], dims=1)...], A=mean(accepts[1:i+1]))
			@info "" varinfo
	
	        # next!(p; showvalues=generate_showvalues(mean(accepts[1:i+1])))
	    end
	
	    samples = samples[Int(n_burn):end, :]
	    accepts = accepts[Int(n_burn):end]
	    @info "Acceptance Ratio" ratio = round(mean(accepts); digits=4)
	    return (samples, accepts)
	end
end

# ╔═╡ 889d4e9d-4a5b-4d06-b847-cc231b285377
begin
	s, a = gibbs_sample2(
	    df,
	    (
	        main.HaRAM(ϵ=1e-3, L=5, γ=2e1),
	        main.HaRAM(ϵ=3e-3, L=7, γ=1e-2),
	    ),
	    n=400,
		m = [2e5; 1e2; 1e-1; 1e+1; 1e+1]
	)
end;

# ╔═╡ 7a37a1eb-1311-4fe3-9285-2f1e9ecbac7e
begin
	θ = fmap(eachcol(s[a, :]), Binv)
	hist = plot(histogram(θ[1], bins=100, label="Δ", normalize=true), histogram(θ[2], bins=100, label="β", normalize=true), layout=(1, 2), size=(800, 250))
	hist
end

# ╔═╡ 1fa7aab6-e87b-4145-83ec-be889c93a685
map(i -> (@pipe Binv[i].(θ[i]) |> (mean(_), median(_))), 1:2)

# ╔═╡ 57be7c01-b5dc-447a-a025-e0065f78e30e
begin
	plt = plot(df.time, df.lcb, ribbon = 3 .* df.seb, label="", lw=2, la=0, ma=0)
	for i in sample(eachindex(θ[1]), 100)
	    plt = plot(plt, df.time .+ θ[1][i], df.lca .- θ[2][i], label="", lw=3, la=0.05, c=:grey)
	end
	
	plt = plot(plt, 
		df.time .+ mean(θ[1]), 
		df.lca .- mean(θ[2]),
		marker=:o, label="x(t-Δ) + β", 
		lw=2, c=:firebrick1
	)
	plt = plot(plt, 
		df.time, 
		df.lcb, 
		marker=:o, label="y(t)", 
		c=:dodgerblue, lw=2
	)
	vline!([minimum(df.time), 4825, 4940, 5185, 5360], label="", legend=:bottomright, ls=:dashdot, c=:black)
	plt
end

# ╔═╡ ba217769-1c79-49b0-9857-80b5f50c02ee
scatter(θ[1], log.(θ[end-1]))

# ╔═╡ 373904bd-62e1-4bcf-8612-d2bccc33f656
begin
	# c = sample(
	# 	π_delay(df),
	# 	Gibbs(
	# 		NUTS(1000, 0.2, :Δ, :β),
	# 		MH(:μ, :τ, :σ2),
	# 	),
	# 	1000
	# )
	# es = hcat([[c[k]...] for k in K]...)
	# θ = fmap(eachcol(es), Binv)
end

# ╔═╡ Cell order:
# ╠═72323508-742e-11ed-14e5-9f08721c77e0
# ╠═254a3c4e-3006-42c9-8d4c-29bb7f008b40
# ╠═221d1eab-6657-4d21-bd0e-6d7f3a92998a
# ╠═0949a1e1-0775-43b7-b9d2-81feabae506f
# ╠═9a9346b0-2d45-4073-99c9-4e1a86f9e190
# ╠═0161373c-ec02-4d78-8206-ec0a2f524229
# ╠═4420364d-a996-44cb-8b5f-4069adff7247
# ╠═889d4e9d-4a5b-4d06-b847-cc231b285377
# ╠═7a37a1eb-1311-4fe3-9285-2f1e9ecbac7e
# ╠═1fa7aab6-e87b-4145-83ec-be889c93a685
# ╠═57be7c01-b5dc-447a-a025-e0065f78e30e
# ╠═ba217769-1c79-49b0-9857-80b5f50c02ee
# ╠═373904bd-62e1-4bcf-8612-d2bccc33f656
