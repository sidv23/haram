### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 90e3e0e2-a880-11ed-3231-7d295e8a086b
begin
	using DrWatson
	import Pkg
	Pkg.activate("../../haram/")
end

# ╔═╡ c19a6c8e-4864-43f1-90f4-26d2c0f5eff3
begin
	using main
	using Plots, Random, Distributions, ProgressMeter
	using Flux, DynamicPPL, Zygote
	
	gr()
	theme(:default)
	cls = palette(:viridis, rev=false)
	default(fmt=:png, levels=7, lw=0.5, msw=0.5, la=0.5)
	ProgressMeter.ijulia_behavior(:clear);
end

# ╔═╡ 2021a423-8aba-4a14-ae96-fbbd4e53f6ff
begin
	Random.seed!(2022)
	s = 1.5
	μ = [[1.0 1.0], [-1.0 -1.0], [1.0 -1.0], [-1.0 1.0],] .* 2
	x0 = [s .* randn(30, 2) .+ μ[1]; s .* randn(100, 2) .+ μ[2]]
	x1 = [s .* randn(10, 2) .+ μ[3]; s .* randn(80, 2) .+ μ[4]]
	xs = [x0; x1]
	ys = [zeros(size(x0, 1)); ones(size(x1, 1))] .|> Int;
end

# ╔═╡ ff6cebe9-3fec-4278-8c26-56fa1e6d728c
begin
	function plot_colorbar(cls)
	    scatter([0, 0], [0, 1],
	        zcolor=[0, 1], clims=(0, 1), xlims=(1, 1.1), c=cls,
	        label="", colorbar_title="", framestyle=:none
	    )
	end
	
	function nn_plot(theta; res=25, c=:viridis)
	    x_range = collect(range(-6; stop=6, length=res))
	    contourf(
	        x_range, x_range,
	        (x, y) -> f(theta)([x, y])[1], c=c
	    )
	    scatter!(
	        Tuple.(eachrow(xs)),
	        c=map(x -> x == 1 ? :firebrick1 : :dodgerblue1, ys),
	        m=map(x -> x == 1 ? :square : :circle, ys),
	        group=ys,
	        legend=:bottomright
	    )
	end
	
	function nn_plot_mean(thetas; c=:viridis, res=25)
	    x_range = collect(range(-6; stop=6, length=res))
	    contourf(
	        x_range, x_range,
	        (x, y) -> mean([f(theta)([x, y])[1] for theta in eachrow(thetas)]), c=c
	    )
	    scatter!(
	        Tuple.(eachrow(xs)),
	        c=map(x -> x == 1 ? :firebrick1 : :dodgerblue1, ys),
	        m=map(x -> x == 1 ? :square : :circle, ys),
	        group=ys,
	        legend=:bottomright
	    )
	end
	
	function make_model(ξ::DynamicPPL.Model, d)
	    U(x) = min(-logjoint(ξ, (; θ=x)), floatmax(Float64))
	    dU(x) = Zygote.gradient(x_ -> Zygote.forwarddiff(U, x_), x)[1]
	    h(x) = max(exp(-U(x)), floatmin(Float64))
	    g(x) = Zygote.gradient(x_ -> Zygote.forwarddiff(f, x_), x)
	    return main.Model(ξ=ξ, d=d, f=h, g=g, U=U, dU=dU)
	end

	F = Chain(
	    Dense(2, 2, tanh),
	    Dense(2, 2, tanh),
	    Dense(2, 1, sigmoid)
	)
	
	θ, f = Flux.destructure(F)
	sigma = sqrt(1.0 / 0.09)
	
	@model function bayes_nn(xs, ys, nθ, f)
	    θ ~ MvNormal(zeros(nθ), sigma .* ones(nθ))
	    F = f(θ)
	    ps = F(xs)
	    for i in eachindex(ys)
	        ys[i] ~ Bernoulli(ps[i])
	    end
	end;
	
	ξ = bayes_nn(xs', ys, length(θ), f)
	model = make_model(ξ, length(θ));
end

# ╔═╡ ec47fd70-f7ca-423b-a2ae-21bd9e16a16a
scatter(Tuple.(eachrow(xs)), group=ys, ratio=1)

# ╔═╡ 8a30c96c-0c3b-477d-ae8c-34188c3c6706
begin
	Random.seed!(2022)
	s1, a1 = mcmc(
	    DualAverage(λ=1, δ=0.5),
	    HMC(),
	    model; n=2e3, n_burn=5e2
	)
	x_hmc = s1[a1, :];
end;

# ╔═╡ 482abb01-885c-4e37-acea-2e5c888993af
begin
	Random.seed!(2022)
	s2, a2 = mcmc(
	    DualAverage(λ=2, δ=0.5),
	    HaRAM(),
	    model; n=2e3, n_burn=5e2
	)
	x_haram = s2[a2, :];
end;

# ╔═╡ cd7b2895-cbe9-4bf6-b563-6bc3d0943139
begin
	st, res = 2, 50
	
	plot(
	    plot(
	        nn_plot_mean(x_haram[1:st:end, :], res=res, c=cls),
	        colorbar=false, title="HaRAM",
	    ),
	    plot(
	        nn_plot_mean(x_hmc[1:st:end, :], res=res, c=cls),
	        colorbar=false, title="HMC",
	    ),
	    plot_colorbar(cls),
	    layout=(@layout [grid(1, 2) a{0.035w}]),
	    link=:all,
	    size=(1000, 400)
	)
end

# ╔═╡ Cell order:
# ╠═90e3e0e2-a880-11ed-3231-7d295e8a086b
# ╠═c19a6c8e-4864-43f1-90f4-26d2c0f5eff3
# ╠═2021a423-8aba-4a14-ae96-fbbd4e53f6ff
# ╠═ff6cebe9-3fec-4278-8c26-56fa1e6d728c
# ╠═ec47fd70-f7ca-423b-a2ae-21bd9e16a16a
# ╠═8a30c96c-0c3b-477d-ae8c-34188c3c6706
# ╠═482abb01-885c-4e37-acea-2e5c888993af
# ╠═cd7b2895-cbe9-4bf6-b563-6bc3d0943139
