using DrWatson
@quickactivate

begin
    using main
    using BlockDiagonals, LinearAlgebra, Plots, ProgressMeter
    using Distributions, MCMCChains, Random

    gr(fmt=:png, levels=5, lw=0.5, msa=0.1, msw=0.5, ma=0.2, msc=:firebrick1, legend=:topright)
    ProgressMeter.ijulia_behavior(:clear);
end

begin
    function generate_model(; m = 5.0 ,s = 0.5, d = 10)
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
    
    function scatterplots(xs, names; baseplt=plot(0,0,label=nothing), l=400, kwargs...)
        L = length(xs)
        ds = sample(1:size(xs[1], 2), 2, replace=false)
        plts = [scatter(baseplt, x[:, ds[1]], x[:, ds[2]], label=names[i]) for (x, i) in zip(xs, eachindex(xs))]
        plot(plts..., axes=false, ticks=false; kwargs...)
    end
    
    function traceplots(xs, names, args...; baseplt=plot(0,0,label=nothing), l=400, kwargs...)
        L = length(xs)
        ds = sample(1:size(xs[1], 2), 1, replace=false)
        plts = [plot(x[:, ds[1]], label=names[i]; args...) for (x, i) in zip(xs, eachindex(xs))]
        plot(plts...; kwargs...)
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
    
    function scatterplot(plt, x; lw=0.1, la=0.2, col=:black, kwargs...)
        scatter(plt, x; kwargs...)
        plot!(x, lw=lw, la=la, c=col, label="")
    end
end



begin
    # Example: $\mathbb R^{2}$
    model = generate_model(d=2);
    cls = palette(:linear_wcmr_100_45_c42_n256, 100, rev=false)
    plt(; lim=(-10, 10)) = contourf(
        repeat([range(lim..., length=100)], 2)...,
        (x, y) -> model.f([x; y])^(1-1e-2),
        c=cls
    )
    plt2 = plt()
end



begin
    Random.seed!(2023)
    @time s1, a1 = mcmc(
        DualAverage(λ=5, δ=0.8),
        HMC(),
        model; n=5e3, n_burn=1e3
    )
end

x_hmc = s1[a1, :]
chain_hmc = Chains(x_hmc)
plt2_hmc = scatterplot(plt(), x_hmc[:, 1:2] |> m2t, label="HMC")

Random.seed!(2023)
@time s2, a2 = mcmc(
    DualAverage(λ=10, δ=0.6),
    HaRAM(),
    model; n=5e3, n_burn=1e3
)
x_haram = s2[a2, :]
chain_haram = Chains(x_haram)
plt2_haram = scatterplot(plt(), x_haram[:, 1:2] |> m2t, label="HMC")
Random.seed!(2023)
@time s3, a3 = mcmc(
    RAM(kernel=MvNormal(model.d, 5), z=randn(model.d)),
    model; n=5e3, n_burn=1e3
)
x_ram = s3[a3, :]
chain_ram = Chains(x_ram)
plt2_ram = scatterplot(plt(), x_ram[:, 1:2] |> m2t, label="RAM")
Random.seed!(2023)
@time s4, w4, a4 = mcmc(
    PEHMC(ϵ=0.8, L=5, N=2),
    model; n=1e3, n_burn=1e2
)
x_pehmc, w_pehmc = [s4[a4, :]...] |> a2m, [w4[a4, :]...]
chain_pehmc = Chains(x_pehmc)
plt2_pehmc = scatterplot(plt(), x_pehmc |> m2t, ms=w_pehmc * 5, label="PE-HMC")
scatterplot(plt(), x_haram[:, 1:2] |> m2t, la=0.1, lw=1, label="RAM")
savefig("./plots/test.pdf")
Random.seed!(2023)
@time s5, a5 = mcmc(
    PT(τ=[1.0, 0.5, 0.1]),
    HMC(ϵ=0.95, L=2),
    model; n=1e4, n_burn=1e3
)
x_pt_hmc = s5[a5, :]
chain_pt_hmc = Chains(x_pt_hmc)
plt2_pt_hmc = scatterplot(plt(), x_pt_hmc[:, 1:2] |> m2t, label="HMC (Parallel Tempered)")
Random.seed!(2023)
@time s6, a6 = mcmc(
    PT(τ=[1.0, 0.5, 0.1]),
    HMC(ϵ=0.5, L=3),
    model; n=1e4, n_burn=1e3
)
x_pt_haram = s6[a6, :]
chain_pt_haram = Chains(x_pt_haram)
plt2_pt_haram = scatterplot(plt(), x_pt_haram[:, 1:2] |> m2t, label="HaRAM (Parallel Tempered)")
names = ("HMC", "HaRAM", "RAM", "PE-HMC")
xs2d = [x_hmc, x_haram, x_ram, x_pehmc]
chains2d = [chain_hmc, chain_haram, chain_ram, chain_pehmc];

plt2_tr = traceplots(xs2d, names, lw=1, layout=(2,2), ylim=(-10,10), l=100, size=(900, 500), title=["" "d=2" "" ""])
plt2_acf = acfplots(chains2d, names, lw=3)
names = ("HMC", "HaRAM", "RAM", "PE-HMC", "PT-HMC", "PT-HaRAM")
xs2d = [x_hmc, x_haram, x_ram, x_pehmc, x_pt_hmc, x_pt_haram]
chains2d = [chain_hmc, chain_haram, chain_ram, chain_pehmc, chain_pt_hmc, chain_pt_haram];


plt2_acf_all = acfplots(chains2d, names, lw=3, title="d=2")
plt2_tr_all = traceplots(xs2d, names, lw=1, layout=(2,3), ylim=(-10,10), l=100, size=(900, 500), title=["" "d=2" "" "" "" "" ""])

# @pipe ["plt2_tr_all", "plt2_acf_all"] .|>
# begin 
#     plot(eval(Meta.parse(_))); savefig.(plotsdir.("gaussian3/d2/" .* _ .* [".pdf", ".svg"] )) 
# end;
# @pipe ["plt2_tr", "plt2_acf", "plt2_hmc", "plt2_haram", "plt2_ram", "plt2_pehmc", "plt2_pt_hmc", "plt2_pt_haram"] .|> 
# begin 
#     plot(eval(Meta.parse(_))); savefig.(plotsdir.("gaussian3/d2/" .* _ .* [".pdf", ".svg"] )) 
# end;
w2_minibatch(xs2d, model)
# Example: $\mathbb R^{10}$
model = generate_model(d=10);
Random.seed!(2023)
@time s1, a1 = mcmc(
    DualAverage(λ=15, δ=0.6),
    HMC(),
    model; n=4e3, n_burn=1e3
)
x_10_hmc = s1[a1, :]
chain_10_hmc = Chains(x_10_hmc)
plt10_hmc = scatterplot(x_10_hmc[:, 1:2] |> m2t, lim=(-10, 10), label="HMC")
Random.seed!(2023)
@time s2, a2 = mcmc(
    DualAverage(λ=20, δ=0.6),
    HaRAM(),
    model; n=4e3, n_burn=1e2
)
x_10_haram = s2[a2, :]
chain_10_haram = Chains(x_10_haram)
plt10_haram = scatterplot(x_10_haram[:, 1:2] |> m2t, lim=(-10, 10), label="HaRAM")
Random.seed!(2023)
@time s3, a3 = mcmc(
    RAM(kernel=MvNormal(model.d, 0.5), z=randn(model.d)),
    model; n=4e3, n_burn=1e3
)
x_10_ram = s3[a3, :]
chain_10_ram = Chains(x_10_ram)
plt10_ram = scatterplot(x_10_ram[:, 1:2] |> m2t, lim=(-10, 10), label="RAM")
Random.seed!(2023)
@time s4, w4, a4 = mcmc(
    PEHMC(ϵ=0.25, L=20, N=2),
    model; n=1e3, n_burn=1e2
)
x_10_pehmc, w_10_pehmc = [s4[a4, :]...] |> a2m, [w4[a4, :]...]
chain_10_pehmc = Chains(x_10_pehmc)
plt10_pehmc = scatterplot(x_10_pehmc[:, 1:2] |> m2t, ms=w_10_pehmc * 5, lim=(-10, 10), label="PEHMC")
Random.seed!(2023)
@time s5, a5 = mcmc(
    PT(τ=[1.0, 0.5, 0.1]),
    HMC(ϵ=0.65, L=20),
    model; n=4e3, n_burn=1e3
)
x_10_pt_hmc = s5[a5, :]
chain_10_pt_hmc = Chains(x_10_pt_hmc)
plt10_pt_hmc = scatterplot(x_10_pt_hmc[:, 1:2] |> m2t, lim=(-10, 10), label="PT-HMC")
Random.seed!(2023)
@time s6, a6 = mcmc(
    PT(τ=[1.0, 0.5, 0.1]),
    HaRAM(ϵ=0.55, L=25, γ=0.25),
    model; n=4e3, n_burn=1e3
)
x_10_pt_haram = s6[a6, :]
chain_10_pt_haram = Chains(x_10_pt_haram)
plt10_pt_haram = scatterplot(x_10_pt_haram[:, 1:2] |> m2t, ma=0.2, lim=(-10, 10), label="PT-HaRAM")
names = ("HMC", "HaRAM", "RAM", "PE-HMC")
xs10d = [x_10_hmc, x_10_haram, x_10_ram, x_10_pehmc]
chains10d = [chain_10_hmc, chain_10_haram, chain_10_ram, chain_10_pehmc];
plt10_acf = acfplots(chains10d, names, lw=3)
plt10_tr = traceplots(xs10d, names, lw=1, layout=(2,2), ylim=(-10,10), l=100, size=(900, 500), title=["d=$(model.d)" "" "" ""])
# @pipe ["plt10_tr", "plt10_acf", "plt10_hmc", "plt10_haram", "plt10_ram", "plt10_pehmc", "plt10_pt_hmc", "plt10_pt_haram"] .|> 
# begin
#     plot(eval(Meta.parse(_))); savefig.(plotsdir.("gaussian3/d10/" .* _ .* [".pdf", ".svg"] )) 
# end;
names = ("HMC", "HaRAM", "RAM", "PE-HMC")
names = ("HMC", "HaRAM", "RAM", "PE-HMC", "pt-hmc", "pt-haram")
xs10d = [x_10_hmc, x_10_haram, x_10_ram, x_10_pehmc, x_10_pt_hmc, x_10_pt_haram]
chains10d = [chain_10_hmc, chain_10_haram, chain_10_ram, chain_10_pehmc, chain_10_pt_hmc, chain_10_pt_haram];

plt10_acf_all = acfplots(chains2d, names, lw=3, title="d=10")
plt10_tr_all = traceplots(xs2d, names, lw=1, layout=(2,3), ylim=(-10,10), l=100, size=(900, 500), title=["" "d=10" "" "" "" "" ""])

# @pipe ["plt10_tr_all", "plt10_acf_all"] .|>
# begin 
#     plot(eval(Meta.parse(_))); savefig.(plotsdir.("gaussian3/d10/" .* _ .* [".pdf", ".svg"] )) 
# end;
w2_minibatch(xs10d, model)
d2 = [
    "plt2_tr", "plt2_acf", "plt2_tr_all", "plt2_acf_all", 
    "plt2_hmc", "plt2_haram", "plt2_ram", "plt2_pehmc", 
    "plt2_pt_hmc", "plt2_pt_haram"
] 
d10= [
    "plt10_tr", "plt10_acf", "plt10_tr_all", "plt10_acf_all",
    "plt10_hmc", "plt10_haram", "plt10_ram", "plt10_pehmc",
    "plt10_pt_hmc", "plt10_pt_haram"
]

using Pipe
@pipe [d2; d10].|> 
begin 
    plot(eval(Meta.parse(_))); savefig.(plotsdir.("gaussian3/d/" .* _ .* [".pdf", ".svg"] )) 
end;
plot(
    (@pipe ["plt2_hmc", "plt2_haram", "plt2_ram", "plt2_pehmc"] .|> plot(eval(Meta.parse(_)), colorbar=false))...,
    layout=(2,2), size=(700, 500)
)
savefig.(plotsdir.("gaussian3/d/" .* "scatter2" .* [".pdf", ".svg"] ))
plot(
    (@pipe ["plt10_hmc", "plt10_haram", "plt10_ram", "plt10_pehmc"] .|> plot(eval(Meta.parse(_)), colorbar=false))...,
    layout=(2,2), size=(700, 500)
)
savefig.(plotsdir.("gaussian3/d/" .* "scatter10" .* [".pdf", ".svg"] ))