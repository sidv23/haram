using DrWatson
@quickactivate "haram"

begin
    using main
    using BlockDiagonals, LinearAlgebra, Plots
    using Distributions, MCMCChains, Random
end


plot_params = Plot_params(plotband=7, quiveralpha=0, colorscheme=palette(:linear_wcmr_100_45_c42_n256, 100, rev=true))

begin
    μ = [[2.18, 5.76], [8.67, 9.59], [4.24, 8.48], [8.41, 1.68], [3.93, 8.82], [3.25, 3.47], [1.70, 0.50], [4.59, 5.60], [6.91, 5.81], [6.87, 5.40], [5.41, 2.65], [2.70, 7.88], [4.98, 3.70], [1.14, 2.39], [8.33, 9.50], [4.93, 1.5], [1.83, 0.09], [2.26, 0.31], [5.54, 6.86], [1.69, 8.11]]

    Σ = [diagm(fill(0.02, 2)) for _ in eachindex(μ)]

    ξ = MixtureModel(
        [MvNormal(x .- 5.0, y) for (x, y) in zip(μ, Σ)],
        fill(1 / length(μ), length(μ))
    )
    model = Model(ξ=ξ)
    plt() = plot_pdf(model, plot_params)
end



begin
    s1, a1 = mcmc(
        HaRAM(ϵ=0.15, L=20, γ=0.8),
        model; n=1e4, n_burn=1e3
    )
    x_haram = s1[a1, :]
    scatter(plt(), x_haram |> m2t, ma=0.2)
end

begin
    s2, a2 = mcmc(
        HMC(ϵ=0.25, L=10),
        model; n=1e4, n_burn=1e3
    )
    x_hmc = s2[a2, :]
    scatter(plt(), x_hmc |> m2t, ma=0.2)
end

begin
    s3, a3 = mcmc(
        RAM(kernel=MvNormal(2, 2.0), z=randn(model.d)),
        model; n=1e4, n_burn=1e3
    )
    x_ram = s3[a3, :]
    scatter(plt(), x_ram |> m2t, ma=0.2)
end

begin
    s4, w4, a4 = mcmc(
        PEHMC(ϵ=0.5, L=10, N=2),
        model; n=1e3, n_burn=1e2
    )
    x_pehmc, w_pehmc = [s4[a4, :]...] |> a2m, [w4[a4, :]...]
    scatter(plt(), x_pehmc |> m2t, ms=w_pehmc * 5)
end

begin
    s5, a5 = mcmc(
        PT(τ=[1.0, 0.5, 0.1]),
        HMC(ϵ=0.5, L=20),
        model; n=1e4, n_burn=1e3
    )
    x_pt_hmc = s5[a5, :]
end

begin
    s5, a5 = mcmc(
        PT(τ=[1.0, 0.5, 0.1]),
        HaRAM(ϵ=0.5, L=20, γ=0.1),
        model; n=1e4, n_burn=1e3
    )
    x_pt_haram = s6[a6, :]
end

contourf(-1:0.1:1, -2:0.1:2, (x, y) -> f([x, y]))