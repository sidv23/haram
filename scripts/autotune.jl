using Revise, DrWatson
@quickactivate

begin
    using main
    using Distributions, BenchmarkTools, Plots, Pipe, ProgressMeter
end

begin
    gr(fmt=:png, levels=3, xguidefontsize=9, msw=0.2)
    theme(:default)
    Plots.gr_cbar_width[] = 0.02
    cls = palette(:linear_wcmr_100_45_c42_n256, 100, rev=false)
    ProgressMeter.ijulia_behavior(:clear)
end

begin
    function plt(d=2; lim=(-5, 5), l=200, bar=false)
        sq = range(lim..., length=l)
        if d == 2
            contourf(sq, sq, (x, y) -> model(d).f([x; y; fill(0, model(d).d - 2)]), c=cls, lw=0.1, colorbar=bar)
        elseif d >= 3
            contourf(sq, sq, (x, y) -> model(d).f([x; y; fill(0, model(d).d - 2)]), c=cls, lw=0.1, colorbar=bar)
        end
    end

    function plot_colorbar(lim)
        scatter([0, 0], [0, 1],
            zcolor=lim, clims=lim, xlims=(1, 1.01), c=cls,
            label="", colorbar_title="", framestyle=:none
        )
    end


    function make_plot(x1, x2, d, lim; kwargs...)
        k = sample(1:d, 2, replace=false)

        l = @layout [
            [grid(1, 2)
            a{0.25h}
            b{0.25h}] c{0.01w}
        ]

        p = plt(d, lim=lim)
        p0 = plot_colorbar((0, 1))
        p1 = scatter(p, x1[:, k[1]], x1[:, k[2]], label="HMC", ratio=1, grid=false; legend=:bottomright)
        p2 = scatter(p, x2[:, k[1]], x2[:, k[2]], label="HaRAM", ratio=1, grid=false; legend=:bottomright)
        p3 = plot(x1[:, 1], ylim=lim, label="HMC", legend=:topright)
        p4 = plot(x2[:, 1], ylim=lim, label="HaRAM", legend=:topright)

        plot(p1, p2, p3, p4, p0, layout=l, title=[i == 5 ? "(d=$d)" : "" for j in 1:1, i in 1:5]; kwargs...)
    end


    model(d) = main.Model(
        ξ=MixtureModel(
            [MvNormal(x, 1.0 / d^0.2) for x in ([-10, +10] .* fill(ones(d)) ./ √d)]
        )
    )
end

begin
    d, lim = 3, (-8, 8)

    @time s, a = mcmc(
        DualAverage(λ=20, δ=0.55),
        HMC(),
        model(d); n=5000, n_burn=1000
    )
    x_hmc = s[a, :]

    @time s, a = mcmc(
        DualAverage(λ=20, δ=0.55),
        HaRAM(),
        model(d); n=5000, n_burn=500
    )
    x_haram = s[a, :]

    x_hmc1, x_haram1 = x_hmc, x_haram
    make_plot(x_hmc, x_haram, d, lim)
end

begin
    d, lim = 10, (-6, 6)

    @time s, a = mcmc(
        DualAverage(λ=20, δ=0.55),
        HMC(),
        model(d); n=5000, n_burn=1000
    )
    x_hmc = s[a, :]

    @time s, a = mcmc(
        DualAverage(λ=20, δ=0.55),
        HaRAM(),
        model(d); n=5000, n_burn=1000
    )
    x_haram = s[a, :]

    x_hmc2, x_haram2 = x_hmc, x_haram

    make_plot(x_hmc, x_haram, d, lim)
end

begin
    d, lim = 50, (-3, 3)

    @time s, a = mcmc(
        DualAverage(λ=20, δ=0.55),
        HMC(),
        model(d); n=5000, n_burn=1000
    )
    x_hmc = s[a, :]

    @time s, a = mcmc(
        DualAverage(λ=20, δ=0.555),
        HaRAM(),
        model(d); n=5000, n_burn=1000
    )
    x_haram = s[a, :]

    x_hmc3, x_haram3 = x_hmc, x_haram
    make_plot(x_hmc, x_haram, d, lim)
end

begin
    d, lim = 100, (-2, 2)

    @time s, a = mcmc(
        DualAverage(λ=30, δ=0.55),
        HMC(),
        model(d); n=5000, n_burn=100
    )
    x_hmc = s[a, :]

    @time s, a = mcmc(
        DualAverage(λ=30, δ=0.55),
        HaRAM(),
        model(d); n=5000, n_burn=100
    )
    x_haram = s[a, :]

    x_hmc4, x_haram4 = x_hmc, x_haram
    make_plot(x_hmc, x_haram, d, lim)
end


gr(fmt=:png, levels=4, xguidefontsize=9, msc=:black, msw=0.1, lw=0.5, tickfontsize=7)

d, lim = 3, (-8, 8)
make_plot(x_hmc1, x_haram1, d, lim)
@pipe "d$d" |> begin
    savefig(plotsdir("autotune/" * _ * ".pdf"))
    savefig(plotsdir("autotune/" * _ * ".svg"))
end

d, lim = 10, (-6, 6)
make_plot(x_hmc2, x_haram2, d, lim)
@pipe "d$d" |> begin
    savefig(plotsdir("autotune/" * _ * ".pdf"))
    savefig(plotsdir("autotune/" * _ * ".svg"))
end

d, lim = 50, (-3, 3)
make_plot(x_hmc3, x_haram3, d, lim)
@pipe "d$d" |> begin
    savefig(plotsdir("autotune/" * _ * ".pdf"))
    savefig(plotsdir("autotune/" * _ * ".svg"))
end

d, lim = 100, (-2, 2)
make_plot(x_hmc4, x_haram4, d, lim)
@pipe "d100" |> begin
    savefig(plotsdir("autotune/" * _ * ".pdf"))
    savefig(plotsdir("autotune/" * _ * ".svg"))
end
