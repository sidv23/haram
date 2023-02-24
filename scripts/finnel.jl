# -*- coding: utf-8 -*-
using DrWatson
@quickactivate "haram"

begin
    using main
    using Plots, Random, Distributions
    using Flux, Turing, ForwardDiff
    gr()
    theme(:default)
end

begin
    function plt(m=model(2), exts=[(-50, 50), (-30, 30)], len=500, cls=palette(:linear_wcmr_100_45_c42_n256, 100); levels=10)
        return contourf(
            range(exts[1]..., length=len),
            range(exts[2]..., length=len),
            # (x, y) -> m.f([x; y])^(1-1e-10), c=cls,
            (x, y) -> -(m.U([x; y])), c=cls,
            fa=0.25,
            levels=levels
        )
    end

    scatterplot(x, baseplt=nothing; kwargs...) = begin
        if !isnothing(baseplt)
            p = plot(baseplt)
        else
            p = plot(x -> 0, 0.01, 0.02, la=0, lw=0, c=:white)
        end
        p = plot(p, Tuple.(eachrow(x)); c=:black, lw=0.1, la=0.5, label="")
        p = scatter(p, Tuple.(eachrow(x)); c=:orange, kwargs...)
    end

    function make_colorbar(m, lims; len=500)
        z = map(
            (x, y) -> -(m.U([x; y])),
            range(lims[1]..., length=len),
            range(lims[2]..., length=len)
        ) |> extrema

        cls = palette(:linear_wcmr_100_45_c42_n256, 100)
        scatter(
            [0, 0], [0, 1],
            zcolor=[0, 3],
            clims=z,
            c=cls,
            xlims=(1, 1.1),
            xshowaxis=false,
            yshowaxis=false,
            label="",
            cbartitle="",
            grid=false
        )
    end

    @model function funnel(d, c=2.0)
        y ~ MixtureModel([Normal(-10, 7), Normal(10, 7)])
        x ~ MixtureModel([MvNormal(x_ * ones(d - 1), 1 * exp(sign(x_) * 0.5 * (y + x_))) for x_ in [-c, +c]])
    end

    function model(d=2; c=10.0)
        ℓ(x) = logjoint(funnel(d, c), (; x=x[1:end-1], y=x[end]))
        U(x) = isfinite(ℓ(x)) ? -ℓ(x) : 1e200
        dU(x) = ForwardDiff.gradient(U, x)
        f(x) = max(exp(-U(x)), 1e-200)
        g(x) = ForwardDiff.gradient(f, x)
        return main.Model(ξ=mod, d=d, f=f, g=g, U=U, dU=dU)
    end

    m = model(c=20.0)
    exts = [(-500, 500), (-25, 25)]
    plt(m, exts, levels=100)
end

begin
    @model function funnel(d, c=2.0)
        μ = 20
        σ = 5.0
        y ~ MixtureModel([Normal(-μ, σ), Normal(μ, σ)])
        x ~ MixtureModel([MvNormal(x_ * ones(d - 1) * 10, 1 * exp(sign(x_) * 0.5 * (y + x_))) for x_ in [-c, +c]])
    end

    m = model(c=25.0)
    exts = [(-500, 500), (-25, 25)]
    plt(m, exts, levels=100)
end

begin
    s1, a1 = mcmc(
        HaRAM(ϵ=0.1, L=200, γ=0.1),
        m; n=5e3, n_burn=1e3
    )
    x_haram = s1[a1, :]
    scatterplot(
        x_haram,
        plt(m, extrema.(eachcol(x_haram)), levels=5);
        msw=0.1, ma=0.5, label="RA-HMC"
    )
end


begin
    s2, a2 = mcmc(
        main.HMC(ϵ=0.6, L=200),
        m; n=5e3, n_burn=1e3
    )
    x_hmc = s2[a2, :]
    scatterplot(
        x_hmc,
        plt(m, extrema.(eachcol(x_haram)), levels=10);
        msw=0.1, ma=0.5, label="HMC"
    )
end


begin
    lims = extrema.(eachcol([x_hmc; x_haram]))
    baseplt = plt(m, lims, levels=100)
    pars = (; msw=0.1, ma=0.5, colorbar=false)

    p1 = scatterplot(x_hmc, baseplt; label="HMC", pars...)
    p2 = scatterplot(x_haram, baseplt; label="RA-HMC", pars...)

    h = make_colorbar(m, lims)

    @show h[1][:xaxis][:showaxis]
    @show h[1][:yaxis][:showaxis]

    l = @layout [grid(1, 2) a{0.01w}]
    plot(p1, p2, h, layout=l, size=(800, 300))
end


begin
    lims = extrema.(eachcol(x_hmc))

    zoom_plt1 = scatterplot(
        x_hmc,
        plt(m, lims, levels=10);
        msw=0.1, ma=0.25, label="HMC",
        xlim=lims[1], ylim=lims[2]
    )

    zoom_plt2 = scatterplot(
        x_haram,
        plt(m, lims, levels=10);
        msw=0.1, ma=0.25, label="RA-HMC",
        xlim=lims[1], ylim=lims[2]
    )

    h = make_colorbar(m, lims)

    l = @layout [grid(1, 2) a{0.01w}]
    plot(zoom_plt1, zoom_plt2, h, layout=l, size=(800, 300))
end