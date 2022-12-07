using DrWatson, Revise;
@quickactivate

begin
    using main
    using DynamicPPL, Bijectors, Turing
    using LinearAlgebra, Plots, Pipe
    using LaTeXStrings
    using Plots.PlotMeasures
    using ForwardDiff
end

gr(levels=7, xguidefontsize=9, msw=0.2, fmt=:svg, display_type=:inline, label="")
cls = palette(:linear_wcmr_100_45_c42_n256, 100, rev=false)

params = (;
    foreground_color_legend=nothing,
    axis=false,
    grid=false,
    border=nothing,
    box=false,
    ticks=false,
    colorbar=false
)

Dy = Normal(0, 1)
Dx = MixtureModel([Normal(-3, 1.0), Normal(3, 1.0)])

L(x, y) = (-(logpdf(Dx, x) + logpdf(Dy, y)))^(-0.99999)

begin
    py = plot(
        x -> pdf(Dy, x),
        fill=0,
        label="",
        legend=:topright,
        permute=(:y, :x),
        fc=:firebrick1,
        lc=:firebrick1,
        fa=0.2, lw=3;
        params...,
        margin=0mm
    )
    hline!([0.0], c=:grey, ls=:dashdot)
    # annotate!(0, 0, L"K(\mathbf{p})")

    px = plot(
        x -> pdf(Dx, x),
        -8, 8,
        fill=0,
        label="",
        legend=:top,
        fc=:dodgerblue,
        lc=:dodgerblue,
        lw=3, fa=0.2;
        params...,
        margin=0mm
    )
    vline!([-3.0, 3.0], c=:grey, ls=:dashdot)
    # annotate!(0, -0.025, L"U(\mathbf{q})")

    Pxy() = begin
        contourf(
            -8:0.1:8, -8:0.1:8,
            (x, y) -> L(x, y);
            params...,
            c=cls,
            padding=0,
            margin=0mm
        )
        vline!([-3.0, 3.0], c=:grey, ls=:dashdot)
        hline!([0.0], c=:grey, ls=:dashdot)
    end

    p0 = plot(0, 0; params...)
    annotate!(0.5, 0.5, L"H(\mathbf{q}, \mathbf{p})")

    l = @layout [
        [b{0.15h} a{0.15w}]
        [d{0.8h} c{0.15w}]
    ]
    plt(p=Pxy(); size=(650, 650)) = plot(px, p0, p, py, layout=l, size=size)
    plt()
end

begin
    model() = main.Model(
        ξ=Dx,
        f=x -> pdf(Dx, x...),
        U=x -> -logpdf(Dx, x...),
        dU=x -> ForwardDiff.derivative(x_ -> -logpdf(Dx, x_), x)
    )

    state = (; q=-4.5, p=0.9, m=2.5)
    p1, _ = OnePath(state, HMC(ϵ=0.25, L=25), model())
    p2, _ = OnePath(state, HaRAM(ϵ=0.25, L=25, γ=0.25), model())
    plt1 = plot(Pxy(), map(x -> (x.q, x.p), p1), marker=:o, c=:dodgerblue, lw=3, ms=5, msw=1.0, lc=:grey, label="HMC")
    plt1 = plot(plt1, map(x -> (x.q, x.p), p2[1:end-24]), marker=:o, c=:orange, lw=3, ms=5, msw=1.0, lc=:grey, label="Repelling")
    plt1 = plot(plt1, map(x -> (x.q, x.p), p2[end-24:end]), marker=:o, c=:chartreuse, lw=3, ms=5, msw=1.0, lc=:grey, label="Attracting")
    plt1 = scatter(plt1, (state.q, state.p), c=:white, ms=9, msw=3.0, legend=:topleft, legendfontsize=9, thickness_scaling=1)
    annotate!(0, 8, L"U(\mathbf{q})")
    annotate!(8.3, 0, L"K(\mathbf{p})")
    # plt(plt1, size=(350, 350))
    plt(plt1, size=(850, 350))
end
# savefig(plt(plt1, size=(550, 550)), plotsdir("illustrations/phase.pdf"))
savefig(plt(plt1, size=(850, 350)), plotsdir("illustrations/phase-wide.pdf"))