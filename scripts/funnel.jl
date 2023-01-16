using DrWatson
@quickactivate "haram"

begin
    using main
    using Plots, Random, Distributions
    using Flux, Turing
    gr()
    theme(:default)
end
begin
    function plt(exts=[(-20, 20), (-10, 10)], len=500, cls=palette(:linear_wcmr_100_45_c42_n256, 100))
        return contourf(
            range(exts[1]..., length=len),
            range(exts[2]..., length=len),
            (x, y) -> model(2).f([x; y])^0.9, c=cls
        )
    end

    @model function funnel(d, a, b, c)
        y ~ Normal(0, a)
        x ~ MixtureModel([MvNormal(x_ * ones(d - 1), exp(sign(x_) * b * (y - x_)) * x_ * 100) for x_ in [-c, +c]])
    end

    function model(d=2; a=sqrt(10.1), b=0.5, c=10.0)
        ℓ(x) = logjoint(funnel(d, a, b, c), (; x=x[1:end-1], y=x[end]))
        U(x) = isfinite(ℓ(x)) ? -ℓ(x) : 1e100
        dU(x) = ForwardDiff.gradient(U, x)
        f(x) = max(exp(-U(x)), 1e-100)
        g(x) = ForwardDiff.gradient(f, x)
        return main.Model(ξ=mod, d=d, f=f, g=g, U=U, dU=dU)
    end
    plt()
end

# @time sample(funnel1(d), Turing.HMC(0.5, 20), Int(1000));
# @time mcmc(main.HaRAM(ϵ=0.5, L=20), model; n=1000, n_burn=1);

s1, a1 = mcmc(
    DualAverage(λ=5, δ=0.55),
    HaRAM(ϵ=0.2, L=5, γ=0.2),
    model(); n=2e3, n_burn=1e3
)
x_haram = s1[a1, :];
scatter(plt(), Tuple.(eachrow(x_haram)))

s2, a2 = mcmc(
    # DualAverage(λ=200, δ=0.95),
    # HMC(ϵ=20, L=5),
    main.HMC(ϵ=0.25, L=5),
    model(); n=5e3, n_burn=5e3
)
x_hmc = s2[a2, :];
scatter(plt(), Tuple.(eachrow(x_hmc)), ms=1.5,  msw=0)






begin
    p1 = scatter(plt(extrema.(eachcol(x_haram)), 500), Tuple.(eachrow(x_haram)), ma=0.2, msw=0, label="", title="HaRAM", colorbar=false)
    xlim, ylim = extrema.(eachcol(x_haram))
    p2 = scatter(plt(extrema.(eachcol(x_haram)), 500), Tuple.(eachrow(x_hmc)), ma=0.2, msw=0, label="", title="HMC", colorbar=false, xlim=xlim, ylim=ylim)
    plot(p1, p2, size=(800, 300))
end