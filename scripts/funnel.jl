using DrWatson
@quickactivate "haram"

begin
    using main
    using Plots, Random, Distributions
    using Flux, Turing
    pyplot()
    theme(:default)
end

function plt(exts=[(-20, 20), (-10, 10)], len=500, cls=palette(:viridis, 100))
    return contourf(
        range(exts[1]..., length=len),
        range(exts[2]..., length=len),
        (x, y) -> model.f([x, y]), c=cls)
end

make_plots(x) = plot(
    [scatter(Tuple.(eachrow(x[:, [i, end]])), label="", ma=0.2, msw=0) for i in 1:1:size(x, 2)-1]...,
    # layout=(1, model.d - 1),
    # size=(150 * model.d, 200)
)

@model function funnel(d, a, b, c)
    y ~ Normal(0.1, a)
    x ~ MixtureModel([MvNormal(x_ * ones(d - 1), exp(b * y)) for x_ in [-c, +c]])
end;

begin
    a, b, c, d = sqrt(3.0), 0.5, 5.0, 20
    U(x) = -logjoint(funnel(d, a, b, c), (; x=x[1:end-1], y=x[end]))
    f(x) = exp(-U(x))
    model = Model(ξ=mod, d=d, f=f, U=U)
end

s1, a1 = mcmc(HaRAM(ϵ=0.1, L=10, γ=0.2),model; n=2e3, n_burn=1e3)
x_haram = s1[a1, :];
make_plots(x_haram[:, [1,2,3,4,end]])

 s2, a2 = mcmc(main.HMC(ϵ=0.2, L=25), model; n=2e3, n_burn=1e3)
x_hmc = s2[a2, :];
make_plots(x_hmc[:,[1,2,3,4,end]])

begin
    p1 = scatter(plt(extrema.(eachcol(x_haram)), 500), Tuple.(eachrow(x_haram)), ma=0.2, msw=0, label="", title="HaRAM", colorbar=false)
    xlim, ylim = extrema.(eachcol(x_haram))
    p2 = scatter(plt(extrema.(eachcol(x_haram)), 500), Tuple.(eachrow(x_hmc)), ma=0.2, msw=0, label="", title="HMC", colorbar=false, xlim=xlim, ylim=ylim)
    plot(p1, p2, size=(800, 300))
end