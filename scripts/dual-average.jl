using DrWatson
@quickactivate

using main
using Distributions, BenchmarkTools, Plots, Pipe

function plt(model; lim=(-5, 5), l=200, kwargs...)
    sq = range(lim..., length=l)
    contourf(sq, sq, (x, y) -> model.f([x; y; fill(0, model.d - 2)]); kwargs...)
end

#################################################

begin
    d = 100
    model = Model(
        ξ=MixtureModel([MvNormal(x, 1.0 / d^0.2) for x in ([-10, +10] .* fill(ones(d)) ./ √d)])
    )
end

#################################################


x, y = mcmc(
    DualAverage(λ=15, δ=0.5),
    HMC(),
    model; n=5000, n_burn=1000
)

newx, newy = mcmc(
    DualAverage(λ=15, δ=0.5),
    HaRAM(),
    model; n=5000, n_burn=1000
)

cls = palette(:linear_wcmr_100_45_c42_n256, rev=false);
@pipe (x[y, :], newx[newy, :]) .|>
      sample(eachcol(_) |> collect, 2, replace=false) .|>
      hcat(_...) .|>
      Tuple.(eachrow(_)) .|>
      scatter(
        plt(
            model, lim=(-3, 3), c=cls, levels=5, lw=0.01
        ), _, label="", ma=0.2, msw=0) |>
      plot(_..., layout=(1, 2), size=(800, 300))


@pipe (x[y, :], newx[newy, :]) .|>
      sample(eachcol(_) |> collect, 1) .|>
      hcat(_...) .|>
      plot(_, label="", ylim=(-3, 3)) |>
      plot(_..., layout=(1, 2), size=(800, 300))
