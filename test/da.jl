using DrWatson
@quickactivate

using main
using Distributions, BenchmarkTools, Plots, Pipe

function plt(d=2; lim=(-5, 5), l=200)
    sq = range(lim..., length=l)
    if d == 2
        contourf(sq, sq, (x, y) -> model.f([x; y; fill(0, model.d - 2)]))
    elseif d >= 3
        contourf(sq, sq, (x, y) -> model.f([x; y; fill(0, model.d - 2)]))
    end
end

#################################################

begin
    d = 50
    model = Model(
        ξ=MixtureModel([MvNormal(x, 1.0 / d^0.2) for x in ([-10, +10] .* fill(ones(d)) ./ √d)])
    )
end

#################################################


x, y = mcmc(
    DualAverage(λ=10, δ=0.5),
    HMC(),
    model; n=5000, n_burn=1000
)

newx, newy = mcmc(
    DualAverage(λ=10, δ=0.5),
    HaRAM(),
    model; n=5000, n_burn=1000
)


@pipe (x[y, :], newx[newy, :]) .|>
      sample(eachcol(_) |> collect, 2, replace=false) .|>
      hcat(_...) .|>
      Tuple.(eachrow(_)) .|>
      scatter(plt(lim=(-3, 3)), _, label="", ma=0.1, msw=0) |>
      plot(_..., layout=(1, 2), size=(800, 300))


@pipe (x[y, :], newx[newy, :]) .|>
      sample(eachcol(_) |> collect, 1) .|>
      hcat(_...) .|>
      plot(_, label="", ylim=(-3, 3)) |>
      plot(_..., layout=(1, 2), size=(800, 300))
