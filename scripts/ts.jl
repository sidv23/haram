using DrWatson;
@quickactivate "haram";
using main
using DynamicPPL, Bijectors
using Zygote
using LinearAlgebra, Plots, Pipe

function make_model(ξ::DynamicPPL.Model, d)
    U(x) = min(-logjoint(ξ, (;
            σ2=invlink(Dists[1], x[1]),
            b0=invlink(Dists[2], x[2]),
            b1=invlink(Dists[3], x[3])
        )), floatmax(Float64))
    dU(x) = Zygote.gradient(x_ -> Zygote.forwarddiff(U, x_), x)[1]
    f(x) = max(exp(-U(x)), floatmin(Float64))
    g(x) = Zygote.gradient(x_ -> Zygote.forwarddiff(f, x_), x)
    return main.Model(ξ=ξ, d=d, f=f, g=g, U=U, dU=dU)
end

Dists = [Uniform(0, 1), Uniform(0, 10), Uniform(0, 10)]

b = stack([
    inverse(bijector(Dists[1])),
    inverse(bijector(Dists[2])),
    inverse(bijector(Dists[3])),
]...)

@model function lr(x, y)
    σ2 ~ transformed(Dists[1])
    b0 ~ transformed(Dists[2])
    b1 ~ transformed(Dists[3])
    y ~ MvNormal(b0 .+ b1 .* x, σ2 * I)
end

x = rand(200)
y = 7 .+ 5.5 .* x .+ rand(Normal(0.0, 0.5), 200)
LR = lr(x, y)
model = make_model(LR, 3)
s, a = mcmc(DualAverage(λ=0.3, δ=0.55), main.HaRAM(), model; n=2e3, n_burn=1e3)


cc = s[a, :] |> eachrow .|> b |> a2m
plot([histogram(cc[:, i]) for i in 1:3]..., layout=(3, 1))


foo(x) = map(z -> z[2] + z[3] * x, eachrow(cc))
predict(x) = map(z -> z[2] + z[3] * x, eachrow(cc)) |> mean
lo(x) = @pipe map(z -> z[2] + z[3] * x, eachrow(cc)) |> quantile(_, 1e-100)
hi(x) = @pipe map(z -> z[2] + z[3] * x, eachrow(cc)) |> quantile(_, 1-1e-100)

begin
    xs = sort(x)
    scatter(x, y, ma=0.5, msw=0, label="")
    plot!(xs, predict.(xs), ribbon=(predict.(xs) .- lo.(xs), hi.(xs) - predict.(xs)), fa=0.5, lw=3, label="")
end