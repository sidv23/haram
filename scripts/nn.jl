using DrWatson
@quickactivate "haram"

begin
    using main
    using Plots, Random, Distributions
    using Flux, Turing
    pyplot()
    theme(:default)
end

begin
    Random.seed!(2022)
    s = 1.25
    μ = [[1.0 1.0], [-1.0 -1.0], [1.0 -1.0], [-1.0 1.0],] .* 2
    x0 = [s .* randn(30, 2) .+ μ[1]; s .* randn(100, 2) .+ μ[2]]
    x1 = [s .* randn(10, 2) .+ μ[3]; s .* randn(80, 2) .+ μ[4]]
    xs = [x0; x1]
    ts = [zeros(size(x0, 1)); ones(size(x1, 1))] .|> Int
end

scatter(Tuple.(eachrow(xs)), group=ts, ratio=1)
savefig(plotsdir("nn/pyplot.pdf"))

model = Chain(
    Dense(2, 3, tanh),
    Dense(3, 2, tanh),
    Dense(2, 1, sigmoid)
)

parameters_initial, reconstruct = Flux.destructure(model)
length(parameters_initial)

alpha = 0.09
sig = sqrt(1.0 / alpha)

npars = length(parameters_initial)
nn = reconstruct(randn(20))
preds = nn(xs')

@model function bayes_nn(xs, ts, nparameters, reconstruct)
    parameters ~ MvNormal(zeros(nparameters), sig .* ones(nparameters))
    nn = reconstruct(parameters)
    preds = nn(xs)
    for i in eachindex(ts)
        ts[i] ~ Bernoulli(preds[i])
    end
end;

U(x) = -1 * logjoint(bayes_nn(xs', ts, length(parameters_initial), reconstruct), (; parameters=x))
f(x) = exp(U(x))
model = Model(ξ=mod, d=20, f=f, U=U)

function plot_colorbar(cls)
    scatter([0, 0], [0, 1],
        zcolor=[0, 1], clims=(0, 1), xlims=(1, 1.1), c=cls,
        label="", colorbar_title="", framestyle=:none
    )
end

function nn_plot(theta; res=25, c=:viridis)
    x_range = collect(range(-6; stop=6, length=res))
    contourf(
        x_range, x_range,
        (x, y) -> reconstruct(theta)([x, y])[1], c=c
    )
    scatter!(
        Tuple.(eachrow(xs)),
        c=map(x -> x == 1 ? :firebrick1 : :dodgerblue1, ts),
        m=map(x -> x == 1 ? :square : :circle, ts),
        group=ts,
        legend=:bottomright
    )
end

function nn_plot_mean(thetas; c=:viridis, res=25)
    x_range = collect(range(-6; stop=6, length=res))
    contourf(
        x_range, x_range,
        (x, y) -> mean([reconstruct(theta)([x, y])[1] for theta in eachrow(thetas)]), c=c
    )
    scatter!(
        Tuple.(eachrow(xs)),
        c=map(x -> x == 1 ? :firebrick1 : :dodgerblue1, ts),
        m=map(x -> x == 1 ? :square : :circle, ts),
        group=ts,
        legend=:bottomright
    )
end


function bayes_plot(; method=1, step=10, cls=palette(:thermal, rev=false), res=50)
    if method == 1
        a = @animate for i in 0:100:size(s1, 1)-50
            plot(
                plot(
                    nn_plot_mean(s1[i+1:i+5, :], res=res, c=cls),
                    colorbar=false, title="HaRAM(i=$i)"
                ),
                plot(
                    nn_plot_mean(s2[i+1:i+5, :], res=res, c=cls),
                    colorbar=false, title="HMC(i=$i)"
                ),
                plot_colorbar(cls),
                layout=(@layout [grid(1, 2) a{0.035w}]),
                link=:all,
                size=(800, 300)
            )
        end
    elseif method == 2
        a = @animate for i in 0:100:size(s1, 1)
            plot(
                plot(
                    nn_plot(s1[i+1, :], res=res, c=cls),
                    colorbar=false, title="HaRAM(i=$i)"
                ),
                plot(
                    nn_plot(s2[i+1, :], res=res, c=cls),
                    colorbar=false, title="HMC(i=$i)"
                ),
                plot_colorbar(cls),
                layout=(@layout [grid(1, 2) a{0.035w}]),
                link=:all,
                size=(800, 300)
            )
        end
    elseif method == 3
        a = plot(
            plot(
                nn_plot_mean(x_haram[1:step:end, :], res=res, c=cls),
                colorbar=false, title="HaRAM"
            ),
            plot(
                nn_plot_mean(x_hmc[1:step:end, :], res=res, c=cls),
                colorbar=false, title="HMC"
            ),
            plot_colorbar(cls),
            layout=(@layout [grid(1, 2) a{0.035w}]),
            link=:all,
            size=(800, 300)
        )
    end
    return a
end


begin
    Random.seed!(2022)
    s1, a1 = mcmc(
        HaRAM(ϵ=0.05, L=4, γ=0.25),
        model; n=5e3, n_burn=1e3
    )
    x_haram = s1[a1, :];
    s2, a2 = mcmc(
        main.HMC(ϵ=0.05, L=4),
        model; n=5e3, n_burn=1e3
    )
    x_hmc = s2[a2, :];
    bayes_plot(method=3)
end

bayes_plot(method=3, cls=palette(:viridis, [0.1:0.1:0.8...], rev=false), step=20, res=100)# |> savefig(plotsdir("nn/contour.pdf"))

gif(bayes_plot(method=1, palette(:viridis, [0.1:0.1:0.8...], rev=false)), plotsdir("gifs/bayes_mean2.gif"), fps=10)
# gif(a, plotsdir("gifs/bayes_nomean.gif"), fps=10)

# savefig(plotsdir("nn/final.pdf"))
# savefig(plotsdir("nn/final.svg"))
