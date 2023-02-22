using DrWatson
@quickactivate "haram"

begin
    using main
    using Plots, Random, Distributions
    using Flux, DynamicPPL, Zygote
    gr()
    theme(:default)
    default(levels=7, lw=0.5, la=0.5, msw=0.5)
end


begin
    
    function make_model(ξ::DynamicPPL.Model, d)
        U(x) = min(-logjoint(ξ, (; θ=x)), floatmax(Float64))
        dU(x) = Zygote.gradient(x_ -> Zygote.forwarddiff(U, x_), x)[1]
        f(x) = max(exp(-U(x)), floatmin(Float64))
        g(x) = Zygote.gradient(x_ -> Zygote.forwarddiff(f, x_), x)
        return main.Model(ξ=ξ, d=d, f=f, g=g, U=U, dU=dU)
    end
    
    
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
        (x, y) -> f(theta)([x, y])[1], c=c
        )
        scatter!(
        Tuple.(eachrow(xs)),
        c=map(x -> x == 1 ? :firebrick1 : :dodgerblue1, ys),
        m=map(x -> x == 1 ? :square : :circle, ys),
        group=ys,
        legend=:bottomright
        )
    end
    
    function nn_plot_mean(thetas; c=:viridis, res=25)
        x_range = collect(range(-6; stop=6, length=res))
        contourf(
        x_range, x_range,
        (x, y) -> mean([f(theta)([x, y])[1] for theta in eachrow(thetas)]), c=c
        )
        scatter!(
        Tuple.(eachrow(xs)),
        c=map(x -> x == 1 ? :firebrick1 : :dodgerblue1, ys),
        m=map(x -> x == 1 ? :square : :circle, ys),
        group=ys,
        legend=:bottomright
        )
    end
    
    
    function bayes_plot(; method=1, step=10, cls=palette(:thermal, rev=false), res=50)
        if method == 1
            a = @animate for i in 0:step:size(s1, 1)-50
                plot(
                plot(
                nn_plot_mean(s1[i+1:i+5, :], res=res, c=cls),
                colorbar=false, title="HMC(i=$i)"
                ),
                plot(
                nn_plot_mean(s2[i+1:i+5, :], res=res, c=cls),
                colorbar=false, title="HaRAM(i=$i)"
                ),
                plot_colorbar(cls),
                layout=(@layout [grid(1, 2) a{0.035w}]),
                link=:all,
                size=(800, 300)
                )
            end
        elseif method == 2
            a = @animate for i in 0:step:size(s2, 1)
                plot(
                plot(
                nn_plot(s1[i+1, :], res=res, c=cls),
                colorbar=false, title="HMC(i=$i)"
                ),
                plot(
                nn_plot(s2[i+1, :], res=res, c=cls),
                colorbar=false, title="HaRAM(i=$i)"
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
            nn_plot_mean(x_hmc[1:step:end, :], res=res, c=cls),
            colorbar=false, title="HMC"
            ),
            plot(
            nn_plot_mean(x_haram[1:step:end, :], res=res, c=cls),
            colorbar=false, title="HaRAM"
            ),
            plot_colorbar(cls),
            layout=(@layout [grid(1, 2) a{0.035w}]),
            link=:all,
            size=(800, 300)
            )
        end
        return a
    end
end


###########################################

begin
    Random.seed!(2022)
    s = 1.25
    μ = [[1.0 1.0], [-1.0 -1.0], [1.0 -1.0], [-1.0 1.0],] .* 2
    x0 = [s .* randn(30, 2) .+ μ[1]; s .* randn(100, 2) .+ μ[2]]
    x1 = [s .* randn(10, 2) .+ μ[3]; s .* randn(80, 2) .+ μ[4]]
    xs = [x0; x1]
    ys = [zeros(size(x0, 1)); ones(size(x1, 1))] .|> Int
end

scatter(Tuple.(eachrow(xs)), group=ys, ratio=1)

F = Chain(
Dense(2, 2, selu),
Dense(2, 2, selu),
Dense(2, 1, sigmoid)
)

θ, f = Flux.destructure(F)

alpha = 0.09
sig = sqrt(1.0 / alpha)

@model function bayes_nn(xs, ys, nθ, f)
    θ ~ MvNormal(zeros(nθ), sig .* ones(nθ))
    F = f(θ)
    ps = F(xs)
    for i in eachindex(ys)
        ys[i] ~ Bernoulli(ps[i])
    end
end;


ξ = bayes_nn(xs', ys, length(θ), f)
model() = make_model(ξ, length(θ))


begin
    Random.seed!(2022)
    s1, a1 = mcmc(
    DualAverage(λ=1, δ=0.5),
    HMC(),
    model(); n=5e3, n_burn=1e3
    )
    x_hmc = s1[a1, :]
    nn_plot_mean(x_hmc, res=50)
end

begin
    Random.seed!(2022)
    s2, a2 = mcmc(
    DualAverage(λ=1, δ=0.5),
    HaRAM(),
    model(); n=5e3, n_burn=1e3
    )
    x_haram = s2[a2, :]
    nn_plot_mean(x_haram, res=50)
end

bayes_plot(method=3, cls=palette(:viridis, rev=false), res=100)# |> 
savefig(plotsdir("nn/contour.pdf"))

gif(bayes_plot(method=1, cls=palette(:viridis, rev=false), step=20), plotsdir("gifs/bayes_mean2.gif"), fps=7)
# gif(a, plotsdir("gifs/bayes_nomean.gif"), fps=10)

# savefig(plotsdir("nn/final.pdf"))
# savefig(plotsdir("nn/final.svg"))
