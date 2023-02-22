using DrWatson
@quickactivate

begin
    using main
    using Distributions, BenchmarkTools, Pipe
    using Plots
end


begin
    path2pts(path) = map(x -> tuple(x.q...), path)

    plotpath(pts, bgplt=plot(0, 0, label=nothing); kwargs...) = begin
        plt = plot(bgplt, pts; kwargs...)
        scatter!(plt, pts; kwargs...)
    end

    animatepath(pts, bgplt=plot(0, 0, label=nothing); kwargs...) = begin
        a = @animate for i in eachindex(pts)
            plt = plot(bgplt, pts[1:i]; kwargs...)
            scatter!(plt, pts[i:i]; kwargs...)
        end
        return a
    end
end


path, accept = OnePath(state, method, model);


# reversibility
begin
    model = Model(ξ=MvNormal([0.0, 0.0], [1.0 0.5; 0.5 1.0]))
    method = HaRAM(ϵ=0.4, L=10, γ=0.9)
    p = randn(2)
    state = (; q=rand(model.ξ), p=p, m=1.0)
    path1 = OnePath(state, method, model)[1]
    @pipe animatepath(path1 |> path2pts, contourplot(model), c=:dodgerblue) |> gif(_, fps=5)
end


# reversibility
begin
    model = Model(ξ=MvNormal([0.0, 0.0], [1.0 0.5; 0.5 1.0]))
    method = HaRAM(ϵ=0.8, L=10, γ=0.33)
    p = [1.0, 1.0]
    state = (; q=rand(model.ξ), p=p, m=1.0)
    path1 = OnePath(state, method, model)[1]
    newstate = path1[end]
    newstate = (; newstate..., p=-newstate.p)
    path2 = OnePath(state, method, model)[1]
    @pipe animatepath([path1; path2] |> path2pts, contourplot(model), c=:dodgerblue) |> gif(_, fps=10)
end