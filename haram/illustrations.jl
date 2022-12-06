using DrWatson, Revise;
@quickactivate

begin
    using main
    using DynamicPPL, Bijectors, Turing
    using LinearAlgebra, Plots, Pipe
end


l = @layout [
    [b{0.1h} a{0.1w}]
    [grid(1, 1) c{0.1w}]
]

Dy = Normal(0, 1)
Dx = MixtureModel([Normal(-2, 1.0), Normal(2, 1.0)])


L(x) = pdf(Dx, x[1]) * pdf(Dy, x[2])

p1 = plot(Normal(), fill=0, label="", axis=false, ticks=false, permute=(:y, :x))
p2 = plot(Normal(), fill=0, label="", axis=false, ticks=false)
p3 = contourf()

plot(, layout=l)