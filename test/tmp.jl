using Distributions, Random, StatsPlots
using Turing
using Random
Random.seed!(12)

p_true = 0.5;
N = 100
data = rand(Bernoulli(p_true), N);
prior_belief = Beta(1, 1);

@model function coinflip(; N::Int)
    p ~ Beta(1, 1) # Prior
    y ~ filldist(Bernoulli(p), N)
    return y
end

coinflip(y::AbstractVector{<:Real}) = coinflip(; N=length(y)) | (; y)
model = coinflip(data)

foo(x) = loglikelihood(model, (; p=x))

loglikelihood(model, (; p=0.5)) * model