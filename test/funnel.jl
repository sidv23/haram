using Revise, DrWatson
@quickactivate

begin
  using main
  using Distributions, BenchmarkTools, Plots, Pipe, ProgressMeter
  using Turing
end

s2v(s::Symbol, v::Any) = @eval (($s) = ($v))
s2v(s::Symbol, v::D) where {D<:Distribution} = @eval (($s) ~ ($v))
v2n(V::Vector{T}, S::NTuple) where {T<:Real} = NamedTuple{S}(V)


@model function funnel(d, c=2.0)
  a ~ Bernoulli(0.5)
  x ~ MvNormal(d, 1.0)
  if a > 0
    Y = 3 .* x[1] .- 2
    X = exp(Y/2) .* x[2:end] .- fill(c, d-1)'
  else
    Y = 3 .* x[1] .+ 2
    X = exp(Y/2) .* x[2:end] .+ fill(c, d)'
  end
  return [Y; X]
end

function model(d=2; c=10.0)
  ℓ(x) = logjoint(funnel(d, c), (; x=x[1:end-1], y=x[end]))
  U(x) = isfinite(ℓ(x)) ? -ℓ(x) : 1e100
  dU(x) = ForwardDiff.gradient(U, x)
  f(x) = max(exp(-U(x)), 1e-100)
  g(x) = ForwardDiff.gradient(f, x)
  return main.Model(ξ=mod, d=d, f=f, g=g, U=U, dU=dU)
end