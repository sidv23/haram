using DrWatson
@quickactivate

using main
using Distributions, BenchmarkTools, MCMCChains, Plots


model = Model(ξ=MvNormal(2, 1.0))
method = HMC(ϵ=0.5, L=20)
method = HaRAM(ϵ=0.5, L=20, γ=0.05)
sample, accept = mcmc(method, model; n=1e4, n_burn=1e4);
x = sample[accept, :]

scatter(contourplot(model), x |> m2t, ma=0.05)


histogram(x[:, 2] |> m2t, bins=100, normalize=true)
plot!(x -> pdf(Normal(0, 1), x), lw=3)


model = Model(ξ=MvNormal(100, 1.0))
s1, a1 = mcmc(HaRAM(ϵ=0.5, L=20, γ=0.05), model; n=1e4, n_burn=1e4);
s2, a2 = mcmc(HMC(ϵ=0.5, L=20), model; n=1e4, n_burn=1e4);
s3, a3 = mcmc(RAM(kernel=MvNormal(model.d, 1.0), z=randn(model.d)), model; n=1e4, n_burn=1e4);

c1 = Chains(s1[a1, :])
c2 = Chains(s2[a2, :])
c3 = Chains(s3[a3, :])

ess_rhat(c1)[:, 2] |> mean
ess_rhat(c2)[:, 2] |> mean
ess_rhat(c3)[:, 2] |> mean

plot( mean(abs.(autocor(c1, lags=0:10)[:, :]), dims=1)', lw=3, label="HaRAM")
plot!(mean(abs.(autocor(c2, lags=0:10)[:, :]), dims=1)', lw=3, label="HMC")
plot!(mean(abs.(autocor(c3, lags=0:10)[:, :]), dims=1)', lw=3, label="RAM")


histogram(s2[a2, rand(1:model.d, 1)] |> m2t, bins=100, normalize=true)
plot!(x -> pdf(Normal(0, 1), x), lw=3)


W2_minibatch(s1[a1, :], randn(size(s1[a1, :])), N=100, k=64)
W2_minibatch(s2[a2, :], randn(size(s2[a2, :])), N=100, k=64)
W2_minibatch(s3[a3, :], randn(size(s3[a3, :])), N=100, k=64)
scatter(contourplot(model), s1[a1, 1:2] |> m2t, ma=0.1, c=:green, msw=0)
contourplot(model)





x = [randn(2) for _ in 1:10]
bs = [randn() for _ in 1:10]
y = [randn(2) for _ in 1:10]
as = [randn() for _ in 1:10]



tmp = Zygote.gradient(x -> (bs'model.U.(x)) .+ NNlib.logsumexp((1 .- bs) .* model.U.(x) .- log(10)), x)[1]

state = (; q=x, b=bs)
dV = state -> Zygote.gradient(state_ -> Zygote.forwarddiff(U, state_), state)[1]
V(state)
dV(state)



Y = [x bs]
foo(x, bs) = (bs'model.U.(x)) .+ NNlib.logsumexp((1 .- bs) .* model.U.(x) .- log(10))
dfoo_x(x, bs) = Zygote.gradient(x_ -> foo(x_, bs), x)[1]
dfoo_bs(x, bs) = Zygote.gradient(bs_ -> foo(x, bs_), bs)[1]
foo(x, bs)
dfoo_x(x, bs)
dfoo_bs(x, bs)




[[randn(2)]; x]

m = PEHMC(ϵ=0.5, L=20, N=3)
Σ = [1.0 0.5; 0.5 1.0]
R = [0.0 -1.0; 1.0 0.0]
model = Model(ξ=MixtureModel([MvNormal([3.0, 3.0], Σ), MvNormal([-3.0, -3.0], R * Σ * R')]))


s, w, a = mcmc(m, model; n=1e3, n_burn=1e2)
s, w = s[a, :], w[a, :]
begin
    X = [s...] |> a2m
    histogram(X[:, 1], weights=[w...], bins=100, normalize=true)
    plot!(x -> 0.5 * (pdf(Normal(-3), x) + pdf(Normal(3), x)), lw=3)
    contourplot(model)
    scatter!(X |> m2t, ms=[w...] .* 5, c=:green)
end

st = InitializeState(randn(2), m, model)
begin
    st = RefreshState(st.q, st, m, model)
    newst, w, a = OneStep(st, m, model)
    println(a)
    contourplot(model)
    if a
        st = newst
        scatter!(map(x -> tuple(x...), newst.q), ms=w .* 10, c=:green)
    else
        scatter!(map(x -> tuple(x...), newst.q), ms=w .* 10, c=:red)
    end
end



plot(X[:, 2])


s1, a1 = mcmc(HaRAM(ϵ=0.5, L=20, γ=0.1), model; n=1e4, n_burn=1e4);
X1 = s1[a1, :]

histogram(X1[:, 1], bins=100, normalize=true)
plot!(x -> 0.5 * (pdf(Normal(-3), x) + pdf(Normal(3), x)), lw=3)
contourplot(model)
scatter!(X1 |> m2t, c=:green, ma=0.05)
plot(X1[:, 1])


foo = (x, y) -> model.U([x; y])
xseq = range(-10, 10; length=1000)
heatmap(xseq, xseq, (x, y) -> foo(x, y) ./ 10000)




s2, a2 = mcmc(
    main.PT(τ=[1.0, 0.1], swap_rate=1.0),
    HaRAM(ϵ=0.5, L=20, γ=0.1), 
    model; n=1e4, n_burn=1e4
);

contourplot(model)
scatter!(s2[a2, :] |> m2t, ma=0.5)
histogram(s2[a2, 2], bins=100, normalize=true)
plot!(x -> 0.5 * (pdf(Normal(-3), x) + pdf(Normal(3), x)), lw=3)
plot(s2[a2, 2])