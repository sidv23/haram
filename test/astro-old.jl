using CSV, DataFrames

function log_likelihood(B, param::AbstractVector{T}, data) where {T}

    delta = param[1]
    c = param[2]
    mu = param[3]
    sigma = param[4]
    sqrttau = param[5]

    tau = sqrttau^2

    time = data[:, 1]
    leng_time = length(time)

    lcA = data[:, 2]
    se_lcA = data[:, 3]
    lcB = data[:, 4]
    se_lcB = data[:, 5]

    # sorting time given delta
    time_d = time .- delta
    time_temp = [time; time_d]
    ord = sortperm(time_temp)
    time_comb = time_temp[ord]
    leng_time_comb = length(time_comb)

    # indicator taking on 1 for X(t - delta) and 0 for X(t)
    ind = [fill(0, leng_time); fill(1, leng_time)][ord]

    lc_temp = [lcA; lcB .- c]
    lc_comb = lc_temp[ord]
    se_lc_temp = [se_lcA; se_lcB]
    se_lc_comb = se_lc_temp[ord]

    # x_star_i, i = 1, 2, ..., 2n
    x = lc_comb .- mu

    # omega_i, i = 1, 2, ..., 2n to be saved
    B = fill(T(1), leng_time_comb) 

    # x_hat_i, i = 1, 2, ..., 2n to be saved
    mu_i = fill(T(-1), leng_time_comb)
    mu_star_i = fill(T(-1), leng_time_comb)

    # a_i, i = 2, ..., 2n
    a_i = exp.(-1 .* (time_comb[2:end] .- time_comb[1:end-1]) ./ tau)

    # omega_i, i = 1, 2, ..., 2n
    var0 = tau * sigma^2 / 2

    B[1] = se_lc_comb[1]^2 / (se_lc_comb[1]^2 + var0)
    for k in 2:leng_time_comb
        B[k] = se_lc_comb[k]^2 / (se_lc_comb[k]^2 +
                                  a_i[k-1]^2 * (1 - B[k-1]) * se_lc_comb[k-1]^2 +
                                  var0 * (1 - a_i[k-1]^2))
    end

    # # x_hat_i, i = 1, 2, ..., 2n
    mu_i[1] = (1 - B[1]) * x[1]
    for k in 2:leng_time_comb
        mu_i[k] = (1 - B[k]) * x[k] + B[k] * a_i[k-1] * mu_i[k-1]
    end

    # mu_star_i = [0.0; a_i .* mu_i[1:leng_time_comb-1]]
    mu_star_i[2:leng_time_comb] = a_i .* mu_i[1:leng_time_comb-1]

    var_star_i = se_lc_comb .^ 2 ./ B

    # log-likelihood
    # sum(dnorm(x, mean = mu_star_i, sd = sqrt(var_star_i), log = TRUE))
    sum([logpdf(Normal(u, v), w) for (u, v, w) in zip(mu_star_i, var_star_i, x)])
end


function log_prior(param)
    delta = param[1] # Unif(-1178.939, 1178.939)
    c = param[2] # Unif(-60, 60)
    mu = param[3] # Unif(-30, 30)
    sigma = param[4] # sigma^2 ~ inverse-Gamma(1, 2 * 10^(-7))
    tau = param[5] # inverse-Gamma(1, 1)

    if (delta < -1178.939 || delta > 1178.939 || c < -60 || c > 60 || mu < -30 || mu > 30)
        -Inf
    else
        -2 * log(sigma^2) - 2 * 10^(-7) / sigma^2 - 2 * log(tau^2) - 1 / tau^2
    end
end

df = CSV.File("data/q0957usno.csv") |> DataFrame

# param = [420, 0.6, -18, 0.01^2, 200.0]
param = [1050, 0.6, -18, 0.01, 200.0]
log_likelihood(param, df)
log_prior(param)
ForwardDiff.gradient(U, [-333420, 0.6, -18, 0.01^2, 200.0])


function U(x)
    log_likelihood(x, df) + log_prior(x)
end

model = Model(
    ξ=Normal(5, 1.0),
    d=5,
    f=x -> exp(-U(x)),
    U=U,
    dU=x -> ForwardDiff.gradient(U, x)
)

st = (; q=param, p=randn(5), m=1.0)
OneStep(st, m, model)

# m = HaRAM(ϵ=0.2, L=10, γ=0.1)
m = RAM(kernel=MvNormal(5, 3), z=randn(5))
samples, accepts = mcmc(m, model; n=100, n_burn=100, init=[0.0, -0.113, 2.658, 0.01, 200])