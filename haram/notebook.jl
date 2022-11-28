### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 8902a846-fbb9-42fc-8742-c9c4a84db52c
begin
    import Pkg
    # Pkg.activate(mktempdir())
  #   Pkg.add([
		# Pkg.PackageSpec(name="BenchmarkTools", version="1.0.0"),
  #       Pkg.PackageSpec(name="CSV", version="0.8.5"),
  #       Pkg.PackageSpec(name="Chain", version="0.4.6"),
  #       Pkg.PackageSpec(name="DataFrames", version="1.1.1"),
  #       Pkg.PackageSpec(name="DifferentialEquations", version="6.17.1"),
  #       Pkg.PackageSpec(name="Distributions", version="0.24.18"),
		# Pkg.PackageSpec(name="LaTeXStrings", version="1.2.1"),
		# Pkg.PackageSpec(name="LazyArrays", version="0.21.5"),
		# Pkg.PackageSpec(name="Plots", version="1.16.2"),
		# Pkg.PackageSpec(name="PlutoUI", version="0.7.9"),
		# Pkg.PackageSpec(name="StatsBase", version="0.33.8"),
		# Pkg.PackageSpec(name="StatsPlots", version="0.14.21"),
		# Pkg.PackageSpec(name="Turing", version="0.16.0")
  #   ])
	using BenchmarkTools
	using CSV
	using DataFrames
	using DifferentialEquations
	using Distributions
	using LaTeXStrings
	using LazyArrays
	using LinearAlgebra
	using Random
	using StatsBase
	using StatsPlots
	using Turing
	using Plots
	using PlutoUI
	using LinearAlgebra: qr
	using Statistics: mean, std
end

# ╔═╡ 31161289-1d4c-46ba-8bd9-e687fb7da29e
begin
	using InteractiveUtils
	with_terminal() do
		versioninfo()
	end
end

# ╔═╡ 4af78efd-d484-4241-9d3c-97cc78e1dbd4
begin
	Turing.setprogress!(false);
	Random.seed!(1);
end

# ╔═╡ 5df4d7d2-c622-11eb-3bbd-bff9668ee5e0
md"""
# Turing Workshop
"""

# ╔═╡ 19c63110-4baa-4aff-ab91-46e5c149f3a2
Resource("https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg", :width => 120, :display => "inline")

# ╔═╡ dceb8312-230f-4e4b-9285-4e23f219b838
Resource("https://github.com/storopoli/Turing-Workshop/blob/master/images/bayes-meme.jpg?raw=true", :width => 250, :display=>"center")

# ╔═╡ cda7dc96-d983-4e31-9298-6148205b54b1
md"""
A little bit about myself:

$(Resource("https://github.com/storopoli/Turing-Workshop/blob/master/images/profile_pic.jpg?raw=true", :width => 100, :align => "right"))

* **Jose Storopoli**, PhD 🌐 [storopoli.io](https://storopoli.io)
* Associate Professor at [**Universidade Nove de Julho** (UNINOVE)](https://uninove.br)
* Teach undergraduates [**Statistics** and **Machine Learning** (using Python](https://storopoli.io/ciencia-de-dados) 😓, but I'm starting ot migrate the content to [**Julia**](https://storopoli.io/Julia-ciencia-de-dados/) 🚀)
* Teach graduate students [**Bayesian Statistics** (using `Stan`)](https://storopoli.io/Estatistica-Bayesiana) and [**Scientific Computing** (using **Julia**](https://storopoli.io/Computacao-Cientifica) 🚀)
* I've made some `Turing` tutorials, you can check them out at [storopoli.io/Bayesian-Julia](https://storopoli.io/Bayesian-Julia)
* You can find me on [Twitter](https://twitter.com/JoseStoropoli) (altough I rarelly use it) or on [LinkedIn](https://www.linkedin.com/in/storopoli/)
"""

# ╔═╡ 2164bf58-75ff-470c-828c-b0165f0d980d
md"""
This workshop can be found in a CreativeCommons [YouTube Video](https://youtu.be/CKSxxJ7RdAU)
"""

# ╔═╡ 55777e4a-b197-4a61-8e57-6ae9792c0564
html"""
<iframe width="560" height="315" src="https://www.youtube.com/embed/CKSxxJ7RdAU" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ 1436305e-37d8-44f1-88d6-4de838580360
md"""
## Bayesian Statistics?!

**Bayesian statistics** is an approach to inferential statistics based on Bayes' theorem, where available knowledge about parameters in a statistical model is updated with the information in observed data. The background knowledge is expressed as a prior distribution and combined with observational data in the form of a likelihood function to determine the posterior distribution. The posterior can also be used for making predictions about future events.

$$\underbrace{P(\theta \mid y)}_{\text{Posterior}} = \frac{\overbrace{P(y \mid  \theta)}^{\text{Likelihood}} \cdot \overbrace{P(\theta)}^{\text{Prior}}}{\underbrace{P(y)}_{\text{Normalizing Costant}}}$$

> No $p$-values! Nobody knows what they are anyway... Not $P(H_0 \mid y)$
"""

# ╔═╡ 08f508c4-233a-4bba-b313-b04c1d6c4a4c
md"""
### Recommended Books
"""

# ╔═╡ 868d8932-b108-41d9-b4e8-d62d31b5465d
md"""
We are not covering Bayesian stuff, but there are some **awesome books**:
"""

# ╔═╡ 653ec420-8de5-407e-91a9-f045e25a6395
md"""
[$(Resource("https://github.com/storopoli/Turing-Workshop/blob/master/images/BDA_book.jpg?raw=true", :width => 100.5*1.5))](https://www.routledge.com/Bayesian-Data-Analysis/Gelman-Carlin-Stern-Dunson-Vehtari-Rubin/p/book/9781439840955)
[$(Resource("https://github.com/storopoli/Turing-Workshop/blob/master/images/SR_book.jpg?raw=true", :width => 104*1.5))](https://www.routledge.com/Statistical-Rethinking-A-Bayesian-Course-with-Examples-in-R-and-STAN/McElreath/p/book/9780367139919)
[$(Resource("https://github.com/storopoli/Turing-Workshop/blob/master/images/ROS_book.jpg?raw=true", :width => 118*1.5))](https://www.cambridge.org/fi/academic/subjects/statistics-probability/statistical-theory-and-methods/regression-and-other-stories)
[$(Resource("https://github.com/storopoli/Turing-Workshop/blob/master/images/Bayes_book.jpg?raw=true", :width => 102*1.5))](https://www.amazon.com/Theory-That-Would-Not-Die/dp/0300188226/)
"""

# ╔═╡ 716cea7d-d771-46e9-ad81-687292004009
md"""
## 1. What is Turing?
"""

# ╔═╡ cb808fd4-6eb2-457e-afa4-58ae1be09aec
md"""
[**`Turing`** (Ge, Xu & Ghahramani, 2018)](http://turing.ml/) is a ecosystem of Julia packages for Bayesian Inference using [probabilistic programming](https://en.wikipedia.org/wiki/Probabilistic_programming). Models specified using `Turing` are easy to read and write -- models work the way you write them. Like everything in Julia, `Turing` is fast [(Tarek, Xu, Trapp, Ge & Ghahramani, 2020)](https://arxiv.org/abs/2002.02702).

Before we dive into how to specify models in Turing. Let's discuss Turing's **ecosystem**.
We have several Julia packages under the Turing's GitHub organization [TuringLang](https://github.com/TuringLang), but I will focus on 5 of those:

* [`Turing.jl`](https://github.com/TuringLang/Turing.jl): main package that we use to **interface with all the Turing ecosystem** of packages and the backbone of the PPL Turing.

* [`MCMCChains.jl`](https://github.com/TuringLang/MCMCChains.jl): is an interface to **summarizing MCMC simulations** and has several utility functions for **diagnostics** and **visualizations**.

* [`DynamicPPL.jl`](https://github.com/TuringLang/DynamicPPL.jl): which specifies a domain-specific language and backend for Turing (which itself is a PPL), modular and written in Julia

* [`AdvancedHMC.jl`](https://github.com/TuringLang/AdvancedHMC.jl): modular and efficient implementation of advanced HMC algorithms. The state-of-the-art HMC algorithm is the **N**o-**U**-**T**urn **S**ampling (NUTS) (Hoffman & Gelman, 2011)

* [`DistributionsAD.jl`](https://github.com/TuringLang/DistributionsAD.jl): defines the necessary functions to enable automatic differentiation (AD) of the `logpdf` function from [`Distributions.jl`](https://github.com/JuliaStats/Distributions.jl) using the packages [`Tracker.jl`](https://github.com/FluxML/Tracker.jl), [`Zygote.jl`](https://github.com/FluxML/Zygote.jl), [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) and [`ReverseDiff.jl`](https://github.com/JuliaDiff/ReverseDiff.jl). The main goal of `DistributionsAD.jl` is to make the output of `logpdf` differentiable with respect to all continuous parameters of a distribution.
"""

# ╔═╡ 0484ae7f-bd8a-4615-a760-5c4b2eef9d3f
md"""
## 2. How to Specify a Model? `@model`
"""

# ╔═╡ 1d467044-bc7d-4df7-bda6-bb8ea6ff0712
md"""
**We specify the model inside a macro** `@model` where we can assign variables in two ways:

* using `~`: which means that a variable follows some probability distribution (Normal, Binomial etc.) and its value is random under that distribution

* using `=`: which means that a variable does not follow a probability distribution and its value is deterministic (like the normal `=` assignment in programming languages)

Turing will perform automatic inference on all variables that you specify using `~`.

Just like you would write in mathematical form:

$$\begin{aligned}
p &\sim \text{Beta}(1,1) \\
\text{coin flip} &\sim \text{Bernoulli}(p)
\end{aligned}$$

> **Example**: Unfair coin with $p$ = 0.7.
"""

# ╔═╡ b1d99482-53f5-4c6b-8c20-c761ff6bdb77
coin_flips = rand(Bernoulli(0.7), 100);

# ╔═╡ 65fa382d-4ef7-432d-8630-27082977185b
@model coin(coin_flips) = begin
	p ~ Beta(1, 1)
	for i ∈ 1:length(coin_flips)
		coin_flips[i] ~ Bernoulli(p)
	end
end;

# ╔═╡ 06f93734-2315-4b36-a39a-09e8167bab1f
begin
	chain_coin = sample(coin(coin_flips), MH(), 100);
	summarystats(chain_coin)
end

# ╔═╡ 9f6b96a7-033d-4c7d-a853-46a0b5af4675
md"""
## 3. How to specify a MCMC sampler (`NUTS`, `HMC`, `MH` etc.)
"""

# ╔═╡ b7667fb4-6e76-4711-b61d-dae5f993531e
md"""
We have [several samplers](https://turing.ml/dev/docs/using-turing/sampler-viz) available:

* `MH()`: **M**etropolis-**H**astings
* `PG()`: **P**article **G**ibbs
* `SMC()`: **S**equential **M**onte **C**arlo
* `HMC()`: **H**amiltonian **M**onte **C**arlo
* `HMCDA()`: **H**amiltonian **M**onte **C**arlo with Nesterov's **D**ual **A**veraging
* `NUTS()`: **N**o-**U**-**T**urn **S**ampling

Just stick your desired `sampler` inside the function `sample(model, sampler, N; kwargs)`.

Play around if you want. Choose your `sampler`:
"""

# ╔═╡ cb168dc1-70e2-450f-b2cf-c8680251ab27
@bind chosen_sampler Radio([
		"MH()",
		"PG(Nₚ) - Number of Particles",
		"SMC()",
		"HMC(ϵ, L) - leaprog step size(ϵ) and number of leaprogs steps (L)",
		"HMCDA(Nₐ, δ, λ) - Number of samples to use for adaptation (Nₐ), target acceptance ratio (δ), and target leapfrog length(λ)",
		"NUTS(Nₐ, δ) - Number of samples to use for adaptation (Nₐ) and target acceptance ratio (δ)"], default = "MH()")

# ╔═╡ 07d408cf-d202-40b2-90c2-5e8630549339
begin
	your_sampler = nothing
	if chosen_sampler == "MH()"
		your_sampler = MH()
	elseif chosen_sampler == "PG(Nₚ) - Number of Particles"
		your_sampler = PG(2)
	elseif chosen_sampler == "SMC()"
		your_sampler = SMC()
	elseif chosen_sampler == "HMC(ϵ, L) - leaprog step size(ϵ) and number of leaprogs steps (L)"
		your_sampler = HMC(0.05, 10)
	elseif chosen_sampler == "HMCDA(Nₐ, δ, λ) - Number of samples to use for adaptation (Nₐ), target acceptance ratio (δ), and target leapfrog length(λ)"
		your_sampler = HMCDA(10, 0.65, 0.3)
	elseif chosen_sampler == "NUTS(Nₐ, δ) - Number of samples to use for adaptation (Nₐ) and target acceptance ratio (δ)"
		your_sampler = NUTS(10, 0.65)
	end
end

# ╔═╡ 744a8a63-647f-4550-adf7-44354fde44be
begin
	chain_coin_2 = sample(coin(coin_flips), your_sampler, 100); # Here is your sampler
	summarystats(chain_coin_2)
end

# ╔═╡ e6365296-cd68-430e-99c5-fb571f39aad5
md"""
### 3.1 MOAH CHAINS!!: `MCMCThreads` and `MCMCDistributed`
"""

# ╔═╡ 927ad0a4-ba68-45a6-9bde-561915503e48
md"""
There is some methods of `Turing`'s `sample()` that accepts either:

* `MCMCThreads()`: uses multithread stuff with [`Threads.jl`](https://docs.julialang.org/en/v1/manual/multi-threading/#man-multithreading)
* `MCMCDistributed()`: uses multiprocesses stuff with [`Distributed.jl`](https://docs.julialang.org/en/v1/manual/distributed-computing/) and uses the [MPI -- Message Passing Interface](https://en.wikipedia.org/wiki/Message_Passing_Interface)


> If you are using `MCMCDistributed()` don't forget the macro `@everywhere` and the `addprocs()` stuff

Just use `sample(model, sampler, MCMCThreads(), N, chains)`

Let's revisit our biased-coin example:
"""

# ╔═╡ ab6c2ba6-4cd8-473a-88c6-b8d61551fb22
begin
	chain_coin_parallel = sample(coin(coin_flips), MH(), MCMCThreads(), 2_000, 2);
	summarystats(chain_coin_parallel)
end

# ╔═╡ 2ab3c34a-1cfc-4d20-becc-5902d08d03e0
md"""
### 3.2 LOOK MUM NO DATA!!: Prior Predictive Checks `Prior()`
"""

# ╔═╡ 924fcad9-75c1-4707-90ef-3e36947d64fe
md"""
It's very important that we check if our **priors make sense**. This is called **Prior Predictive Check** (Gelman et al., 2020b). Obs: I will not cover **Posterior Predictive Check** because is mostly the same procedure in `Turing`.
"""

# ╔═╡ fc8e40c3-34a1-4b2e-bd1b-893d7998d359
md"""
$(Resource("https://github.com/storopoli/Turing-Workshop/blob/master/images/bayesian_workflow.png?raw=true", :width => 700))

Based on Gelman et al. (2020b)
"""

# ╔═╡ fb366eb1-4ab0-4e7a-83ed-d531978c06a0
md"""
Predictive checks are a great way to **validate a model**. The idea is to **generate data from the model** using **parameters from draws from the prior or posterior**. *Prior predictive check* is when we simulate data using model's parameters values drawn fom the *prior* distribution, and *posterior* predictive check is is when we simulate data using model's parameters values drawn fom the *posterior* distribution.

The workflow we do when specifying and sampling Bayesian models is not linear or acyclic (Gelman et al., 2020b). This means that we need to iterate several times between the different stages in order to find a model that captures best the data generating process with the desired assumptions.

This is quite easy in `Turing`. We need to create a *prior* distribution for our model. To accomplish this, instead of supplying a MCMC sampler like `NUTS()` or `MH()`, we supply the "sampler" `Prior()` inside `Turing`'s `sample()` function:
"""

# ╔═╡ 0fe83f55-a379-49ea-ab23-9defaab05890
begin
	prior_chain_coin = sample(coin(coin_flips), Prior(), 1_000)
	summarystats(prior_chain_coin)
end

# ╔═╡ 3aa95b4b-aaf8-45cf-8bc5-05b65b4bcccf
md"""
Now we can perform predictive checks using both the *prior* (`prior_chain_coin`) or *posterior* (`chain_coin`) distributions. To draw from the prior and posterior predictive distributions we instantiate a "predictive model", i.e. a `Turing` model but with the observations set to `missing`, and then calling `predict()` on the predictive model and the previously drawn samples.

Let's do the *prior* predictive check:
"""

# ╔═╡ dd27ee5f-e442-42d7-a39b-d76328d2e59f
begin
	missing_data = Vector{Missing}(missing, length(coin_flips)); # vector of `missing`
	model_predict = coin(missing_data); # instantiate the "predictive model"
	prior_check = predict(model_predict, prior_chain_coin);
	describe(DataFrame(summarystats(prior_check)))
end

# ╔═╡ c4808b43-bc0f-4254-abf1-1adc19135dc7
md"""
### 3.3 Posterior Predictive Checks

The *posterior* predictive check is trivial, just do the same but with the posterior `chain_coin`:
"""

# ╔═╡ 1773d8c3-4651-4128-9442-e7c858bc4a43
begin
	posterior_check = predict(model_predict, chain_coin);
	describe(DataFrame(summarystats(posterior_check)))
end

# ╔═╡ 5674f7aa-3205-47c7-8367-244c6419ce69
md"""
## 4. How to inspect chains and plot stuff with `MCMCChains.jl`
"""

# ╔═╡ 83cc80c1-d97e-4b82-872e-e5493d2b62ab
md"""
We can inspect and plot our model's chains and its underlying parameters with [`MCMCChains.jl`](https://turinglang.github.io/MCMCChains.jl/stable/)

1. **Inspecting Chains**
   * **Summary Statistics**: just do `summarystats(chain)`
   * **Quantiles** (Median, etc.): just do `quantile(chain)`
   * What if I just want a **subset** of parameters?: just do `group(chain, :parameter)` or index with `chain[:, 1:6, :]` or `chain[[:parameters,...]]`
"""

# ╔═╡ 475be60f-1876-4086-9725-3bf5f52a3e43
summarystats(chain_coin_parallel)

# ╔═╡ f6bc0cfd-a1d9-48e5-833c-f33bf1b89d45
quantile(chain_coin_parallel)

# ╔═╡ ed640696-cae6-47e1-a4df-0655192e0855
quantile(group(chain_coin_parallel, :p))

# ╔═╡ bc9fa101-8854-4af5-904a-f0b683fb63b1
summarystats(chain_coin_parallel[:, 1:1, :])

# ╔═╡ c82687d1-89d0-4ecd-bed7-1708ba8b2662
md"""
2. **Plotting Chains**: Now we have several options. The default `plot()` recipe will plot a `traceplot()` side-by-side with a `mixeddensity()`.

   First, we have to choose either to plot **parameters**(`:parameter`) or **chains**(`:chain`) with the keyword `colordim`.
"""

# ╔═╡ 270c0b90-cce1-4092-9e29-5f9deda2cb7d
plot(chain_coin_parallel; colordim=:chain, dpi=300)

# ╔═╡ c4146b8b-9d11-446e-9765-8d5283a6d445
plot(chain_coin_parallel; colordim=:parameter, dpi=300)

# ╔═╡ 3d09c8c3-ce95-4f26-9136-fedd601e2a70
md"""
Second, we have several plots to choose from:
* `traceplot()`: used for inspecting Markov chain **convergence**
* `meanplot()`: running average plots per interaction
* `density()`: **density** plots
* `histogram()`: **histogram** plots
* `mixeddensity()`: **mixed density** plots
* `autcorplot()`: **autocorrelation** plots
"""

# ╔═╡ 8d9bdae2-658d-45bf-9b25-50b6efbe0cdf
plot(
	traceplot(chain_coin_parallel, title="traceplot"),
	meanplot(chain_coin_parallel, title="meanplot"),
	density(chain_coin_parallel, title="density"),
	histogram(chain_coin_parallel, title="histogram"),
	mixeddensity(chain_coin_parallel, title="mixeddensity"),
	autocorplot(chain_coin_parallel, title="autocorplot"),
	dpi=300, size=(840, 600)
)

# ╔═╡ 41b014c2-7b49-4d03-8741-51c91b95f64c
md"""
There is also the option to **construct your own plot** with `plot()` and the keyword `seriestype`:
"""

# ╔═╡ 2f08c6e4-fa7c-471c-ad9f-9d036e3027d5
plot(chain_coin_parallel, seriestype = (:meanplot, :autocorplot), dpi=300)

# ╔═╡ 5f639d2d-bb96-4a33-a78e-d5b9f0e8d274
md"""
Finally there is one special plot that makes a **cornerplot** (requires `StatPlots`) of parameters in a chain:

> Obs: I will hijack a multi-parameter model from *below* to show the cornerplot
"""

# ╔═╡ c70ebb70-bd96-44a5-85e9-871b0e478b1a
md"""
## 5. Better tricks to avoid `for`-loops inside `@model` (`lazyarrays` and `filldist`)
"""

# ╔═╡ 36258bdd-f617-48f6-91c9-e8bbff78ebd8
md"""
**Using Logistic Regression**
"""

# ╔═╡ 6630eb47-77f6-48e9-aafe-55bda275449c
md"""
First the Naïve model *with* `for`-loops:
"""

# ╔═╡ 37e751c7-8b6c-47d9-8013-97015d1e1fb2
@model logreg(X,  y; predictors=size(X, 2)) = begin
	#priors
	α ~ Normal(0, 2.5)
	β = Vector{Float64}(undef, predictors)
	for i ∈ 1:predictors
		β[i] ~ Normal()
	end

	#likelihood
	for i ∈ 1:length(y)
		y[i] ~ BernoulliLogit(α +  X[i, :] ⋅ β)
	end
end;

# ╔═╡ 7a21e7a0-322b-4f8e-9d8b-a2f452f7e092
md"""
* `Turing`'s `BernoulliLogit()` is a logit-parameterised Bernoulli distribution that convert logodds to probability.
"""

# ╔═╡ f8f59ebb-bb1e-401f-97b5-507634badb3f
md"""
Now a model *without* `for`-loops
"""

# ╔═╡ 15795f79-7d7b-43d2-a4b4-99ad968a7f72
@model logreg_vectorized(X,  y; predictors=size(X, 2)) = begin
	#priors
	α ~ Normal(0, 2.5)
	β ~ filldist(Normal(), predictors)

	#likelihood
	# y .~ BernoulliLogit.(α .+ X * β)
	y ~ arraydist(LazyArray(@~ BernoulliLogit.(α .+ X * β)))
end;

# ╔═╡ dd5fbb2a-4220-4e47-945a-6870b799c50d
md"""
* `Turing`'s `arraydist()` function wraps an array of distributions returning a new distribution sampling from the individual distributions.

* `LazyArrays`' `LazyArray()` constructor wrap a lazy object that wraps a computation producing an `array` to an `array`. Last, but not least, the macro `@~` creates a broadcast and is a nice short hand for the familiar dot `.` broadcasting operator in Julia. This is an efficient way to tell Turing that our `y` vector is distributed lazily as a `BernoulliLogit` broadcasted to `α` added to the product of the data matrix `X` and `β` coefficient vector.
"""

# ╔═╡ 0cc8e12c-9b72-41ec-9c13-d9ae0bdc6100
md"""
For our example, I will use a famous dataset called `wells` (Gelman & Hill, 2007), which is data from a survey of 3,200 residents in a small area of Bangladesh suffering from arsenic contamination of groundwater. Respondents with elevated arsenic levels in their wells had been encouraged to switch their water source to a safe public or private well in the nearby area and the survey was conducted several years later to learn which of the affected residents had switched wells. It has 3,200 observations and the following variables:

* `switch` – binary/dummy (0 or 1) for well-switching.

* `arsenic` – arsenic level in respondent's well.

* `dist` – distance (meters) from the respondent's house to the nearest well with safe drinking water.

* `association` – binary/dummy (0 or 1) if member(s) of household participate in community organizations.

* `educ` – years of education (head of household).
"""

# ╔═╡ fce0f511-3b00-4079-85c6-9b2d2d7c04cb
begin
	# Logistic Regression
	wells = CSV.read(download("https://github.com/storopoli/Turing-Workshop/blob/master/data/wells.csv?raw=true"), DataFrame);
	X_wells = Matrix(select(wells, Not(:switch)));
	y_wells = wells[:, :switch];
end

# ╔═╡ 5ba6b247-8277-4100-abe7-8d06af04a011
md"""
Why do that?

1. Well, you'll have nice performance gains
"""

# ╔═╡ 0f000fc4-1a7b-4522-8355-8df572ee8800
with_terminal() do
	@btime sample(logreg($X_wells, $y_wells), MH(), 100);
end

# ╔═╡ 8a87e324-f3d9-4162-88ab-3833a6d1fc2e
with_terminal() do
	@btime sample(logreg_vectorized($X_wells, $y_wells), MH(), 100);
end

# ╔═╡ 3c954cbc-aed7-4d22-b578-a80ce62ebb49
md"""
2. Some [autodiff backends only works without `for`-loops inside the `@model`](https://turing.ml/dev/docs/using-turing/performancetips#special-care-for-codetrackercode-and-codezygotecode):
   * [`Tracker.jl`](https://github.com/FluxML/Tracker.jl)
   * [`Zygote.jl`](https://github.com/FluxML/Zygote.jl)
"""

# ╔═╡ 521e2473-1aba-43be-951a-25537062891e
md"""
### 5.1 Which [autodiff backend](https://turing.ml/dev/docs/using-turing/autodiff) to use?
"""

# ╔═╡ bafc91d2-8cae-4af8-b5ed-8199eef40c4d
md"""
We have mainly two [types of autodiff](https://en.wikipedia.org/wiki/Automatic_differentiation) (both uses the chain rule $\mathbb{R}^N \to \mathbb{R}^M$)

* **Forward Autodiff**: The **independent** variable is fixed and differentiation is performed in a *forward* manner. Preffered when $N < M$
   * [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl): current (version 0.16) Turing's default, `:forwarddiff`

* **Reverse Autodiff**: The **dependent** variable is fixed and differentiation is performed in a *backward* manner. Preffered when $N > M$
   * [`Tracker.jl`](https://github.com/FluxML/Tracker.jl): `:tracker`
   * [`Zygote.jl`](https://github.com/FluxML/Zygote.jl): `:zygote`
   * [`ReverseDiff.jl`](https://github.com/JuliaDiff/ReverseDiff.jl): `:reversediff`

Checkout this video is awesome to learn what Automatic Differentiation is!
"""

# ╔═╡ a2292bc1-3379-450d-beb5-ae8f41b69be8
html"""<iframe width="560" height="315" src="https://www.youtube.com/embed/wG_nF1awSSY" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>"""

# ╔═╡ 38055b57-f983-4440-bef5-0ab6d180ff1e
md"""
To change `Turing`'s autodiff backend just type:

```julia
Turing.setadbackend(:zygote)
```

or 

```julia
Turing.setadbackend(:tracker)
```

Note that you need to import the backend:

```julia
using Zygote
```
"""

# ╔═╡ 7d4d06ca-f96d-4b1e-860f-d9e0d6eb6723
md"""
## 6. Take me up! Let's get Hierarchical (Hierarchical Models)
"""

# ╔═╡ c64d355f-f5a2-46a5-86f3-2d02da98f305
md"""
Bayesian **hierarchical** models (also called **multilevel** models) are a statistical model written at **multiple levels** (hierarchical form) that estimates the parameters of the posterior distribution using the Bayesian approach. The sub-models combine to form the hierarchical model, and **Bayes' theorem is used to integrate them with the observed data** and to account for all the **uncertainty** that is present.

Hierarchical modeling is used when **information is available at several different levels of observation units**. The hierarchical form of analysis and organization helps to understand multiparameter problems and also plays an important role in the development of computational strategies.
"""

# ╔═╡ 262cb245-0bc1-4a36-b0bc-de52c08ccde0
md"""
"Even though observations directly inform only a single set of parameters, the latent population model couples the individual parameters and provides a backdoor for observations to inform all of the contexts. For example the observations from the $k$th context, $y_k$, directly inform the parameters that quantify the behavior of that context, $\theta_k$. Those parameters, however, directly inform the population parameters $\phi$ which then inform all of the other contexts through the population model. Similarly observations that directly inform the other contexts indirectly inform the population parameters which then feeds back into the $k$th context."

[Betancourt (2020)](https://betanalpha.github.io/assets/case_studies/hierarchical_modeling.html)
"""

# ╔═╡ 3ecc92b8-6a10-4f51-93d7-72449e248dc2
Resource("https://github.com/storopoli/Turing-Workshop/blob/master/images/multilevel_models.png?raw=true", :width => 1_000)

# ╔═╡ a0c7ca50-3a3f-483c-ae01-fd774e0c072d
md"""
> figure adapted from [Michael Betancourt (CC-BY-SA-4.0)](https://betanalpha.github.io/assets/case_studies/hierarchical_modeling.html)
"""

# ╔═╡ cb3dd785-11ff-42fe-ab85-0dd03e45209e
md"""
### 6.1 Hyperprior

As the priors of the parameters are sampled from another prior of the hyperparameter (upper-level's parameter), which are called hyperpriors. This makes one group's estimates help the model to better estimate the other groups by providing more **robust and stable estimates**.

We call the global parameters as **population effects** (or population-level effects, also sometimes called fixed effects) and the parameters of each group as **group effects** (or group-level effects, also sometimes called random effects). That is why multilevel models are also known as mixed models in which we have both fixed effects and random effects.
"""

# ╔═╡ 4812f80e-79a9-4519-9e4d-a45127ca6a49
md"""
### 6.2 Three Approaches to Multilevel Models

Multilevel models generally fall into three approaches:

1. **Random-intercept model**: each group receives a **different intercept** in addition to the global intercept.

2. **Random-slope model**: each group receives **different coefficients** for each (or a subset of) independent variable(s) in addition to a global intercept.

3. **Random-intercept-slope model**: each group receives **both a different intercept and different coefficients** for each independent variable in addition to a global intercept.
"""

# ╔═╡ 318697fe-1fbc-4ac3-a2aa-5ecf775072d4
md"""
#### Random-Intercept Model

The first approach is the **random-intercept model** in which we specify a different intercept for each group,
in addition to the global intercept. These group-level intercepts are sampled from a hyperprior.

To illustrate a multilevel model, I will use a linear regression example with a Gaussian/normal likelihood function.
Mathematically a Bayesian multilevel random-slope linear regression model is:

$$\begin{aligned}
\mathbf{y} &\sim \text{Normal}\left( \alpha + \alpha_j + \mathbf{X} \cdot \boldsymbol{\beta}, \sigma \right) \\
\alpha &\sim \text{Normal}(\mu_\alpha, \sigma_\alpha) \\
\alpha_j &\sim \text{Normal}(0, \tau) \\
\boldsymbol{\beta} &\sim \text{Normal}(\mu_{\boldsymbol{\beta}}, \sigma_{\boldsymbol{\beta}}) \\
\tau &\sim \text{Cauchy}^+(0, \psi_{\alpha})\\
\sigma &\sim \text{Exponential}(\lambda_\sigma)
\end{aligned}$$
"""

# ╔═╡ 9acc7a1c-f638-4a2e-ad67-c16cff125c86
@model varying_intercept(X, idx, y; n_gr=length(unique(idx)), predictors=size(X, 2)) = begin
    # priors
    α ~ Normal(mean(y), 2.5 * std(y))       # population-level intercept
    β ~ filldist(Normal(0, 2), predictors)  # population-level coefficients
    σ ~ Exponential(1 / std(y))             # residual SD
    
	# prior for variance of random intercepts
    # usually requires thoughtful specification
    τ ~ truncated(Cauchy(0, 2), 0, Inf)     # group-level SDs intercepts
    αⱼ ~ filldist(Normal(0, τ), n_gr)       # group-level intercepts

    # likelihood
    ŷ = α .+ X * β .+ αⱼ[idx]
    y ~ MvNormal(ŷ, σ)
end;

# ╔═╡ 885fbe97-edd6-44d2-808d-8eeb1e9cb2b4
md"""
#### Random-Slope Model

The second approach is the **random-slope model** in which we specify a different slope for each group,
in addition to the global intercept. These group-level slopes are sampled from a hyperprior.

To illustrate a multilevel model, I will use a linear regression example with a Gaussian/normal likelihood function.
Mathematically a Bayesian multilevel random-slope linear regression model is:

$$\begin{aligned}
\mathbf{y} &\sim \text{Normal}\left( \alpha + \mathbf{X} \cdot \boldsymbol{\beta}_j \cdot \boldsymbol{\tau}, \sigma \right) \\
\alpha &\sim \text{Normal}(\mu_\alpha, \sigma_\alpha) \\
\boldsymbol{\beta}_j &\sim \text{Normal}(0, 1) \\
\boldsymbol{\tau} &\sim \text{Cauchy}^+(0, \psi_{\boldsymbol{\beta}})\\
\sigma &\sim \text{Exponential}(\lambda_\sigma)
\end{aligned}$$
"""

# ╔═╡ 7f526d1f-bd56-4e51-9f7b-ce6b5a2a1853
@model varying_slope(X, idx, y; n_gr=length(unique(idx)), predictors=size(X, 2)) = begin
    # priors
    α ~ Normal(mean(y), 2.5 * std(y))                   # population-level intercept
    σ ~ Exponential(1 / std(y))                         # residual SD
    
	# prior for variance of random slopes
    # usually requires thoughtful specification
    τ ~ filldist(truncated(Cauchy(0, 2), 0, Inf), n_gr) # group-level slopes SDs
    βⱼ ~ filldist(Normal(0, 1), predictors, n_gr)       # group-level standard normal slopes

    # likelihood
    ŷ = α .+ X * βⱼ * τ
    y ~ MvNormal(ŷ, σ)
end;

# ╔═╡ f7971da6-ead8-4679-b8cf-e3c35c93e6cf
md"""
#### Random-intercept-slope Model

The third approach is the **random-intercept-slope model** in which we specify a different intercept
and  slope for each group, in addition to the global intercept.
These group-level intercepts and slopes are sampled from hyperpriors.

To illustrate a multilevel model, I will use a linear regression example with a Gaussian/normal likelihood function.
Mathematically a Bayesian multilevel random-intercept-slope linear regression model is:

$$\begin{aligned}
\mathbf{y} &\sim \text{Normal}\left( \alpha + \alpha_j + \mathbf{X} \cdot \boldsymbol{\beta}_j \cdot \boldsymbol{\tau}_{\boldsymbol{\beta}}, \sigma \right) \\
\alpha &\sim \text{Normal}(\mu_\alpha, \sigma_\alpha) \\
\alpha_j &\sim \text{Normal}(0, \tau_{\alpha}) \\
\boldsymbol{\beta}_j &\sim \text{Normal}(0, 1) \\
\tau_{\alpha} &\sim \text{Cauchy}^+(0, \psi_{\alpha})\\
\boldsymbol{\tau}_{\boldsymbol{\beta}} &\sim \text{Cauchy}^+(0, \psi_{\boldsymbol{\beta}})\\
\sigma &\sim \text{Exponential}(\lambda_\sigma)
\end{aligned}$$
"""

# ╔═╡ 546726af-5420-4a4f-8c0c-fe96a2ba43bc
@model varying_intercept_slope(X, idx, y; n_gr=length(unique(idx)), predictors=size(X, 2)) = begin
    # priors
    α ~ Normal(mean(y), 2.5 * std(y))                    # population-level intercept
    σ ~ Exponential(1 / std(y))                          # residual SD
    
	# prior for variance of random intercepts and slopes
    # usually requires thoughtful specification
    τₐ ~ truncated(Cauchy(0, 2), 0, Inf)                 # group-level SDs intercepts
    τᵦ ~ filldist(truncated(Cauchy(0, 2), 0, Inf), n_gr) # group-level slopes SDs
    αⱼ ~ filldist(Normal(0, τₐ), n_gr)                   # group-level intercepts
    βⱼ ~ filldist(Normal(0, 1), predictors, n_gr)        # group-level standard normal slopes

    # likelihood
    ŷ = α .+ αⱼ[idx] .+ X * βⱼ * τᵦ
    y ~ MvNormal(ŷ, σ)
end;

# ╔═╡ 9ebac6ba-d213-4ed8-a1d5-66b841fafa00
md"""
## 7. Crazy Stuff
"""

# ╔═╡ 45c342fd-b893-46aa-b2ee-7c93e7a1d207
md"""
There is a **lot** of *crazy* stuff you can do with `Turing` and Bayesian models.

Here I will cover:

1. **Discrete Parameters (HMM)**

2. **Models with ODEs**
"""

# ╔═╡ d44c7baa-80d2-4fdb-a2de-35806477dd58
md"""
### 7.1 Discrete Parameters (HMM)
"""

# ╔═╡ c1b2d007-1004-42f5-b65c-b4e2e7ff7d8e
Resource("https://github.com/storopoli/Turing-Workshop/blob/master/images/HMM.png?raw=true", :width => 400)

# ╔═╡ c1dcfd47-9e25-470b-a1b3-ab66bfac59d6
md"""
 $\mu_1$ = $(@bind μ₁_sim Slider(1:1:10, default = 1, show_value=true))

 $\mu_2$ = $(@bind μ₂_sim Slider(1:1:10, default = 5, show_value=true))
"""

# ╔═╡ f1153918-0748-4400-ae8b-3b59f8c5d755
md"""
I **love** [`Stan`](https://mc-stan.org), use it on a daily basis. But `Stan` has some quirks. Particularly, NUTS and HMC samplers **cannot tolerate discrete parameters**.

Solution? We have to **marginalize** them.

First, I will show the `Stan` example of a Hidden Markov Model (HMM) with marginalization. And then let's see how `Turing` fare with the same problem.

"""

# ╔═╡ ad6c4533-cd56-4f6f-b10d-d7bc3145ba16
md"""
We have several ways to marginalize discrete parameters in HMM:

1. **Filtering** (a.k.a [Forward Algorithm](https://en.wikipedia.org/wiki/Forward_algorithm)) <---- we'll cover this one
2. **Smoothing** (a.k.a [Forward-Backward Algorithm](https://en.wikipedia.org/wiki/Forward%E2%80%93backward_algorithm))
3. **MAP Estimation** (a.k.a [Viterbi Algorithm](https://en.wikipedia.org/wiki/Viterbi_algorithm))

A very good reference is [Damiano, Peterson & Weylandt (2017)](https://github.com/luisdamiano/stancon18)

"""

# ╔═╡ 2ef397a6-f7fb-4fc2-b918-40ab545ce19f
md"""
#### Forward Algorithm
"""

# ╔═╡ 8d347172-2d26-4d9e-954d-b8924ed4c9e2
md"""
We define the forward variables $\boldsymbol{\alpha}_t$, beginning at time $t = 1$, as follows

$$\boldsymbol{\alpha}_1 = \boldsymbol{\delta}^{(1)} \boldsymbol{\Gamma}(y_1), \qquad \boldsymbol{\alpha}_{t} = \boldsymbol{\alpha}_{t-1} \boldsymbol{\Gamma} \boldsymbol{\delta}^{t-1}(y_{t})$$

where:

* Observed variable of interest at time $t$: $y_t$

* Unobserved transition probability matrix at time $t$: $\boldsymbol{\Gamma}^{(t)} \in \mathbb{R}^{j \times j}$

* Unobserved state distribution at time $t=1$: $\boldsymbol{\delta}^{(1)} \sim \text{Dirichlet}(\boldsymbol{\nu})$

Then, the marginal likelihood is obtained by summing over:

$$\sum^J_{j=1} \alpha_T (j) = \boldsymbol{\alpha}_T \mathbf{1}^{\top}$$

(dot product of the vector of $\alpha$s with a row vector of $1$s)

Note that one of the assumptions is $\boldsymbol{\Gamma}$ is **fullrank** (no linear dependence) and **ergodic** (it converges to a unique stationary distribution in $\lim_{t \to \infty}$)
"""

# ╔═╡ ca962c0e-4620-4888-b7c3-aa7f6d7899e9
md"""
As an example, I will use [Leos-Barajas & Michelot, 2018](http://arxiv.org/abs/1806.10639)'s.

It's a $2$-state HMM with Gaussian state-dependent distributions for the observation process $X_t$. That is, at each time step $(t=1,2,\dots)$, we have

$$Y_t \mid S_t = j \sim N(\mu_j,\sigma^2)$$

for $j \in \{1, 2\}$

The marginal likelihood of the model can be written with the forward algorithm using a transition matrix $\boldsymbol{\Gamma}$:

$$\boldsymbol{\Gamma}(x_t) =
\begin{pmatrix}
\phi(x_t \vert \mu_1, \sigma^2) & 0 \\
0 & \phi(x_t \vert \mu_2, \sigma^2) \\
\end{pmatrix}$$

where $\phi$ is the Gaussian PDF.

We are interested in knowing $\boldsymbol{\Gamma}$ and also $\boldsymbol{\mu}$!
"""

# ╔═╡ 6fd49295-d0e3-4b54-aeae-e9cd07a5281c
md"""
#### Random Data
"""

# ╔═╡ 58c5460f-c7f4-4a0a-9e18-71b9580e9148
begin
	const T = 2 # Number of States
	
	# Transition Probabilities
	const Γ = Matrix([0.9 0.1; 0.1 0.9])
	# initial distribution set to the stationary distribution
	const δ = (Diagonal(ones(T)) - Γ .+ 1) \ ones(T)
	# State-Dependent Gaussian means
	const μ = [1, 5]
	
	const n_obs = 1_000
	S = Vector{Int64}(undef, n_obs)
	y = Vector{Float64}(undef, n_obs)
	
	# initialise state and observation
	S[1] = sample(1:T, aweights(δ))
	y[1] = rand(Normal(μ[S[1]], 2))
	
	# simulate state and observation processes forward
	for t in 2:n_obs
	    S[t] = sample(1:T, aweights(Γ[S[t - 1], :]))
	    y[t] = rand(Normal(μ[S[t]], 2))
	end
end

# ╔═╡ 46ba21ab-bce5-4eed-bd63-aae7340c8180
begin
	# State-Dependent Gaussian means
	μ_sim = [μ₁_sim, μ₂_sim]
	
	S_sim = Vector{Int64}(undef, n_obs)
	y_sim = Vector{Float64}(undef, n_obs)
	
	# initialise state and observation
	S_sim[1] = sample(1:T, aweights(δ))
	y_sim[1] = rand(Normal(μ[S[1]], 2))
	
	# simulate state and observation processes forward
	for t in 2:n_obs
	    S_sim[t] = sample(1:T, aweights(Γ[S_sim[t - 1], :]))
	    y_sim[t] = rand(Normal(μ_sim[S_sim[t]], 2))
	end
	Plots.gr(dpi=300)
	scatter(y_sim, mc= S_sim, xlabel=L"t", ylabel=L"y", label=false, ylim=(-5,13), yticks=(vcat(0, μ_sim, 10), vcat("0", "μ₁", "μ₂", "10")))
	hline!([μ₁_sim,μ₂_sim], lw=4, label=false, c=:black, style=:dash)
end

# ╔═╡ 5d3d2abb-85e3-4371-926e-61ff236253f1
md"""
Here is the `Stan` code (I've simplified from Leos-Barajas & Michelot's original code) :

> Note that we are using the `log_sum_exp()` trick
"""

# ╔═╡ 247a02e5-8599-43fd-9ee5-32ba8b827477
md"""
```cpp
data {
  int<lower=1> K; // number of states
  int<lower=1> T; // length of data set
  real y[T]; // observations
}
parameters {
  positive_ordered[K] mu; // state-dependent parameters
  simplex[K] theta[K]; // N x N tpm
}
model{
  // priors
  mu ~ student_t(3, 0, 1);
  for (k in 1:K)
    theta[k] ~ dirichlet([0.5, 0.5]);

  // Compute the marginal probability over possible sequences
  vector[K] acc;
  vector[K] lp;

  // forward algorithm implementation
  for(k in 1:K) // first observation
    lp[k] = normal_lpdf(y[1] | mu[k], 2);
  for (t in 2:T) {     // looping over observations
      for (k in 1:K){   // looping over states
          acc[k] = log_sum_exp(log(theta[k]) + lp) +
            normal_lpdf(y[t] | mu[k], 2);
      }
      lp = acc;
    }
  target += log_sum_exp(lp);
}
```

Obs: `log_sum_exp(a, b) = log(exp(a) + exp(b))`
"""

# ╔═╡ 6db0245b-0461-4db0-9462-7a5f80f7d589
md"""
Here's how we would do in `Turing`

> Note the Composite MCMC Sampler 
"""

# ╔═╡ b5a79826-151e-416e-b0a2-1a58eec9196c
begin
	@model hmm(y, K::Int64; T=length(y)) = begin
		# state sequence in a Libtask.TArray
		s = tzeros(Int, T)

		# Transition Probability Matrix.
		θ = Vector{Vector}(undef, K)

		# Priors
		μ ~ filldist(truncated(TDist(3), 1, 6), K)
			
		for i = 1:K
			θ[i] ~ Dirichlet([0.5, 0.5])
		end
		
		# Positive Ordered
		if any(μ[i] > μ[i+1] for i in 1:(length(μ) - 1))
        	# update the joint log probability of the samples
        	# we set it to -Inf such that the samples are rejected
        	Turing.@addlogprob!(-Inf)
		end

		# first observation
		s[1] ~ Categorical(K)
		y[1] ~ Normal(μ[s[1]], 2)

		# looping over observations
		for i = 2:T
			s[i] ~ Categorical(vec(θ[s[i - 1]]))
			y[i] ~ Normal(μ[s[i]], 2)
		end
	end;

	composite_sampler = Gibbs(NUTS(10, 0.65, :μ, :θ),
					PG(1, :s));

	hmm_chain = sample(hmm(y, 2), composite_sampler, 20);
	summarystats(hmm_chain[:, 1:6, :]) #only μ and θ
end

# ╔═╡ cd410368-9022-4030-86a0-1d125e76bc62
md"""
> Obs: probably in the future we'll have better implementation for positive ordered constraints in `Turing`. It will reside in the [`Bijectors.jl`](https://github.com/TuringLang/Bijectors.jl) package. Actually check this [PR](https://github.com/TuringLang/Bijectors.jl/pull/186), it seems positive ordered is coming to `Turing`.
"""

# ╔═╡ 9b0b62cb-2c61-4d47-a6c7-09c0c1a75a24
md"""
### 7.2 ODEs in `Turing` (SIR Model)
"""

# ╔═╡ 9b020402-ea15-4f52-9fff-c70d275b97ac
Resource("https://github.com/storopoli/Turing-Workshop/blob/master/images/SIR.png?raw=true", :width => 400)

# ╔═╡ c81f4877-024f-4dc8-b7ce-e781ab6101f3
md"""
The Susceptible-Infected-Recovered (SIR) model splits the population in three time-dependent compartments: the susceptible, the infected (and infectious), and the recovered (and not infectious) compartments. When a susceptible individual comes into contact with an infectious individual, the former can become infected for some time, and then recover and become immune. The dynamics can be summarized in a system ODEs:
"""

# ╔═╡ f2272fd5-5132-4a6e-b2ff-136dc2fb2903
md"""
$$\begin{aligned}
 \frac{dS}{dt} &= -\beta  S \frac{I}{N}\\
 \frac{dI}{dt} &= \beta  S  \frac{I}{N} - \gamma  I \\
 \frac{dR}{dt} &= \gamma I
\end{aligned}$$

where

*  $S(t)$ is the number of people susceptible to becoming infected (no immunity),

*  $I(t)$ is the number of people currently infected (and infectious),

*  $R(t)$ is the number of recovered people (we assume they remain immune indefinitely),

*  $\beta$ is the constant rate of infectious contact between people,

*  $\gamma$ the constant recovery rate of infected individuals.

"""

# ╔═╡ 2d230fea-dcf2-41e6-a477-2a2334f56990
md"""
#### How to code and ODE in Julia?

It's very easy:

1. Use [`DifferentialEquations.jl`](https://diffeq.sciml.ai/)
2. Create a ODE function
3. Choose:
   * Initial Conditions: $u_0$
   * Parameters: $p$
   * Time Span: $t$
   * *Optional*: [Solver](https://diffeq.sciml.ai/stable/solvers/ode_solve/) or leave blank for auto

PS: If you like SIR models checkout [`epirecipes/sir-julia`](https://github.com/epirecipes/sir-julia)
"""

# ╔═╡ 44f9935f-c5a5-4f08-a94b-7f6ee70df358
function sir_ode!(du, u, p, t)
	    (S, I, R) = u
	    (β, γ) = p
	    N = S + I + R
	    infection = β * I / N * S
	    recovery = γ * I
	    @inbounds begin
	        du[1] = -infection
	        du[2] = infection - recovery
	        du[3] = recovery
	    end
	    nothing
	end;

# ╔═╡ 92e17d42-c6d1-4891-99a9-4a3be9e2decf
md"""
 $I_0$ = $(@bind I₀ Slider(1:1:20, default = 1, show_value=true))

 $\beta$ = $(@bind sim_β Slider(0.1:0.2:3, default = 1.9, show_value=true))

 $\gamma$ = $(@bind sim_γ Slider(0.1:0.1:1.5, default = 0.9, show_value=true))
"""

# ╔═╡ 39902541-5243-4fa9-896c-36db93d9fcea
begin
	u = [763, I₀, 0];
	p_sim = [sim_β, sim_γ];
	tspan_sim = (0.0, 15.0);
end

# ╔═╡ 646ab8dc-db5a-4eb8-a08b-217c2f6d86be
begin
	Plots.gr(dpi=300)
	problem = ODEProblem(sir_ode!, [763, I₀, 0], tspan_sim, p_sim)
	solution = solve(problem, Tsit5(), saveat=1.0)
	plot(solution, label=[L"S" L"I" L"R"], lw=3)
	xlabel!("days")
	ylabel!("N")
end

# ╔═╡ 5c017766-445d-4f4b-98f1-ae63e78ec34b
md"""
As an example, I will use [Grinsztajn, Semenova, Margossian & Riou. 2021)](https://arxiv.org/abs/2006.02985)'s.

It's a boarding school:

> Outbreak of **influenza A (H1N1)** in 1978 at a British boarding school. The data consists of the daily number of students in bed, spanning over a time interval of 14 days. There were **763 male students** who were mostly full boarders and 512 of them became ill.  The outbreak lasted from the 22nd of January to the 4th of February. It is reported that **one infected boy started the epidemic**, which spread rapidly in the relatively closed community of the boarding school.

The data are freely available in the R package `{outbreaks}`, maintained as part of the [R Epidemics Consortium](http://www.repidemicsconsortium.org).
"""

# ╔═╡ 0a76f019-4853-4ba3-9af8-9f33e1d4c956
begin
	# Boarding School SIR
	boarding_school = CSV.read(download("https://github.com/storopoli/Turing-Workshop/blob/master/data/influenza_england_1978_school.csv?raw=true"), DataFrame);
	cases = boarding_school.in_bed;
end

# ╔═╡ b0cc8694-b7ab-4d23-a208-055299840334
plot(boarding_school.date, cases, markershape=:o, dpi=300, xlab=L"t", ylab="cases", label=false, title="Boarding School H1N1 Outbreak")

# ╔═╡ 680f104e-80b4-443f-b4bc-532df758c162
md"""
Here's how we would do in `Turing`:

> Note the ODE system inside `@model`
"""

# ╔═╡ ddfc38fc-b47d-4ea5-847a-e9cbee3aa0a1
@model sir(cases, I₀) = begin
  # Calculate number of timepoints
  l = length(cases)
  N = 763
  S₀ = N - I₀
  R₀ = 0

  # Priors
  β ~ TruncatedNormal(2, 1,  1e-6, 10)     # using 10 instead of `Inf` because numerical issues arose
  γ ~ TruncatedNormal(0.4, 0.5,  1e-6, 10) # using 10 instead of `Inf` because numerical issues arose
  ϕ⁻ ~ truncated(Exponential(5), 1, 20)
  ϕ = 1.0 / ϕ⁻

  # ODE Stuff
  u = float.([S₀, I₀, R₀])
  p = [β, γ]
  tspan = (0.0, float(l))
  prob = ODEProblem(sir_ode!,
          u,
          tspan,
          p)
  sol = solve(prob,
              Tsit5(), # You can change the solver (similar to RK45)
              saveat=1.0)
  solᵢ = Array(sol)[2, 2:end] # Infected

  # Likelihood
  for i in 1:l
    solᵢ[i] = max(1e-6, solᵢ[i]) # numerical issues arose
    cases[i] ~ NegativeBinomial(solᵢ[i], ϕ)
  end
end;

# ╔═╡ ee2616ca-2602-4823-9cfb-123b958701c4
begin
	sir_chain = sample(sir(cases, 1), NUTS(1_000, 0.65), MCMCThreads(), 2_000, 2);
	summarystats(sir_chain[:, 1:2, :]) # only β and γ
end

# ╔═╡ 3f7c469a-c366-49dd-b09c-ae9b2b5db3fd
corner(sir_chain, dpi=300)

# ╔═╡ 7a62c034-3709-483a-a663-7fe5e09cb773
begin
	Plots.gr(dpi=300)
	plot(sir_chain[:, 1:2, :]) # only β and γ
end

# ╔═╡ 7f1fd9b4-517a-4fec-89bb-4d696dadbc3d
md"""
## 8.1 Computational Tricks
"""

# ╔═╡ 81e29fc7-b5d3-46d8-aeac-fb8e6dc11b16
md"""
### 8.1 Non-centered parametrization (Funnel of Death)
"""

# ╔═╡ 5291b260-9a68-4c8b-aff4-7797804ccc95
md"""
Sometimes our posterior has **crazy geometries** that makes our MCMC sampler (including NUTS and HMC) to have a hard time to sample from it.

This example is from Neal (2003) and is called Neal's Funnel (altough some call it Funnel of Death). It exemplifies the difficulties of sampling from some hierarchical models. Here I will show a 2-D example with $x$ and $y$:

$$p(y,x) = \text{Normal}(y \mid 0,3) \times
\text{normal}\left(x \mid 0,\exp\left(\frac{y}{2}\right)\right)$$
"""

# ╔═╡ 08dbe330-670d-48d5-b704-2421e687bff1
begin
	funnel_y = rand(Normal(0, 3), 10_000)
	funnel_x = rand(Normal(), 10_000) .* exp.(funnel_y / 2)
	Plots.gr(dpi=300)
	scatter((funnel_x, funnel_y),
	        label=false, ma=0.3,
	        xlabel=L"x", ylabel=L"y",
	        xlims=(-100, 100))
end

# ╔═╡ c109b759-7b73-4593-b9ea-8cc97b61d6fe
md"""
#### Whats the problem here?

* At the *bottom* of the funnel: **low** $\epsilon$ and **high** $L$
* At the *top* of the funnel: **high** $\epsilon$ and **low** $L$

HMC you have to set your $\epsilon$ and $L$ so it's fixed.

NUTS can automatically set $\epsilon$ and $L$ during warmup (it can vary) but it's fixed during sampling.

So basically you are screwed if you do not reparametrize!
"""

# ╔═╡ fe0fefb6-2755-4319-a944-bbbc7843aead
begin
	Plots.plotly(dpi=300)
	x = -2:0.01:2;
	kernel(x, y) = logpdf(Normal(0, exp(y / 2)), x)
	surface(x, x, kernel, xlab="x", ylab="y", zlab="log(PDF)")
end

# ╔═╡ 60494b7c-1a08-4846-8a80-12533552a697
md"""
#### Reparametrization
What if we reparameterize so that we can express $y$ and $x$ as standard normal distributions, by using a reparameterization trick:

$$\begin{aligned}
x^* &\sim \text{Normal}(0, 1)\\
x &= x^* \cdot \sigma_x + \mu_x
\end{aligned}$$

This also works for multivariate stuff
"""

# ╔═╡ b57195f9-c2a1-4676-96f9-faee84f7fc26
md"""
#### Non-Centered Reparametrization of the Funnel of Death

We can provide the MCMC sampler a better-behaved posterior geometry to explore:


$$\begin{aligned}
p(y^*,x^*) &= \text{Normal}(y^* \mid 0,0) \times
\text{Normal}(x^* \mid 0,0)\\
y &= 3y^* + 0\\
x &= \exp \left( \frac{y}{2} \right) x^* + 0
\end{aligned}$$

Below there is is the Neal's Funnel reparameterized as standard normal:
"""

# ╔═╡ 438d437e-7b00-4a13-8f8a-87fdc332a190
begin
	Plots.plotly(dpi=300)
	kernel_reparameterized(x, y) = logpdf(Normal(), x)
	surface(x, x,  kernel_reparameterized, xlab="x", ylab="y", zlab="log(PDF)")
end

# ╔═╡ 800fe4ba-e1a4-4e94-929f-7d66516e7bd6
md"""
#### Non-Centered Reparametrization of a Hierarchical Model
"""

# ╔═╡ 7d5a29c6-e71e-4ccb-a1c2-7fba663f038c
@model varying_intercept_ncp(X, idx, y; n_gr=length(unique(idx)), predictors=size(X, 2)) = begin
    # priors
    α ~ Normal(mean(y), 2.5 * std(y))       # population-level intercept
    β ~ filldist(Normal(0, 2), predictors)  # population-level coefficients
    σ ~ Exponential(1 / std(y))             # residual SD

	# prior for variance of random intercepts
    # usually requires thoughtful specification
    τ ~ truncated(Cauchy(0, 2), 0, Inf)    # group-level SDs intercepts
    zⱼ ~ filldist(Normal(0, 1), n_gr)      # NCP group-level intercepts

    # likelihood
    ŷ = α .+ X * β .+ zⱼ[idx] .* τ
    y ~ MvNormal(ŷ, σ)
end;

# ╔═╡ 3364e9f6-b2af-45f1-b1f3-ef7b0cd4910a
md"""
To reconstruct the original `αⱼ`s just multiply `zⱼ[idx] .* τ`:

```julia
τ = summarystats(chain_ncp)[:τ, :mean]
αⱼ = mapslices(x -> x * τ, chain_ncp[:,namesingroup(chain_ncp, :zⱼ),:].value.data, dims=[2])
chain_ncp_reconstructed = hcat(Chains(αⱼ, ["αⱼ[$i]" for i in 1:length(unique(idx))]), chain_ncp)
```

> [`mapslices()`](https://docs.julialang.org/en/v1/base/arrays/#Base.mapslices) is a `Base` Julia function that maps a function `f` to each `slice` (column) of an `Array`.
"""

# ╔═╡ 26265a91-2c8e-46d8-9a87-a2d097e7433a
md"""
### 8.2 $\mathbf{QR}$ decomposition
"""

# ╔═╡ 2eeb402e-c5f9-449c-af19-ff8f2e6c7246
md"""

Back in "Linear Algebra 101" we've learned that any matrix (even retangular ones) can be factored into the product of two matrices:

*  $\mathbf{Q}$: an orthogonal matrix (its columns are orthogonal unit vectors meaning $\mathbf{Q}^T = \mathbf{Q}^{-1})$.
*  $\mathbf{R}$: an upper triangular matrix.

This is commonly known as the [**QR Decomposition**](https://en.wikipedia.org/wiki/QR_decomposition):

$$\mathbf{A} = \mathbf{Q} \cdot \mathbf{R}$$

But what can we do with QR decomposition? It can speed up `Turing`'s sampling by a huge factor while also **decorrelating** the columns of $\mathbf{X}$, *i.e.* the independent variables.

The orthogonal nature of QR decomposition alters the posterior's topology and makes it easier for HMC or other MCMC samplers to explore it.

Now let's us incorporate QR decomposition in the logistic regression model.
Here, I will use the "thin" instead of the "fat" QR, which scales the $\mathbf{Q}$ and $\mathbf{R}$ matrices by a factor of $\sqrt{n-1}$ where $n$ is the number of rows of $\mathbf{X}$. In practice it is better implement the thin QR decomposition, which is to be preferred to the fat QR decomposition. It is numerically more stable. Mathematically, the thin QR decomposition is:

$$\begin{aligned}
x &= \mathbf{Q}^* \mathbf{R}^* \\
\mathbf{Q}^* &= \mathbf{Q} \cdot \sqrt{n - 1} \\
\mathbf{R}^* &= \frac{1}{\sqrt{n - 1}} \cdot \mathbf{R}\\
\boldsymbol{\mu}
&= \alpha + \mathbf{X} \cdot \boldsymbol{\beta} \\
&= \alpha + \mathbf{Q}^* \cdot \mathbf{R}^* \cdot \boldsymbol{\beta} \\
&= \alpha + \mathbf{Q}^* \cdot (\mathbf{R}^* \cdot \boldsymbol{\beta}) \\
&= \alpha + \mathbf{Q}^* \cdot \widetilde{\boldsymbol{\beta}} \\
\end{aligned}$$

Then we can recover original $\boldsymbol{\beta}$ with:

$$\boldsymbol{\beta} = \mathbf{R}^{*-1} \cdot \widetilde{\boldsymbol{\beta}}$$

Here's applied to our Logistic Regression example:

> Look at the `ess` in both examples
"""

# ╔═╡ 6870ca6d-256d-4a38-970e-1c26ceba9fa4
begin
	Q, R = qr(X_wells);
	Q_ast = Matrix(Q) * sqrt(size(X_wells, 1) - 1);
	R_ast = R / sqrt(size(X_wells, 1) - 1);
end

# ╔═╡ e5dac5c5-4644-443f-aa79-e43b399712c0
begin
	chain_log_reg = sample(logreg_vectorized(X_wells, y_wells), NUTS(1_000, 0.65), 2_000);
	summarystats(chain_log_reg)
end

# ╔═╡ 85f98ea6-9351-4527-8b8e-b2827a7735ff
begin
	chain_qr = sample(logreg_vectorized(Q_ast, y_wells), NUTS(1_000, 0.65), 2_000);
	summarystats(chain_qr)
end

# ╔═╡ 859ce60b-2f32-44d1-919a-dbdaf1be38fb
md"""
Now we have to reconstruct our $\boldsymbol{\beta}$s:

> [`mapslices()`](https://docs.julialang.org/en/v1/base/arrays/#Base.mapslices) is a `Base` Julia function that maps a function `f` to each `slice` (column) of an `Array`.
"""

# ╔═╡ 0377939c-00ac-42ae-b981-cdc897421588
begin
	betas = mapslices(x -> R_ast^-1 * x, chain_qr[:, namesingroup(chain_qr, :β),:].value.data, dims=[2]);
	
	chain_qr_reconstructed = hcat(Chains(betas, ["real_β[$i]" for i in 1:size(Q_ast, 2)]), chain_qr);
	
	summarystats(chain_qr_reconstructed)
end

# ╔═╡ 2f907e0d-171e-44c3-a531-5f11da08b3cf
md"""
## Pluto Stuff
"""

# ╔═╡ 31b6d4ec-d057-44ca-875b-0c3257895dd3
PlutoUI.TableOfContents(aside=true)

# ╔═╡ 98ece9fe-dfcc-4dd8-bd47-049217d2afcf
md"""
## References

Betancourt, M. (2020). Hierarchical Modeling. Available on: https://betanalpha.github.io/assets/case_studies/hierarchical_modeling.html

Damiano, L., Peterson, B., & Weylandt, M. (2017). A Tutorial on Hidden Markov Models using Stan. https://github.com/luisdamiano/stancon18 (Original work published 2017)

Ge, H., Xu, K., & Ghahramani, Z. (2018). Turing: A Language for Flexible Probabilistic Inference. International Conference on Artificial Intelligence and Statistics, 1682–1690. http://proceedings.mlr.press/v84/ge18b.html

Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). *Bayesian Data Analysis*. Chapman and Hall/CRC.

Gelman, A., Hill, J., & Vehtari, A. (2020a). *Regression and other stories*. Cambridge University Press.

Gelman, A., Vehtari, A., Simpson, D., Margossian, C. C., Carpenter, B., Yao, Y., Kennedy, L., Gabry, J., Bürkner, P.-C., & Modrák, M. (2020b). Bayesian Workflow. ArXiv:2011.01808 [Stat]. http://arxiv.org/abs/2011.01808

Grinsztajn, L., Semenova, E., Margossian, C. C., & Riou, J. (2021). Bayesian workflow for disease transmission modeling in Stan. ArXiv:2006.02985 [q-Bio, Stat]. http://arxiv.org/abs/2006.02985

Hoffman, M. D., & Gelman, A. (2011). The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo. Journal of Machine Learning Research, 15(1), 1593–1623.

Leos-Barajas, V., & Michelot, T. (2018). An Introduction to Animal Movement Modeling with Hidden Markov Models using Stan for Bayesian Inference. ArXiv:1806.10639 [q-Bio, Stat]. http://arxiv.org/abs/1806.10639

McElreath, R. (2020). *Statistical rethinking: A Bayesian course with examples in R and Stan*. CRC press.

McGrayne, S.B (2012). *The Theory That Would Not Die: How Bayes' Rule Cracked the Enigma Code, Hunted Down Russian Submarines, and Emerged Triumphant from Two Centuries of Controversy* Yale University Press.

Neal, R. M. (2003). Slice Sampling. The Annals of Statistics, 31(3), 705–741.

Tarek, M., Xu, K., Trapp, M., Ge, H., & Ghahramani, Z. (2020). DynamicPPL: Stan-like Speed for Dynamic Probabilistic Models. ArXiv:2002.02702 [Cs, Stat]. http://arxiv.org/abs/2002.02702
"""

# ╔═╡ e66e67e8-8ac2-41a3-9926-3f0ac3b9c47d
md"""
## License

This content is licensed under [Creative Commons Attribution-ShareAlike 4.0 Internacional](http://creativecommons.org/licenses/by-sa/4.0/).

[![CC BY-SA 4.0](https://licensebuttons.net/l/by-sa/4.0/88x31.png)](http://creativecommons.org/licenses/by-sa/4.0/)
"""

# ╔═╡ 634c9cc1-5a93-42b4-bf51-17dadfe488d6
md"""
## Environment
"""

# ╔═╡ 50e01181-1911-426b-9228-4663a1297619
with_terminal() do
	deps = [pair.second for pair in Pkg.dependencies()]
	deps = filter(p -> p.is_direct_dep, deps)
	deps = filter(p -> !isnothing(p.version), deps)
	list = ["$(p.name) $(p.version)" for p in deps]
	sort!(list)
	println(join(list, '\n'))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
InteractiveUtils = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"
Turing = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"

[compat]
BenchmarkTools = "~1.3.2"
CSV = "~0.10.7"
DataFrames = "~1.4.3"
DifferentialEquations = "~7.6.0"
Distributions = "~0.25.79"
LaTeXStrings = "~1.3.0"
LazyArrays = "~0.22.12"
Plots = "~1.36.4"
PlutoUI = "~0.7.48"
StatsBase = "~0.33.21"
StatsPlots = "~0.15.4"
Turing = "~0.22.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.3"
manifest_format = "2.0"
project_hash = "f25993387a85d360f04d72f5a2c883604b5e3032"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[deps.AbstractMCMC]]
deps = ["BangBang", "ConsoleProgressMonitor", "Distributed", "Logging", "LoggingExtras", "ProgressLogging", "Random", "StatsBase", "TerminalLoggers", "Transducers"]
git-tree-sha1 = "5c26c7759412ffcaf0dd6e3172e55d783dd7610b"
uuid = "80f14c24-f653-4e6a-9b94-39d6b0f70001"
version = "4.1.3"

[[deps.AbstractPPL]]
deps = ["AbstractMCMC", "DensityInterface", "Setfield", "SparseArrays"]
git-tree-sha1 = "6320752437e9fbf49639a410017d862ad64415a5"
uuid = "7a57a42e-76ec-4ea3-a279-07e840d6d9cf"
version = "0.5.2"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "52b3b436f8f73133d7bc3a6c71ee7ed6ab2ab754"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.3"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.AdvancedHMC]]
deps = ["AbstractMCMC", "ArgCheck", "DocStringExtensions", "InplaceOps", "LinearAlgebra", "ProgressMeter", "Random", "Requires", "Setfield", "Statistics", "StatsBase", "StatsFuns", "UnPack"]
git-tree-sha1 = "0091e2e4d0a7125da0e3ad8c7dbff9171a921461"
uuid = "0bf59076-c3b1-5ca4-86bd-e02cd72cde3d"
version = "0.3.6"

[[deps.AdvancedMH]]
deps = ["AbstractMCMC", "Distributions", "Random", "Requires"]
git-tree-sha1 = "d7a7dabeaef34e5106cdf6c2ac956e9e3f97f666"
uuid = "5b7e9947-ddc0-4b3f-9b55-0d8042f74170"
version = "0.6.8"

[[deps.AdvancedPS]]
deps = ["AbstractMCMC", "Distributions", "Libtask", "Random", "StatsFuns"]
git-tree-sha1 = "9ff1247be1e2aa2e740e84e8c18652bd9d55df22"
uuid = "576499cb-2369-40b2-a588-c64705576edc"
version = "0.3.8"

[[deps.AdvancedVI]]
deps = ["Bijectors", "Distributions", "DistributionsAD", "DocStringExtensions", "ForwardDiff", "LinearAlgebra", "ProgressMeter", "Random", "Requires", "StatsBase", "StatsFuns", "Tracker"]
git-tree-sha1 = "67fcc7d46c26250e89fc62798fbe07b5ee264c6f"
uuid = "b5ca4192-6429-45e5-a2d9-87aec30a685c"
version = "0.1.6"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.ArrayInterface]]
deps = ["ArrayInterfaceCore", "Compat", "IfElse", "LinearAlgebra", "Static"]
git-tree-sha1 = "6d0918cb9c0d3db7fe56bea2bc8638fc4014ac35"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "6.0.24"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "c46fb7dd1d8ca1d213ba25848a5ec4e47a1a1b08"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.26"

[[deps.ArrayInterfaceGPUArrays]]
deps = ["Adapt", "ArrayInterfaceCore", "GPUArraysCore", "LinearAlgebra"]
git-tree-sha1 = "fc114f550b93d4c79632c2ada2924635aabfa5ed"
uuid = "6ba088a2-8465-4c0a-af30-387133b534db"
version = "0.2.2"

[[deps.ArrayInterfaceOffsetArrays]]
deps = ["ArrayInterface", "OffsetArrays", "Static"]
git-tree-sha1 = "3d1a9a01976971063b3930d1aed1d9c4af0817f8"
uuid = "015c0d05-e682-4f19-8f0a-679ce4c54826"
version = "0.1.7"

[[deps.ArrayInterfaceStaticArrays]]
deps = ["Adapt", "ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "LinearAlgebra", "Static", "StaticArrays"]
git-tree-sha1 = "f12dc65aef03d0a49650b20b2fdaf184928fd886"
uuid = "b0d46f97-bff5-4637-a19a-dd75974142cd"
version = "0.1.5"

[[deps.ArrayInterfaceStaticArraysCore]]
deps = ["Adapt", "ArrayInterfaceCore", "LinearAlgebra", "StaticArraysCore"]
git-tree-sha1 = "93c8ba53d8d26e124a5a8d4ec914c3a16e6a0970"
uuid = "dd5226c6-a4d4-4bc7-8575-46859f9c95b9"
version = "0.1.3"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "7c3a0ce8d64ad84438b1bf60eca4f003a84e33e1"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.8.14"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "1dd4d9f5beebac0c03446918741b1a03dc5e5788"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.6"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "45fc0722d59676040e356e7c059d6e51c2b24382"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.17.8"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "7fe6d92c4f281cf4ca6f2fba0ce7b299742da7ca"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.37"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[deps.Bijectors]]
deps = ["ArgCheck", "ChainRulesCore", "ChangesOfVariables", "Compat", "Distributions", "Functors", "InverseFunctions", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "MappedArrays", "Random", "Reexport", "Requires", "Roots", "SparseArrays", "Statistics"]
git-tree-sha1 = "a3704b8e5170f9339dff4e6cb286ad49464d3646"
uuid = "76274a88-744f-5084-9051-94815aaf08c4"
version = "0.10.6"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.BoundaryValueDiffEq]]
deps = ["BandedMatrices", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase", "SparseArrays"]
git-tree-sha1 = "d4d28a2ac12b6252e6b09a6eebf71f6f8795de49"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "2.10.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "Static"]
git-tree-sha1 = "a7157ab6bcda173f533db4c93fc8a27a48843757"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.30"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "c5fd7cd27ac4aed0acf4b73948f0110ff2a854b2"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.7"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRules]]
deps = ["Adapt", "ChainRulesCore", "Compat", "Distributed", "GPUArraysCore", "IrrationalConstants", "LinearAlgebra", "Random", "RealDot", "SparseArrays", "Statistics", "StructArrays"]
git-tree-sha1 = "0c8c8887763f42583e1206ee35413a43c91e2623"
uuid = "082447d4-558c-5d27-93f4-14fc19e9eca2"
version = "1.45.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "d61300b9895f129f4bd684b2aff97cf319b6c493"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.11"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "64df3da1d2a26f4de23871cd1b6482bb68092bd5"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.14.3"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "aaabba4ce1b7f8a9b34c015053d3b1edf60fa49c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.4.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.ConsoleProgressMonitor]]
deps = ["Logging", "ProgressMeter"]
git-tree-sha1 = "3ab7b2136722890b9af903859afcf457fa3059e8"
uuid = "88cd18e8-d9cc-4ea6-8889-5259c0d15c8b"
version = "0.1.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "e08915633fcb3ea83bf9d6126292e5bc5c739922"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.13.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "0f44494fe4271cc966ac4fea524111bef63ba86c"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.4.3"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "UnPack"]
git-tree-sha1 = "90f6f571a5aa9e4f79b9c0f246cdbf507506f5a2"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.40.1"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterfaceCore", "ChainRulesCore", "DataStructures", "Distributions", "DocStringExtensions", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "MuladdMacro", "Parameters", "Printf", "RecursiveArrayTools", "Reexport", "Requires", "SciMLBase", "Setfield", "SimpleNonlinearSolve", "SparseArrays", "Static", "StaticArrays", "Statistics", "Tricks", "ZygoteRules"]
git-tree-sha1 = "65805bb205e8d011fc91da87d41d14394db5d791"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.108.0"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "Markdown", "NLsolve", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "485503846a90b59f3b79b39c2d818496bf50d197"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.24.3"

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "ResettableStacks", "SciMLBase", "StaticArrays", "Statistics"]
git-tree-sha1 = "b520f8b55f41c3a86fd181d342f23cccd047846e"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.14.2"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "c5b6685d53f933c11404a3ae9822afe30d522494"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.12.2"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "JumpProcesses", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "418dad2cfd7377a474326bc86a23cb645fac6527"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.6.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "a7756d098cbabec6b3ac44f369f74915e8cfd70a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.79"

[[deps.DistributionsAD]]
deps = ["Adapt", "ChainRules", "ChainRulesCore", "Compat", "DiffRules", "Distributions", "FillArrays", "LinearAlgebra", "NaNMath", "PDMats", "Random", "Requires", "SpecialFunctions", "StaticArrays", "StatsBase", "StatsFuns", "ZygoteRules"]
git-tree-sha1 = "0c139e48a8cea06c6ecbbec19d3ebc5dcbd7870d"
uuid = "ced4e74d-a319-5a8a-b0ac-84af2272839c"
version = "0.6.43"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "c36550cb29cbe373e95b3f40486b9a4148f89ffd"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.2"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPPL]]
deps = ["AbstractMCMC", "AbstractPPL", "BangBang", "Bijectors", "ChainRulesCore", "ConstructionBase", "Distributions", "DocStringExtensions", "LinearAlgebra", "MacroTools", "OrderedCollections", "Random", "Setfield", "Test", "ZygoteRules"]
git-tree-sha1 = "9a795bb2fe860e2b0a19136429a173de9f8c2774"
uuid = "366bfd00-2699-11ea-058f-f148b4cae6d8"
version = "0.21.3"

[[deps.EllipticalSliceSampling]]
deps = ["AbstractMCMC", "ArrayInterfaceCore", "Distributions", "Random", "Statistics"]
git-tree-sha1 = "4cda4527e990c0cc201286e0a0bfbbce00abcfc2"
uuid = "cad2338a-1db2-11e9-3401-43bc07c9ede2"
version = "1.0.0"

[[deps.EnumX]]
git-tree-sha1 = "e5333cd1e1c713ee21d07b6ed8b0d8853fabe650"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.3"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceGPUArrays", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "9837d3f3a904c7a7ab9337759c0093d3abea1d81"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.22.0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "LinearAlgebra", "Polyester", "Static", "StrideArraysCore"]
git-tree-sha1 = "24db26ecc4c8a00584672d3b4c6cb0bb3dad9d51"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.2.3"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "14a6f7a21125f715d935fe8f83560ee833f7d79d"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "1.2.7"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "802bfc139833d2ba893dd9e62ba1767c88d708ae"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.5"

[[deps.FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "04ed1f0029b6b3af88343e439b995141cb0d0b8d"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.17.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "10fa12fe96e4d76acfa738f4df2126589a67374f"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.33"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "a5e6e7f12607e90d71b09e6ce2c965e41b337968"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.1"

[[deps.Functors]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "a2657dd0f3e8a61dbe70fc7c122038bd33790af5"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.3.0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "6872f5ec8fd1a38880f027a26739d42dcda6691f"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.2"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "051072ff2accc6e0e87b708ddee39b18aa04a0bc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.71.1"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "501a4bf76fd679e7fcd678725d5072177392e756"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.71.1+0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "fb69b2a645fa69ba5f474af09221b9308b160ce6"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.3"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "fb83fbe02fe57f2c068013aa94bcdf6760d3a7a7"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "ba2d094a88b6b287bd25cfa86f301e7693ffae2f"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.7.4"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "e1acc37ed078d99a714ed8376446f92a5535ca65"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.5.5"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "f64b890b2efa4de81520d2b0fbdc9aadb65bdf53"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.13"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "0cf92ec945125946352f3d46c96976ab972bde6f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.3.2"

[[deps.InplaceOps]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "50b41d59e7164ab6fda65e71049fee9d890731ff"
uuid = "505f98c9-085e-5b2c-8e89-488be7bf1f34"
version = "0.3.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "842dd89a6cb75e02e85fdd75c760cdc43f5d6863"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.6"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "3f91cd3f56ea48d4d2a75c2a65455c5fc74fa347"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.JumpProcesses]]
deps = ["ArrayInterfaceCore", "DataStructures", "DiffEqBase", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "Markdown", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "2cc453fad790410a40a6efe38e46c9c5a9c6fa41"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.2.3"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "764164ed65c30738750965d55652db9c94c59bfe"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.4.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "9816b296736292a80b9a3200eb7fbb57aaa3917a"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.5"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "92256444f81fb094ff5aa742ed10835a621aef75"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.8.4"

[[deps.KrylovKit]]
deps = ["ChainRulesCore", "GPUArraysCore", "LinearAlgebra", "Printf"]
git-tree-sha1 = "1a5e1d9941c783b0119897d29f2eb665d876ecf3"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.6.0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LRUCache]]
git-tree-sha1 = "d862633ef6097461037a00a13f709a62ae4bdfdd"
uuid = "8ac3fa9e-de4c-5943-b1dc-09c6b5f20637"
version = "1.4.0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "7e34177793212f6d64d045ee47d2883f09fffacc"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.12"

[[deps.LazyArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "MacroTools", "MatrixFactorizations", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8d9837955df02124bc74f0deb38adbd3a1437846"
uuid = "5078a376-72f3-5289-bfd5-ec5146d43c02"
version = "0.22.12"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LeftChildRightSiblingTrees]]
deps = ["AbstractTrees"]
git-tree-sha1 = "fb6803dafae4a5d62ea5cab204b1e657d9737e7f"
uuid = "1d6d02ad-be62-4b6b-8a6d-2f90e265016e"
version = "0.2.0"

[[deps.LevyArea]]
deps = ["LinearAlgebra", "Random", "SpecialFunctions"]
git-tree-sha1 = "56513a09b8e0ae6485f34401ea9e2f31357958ec"
uuid = "2d8b4e74-eb68-11e8-0fb9-d5eb67b50637"
version = "1.0.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtask]]
deps = ["FunctionWrappers", "LRUCache", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "dfa6c5f2d5a8918dd97c7f1a9ea0de68c2365426"
uuid = "6f1fad26-d15e-5dc8-ae53-837a1d7b8c9f"
version = "0.7.5"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterfaceCore", "DocStringExtensions", "FastLapackInterface", "GPUArraysCore", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "SnoopPrecompile", "SparseArrays", "SuiteSparse", "UnPack"]
git-tree-sha1 = "5b179a11483a789650a7d1ce7131aa0c416991db"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.29.1"

[[deps.LogDensityProblems]]
deps = ["ArgCheck", "DocStringExtensions", "Random", "Requires", "UnPack"]
git-tree-sha1 = "c3e1189191e4528b605070972d7d4e9cd91dd96b"
uuid = "6fdf6af0-433a-55f7-b3ed-c6c6e0b8df7c"
version = "1.0.3"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "946607f84feb96220f480e0422d3484c49c00239"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.19"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDDualNumbers", "SIMDTypes", "SLEEFPirates", "SnoopPrecompile", "SpecialFunctions", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "da5317a78e2a9f692730345cf3ff820109f406d3"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.141"

[[deps.MCMCChains]]
deps = ["AbstractMCMC", "AxisArrays", "Compat", "Dates", "Distributions", "Formatting", "IteratorInterfaceExtensions", "KernelDensity", "LinearAlgebra", "MCMCDiagnosticTools", "MLJModelInterface", "NaturalSort", "OrderedCollections", "PrettyTables", "Random", "RecipesBase", "Serialization", "Statistics", "StatsBase", "StatsFuns", "TableTraits", "Tables"]
git-tree-sha1 = "f5f347b828fd95ece7398f412c81569789361697"
uuid = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
version = "5.5.0"

[[deps.MCMCDiagnosticTools]]
deps = ["AbstractFFTs", "DataAPI", "Distributions", "LinearAlgebra", "MLJModelInterface", "Random", "SpecialFunctions", "Statistics", "StatsBase", "Tables"]
git-tree-sha1 = "dab8e9e9cd714fd3adc2344a97504e7e64e78546"
uuid = "be115224-59cd-429b-ad48-344e309966f0"
version = "0.1.5"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MLJModelInterface]]
deps = ["Random", "ScientificTypesBase", "StatisticalTraits"]
git-tree-sha1 = "c8b7e632d6754a5e36c0d94a4b466a5ba3a30128"
uuid = "e80e1ace-859a-464e-9ed9-23947d8ae3ea"
version = "1.8.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "4a02641a5b58e09bd04123ecb6eda64bbaee15b0"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "0.9.3"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "4d5917a26ca33c66c8e5ca3247bd163624d35493"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "efe9c8ecab7a6311d4b91568bd6c88897822fabe"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.0"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NNlib]]
deps = ["Adapt", "ChainRulesCore", "LinearAlgebra", "Pkg", "Requires", "Statistics"]
git-tree-sha1 = "00bcfcea7b2063807fdcab2e0ce86ef00b8b8000"
uuid = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
version = "0.8.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NamedArrays]]
deps = ["Combinatorics", "DataStructures", "DelimitedFiles", "InvertedIndices", "LinearAlgebra", "Random", "Requires", "SparseArrays", "Statistics"]
git-tree-sha1 = "2fd5787125d1a93fbe30961bd841707b8a80d75b"
uuid = "86f7a689-2022-50b4-a561-43c23ac3c673"
version = "0.9.6"

[[deps.NaturalSort]]
git-tree-sha1 = "eda490d06b9f7c00752ee81cfa451efe55521e21"
uuid = "c020b1a1-e9b0-503a-9c33-f039bfc54a85"
version = "1.0.0"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "440165bf08bc500b8fe4a7be2dc83271a00c0716"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "f71d8950b724e9ff6110fc948dff5a329f901d64"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "df6830e37943c7aaa10023471ca47fb3065cc3c4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "1903afc76b7d01719d9c30d3c7d501b61db96721"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.4"

[[deps.Optimisers]]
deps = ["ChainRulesCore", "Functors", "LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "c8884e88057968a715ddfd39f5c18b4c03aa3b64"
uuid = "3bd65402-5787-11e9-1adc-39752487f4e2"
version = "0.2.11"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceGPUArrays", "ArrayInterfaceStaticArrays", "ArrayInterfaceStaticArraysCore", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "Polyester", "PreallocationTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SnoopPrecompile", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "c3177fec02e122c412558eb190dc2cae05e2736b"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.33.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "b64719e8b4504983c7fca6cc9db3ebc8acc2a4d6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "21303256d239f6b484977314674aef4bb1fe4420"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "c411cd689fbcdd2c01b79db2817114e4321c6fa2"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.36.4"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "efc140104e6d0ae3e7e30d56c98c4a927154d684"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.48"

[[deps.PoissonRandom]]
deps = ["Random"]
git-tree-sha1 = "45f9da1ceee5078267eb273d065e8aa2f2515790"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.3"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "7446b311ce8f2c3f48be75a2c4b53ff02dd0c1df"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.6.18"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "050ca4aa2ca31484b51b849d8180caf8e4449c49"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.1.11"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ForwardDiff"]
git-tree-sha1 = "3953d18698157e1d27a51678c89c88d53e071a42"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "d8ed354439950b34ab04ff8f3dfd49e11bc6c94b"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "97aa253e65b784fd13e83774cadc95b38011d734"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.6.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "7a1a306b72cfa60634f03a911405f4e64d1b718b"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.6.0"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "d12e612bba40d189cead6ff857ddb67bd2e6a387"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "e974477be88cb5e3040009f3767611bc6357846f"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.11"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "DocStringExtensions", "FillArrays", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "Tables", "ZygoteRules"]
git-tree-sha1 = "4c7a6462350942da60ea5749afd7cea58017301f"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.32.2"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "SnoopPrecompile", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "0a2dfb3358fcde3676beb75405e782faa8c9aded"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.Roots]]
deps = ["ChainRulesCore", "CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "a3db467ce768343235032a1ca0830fc64158dadf"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.8"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "50314d2ef65fce648975a8e80ae6d8409ebbf835"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.5"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDDualNumbers]]
deps = ["ForwardDiff", "IfElse", "SLEEFPirates", "VectorizationBase"]
git-tree-sha1 = "dd4195d308df24f33fb10dde7c22103ba88887fa"
uuid = "3cdde19b-5bb0-4aaf-8931-af3e248e098b"
version = "0.1.1"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "c8679919df2d3c71f74451321f1efea6433536cc"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.37"

[[deps.SciMLBase]]
deps = ["ArrayInterfaceCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Preferences", "RecipesBase", "RecursiveArrayTools", "RuntimeGeneratedFunctions", "StaticArraysCore", "Statistics", "Tables"]
git-tree-sha1 = "8ee3f9c8f20744e336fd89617de9835cd4a3f6e4"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.74.1"

[[deps.ScientificTypesBase]]
git-tree-sha1 = "a8e18eb383b5ecf1b5e6fc237eb39255044fd92b"
uuid = "30f210dd-8aff-4c5f-94ba-8e64358c1161"
version = "3.0.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "efd23b378ea5f2db53a55ae53d3133de4e080aa9"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.16"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ArrayInterfaceCore", "FiniteDiff", "ForwardDiff", "Reexport", "SciMLBase", "SnoopPrecompile", "StaticArraysCore"]
git-tree-sha1 = "7f3f0086731cb2eeb50b772389231dbc0b630ddb"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "0.1.2"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "472216c5af9f2f1fce02b760651fe024c75187bd"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.29.0"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "0559586098f3cbd2e835461254ea2fcffa4a61ba"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "4e051b85454b4e4f66e6a6b7bdc452ad9da3dcf6"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.10"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.StatisticalTraits]]
deps = ["ScientificTypesBase"]
git-tree-sha1 = "30b9236691858e13f167ce829490a68e1a597782"
uuid = "64bff920-2084-43da-a3e6-9bb72801c0c9"
version = "3.2.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5783b877201a82fc0014cbf381e7e6eb130473a4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.0.1"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "DataValues", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "e0d5bc26226ab1b7648278169858adcfbd861780"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.4"

[[deps.SteadyStateDiffEq]]
deps = ["DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "d06cb5c8f02490c85d46ec0ed303fd388941bd2e"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "1.9.1"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqNoiseProcess", "DocStringExtensions", "FillArrays", "FiniteDiff", "ForwardDiff", "JumpProcesses", "LevyArea", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SciMLBase", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "df694865e88228da8308ab618cace615305ce9e4"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.57.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "0b1ea3e3fdf93b42e9c0f58347168618b6acc259"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.3.17"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArraysCore", "Tables"]
git-tree-sha1 = "13237798b407150a6d2e2bce5d793d7d9576e99e"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.13"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+0"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "Reexport", "SciMLBase", "SnoopPrecompile", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "ec007bb878acd31cc25d74d139dff5a9693127a3"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.11.1"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg", "SuiteSparse_jll"]
git-tree-sha1 = "04777432d74ec5bc91ca047c9e0e0fd7f81acdb6"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.1+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TerminalLoggers]]
deps = ["LeftChildRightSiblingTrees", "Logging", "Markdown", "Printf", "ProgressLogging", "UUIDs"]
git-tree-sha1 = "f53e34e784ae771eb9ccde4d72e578aa453d0554"
uuid = "5d786b92-1e48-4d6f-9151-6b4477ca9bed"
version = "0.1.6"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "f8629df51cab659d70d2e5618a430b4d3f37f2c3"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.0"

[[deps.Tracker]]
deps = ["Adapt", "DiffRules", "ForwardDiff", "Functors", "LinearAlgebra", "LogExpFunctions", "MacroTools", "NNlib", "NaNMath", "Optimisers", "Printf", "Random", "Requires", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "d963aad627fd7af56fbbfee67703c2f7bfee9dd7"
uuid = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
version = "0.2.22"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "8a75929dcd3c38611db2f8d08546decb514fcadf"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.9"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "77fea79baa5b22aeda896a8d9c6445a74500a2c2"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.74"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "SnoopPrecompile", "Static", "VectorizationBase"]
git-tree-sha1 = "9126e77b5d1afee9f94e8a66165e3716603d6054"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.15"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.Turing]]
deps = ["AbstractMCMC", "AdvancedHMC", "AdvancedMH", "AdvancedPS", "AdvancedVI", "BangBang", "Bijectors", "DataStructures", "Distributions", "DistributionsAD", "DocStringExtensions", "DynamicPPL", "EllipticalSliceSampling", "ForwardDiff", "Libtask", "LinearAlgebra", "LogDensityProblems", "MCMCChains", "NamedArrays", "Printf", "Random", "Reexport", "Requires", "SciMLBase", "Setfield", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Tracker"]
git-tree-sha1 = "8a40377bcc4b054ebdc8f680e96cd73a4a6fe2e6"
uuid = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"
version = "0.22.0"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "fc79d0f926592ecaeaee164f6a4ca81b51115c3b"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.56"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╠═4af78efd-d484-4241-9d3c-97cc78e1dbd4
# ╟─5df4d7d2-c622-11eb-3bbd-bff9668ee5e0
# ╟─19c63110-4baa-4aff-ab91-46e5c149f3a2
# ╟─dceb8312-230f-4e4b-9285-4e23f219b838
# ╟─cda7dc96-d983-4e31-9298-6148205b54b1
# ╟─2164bf58-75ff-470c-828c-b0165f0d980d
# ╟─55777e4a-b197-4a61-8e57-6ae9792c0564
# ╟─1436305e-37d8-44f1-88d6-4de838580360
# ╟─08f508c4-233a-4bba-b313-b04c1d6c4a4c
# ╟─868d8932-b108-41d9-b4e8-d62d31b5465d
# ╟─653ec420-8de5-407e-91a9-f045e25a6395
# ╟─716cea7d-d771-46e9-ad81-687292004009
# ╟─cb808fd4-6eb2-457e-afa4-58ae1be09aec
# ╟─0484ae7f-bd8a-4615-a760-5c4b2eef9d3f
# ╟─1d467044-bc7d-4df7-bda6-bb8ea6ff0712
# ╠═b1d99482-53f5-4c6b-8c20-c761ff6bdb77
# ╠═65fa382d-4ef7-432d-8630-27082977185b
# ╠═06f93734-2315-4b36-a39a-09e8167bab1f
# ╟─9f6b96a7-033d-4c7d-a853-46a0b5af4675
# ╟─b7667fb4-6e76-4711-b61d-dae5f993531e
# ╟─cb168dc1-70e2-450f-b2cf-c8680251ab27
# ╟─07d408cf-d202-40b2-90c2-5e8630549339
# ╠═744a8a63-647f-4550-adf7-44354fde44be
# ╟─e6365296-cd68-430e-99c5-fb571f39aad5
# ╟─927ad0a4-ba68-45a6-9bde-561915503e48
# ╠═ab6c2ba6-4cd8-473a-88c6-b8d61551fb22
# ╟─2ab3c34a-1cfc-4d20-becc-5902d08d03e0
# ╟─924fcad9-75c1-4707-90ef-3e36947d64fe
# ╟─fc8e40c3-34a1-4b2e-bd1b-893d7998d359
# ╟─fb366eb1-4ab0-4e7a-83ed-d531978c06a0
# ╠═0fe83f55-a379-49ea-ab23-9defaab05890
# ╟─3aa95b4b-aaf8-45cf-8bc5-05b65b4bcccf
# ╠═dd27ee5f-e442-42d7-a39b-d76328d2e59f
# ╟─c4808b43-bc0f-4254-abf1-1adc19135dc7
# ╠═1773d8c3-4651-4128-9442-e7c858bc4a43
# ╟─5674f7aa-3205-47c7-8367-244c6419ce69
# ╟─83cc80c1-d97e-4b82-872e-e5493d2b62ab
# ╠═475be60f-1876-4086-9725-3bf5f52a3e43
# ╠═f6bc0cfd-a1d9-48e5-833c-f33bf1b89d45
# ╠═ed640696-cae6-47e1-a4df-0655192e0855
# ╠═bc9fa101-8854-4af5-904a-f0b683fb63b1
# ╟─c82687d1-89d0-4ecd-bed7-1708ba8b2662
# ╠═270c0b90-cce1-4092-9e29-5f9deda2cb7d
# ╠═c4146b8b-9d11-446e-9765-8d5283a6d445
# ╟─3d09c8c3-ce95-4f26-9136-fedd601e2a70
# ╠═8d9bdae2-658d-45bf-9b25-50b6efbe0cdf
# ╟─41b014c2-7b49-4d03-8741-51c91b95f64c
# ╠═2f08c6e4-fa7c-471c-ad9f-9d036e3027d5
# ╟─5f639d2d-bb96-4a33-a78e-d5b9f0e8d274
# ╠═3f7c469a-c366-49dd-b09c-ae9b2b5db3fd
# ╟─c70ebb70-bd96-44a5-85e9-871b0e478b1a
# ╟─36258bdd-f617-48f6-91c9-e8bbff78ebd8
# ╟─6630eb47-77f6-48e9-aafe-55bda275449c
# ╠═37e751c7-8b6c-47d9-8013-97015d1e1fb2
# ╟─7a21e7a0-322b-4f8e-9d8b-a2f452f7e092
# ╟─f8f59ebb-bb1e-401f-97b5-507634badb3f
# ╠═15795f79-7d7b-43d2-a4b4-99ad968a7f72
# ╟─dd5fbb2a-4220-4e47-945a-6870b799c50d
# ╟─0cc8e12c-9b72-41ec-9c13-d9ae0bdc6100
# ╠═fce0f511-3b00-4079-85c6-9b2d2d7c04cb
# ╟─5ba6b247-8277-4100-abe7-8d06af04a011
# ╠═0f000fc4-1a7b-4522-8355-8df572ee8800
# ╠═8a87e324-f3d9-4162-88ab-3833a6d1fc2e
# ╟─3c954cbc-aed7-4d22-b578-a80ce62ebb49
# ╟─521e2473-1aba-43be-951a-25537062891e
# ╟─bafc91d2-8cae-4af8-b5ed-8199eef40c4d
# ╟─a2292bc1-3379-450d-beb5-ae8f41b69be8
# ╟─38055b57-f983-4440-bef5-0ab6d180ff1e
# ╟─7d4d06ca-f96d-4b1e-860f-d9e0d6eb6723
# ╟─c64d355f-f5a2-46a5-86f3-2d02da98f305
# ╟─262cb245-0bc1-4a36-b0bc-de52c08ccde0
# ╟─3ecc92b8-6a10-4f51-93d7-72449e248dc2
# ╟─a0c7ca50-3a3f-483c-ae01-fd774e0c072d
# ╟─cb3dd785-11ff-42fe-ab85-0dd03e45209e
# ╟─4812f80e-79a9-4519-9e4d-a45127ca6a49
# ╟─318697fe-1fbc-4ac3-a2aa-5ecf775072d4
# ╠═9acc7a1c-f638-4a2e-ad67-c16cff125c86
# ╟─885fbe97-edd6-44d2-808d-8eeb1e9cb2b4
# ╠═7f526d1f-bd56-4e51-9f7b-ce6b5a2a1853
# ╟─f7971da6-ead8-4679-b8cf-e3c35c93e6cf
# ╠═546726af-5420-4a4f-8c0c-fe96a2ba43bc
# ╟─9ebac6ba-d213-4ed8-a1d5-66b841fafa00
# ╟─45c342fd-b893-46aa-b2ee-7c93e7a1d207
# ╟─d44c7baa-80d2-4fdb-a2de-35806477dd58
# ╟─c1b2d007-1004-42f5-b65c-b4e2e7ff7d8e
# ╟─c1dcfd47-9e25-470b-a1b3-ab66bfac59d6
# ╟─46ba21ab-bce5-4eed-bd63-aae7340c8180
# ╟─f1153918-0748-4400-ae8b-3b59f8c5d755
# ╟─ad6c4533-cd56-4f6f-b10d-d7bc3145ba16
# ╟─2ef397a6-f7fb-4fc2-b918-40ab545ce19f
# ╟─8d347172-2d26-4d9e-954d-b8924ed4c9e2
# ╟─ca962c0e-4620-4888-b7c3-aa7f6d7899e9
# ╟─6fd49295-d0e3-4b54-aeae-e9cd07a5281c
# ╠═58c5460f-c7f4-4a0a-9e18-71b9580e9148
# ╟─5d3d2abb-85e3-4371-926e-61ff236253f1
# ╟─247a02e5-8599-43fd-9ee5-32ba8b827477
# ╟─6db0245b-0461-4db0-9462-7a5f80f7d589
# ╠═b5a79826-151e-416e-b0a2-1a58eec9196c
# ╟─cd410368-9022-4030-86a0-1d125e76bc62
# ╟─9b0b62cb-2c61-4d47-a6c7-09c0c1a75a24
# ╟─9b020402-ea15-4f52-9fff-c70d275b97ac
# ╟─c81f4877-024f-4dc8-b7ce-e781ab6101f3
# ╟─f2272fd5-5132-4a6e-b2ff-136dc2fb2903
# ╟─2d230fea-dcf2-41e6-a477-2a2334f56990
# ╠═44f9935f-c5a5-4f08-a94b-7f6ee70df358
# ╟─39902541-5243-4fa9-896c-36db93d9fcea
# ╟─92e17d42-c6d1-4891-99a9-4a3be9e2decf
# ╟─646ab8dc-db5a-4eb8-a08b-217c2f6d86be
# ╟─5c017766-445d-4f4b-98f1-ae63e78ec34b
# ╠═0a76f019-4853-4ba3-9af8-9f33e1d4c956
# ╟─b0cc8694-b7ab-4d23-a208-055299840334
# ╟─680f104e-80b4-443f-b4bc-532df758c162
# ╠═ddfc38fc-b47d-4ea5-847a-e9cbee3aa0a1
# ╠═ee2616ca-2602-4823-9cfb-123b958701c4
# ╠═7a62c034-3709-483a-a663-7fe5e09cb773
# ╟─7f1fd9b4-517a-4fec-89bb-4d696dadbc3d
# ╟─81e29fc7-b5d3-46d8-aeac-fb8e6dc11b16
# ╟─5291b260-9a68-4c8b-aff4-7797804ccc95
# ╟─08dbe330-670d-48d5-b704-2421e687bff1
# ╟─c109b759-7b73-4593-b9ea-8cc97b61d6fe
# ╟─fe0fefb6-2755-4319-a944-bbbc7843aead
# ╟─60494b7c-1a08-4846-8a80-12533552a697
# ╟─b57195f9-c2a1-4676-96f9-faee84f7fc26
# ╟─438d437e-7b00-4a13-8f8a-87fdc332a190
# ╟─800fe4ba-e1a4-4e94-929f-7d66516e7bd6
# ╠═7d5a29c6-e71e-4ccb-a1c2-7fba663f038c
# ╟─3364e9f6-b2af-45f1-b1f3-ef7b0cd4910a
# ╟─26265a91-2c8e-46d8-9a87-a2d097e7433a
# ╟─2eeb402e-c5f9-449c-af19-ff8f2e6c7246
# ╠═6870ca6d-256d-4a38-970e-1c26ceba9fa4
# ╠═e5dac5c5-4644-443f-aa79-e43b399712c0
# ╠═85f98ea6-9351-4527-8b8e-b2827a7735ff
# ╟─859ce60b-2f32-44d1-919a-dbdaf1be38fb
# ╠═0377939c-00ac-42ae-b981-cdc897421588
# ╟─2f907e0d-171e-44c3-a531-5f11da08b3cf
# ╠═31b6d4ec-d057-44ca-875b-0c3257895dd3
# ╠═8902a846-fbb9-42fc-8742-c9c4a84db52c
# ╟─98ece9fe-dfcc-4dd8-bd47-049217d2afcf
# ╟─e66e67e8-8ac2-41a3-9926-3f0ac3b9c47d
# ╟─634c9cc1-5a93-42b4-bf51-17dadfe488d6
# ╟─31161289-1d4c-46ba-8bd9-e687fb7da29e
# ╟─50e01181-1911-426b-9228-4663a1297619
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
