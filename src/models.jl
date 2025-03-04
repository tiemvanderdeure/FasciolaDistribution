#### Write all functions used as callable structs
# Growing-degree day curve
struct FixN{F, A}
    f::F
    args::A
end
FixN(f, args...) = FixN(f, args)
(f::FixN)(x) = f.f(f.args..., x)

function gdd_time(days, Tmin, x; maxval = Inf)
    if x > Tmin
        days / (x - Tmin)
    else
        maxval
    end
end
gaussian(rmax, Topt, σ, x) = rmax * exp(-(x - Topt)^2 / (2 * σ^2))
function quadratic(q, Tmin, Tmax, x; minval = zero(x))
    if Tmin < x < Tmax
        q * (x - Tmin) * (Tmax - x)
    else
        minval
    end
end
function linear(a, b, x)
    a + b * x
end
# linear between tmin and tmax, and then sigmoidally declining
function linear_sigmoid(r_Tmax, k, Tmin, Tmax, x; minval = zero(x))
    if x < Tmin
        minval
    elseif x > Tmax
        2 * r_Tmax * GLM.logistic(k * (Tmax - x))
    else
        β = (r_Tmax - 1) / (Tmax - Tmin)
        linear(1, β, x - Tmin)
        # 1 - (1 - rate_Tmax) * (x - Tmin) / (Tmax - Tmin)
    end
end

### Define the statistical models

# This is used for prepatent period and hatching time
@model function censored_gdd_hierarchical(
    T, # temperature of observatoins
    upper, # upper bound of the observations (experiment length)
    ids, # group ids
    y, # observed values
    ngroups = maximum(ids); 
    Tmin_prior = Normal(10, 5), 
    gdd_prior = LogNormal(6, 2)
)
    σ ~ Exponential(10)
    Tmin ~ Tmin_prior
    gdd_mu ~ gdd_prior
    σ_gdd ~ LogNormal(3,2)
    gdd ~ filldist(Normal(gdd_mu, σ_gdd), ngroups)
    y_hat = gdd_time.(gdd[ids], Tmin, T; maxval = floatmax())
    for i in eachindex(y)
        if ismissing(upper[i])
            y[i] ~ Normal(y_hat[i], σ)
        else
            y[i] ~ censored(Normal(y_hat[i], σ), nothing, upper[i])
        end
    end
    return FixN(gdd_time, gdd_mu, Tmin)
end

@model function hatching_success(
    T, n, k, ids, n_groups = maximum(ids); 
    priors = (Tmin = Normal(10, 5), Tmax = Normal(30, 5), k = Exponential(0.5))
    Tmin ~ Tmin_prior
    Tmax ~ Normal(30, 5)
    k ~ Exponential(0.5)
    rTmax ~ Beta(1,1) # 
    rmax ~ filldist(Beta(1,1), n_groups) # this varies per group
    y_hat = linear_sigmoid.(rTmax, k, Tmin, Tmax, T; minval = eps()) .* rmax[ids]
    k .~ Binomial.(n, y_hat)
    return FixN(linear_sigmoid, rTmax, k, Tmin, Tmax)
end
@model function gaussian_rmax_varying(T, ids, n_group = maximum(ids))
    Topt ~ Normal(20, 5)
    rmax_mu ~ Normal(6,2)
    breadth ~ Exponential(5)
    σ_rmax ~ Exponential(1)
    σ ~ Exponential(1)
    rmax ~ filldist(LogNormal(rmax_mu, σ_rmax), n_group)
    y_hat = gaussian.(rmax[ids], Topt, breadth, T)
    y ~ MvNormal(y_hat, σ)
    return FixN(gaussian, exp(rmax_mu), Topt, breadth)
end

@model function infectionefficiency(
    T_bin, n, k, T_rad, ids, y, ngroups = maximum(ids; init = 0); 
    Topt_prior = Normal(20,5), width_prior = Exponential(10)
)
    Topt ~ Topt_prior
    width ~ width_prior
    Tmax = Topt + width
    Tmin = Topt - width
    rmax ~ Uniform(0,1) # infection efficiency must be between 0 and 1
    if ngroups > 0
        rmax_radioactive ~ filldist(Uniform(0, 1000), ngroups)
        σ ~ Exponential(50)
        y_hat_rad = gaussian.(rmax_radioactive[ids], Topt, width, T_rad)
        y ~ MvNormal(y_hat_rad, σ)
    end

    #q = rmax / ((Tmax-Tmin)^2 / 4.0)
#    y_hat = quadratic.(q, Tmin, Tmax, T_bin; minval = eps())
    y_hat = gaussian.(rmax, Topt, width, T_bin)
    for i in eachindex(k)
        k[i] ~ Binomial(n[i], y_hat[i])
    end
    return FixN(gaussian, rmax, Topt, width)
end


#=

@model function infectionefficiency_gig(T_bin, n, k; Topt_prior = Normal(20,5), width_prior = Exponential(5))
    Topt ~ Topt_prior
    width ~ width_prior
    Tmax = Topt + width
    Tmin = Topt - width
    rmax ~ Uniform(0,1) # infection efficiency must be between 0 and 1
    q = rmax / ((Tmax-Tmin)^2 / 4.0)
    y_hat = quadratic.(q, Tmin, Tmax, T_bin; minval = eps())
    for i in eachindex(k)
        k[i] ~ Binomial(n[i], y_hat[i])
    end
    return FixN(quadratic, q, Tmin, Tmax)
end

@model function censored_gdd(T, upper, y; Tmin_prior = Normal(10, 5), gdd_prior = LogNormal(4, 2))
    σ ~ Exponential(10)
    Tmin ~ Tmin_prior
    gdd ~ gdd_prior
    gdd_time = Gdd(gdd, Tmin)
    for i in eachindex(y)
        y_hat = gdd_time(T[i])
        if ismissing(upper[i])
            y[i] ~ Normal(y_hat, σ)
        else
            y[i] ~ censored(Normal(y_hat, σ), nothing, upper[i])
        end
    end
    return gdd_time
end
=#
#=
@model function linear_model_binomial(T, n, k; intercept_prior = Normal(0,1), a_prior = Normal(0,1))
    intercept ~ intercept_prior
    a ~ a_prior
    y_hat = GLM.logistic.(intercept .+ a .* T) 
    k .~ Binomial.(n, y_hat)
    return Linear(intercept, a)
end
=#


#=
@model function survival(T, expl, n, k; intercept_prior = Normal(-5,1), a_prior = Normal(0,1))
    intercept ~ intercept_prior
    a ~ a_prior
    dr = exp.(intercept .+ a .* T)
    y_hat = exp.(-dr .* expl)
    k .~ Binomial.(n, y_hat)
    return (; y_hat)
end
=#


