
function objectiveFunction(param::parameters)
    p_join = (param, param.LockDownPolicy, param.TestingPolicy)

    return objectiveFunction(p_join)
end

function objectiveFunction(p_join)
    param = p_join[1]
    sol = solveScenario(p_join)
    if sol.retcode != :Success
        printstyled("LockDownPolicy = $(param.LockDownPolicy)\n",color=:red)
        return -Inf
    end
    return objectiveFunction(sol,p_join)
end

function objectiveFunction(sol,p_join)
    param = p_join[1]
    result = objectiveProduction(sol,p_join) + objectiveDeaths(sol,p_join)

    return result
end

function objectiveProduction(p_join)
    param = p_join[1]
    sol = solveScenario(p_join)
    if sol.retcode != :Success
        return -Inf
    end
    return objectiveProduction(sol,p_join)
end

function objectiveProduction(sol,p_join)
    param = p_join[1]
    
    result, err = quadgk( t-> objInstProduction(p_join,sol,t),0.0,sol.t[end],atol=0.1)
    return result
end

function objectiveDeaths(p_join)
    param = p_join[1]
    sol = solveScenario(p_join)
    if sol.retcode != :Success
        return -Inf
    end
    return objectiveDeaths(sol,p_join)
end

function objectiveDeaths(sol,p_join)
    param = p_join[1]
    result, err = quadgk( t-> - objInstDeaths(p_join,sol,t),0.0,sol.t[end],atol=0.1)
    return result
end

function objInstProduction(p_join,sol,t)
    param = p_join[1]
    if param.MAX_functional == "totalSocialWelfare"
        return U(p_join,sol,t)
    elseif param.MAX_functional == "VSL"
        return Y(p_join,sol,t)
    else
        error("MAX_functional unknown")
    end
end

function objInstDeaths(p_join,sol,t)
    param = p_join[1]
    if param.MAX_functional == "totalSocialWelfare"
        return param.θ * V(p_join,sol,t)
    elseif param.MAX_functional == "VSL"
        return param.VSL * totalDeaths(p_join,sol,t)
    else
        error("MAX_functional unknown")
    end
end

function totalDeaths(p_join,sol,t)
    param = p_join[1]
    S,C,E1,E2,I1,I2,R,V1,Vw,D,Dc  = sol(t)
    return param.N₀ * Dc
end

function Φ(p_join,x)
    param = p_join[1]
    return param.ρ₀*x + exp(param.ρ₁/(param.N₀-x)) - exp(param.ρ₁/param.N₀)
end

function x(p_join,sol,t)
    param = p_join[1]
    pₜ,τₜ = extract(p_join,t)
    S,C,E1,E2,I1,I2,R,V1,Vw,D,Dc = sol(t)
    return  param.N₀ * τₜ*(C + E1 + E2 + param.r*S)
end

function A(p_join,t)
    param = p_join[1]
    pₜ,τₜ = extract(p_join,t)
    return param.A₀*abs(1-pₜ)^param.Δ
end

function L(p_join,sol,t)
    param = p_join[1]
    pₜ,τₜ = extract(p_join,t)
    S,C,E1,E2,I1,I2,R,V1,Vw,D,Dc = sol(t)
    return  param.N₀*abs(S + param.ϵ_C*(C + E1 + E2) + R)
end

function Y(p_join,sol,t)
    param = p_join[1]
    pₜ,τₜ = extract(p_join,t)
    return A(p_join,t)*( abs(1-pₜ)*L(p_join,sol,t) )^param.α - 
                param.testing_cost*Φ(p_join,x(p_join,sol,t))
end

function U(p_join,sol,t)
    param = p_join[1]
    y = Y(p_join,sol,t)
    if y >= 0
        return y^(1-param.σ)/(1-param.σ)
        # return log(y)
    else
        return 0
    end
end

function V(p_join,sol,t)
    param = p_join[1]
    S,C,E1,E2,I1,I2,R,V1,Vw,D,Dc = sol(t)
    return (param.N₀ * Dc)^param.ω/param.ω
end

function TSW(p_join,sol,t)
    param = p_join[1]
    return U(p_join,sol,t) - param.θ * V(p_join,sol,t)
end


function optimizeScenario!(param; Silent::Bool = false)

    printstyled("optimizing scenario: " * param.plotComment *"\n") 

    @unpack NLockDown, NTestingPhases, T_horizon, LockDelay, TestingDelay = param
    @unpack τₗ,τᵤ = param

    # NLOPT
    model  = Model(NLopt.Optimizer);
    if param.optimizationType[1] == "global"
        set_optimizer_attribute(model, "algorithm", :GN_ORIG_DIRECT)
    elseif param.optimizationType[1] == "local"
        set_optimizer_attribute(model, "algorithm", :LN_COBYLA)
    else
        errror("optimizationType can be only global or local")
    end

    set_optimizer_attribute(model,"maxtime",param.OPT_MAX_TIME[1])


    JuMP.@variables(model,begin
        0 <= LockDownPolicy_optim_duration[1:NLockDown] <= T_horizon - LockDelay
        0 <= LockDownPolicy_optim_intensity[1:NLockDown] <= 1.0
        0 <= TestingPolicy_optim_duration[1:NTestingPhases] <= T_horizon - TestingDelay
        0 <= TestingPolicy_optim_intensity[1:NTestingPhases] <= τᵤ
    end)

    set_start_value.(LockDownPolicy_optim_duration, param.LockDownPolicy[1:NLockDown])
    set_start_value.(LockDownPolicy_optim_intensity, param.LockDownPolicy[(NLockDown+1):(2*NLockDown)])
    set_start_value.(TestingPolicy_optim_duration, param.TestingPolicy[1:NTestingPhases])
    set_start_value.(TestingPolicy_optim_intensity, param.TestingPolicy[(NTestingPhases+1):(2*NTestingPhases)])

    if Silent == false
        printstyled("starting with LockDown Duration = $(param.LockDownPolicy[1:NLockDown])\n",color=:yellow)
        printstyled("starting with LockDown Intensity = $(param.LockDownPolicy[(NLockDown+1):(2*NLockDown)])\n",color=:yellow)
        printstyled("starting with Testing Duration = $(param.TestingPolicy[1:NTestingPhases])\n",color=:yellow)
        printstyled("starting with Testing Intensity = $(param.TestingPolicy[(NTestingPhases+1):(2*NTestingPhases)])\n",color=:yellow)
    end

    JuMP.@constraint(model, LockDuration, 
        sum(LockDownPolicy_optim_duration[1:NLockDown]) <= (T_horizon - LockDelay))
    JuMP.@constraint(model, TestingDuration,
        sum(TestingPolicy_optim_duration[1:NTestingPhases]) <= (T_horizon - TestingDelay))

    function objectiveFunctionJuMP(args...)
        LDP = collect(args[1:(2*NLockDown)])       # <- converted to array 
        TP = collect(args[(2*NLockDown+1):end])
        p_join = (param,LDP,TP)
        val = objectiveFunction(p_join)/param.N₀
        if Silent == false
            println("LockDown = $(ForwardDiff.value.(LDP)), Testing = $(ForwardDiff.value.(TP)), value = $(ForwardDiff.value(val))")
        end
        return val
    end
    JuMP.register(model, :objectiveFunctionJuMP, 2*(NLockDown+NTestingPhases), 
                                        objectiveFunctionJuMP, autodiff=true)
    @NLobjective(model,Max,objectiveFunctionJuMP(LockDownPolicy_optim_duration...,
                                                LockDownPolicy_optim_intensity...,
                                                TestingPolicy_optim_duration...,
                                                TestingPolicy_optim_intensity...))
    
    if Silent == false
        JuMP.unset_silent(model)
        JuMP.optimize!(model)
    else 
        JuMP.set_silent(model)
        JuMP.optimize!(model)
        # redirect_stdout(( ()-> JuMP.optimize!(model) ),open("/dev/null", "w"))
    end

    param.LockDownPolicy .= [value.(model[:LockDownPolicy_optim_duration])...,value.(model[:LockDownPolicy_optim_intensity])...]
    param.TestingPolicy .= [value.(model[:TestingPolicy_optim_duration])...,value.(model[:TestingPolicy_optim_intensity])...]

    timeSpentH = solve_time(model)/3600
    timeSpentHRound = round(timeSpentH,digits=1)

    printstyled("scenario '" * param.plotComment * "' finished optimization with terminations status " * string(termination_status(model)) * " in $timeSpentHRound hours\n", color=:green) 
    
    # multithreading print info
    if @isdefined GLOBAL_MAX_TIME
        NrunnigThreadsNow = sum(GLOBAL_RUNNING_THREADS)
        global GLOBAL_MAX_TIME -= NrunnigThreadsNow*timeSpentH
        global GLOBAL_RUNNING_THREADS[Threads.threadid()] = 0 
        if NrunnigThreadsNow > 1
            remainingMeanTime = GLOBAL_MAX_TIME/(NrunnigThreadsNow-1)
        else 
            remainingMeanTime = 1
        end
        printstyled("remaining time = $(round(remainingMeanTime,digits=1)) h\n",color=:yellow)
        printstyled("threads remaining = $(Int(NrunnigThreadsNow-1)) \n",color=:yellow)
    end

    return nothing
end
