
function plotScenario(param; 
    comment::String = "",T_plot::Float64 = 0.85, Close::Bool = false,
    Show::Bool = true, Save::Bool = false, fName::String = "noname.eps", LF::Bool = false)
    
    sol = solveScenario(param)
    
    tᵢ_final = findfirst(sol.t .>= param.T_horizon * T_plot)
    times = sol.t[1:tᵢ_final]
    
    p_join = (param, param.LockDownPolicy, param.TestingPolicy)
    pτ = [extract(p_join,t) for t in times]
    pₜ = [el[1] for el in pτ]
    τₜ = [el[2] for el in pτ]
    
    S,C,E1,E2,I1,I2,R,V1,Vw,D,Dc = [row[1:tᵢ_final] for row in eachrow(sol)]
    if LF
        param_LF = copy(param)
        param_LF.LockDownPolicy .= zeros(2*param_LF.NLockDown)
        param_LF.TestingPolicy .= zeros(2*param_LF.NTestingPhases)
        sol_LF = solveScenario(param_LF)
        p_join_LF = (param_LF, param_LF.LockDownPolicy, param_LF.TestingPolicy)
        pₜ_LF = zeros(length(times))
        τₜ_LF = zeros(length(times))
        S_LF,C_LF,E1_LF,E2_LF,I1_LF,I2_LF,R_LF,V1_LF,Vw_LF,D_LF,Dc_LF = [row[1:tᵢ_final] for row in eachrow(sol_LF)]
    end

    @suppress begin
    # hide plot and set window size
    ioff()
    fig = gcf()
    fig.set_figheight(14)
    fig.set_figwidth(25)

    # LockDownPolicy
    subplot(3,4,1)
    grid(true)
    PyPlot.title("LockDownPolicy",fontsize=10)
    plot(times,pₜ)
    ylim(-0.1,1.1)
    if LF
        LFp = twinx()
        LFp.tick_params(axis="y",color="red",labelcolor="red")
        plot(times,pₜ_LF,color="black")
        ylim(-0.1,1.1)
    end


    # Production in % w.r.t initial production
    Yₜ = [Y(p_join,sol,t) for t in times]
    subplot(3,4,2)
    grid(true)
    PyPlot.title("Production in % w.r.t initial production",fontsize=10)
    plot(times,100 * Yₜ/Yₜ[1])
    if LF
        LFp = twinx()
        LFp.tick_params(axis="y",color="red",labelcolor="red")
        Yₜ_LF = [Y(p_join_LF,sol_LF,t) for t in times]
        plot(times,100 * Yₜ_LF/Yₜ_LF[1],color="black")
    end
    
    # S - Susceptible hosts
    subplot(3,4,3)
    grid(true)
    PyPlot.title("S - Susceptible hosts",fontsize=10)
    ylabel("% of population")
    plot(times, 100 * S)
    if LF
        LFp = twinx()
        LFp.tick_params(axis="y",color="red",labelcolor="red")
        plot(times, 100 * S_LF,color="black")
    end

    # C - Number of contacts 
    subplot(3,4,4)
    grid(true)
    PyPlot.title("C - Number of contacts",fontsize=10)
    ylabel("% of population")
    plot(times,100 * C)
    if LF
        LFp = twinx()
        LFp.tick_params(axis="y",color="red",labelcolor="red")
        plot(times, 100 * C_LF,color="black")
    end

    # E1 - Exposed hosts, no symptoms
    subplot(3,4,5)
    grid(true)
    PyPlot.title("E1 - Exposed hosts, no symptoms",fontsize=10)
    ylabel("% of population")
    plot(times,100 * E1)
    if LF
        LFp = twinx()
        LFp.tick_params(axis="y",color="red",labelcolor="red")
        plot(times, 100 * E1_LF,color="black")
    end

    # E2 - Exposed hosts, no symptoms,   transmitting
    subplot(3,4,6)
    grid(true)
    PyPlot.title("E2 - Exposed hosts, no symptoms,   transmitting",fontsize=10)
    ylabel("% of population")
    plot(times,100 * E2)
    if LF
        LFp = twinx()
        LFp.tick_params(axis="y",color="red",labelcolor="red")
        plot(times, 100 * E2_LF,color="black")
    end

    # I1 - non severe infected,  partially-isolated hosts
    subplot(3,4,7)
    grid(true)
    PyPlot.title("I1 - non severe infected,  partially-isolated hosts",fontsize=10)
    ylabel("% of population")
    plot(times,100 * I1)
    if LF
        LFp = twinx()
        LFp.tick_params(axis="y",color="red",labelcolor="red")
        plot(times, 100 * I1_LF,color="black")
    end

    # I2 - hospitalized infected,    partially-isolated hosts
    subplot(3,4,8)
    grid(true)
    PyPlot.title("I2 - hospitalized infected,    partially-isolated hosts",fontsize=10)
    ylabel("% of population")
    plot(times,100 * I2)
    if LF
        LFp = twinx()
        LFp.tick_params(axis="y",color="red",labelcolor="red")
        plot(times, 100 * I2_LF,color="black")
    end

    # R - Recovered and immune hosts
    subplot(3,4,9)
    grid(true)
    PyPlot.title("R - Recovered and immune hosts",fontsize=10)
    ylabel("% of population")
    plot(times,100 * R)
    if LF
        LFp = twinx()
        LFp.tick_params(axis="y",color="red",labelcolor="red")
        plot(times, 100 * R_LF,color="black")
    end

    # τ - Testing Rate
    subplot(3,4,10)
    grid(true)
    PyPlot.title("τ - Testing Rate",fontsize=10)
    plot(times,τₜ)
    if LF
        LFp = twinx()
        LFp.tick_params(axis="y",color="red",labelcolor="red")
        plot(times,τₜ_LF,color="black")
    end
    

    # Dc - Pandemia deaths
    subplot(3,4,11)
    grid(true)
    PyPlot.title("Dc - Pandemia deaths",fontsize=10)
    ylabel("% of population")
    plot(times,100 * Dc)
    if LF
        LFp = twinx()
        LFp.tick_params(axis="y",color="red",labelcolor="red")
        plot(times, 100 * Dc_LF,color="black")
    end

    # Comments
    obj = round(objectiveFunction(p_join),digits=1)
    obj_prod = round(objectiveProduction(p_join),digits=1)
    obj_deaths = round(objectiveDeaths(p_join),digits=1)
    if comment == ""
        comment = param.plotComment
    end
    label = comment * ",\nObjective = " * string(obj) * ",\nProd = " * string(obj_prod) * ",\nDeaths = " * string(obj_deaths)
    subplot(3,4,12)
    grid(true)
    plot(0,0,label=label)
    legend()

    tight_layout()
    # resolve options
    if Save == true
        savefig(fName)
    end
    if Show == true
        ion()
        show()
    end
    if Close == true
        close("all")
    end

    end #suppress
end

function plotAndSave(paramList::Vector{parameters}; fName::String = "noname.eps",LF::Bool = true)
    lastPlot = paramList[end]
    for param in paramList[1:end-1]
        plotScenario(param,Show=false,comment = param.plotComment,LF = LF)
    end 
    plotScenario(lastPlot,Show=false,Save=true,fName=fName,
                    comment = lastPlot.plotComment,Close=true, LF = LF)
    return nothing
end

function plotSolution(sol,solPreLockDown,pₜ,param; T_plot::Float64 = 0.85)

    @unpack_parameters param

    pₜInterp = LinearInterpolation(sol.t, pₜ)
    tᵢ_final = findfirst(sol.t .>= param.T_horizon * T_plot)
    times = [solPreLockDown.t[1:end-1]...,sol.t[1:tᵢ_final]...]

    pₜ = [zeros(length(solPreLockDown.t)-1)...,pₜInterp.(sol.t[1:tᵢ_final])...]

    S,C,E1,E2,I1,I2,R,V1,Vw,D,Dc = vcat.(N₀*[row[1:end-1] for row in eachrow(solPreLockDown)],[row[1:tᵢ_final] for row in eachrow(sol)])


    for t in 1:length(times)
        state = [S[t],C[t],E1[t],E2[t],I1[t],I2[t],R[t],V1[t],Vw[t],D[t],Dc[t]]
        E1[t] = E1[t]*(1-Iₒ(state,param))
        E2[t] = E2[t]*(1-Iₒ(state,param))
        I1[t] = I1[t]*(1-Iₒ(state,param))
        I2[t] = I2[t]*(1-Iₒ(state,param))
    end

    SI = LinearInterpolation(times,S)
    CI = LinearInterpolation(times,C)
    E1I = LinearInterpolation(times,E1)
    E2I = LinearInterpolation(times,E2)
    I1I = LinearInterpolation(times,I1)
    I2I = LinearInterpolation(times,I2)
    RI = LinearInterpolation(times,R)
    V1I = LinearInterpolation(times,V1)
    VwI = LinearInterpolation(times,Vw)
    DI = LinearInterpolation(times,D)
    DcI = LinearInterpolation(times,Dc)
    pₜI = LinearInterpolation(times,pₜ)

    State_t(t) = [SI(t),CI(t),E1I(t),E2I(t),I1I(t),I2I(t),RI(t),V1I(t),VwI(t),DI(t),DcI(t)]
    A(t) =  A₀ * (1-pₜI(t))^Δ
    L(t) =  ( SI(t) + ϵ_C*(CI(t) + 
                (E1I(t) + E2I(t)))*(1-Iₒ(State_t(t),param)) + RI(t) )
    Y(t) =  A(t)*((1-pₜI(t))*L(t))^α 
    U(Y) =  (Y)^(1-σ)/(1-σ)
    V(Dc) =  (Dc)^ω/ω	
    TSW_tProd(t) =  U(Y(t)) 
    TSW_tDeath(t) =  - θ*V(DcI(t))

    TSWProd, err = quadgk( t-> TSW_tProd(t),0.0,times[end],atol=0.1)
    TSWDeath, err = quadgk( t-> TSW_tDeath(t),0.0,times[end],atol=0.1)
    TSWAll = TSWProd + TSWDeath

    @suppress begin
        # hide plot and set window size
    ioff()
    fig = gcf()
    fig.set_figheight(14)
    fig.set_figwidth(25)

    # LockDownPolicy
    subplot(3,4,1)
    grid(true)
    PyPlot.title("LockDownPolicy",fontsize=10)
    plot(times,pₜI.(times))
    ylim(-0.1,1.1)
    
    # Production in % w.r.t initial production
    subplot(3,4,2)
    grid(true)
    PyPlot.title("Production in % w.r.t initial production",fontsize=10)
    plot(times,100 * Y.(times)/Y(0))
    
    # S - Susceptible hosts
    subplot(3,4,3)
    grid(true)
    PyPlot.title("S - Susceptible hosts",fontsize=10)
    ylabel("% of population")
    plot(times, 100 * SI.(times)/N₀)

    # C - Number of contacts 
    subplot(3,4,4)
    grid(true)
    PyPlot.title("C - Number of contacts",fontsize=10)
    ylabel("% of population")
    plot(times,100 * CI.(times)/N₀)


    # E1 - Exposed hosts, no symptoms
    subplot(3,4,5)
    grid(true)
    PyPlot.title("E1 - Exposed hosts, no symptoms",fontsize=10)
    ylabel("% of population")
    plot(times,100 * E1I.(times)/N₀)

    # E2 - Exposed hosts, no symptoms,   transmitting
    subplot(3,4,6)
    grid(true)
    PyPlot.title("E2 - Exposed hosts, no symptoms,   transmitting",fontsize=10)
    ylabel("% of population")
    plot(times,100 * E2I.(times)/N₀)

    # I1 - non severe infected,  partially-isolated hosts
    subplot(3,4,7)
    grid(true)
    PyPlot.title("I1 - non severe infected,  partially-isolated hosts",fontsize=10)
    ylabel("% of population")
    plot(times,100 * I1I.(times)/N₀)


    # I2 - hospitalized infected,    partially-isolated hosts
    subplot(3,4,8)
    grid(true)
    PyPlot.title("I2 - hospitalized infected,    partially-isolated hosts",fontsize=10)
    ylabel("% of population")
    plot(times,100 * I2I.(times)/N₀)


    # R - Recovered and immune hosts
    subplot(3,4,9)
    grid(true)
    PyPlot.title("R - Recovered and immune hosts",fontsize=10)
    ylabel("% of population")
    plot(times,100 * RI.(times)/N₀)

    # τ - Testing Rate
    subplot(3,4,10)
    grid(true)
    PyPlot.title("τ - Testing Rate",fontsize=10)
    plot(times,zeros(length(times)))
    

    # Dc - Pandemia deaths
    subplot(3,4,11)
    grid(true)
    PyPlot.title("Dc - Pandemia deaths",fontsize=10)
    ylabel("% of population")
    plot(times,100 * DcI.(times)/N₀)

    # Comments
    label = "continuous" * ",\nObjective = " * string(TSWAll) * ",\nProd = " * string(TSWProd) * ",\nDeaths = " * string(TSWDeath)
    subplot(3,4,12)
    grid(true)
    plot(0,0,label=label)
    legend()

    tight_layout()
    # resolve options
    ion()
    show()

    end #suppress
end 