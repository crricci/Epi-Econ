
@with_kw struct parameters

    #INITIAL CONDITION
    N₀::Float64 = 1e5                       # initial number of inidividuals 
    N₀_infect::Float64 = 10.0 		        # initial infects (in E2)

    #EPIDEMIC
    μ::Float64 = 1/70 	   # natural death rate
    μᵢ₁::Float64 = 0.02    # death rate of I₁
    μᵣ::Float64 = 0.0      # death rate of R
    μₑ::Float64 = 0.01     # death rate of E₂
    pₘ::Float64 = 1/3      # death rate of severe infected individuals
    γ::Float64 = 0.05      # share of infection between E₂ and I
    ϕ::Float64 = 0.5  	   # transition from R to S
    φ₁::Float64 = 365/5    # transition from E₁ to E₂
    φ₂::Float64 = 365/5    # transition from E₂ to I
    δ::Float64 = 365/18    # transition from C to S
    δ₁::Float64 = 365/18   # transition from I₁ to R
    δ₂::Float64 = 365/13   # transition from I₂ to r
    β::Float64 = 150.0     # infectious rate 

    #HOSPITALIZATION
    ωₕ::Float64 = 365/11   # hospitalization rate
    η::Float64 = 0.15      # hospitalization (fraction, in (0,1))

    #VACCINATION
    TVac::Float64 = 0.4                     # initial time of vax
    ϵ_V::Float64 = 1.0                      # (no) chance for V₁ to be infected 
    w::Float64 = 365/35                     # waiting time for second dose
    ν::Float64 = 0.0                        # vaccination rate 

    #PRODUCTION                 
    A₀::Float64 = 1e5                       # basic tecnological value
    α::Float64 = 2/3                        # cobb douglas
    Δ::Float64 = 1/5                        #   
    ϵ_C::Float64 = 0.5                      # work-efficiency of people in C
    yearly_US_GDP::Float64 = 21.43*1e12	    # yearly US total GDP in 2019 in $
    Y₀::Float64 = A₀*N₀^α                   # yearly model total GDP with no pandemic in model unit
    VSL_USD::Float64 = 9.7*1e6              # VSL in USDollars
    VSL::Float64 = VSL_USD/yearly_US_GDP*Y₀ # VSL in model Unit

    #SATISFACTION                   
    σ::Float64 = 1/2                        # 
    ω::Float64 = 3.0                        # 
    θ::Float64 = 1e-4                       # relevance of deaths

    #TESTING
    testing_cost::Float64 = 1.0 	        # if 0 no testing cost
    c::Float64 = 1.0       		            # contact information(0.0 = perfect)
    sₚ::Float64 = 0.9 				        # specificity of testing
    sₑ::Float64 = 0.9 				        # sensitivity of testing
    r::Float64 = 0.0     			        # probability of random test
    testCost::Float64 = 50.0 		        # cost for one test in $ 

    ρ₀::Float64 = Y₀*testCost/yearly_US_GDP # cost for one test in model unit
    ρ₁::Float64 = 0.1                       # 
    τₗ::Float64 = 365/20                    # minimum testing rate
    τᵤ::Float64 = 365/2                     # maximum testing rate
    
    #SCENARIO  
    T_horizon::Float64 = 3.5   	            # N of Years 
    
    #LOCKDOWN POLICY
    LockDown_default::Float64 = 0.0         # defautl LockDown Intensity    
    NLockDown::Int = 2                      # N of LockDown phases
    LockDelay::Float64 = 2*7/365 	        # time of delay of LockDown
    LPhaseDelay::Float64 = 5/365            # time of smoothing between phases

    #TESTING POLICY
    τ_default::Float64 = 0.0 		        # default testing rate
    NTestingPhases::Int = 1                 # N of Testing phases, 0 -> no testing 
    TestingDelay::Float64 = 7/365           # time of delay of Testing
    TPhaseDelay::Float64 = 5/365            # time of smoothing between phases
    optimizationType::Vector{String} = ["global"]     # "global" or "local"

    #EXTRA
    LockDownPolicy::Vector{Float64} = zeros(2*NLockDown)      # only for output
    TestingPolicy::Vector{Float64} = zeros(2*NTestingPhases)  # only for output
    continuity_functional::String = "policy_constant" 
    MAX_functional::String = "totalSocialWelfare"       #["totalSocialWelfare,VSL"] 
    labels::Vector{String} = ["S - Susceptible hosts",
                    "C - Number of contacts",
                    "E1 - Exposed hosts, no symptoms",
                    "E2 - Exposed hosts, no symptoms,   transmitting",
                    "I1 - non severe infected,  partially-isolated hosts",
                    "I2 - hospitalized infected,    partially-isolated hosts",
                    "R - Recovered and immune hosts",
                    "V1 - People vaccinated for the first time",
                    "Vw - People waiting for their second dose",
                    "D - Natural deaths",
                    "Dc - Pandemia deaths"]

    #ODE
    u0 = initialCondition(N₀,N₀_infect)
    ODEalg = AutoVern9
    ODEalgStiff = RadauIIA5
    ODEcb = DiscreteCallback(condition_cut,affect_cut!,save_positions = (false,false))
    ODEinfectIndex = [3,4,5,6]
    ODEprob = function ODEmodel(p_join::Tuple{parameters, AbstractVector{T}, AbstractVector{S}}) where {T,S}
        param = p_join[1]
        if T != Union{}
            return ODEProblem(df!,T.(param.u0),(0.0,T_horizon),p_join)
        elseif S != Union{}
            return ODEProblem(df!,S.(param.u0),(0.0,T_horizon),p_join)
        else 
            return ODEProblem(df!,param.u0,(0.0,T_horizon),p_join)
        end
    end
    ODEt_save::StepRangeLen{Float64, Base.TwicePrecision{Float64}, 
                            Base.TwicePrecision{Float64}} = 0.0:1/365:T_horizon
    ODEabstol::Float64 = 1e-16
    ODEreltol::Float64 = 1e-14
    ODEmaxiter::Int = Int(1e6)
    ODEsolAtDelay::Vector{Float64} = zeros(11)

    σGuss::Float64 = 1e-5
    XGauss::Normal{Float64} = Normal(0,σGuss)

    #PLOTS
    plotComment::String = ""

    #JuMP
    OPT_MAX_TIME::Vector{Float64} = [3600 * 8.0]


    #CHECKS 
    @assert plotComment != ""
end

Base.copy(x::T) where T = T([deepcopy(getfield(x, k)) for k ∈ fieldnames(T)]...)
