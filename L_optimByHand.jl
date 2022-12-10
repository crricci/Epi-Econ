
using ProgressMeter

function optimize1PhaseLock(param; accuracy = 0.01)
    printstyled("optimizing scenario: " * param.plotComment *"\n") 

    @assert param.NLockDown == 1
    @assert param.NTestingPhases == 0

    durations = 0.0:accuracy:param.T_horizon-param.LockDelay
    # intensities = 0.0:accuracy:1.0
    intensities = 0.0:accuracy:(1.0-accuracy) # for σ > 1
    allPolicies = [[duration,intensity] for duration in durations for intensity in intensities]

    param_LF = copy(param)
    param_LF.LockDownPolicy .= zeros(2*param_LF.NLockDown)
    param_LF.TestingPolicy .= zeros(2*param_LF.NTestingPhases)
    TSW_LF = objectiveFunction(param_LF)

    TSW_best = TSW_LF
    policyBest = [0.0,0.0]
    @showprogress for policy in allPolicies
        param.LockDownPolicy .= policy
        TSW_policy = objectiveFunction(param)
        if TSW_policy > TSW_best 
            TSW_best = TSW_policy
            policyBest = policy
        end
    end

    param.LockDownPolicy .= policyBest
    # plotScenario(param)
end

function optimize2PhaseLock(param; accuracy = 0.1)
    printstyled("optimizing scenario: " * param.plotComment *"\n") 

    @assert param.NLockDown == 2
    @assert param.NTestingPhases == 0

    intensities1 = 0.0:accuracy:1.0
    intensities2 = 0.0:accuracy:1.0
    totalDurations = 0.0:accuracy:param.T_horizon-param.LockDelay
    allPolicies = [[duration1,totalDuration-duration1,intensity1,intensity2] 
    for totalDuration in totalDurations 
    for duration1 in 0.0:accuracy:totalDuration
    for intensity1 in intensities1 
    for intensity2 in intensities2 ] 

    param_LF = copy(param)
    param_LF.LockDownPolicy .= zeros(2*param_LF.NLockDown)
    param_LF.TestingPolicy .= zeros(2*param_LF.NTestingPhases)
    TSW_LF = objectiveFunction(param_LF)

    TSW_best = TSW_LF
    policyBest = zeros(2*param.NLockDown)
    @showprogress for policy in allPolicies
        param.LockDownPolicy .= policy
        TSW_policy = objectiveFunction(param)
        if TSW_policy > TSW_best 
            TSW_best = TSW_policy
            policyBest = policy
        end
    end

    param.LockDownPolicy .= policyBest
    # plotScenario(param)

end

function optimize1PhaseTesting(param; accuracy = 0.01)
    printstyled("optimizing scenario: " * param.plotComment *"\n") 

    @assert param.NTestingPhases == 1
    @assert param.NLockDown == 0

    accuracySteps = Int(1/accuracy)
    intensities = LinRange(0.0,param.τᵤ,accuracySteps)
    
    durations = 0.0:accuracy:param.T_horizon-param.TestingDelay
    allPolicies = [[duration,intensity] for duration in durations for intensity in intensities]

    param_LF = copy(param)
    param_LF.LockDownPolicy .= zeros(2*param_LF.NLockDown)
    param_LF.TestingPolicy .= zeros(2*param_LF.NTestingPhases)
    TSW_LF = objectiveFunction(param_LF)

    TSW_best = TSW_LF
    policyBest = zeros(2*param.NTestingPhases)
    @showprogress for policy in allPolicies
        param.TestingPolicy .= policy
        TSW_policy = objectiveFunction(param)
        if TSW_policy > TSW_best 
            TSW_best = TSW_policy
            policyBest = policy
        end
    end

    param.TestingPolicy .= policyBest
    # plotScenario(param)
end


function optimize2PhaseTesting(param; accuracy = 0.1)
    printstyled("optimizing scenario: " * param.plotComment *"\n") 

    @assert param.NTestingPhases == 2
    @assert param.NLockDown == 0

    accuracySteps = Int(1/accuracy)
    intensities1 = LinRange(0.0,param.τᵤ,accuracySteps)
    intensities2 = LinRange(0.0,param.τᵤ,accuracySteps)

    totalDurations = 0.0:accuracy:param.T_horizon-param.TestingDelay
    allPolicies = [[duration1,totalDuration-duration1,intensity1,intensity2] 
    for totalDuration in totalDurations 
    for duration1 in 0.0:accuracy:totalDuration
    for intensity1 in intensities1 
    for intensity2 in intensities2 ] 

    param_LF = copy(param)
    param_LF.LockDownPolicy .= zeros(2*param_LF.NLockDown)
    param_LF.TestingPolicy .= zeros(2*param_LF.NTestingPhases)
    TSW_LF = objectiveFunction(param_LF)

    TSW_best = TSW_LF
    policyBest = zeros(2*param.NTestingPhases)
    @showprogress for policy in allPolicies
        param.TestingPolicy .= policy
        TSW_policy = objectiveFunction(param)
        if TSW_policy > TSW_best 
            TSW_best = TSW_policy
            policyBest = policy
        end
    end

    param.TestingPolicy .= policyBest
    # plotScenario(param)

end

function optimize1PhaseLock1PhaseTest(param; accuracy = 0.1)
    printstyled("optimizing scenario: " * param.plotComment *"\n") 

    @assert param.NLockDown == 1
    @assert param.NTestingPhases == 1

    accuracySteps = Int(1/accuracy)
    intensitiesLock = 0.0:accuracy:1.0
    intensitiesTest = LinRange(0.0,param.τᵤ,accuracySteps)
    durationsLock = 0.0:accuracy:param.T_horizon-param.LockDelay
    durationsTest = 0.0:accuracy:param.T_horizon-param.TestingDelay

    allPolicies = [[durationLock,intensityLock,durationTest,intensityTest] 
                    for durationLock in durationsLock
                    for intensityLock in intensitiesLock
                    for durationTest in durationsTest
                    for intensityTest in intensitiesTest]

    param_LF = copy(param)
    param_LF.LockDownPolicy .= zeros(2*param_LF.NLockDown)
    param_LF.TestingPolicy .= zeros(2*param_LF.NTestingPhases)
    TSW_LF = objectiveFunction(param_LF)

    TSW_best = TSW_LF
    policyBest = [0.0,0.0,0.0,0.0]
    @showprogress for policy in allPolicies
        param.LockDownPolicy .= policy[1:2]
        param.TestingPolicy .= policy[3:4]

        TSW_policy = objectiveFunction(param)
        if TSW_policy > TSW_best 
            TSW_best = TSW_policy
            policyBest = policy
        end
    end

    param.LockDownPolicy .= policyBest[1:2]
    param.TestingPolicy .= policyBest[3:4]
    # plotScenario(param)
end


function optimize1PhaseLockGrid(param; accuracy = 0.01)
    printstyled("optimizing scenario: " * param.plotComment *"\n") 

    @assert param.NLockDown == 1
    @assert param.NTestingPhases == 0

    durations = 0.0:accuracy:param.T_horizon-param.LockDelay
    intensities = 0.0:accuracy:(1.0-accuracy) # for σ > 1
    allPolicies = [[duration,intensity] for duration in durations for intensity in intensities]
    val = zeros(length(allPolicies))

    param_LF = copy(param)
    param_LF.LockDownPolicy .= zeros(2*param_LF.NLockDown)
    param_LF.TestingPolicy .= zeros(2*param_LF.NTestingPhases)
    TSW_LF = objectiveFunction(param_LF)

    TSW_best = TSW_LF
    policyBest = [0.0,0.0]
    @showprogress for (i,policy) in enumerate(allPolicies)
        param.LockDownPolicy .= policy
        TSW_policy = objectiveFunction(param)
        val[i] = TSW_policy
        if TSW_policy > TSW_best 
            TSW_best = TSW_policy
            policyBest = policy
        end
    end

    param.LockDownPolicy .= policyBest
    # plotScenario(param)
    return val
end


