include("main.jl")

# define scenarios
pL4T16r0c0tau60 = parameters(τᵤ = 60, NLockDown = 1, LockDelay = 4*7/365, NTestingPhases = 1, TestingDelay = 16*7/365, r = 0, c = 0, plotComment = "LockDelay = 4 weeks, TestDelay = 16 Weeks, r = 0, c = 0, τmax = 60")
pL4T16r0c1tau60 = parameters(τᵤ = 60, NLockDown = 1, LockDelay = 4*7/365, NTestingPhases = 1, TestingDelay = 16*7/365, r = 0, c = 1, plotComment = "LockDelay = 4 weeks, TestDelay = 16 Weeks, r = 0, c = 1, τmax = 60")
pL4T16r0c2tau60 = parameters(τᵤ = 60, NLockDown = 1, LockDelay = 4*7/365, NTestingPhases = 1, TestingDelay = 16*7/365, r = 0, c = 2, plotComment = "LockDelay = 4 weeks, TestDelay = 16 Weeks, r = 0, c = 2, τmax = 60")
pL4T16r1c0tau60 = parameters(τᵤ = 60, NLockDown = 1, LockDelay = 4*7/365, NTestingPhases = 1, TestingDelay = 16*7/365, r = 1, c = 0, plotComment = "LockDelay = 4 weeks, TestDelay = 16 Weeks, r = 1, c = 0, τmax = 60")

pL4T16r0c0eta005 = parameters(η = 0.05, NLockDown = 1, LockDelay = 4*7/365, NTestingPhases = 1, TestingDelay = 16*7/365, r = 0, c = 0, plotComment = "LockDelay = 4 weeks, TestDelay = 16 Weeks, r = 0, c = 0, eta = 0.05")
pL4T16r0c1eta005 = parameters(η = 0.05, NLockDown = 1, LockDelay = 4*7/365, NTestingPhases = 1, TestingDelay = 16*7/365, r = 0, c = 1, plotComment = "LockDelay = 4 weeks, TestDelay = 16 Weeks, r = 0, c = 1, eta = 0.05")
pL4T16r0c2eta005 = parameters(η = 0.05, NLockDown = 1, LockDelay = 4*7/365, NTestingPhases = 1, TestingDelay = 16*7/365, r = 0, c = 2, plotComment = "LockDelay = 4 weeks, TestDelay = 16 Weeks, r = 0, c = 2, eta = 0.05")
pL4T16r1c0eta005 = parameters(η = 0.05, NLockDown = 1, LockDelay = 4*7/365, NTestingPhases = 1, TestingDelay = 16*7/365, r = 1, c = 0, plotComment = "LockDelay = 4 weeks, TestDelay = 16 Weeks, r = 1, c = 0, eta = 0.05")

paramList = [
pL4T16r0c0tau60
pL4T16r0c1tau60
pL4T16r0c2tau60
pL4T16r1c0tau60
pL4T16r0c0eta005
pL4T16r0c1eta005
pL4T16r0c2eta005
pL4T16r1c0eta005
]
NScenarios = length(paramList)

# fix the execution time for multithreading
for (i,param) in enumerate(paramList)
    param.OPT_MAX_TIME[1] = NScenarios * param.OPT_MAX_TIME[1]
end

# define global variables only for debug
global GLOBAL_MAX_TIME = mean([p.OPT_MAX_TIME[1] for p in paramList])/3600
global GLOBAL_RUNNING_THREADS = ones(NScenarios)

# run all scenarios in parallel (modify each), global optimization
printstyled("starting optimization of $NScenarios scenarios\n", color=:yellow) 
for i in 1:NScenarios
    optimizeScenario!(paramList[i],Silent=true)
end
printstyled("done optimizing global\n", color=:yellow) 


# set optimization to local
for (i,param) in enumerate(paramList)
    param.optimizationType[1] = "local"
end

# run all scenarios in parallel (modify each), local optimization
Threads.@threads for i in 1:NScenarios
    optimizeScenario!(paramList[i],Silent=true)
end

# plot and save multiple scenarios all together 
plotAndSave([pL4T16r0c0tau60,pL4T16r0c1tau60,pL4T16r0c2tau60,pL4T16r1c0tau60], fName = "pL4T16tau60.eps")
plotAndSave([pL4T16r0c0eta005,pL4T16r0c1eta005,pL4T16r0c2eta005,pL4T16r1c0eta005], fName = "pL4T16eta005.eps")


# save results 
LockDownAll = [p.LockDownPolicy for p in paramList]
TestingAll = [p.TestingPolicy for p in paramList]
@save "LockDownAllHAND.jld2" LockDownAll
@save "TestingAllHAND.jld2" TestingAll