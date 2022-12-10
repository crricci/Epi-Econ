
function df!(du,u,p_join,t)
	S,C,E1,E2,I1,I2,R,V1,Vw,D,Dc = u
	
	param = p_join[1]
	pₜ,τₜ = extract(p_join,t)
	@unpack_parameters param
	
	du[1] = dS = μ-μ*S-ν*S-(1-pₜ)^2*β*S*(E2*(1-γ)+γ*I1)*(1+c) + ϕ*R + 
		τₜ*(r+(1-r)/(1+c))*C*sₚ + δ*C - τₜ*r*S*(1-sₚ)
	du[2] = dC = (1-pₜ)^2*β*(E2*(1-γ)+γ*I1)*(c*S-(1-ϵ_C)*C) + τₜ*r*S*(1-sₚ) -
		μ*C -τₜ*(r+(1-r)/(1+c))C*sₚ- δ*C
	du[3] = dE1 = (1-pₜ)^2 * (S + (1-ϵ_C)*C ) * β * (E2*(1-γ)+γ*I1) -
		(μ + φ₁)*E1 - τₜ*(r+(1-r)/(1+c))sₑ*E1
	du[4] = dE2 = φ₁*E1 - (μ + φ₂)*E2 - τₜ*(r+(1-r)/(1+c))sₑ*E2
	du[5] = dI1 = φ₂*E2 + τₜ*(r+(1-r)/(1+c))sₑ*(E1+E2) - (μ + μᵢ₁ + (1-η)*δ₁ + η*δ₁)*I1
	du[6] = dI2 = η*δ₁*I1 - (μ + pₘ*δ₂ + (1-pₘ)*δ₂)*I2
	du[7] = dR = (1-η)*δ₁*I1 +  (1-pₘ)*δ₂*I2 - (μ + μᵣ + ϕ)*R + ν*Vw
	du[8] = dV1 = ν*S - (1-ϵ_V)*V1*β*(E2*(1-γ)+γ*I1) - (μ + w)*V1
	du[9] = dVw = w*V1 - (ν + μ)*Vw
	du[10] = dD = μ
	du[11] = dDc = μᵢ₁*I1 + pₘ*δ₂*I2 + μₑ*E2 + μᵣ*R

	return nothing
end

function extract(p_join,t)
	param, LockDownPolicy, TestingPolicy = p_join

	@unpack continuity_functional, T_horizon = param
	@unpack NLockDown, LockDelay, LPhaseDelay = param
	@unpack NTestingPhases, TestingDelay, TPhaseDelay = param

	if NLockDown == 0
		pₜ = param.LockDown_default
	else
		if param.continuity_functional == "policy_constant"
			pₜ = policy_constant(LockDownPolicy,NLockDown,LockDelay,LPhaseDelay,T_horizon,t)
		elseif param.continuity_functional == "policy_linear"
			pₜ = policy_linear(LockDownPolicy,NLockDown,LockDelay,LPhaseDelay,T_horizon,t)
		elseif param.continuity_functional == "policy_cubic"
			pₜ = policy_cubic(LockDownPolicy,NLockDown,LockDelay,LPhaseDelay,T_horizon,t)
		else 
			printstyled("continuity_functional of wrong type\n",color=:red)
		end
	end

	if NTestingPhases == 0 
		τₜ = param.τ_default
	else
		if param.continuity_functional == "policy_constant"
			τₜ = policy_constant(TestingPolicy,NTestingPhases,TestingDelay,TPhaseDelay,T_horizon,t)
		elseif param.continuity_functional == "policy_linear"
			τₜ = policy_linear(TestingPolicy,NTestingPhases,TestingDelay,TPhaseDelay,T_horizon,t)
		elseif param.continuity_functional == "policy_cubic"
			τₜ = policy_cubic(TestingPolicy,NTestingPhases,TestingDelay,TPhaseDelay,T_horizon,t)
		else 
			printstyled("continuity_functional of wrong type\n",color=:red)
		end
	end
	
	return pₜ,τₜ
end


function policy_constant(p,NPhases,delay,phaseDelay,T_horizon,t)
	""" durations are d1,d2... 
		times are t1 = LockDelay, t2 = LockDelay+d1,	t3 = t1+d2... 
		intervals are [0,t1 = LockDelay),[t1 = LockDelay,t2),[t2,t3) ... 
		intensities are i1 = 0.0, i2, i3 ...
	"""
	durations = [delay; p[1:NPhases]; T_horizon-delay-sum(p[1:NPhases])]				
	intensities = [0.0; p[NPhases+1:2*NPhases]; 0.0]		
	ti = cumsum(durations); ti[end] = T_horizon	#for rounding error		

	if t <= ti[1] || t >= ti[end]
		return 0.0
	end
	i = findfirst(t .<= ti)
	if isnothing(i)
		printstyled("ERROR IN policy, t = $t, p = $p\n",color=:red)
	end
	return intensities[i]
end

function policy_linear(p,NPhases,delay,phaseDelay,T_horizon,t)
	""" durations are d1,d2... 
		times are t1 = delay, t2 = delay+d1,	t3 = t1+d2... 
		intervals are [0,t1 = delay),[t1 = delay,t2),[t2,t3) ... 
		intensities are i1 = 0.0, i2, i3 ...
	"""
	durations = [delay; p[1:NPhases]; T_horizon-delay-sum(p[1:NPhases])]				
	intensities = [0.0; p[NPhases+1:2*NPhases]; 0.0]		
	ti = cumsum(durations); ti[end] = T_horizon	#for rounding error

	# if before delay or after last phase -> no LockDown
	if t <= (ti[1] - phaseDelay) || t >= (ti[end] + phaseDelay)
		return 0.0
	end

	#check if we are close to any of the phase-change points ti
	m(j) = (intensities[j+1]-intensities[j]) / (2*phaseDelay)
	q(j) = intensities[j] - m(j)*(ti[j] - phaseDelay)
	for j in 1:NPhases+1
		if abs(t-ti[j]) < phaseDelay
			return m(j)*t + q(j)
		end
	end

	#if none of the above then return the intensity of the phase we are in 
	i = findfirst(t .<= ti)			
	if isnothing(i)
		printstyled("ERROR IN policy, t = $t, p = $p\n",color=:red)
	end
	return intensities[i]
end

function policy_cubic(p,NPhases,delay,phaseDelay,T_horizon,t)
	""" durations are d1,d2... 
		times are t1 = delay, t2 = delay+d1,	t3 = t1+d2... 
		intervals are [0,t1 = delay),[t1 = delay,t2),[t2,t3) ... 
		intensities are i1 = 0.0, i2, i3 ...
	"""
	durations = [delay; p[1:NPhases]; T_horizon-delay-sum(p[1:NPhases])]				
	intensities = [0.0; p[NPhases+1:2*NPhases]; 0.0]		
	ti = cumsum(durations); ti[end] = T_horizon	#for rounding error

	# if before delay or after last phase -> no LockDown
	if t <= (ti[1] - phaseDelay) || t >= (ti[end] + phaseDelay)
		return 0.0
	end

	#check if we are close to any of the phase-change points ti
	for j in 1:NPhases+1
		if abs(t-ti[j]) < phaseDelay
			A = [(ti[j]-phaseDelay)^3 (ti[j]-phaseDelay)^2 (ti[j]-phaseDelay)^1 1;
						(ti[j]+phaseDelay)^3 (ti[j]+phaseDelay)^2 (ti[j]+phaseDelay)^1 1;
		  				3*(ti[j]-phaseDelay)^2 2*(ti[j]-phaseDelay)^1 1 0;
		  				3*(ti[j]+phaseDelay)^2 2*(ti[j]+phaseDelay)^1 1 0]
			a,b,c,d = A \ [intensities[j], intensities[j+1],0,0]
			return a*t^3 + b*t^2 + c*t + d
		end
	end

	#if none of the above then return the intensity of the phase we are in 
	i = findfirst(t .<= ti)		
	if isnothing(i)
		printstyled("ERROR IN policy, t = $t, p = $p\n",color=:red)
	end
	return intensities[i]
end

function solveScenario(p_join)
	param = p_join[1]
    sol = DifferentialEquations.solve(param.ODEprob(p_join),
								param.ODEalg(param.ODEalgStiff()),
								saveat=param.ODEt_save,
								maxiters=1e6,
                                isoutofdomain=(u,p,t) -> any(x -> x < 0, u),
								callback = param.ODEcb,
                                abstol=param.ODEabstol,reltol=param.ODEreltol)
	return sol
end

function solveScenario(param::parameters)
	p_join = (param, param.LockDownPolicy, param.TestingPolicy)
	sol = solveScenario(p_join)
	return sol
end	


function initialCondition(N₀,N₀_infect)
	S = (N₀-N₀_infect) / N₀
	C = 0.0
	E1 = 0.0
	E2 = N₀_infect / N₀
	I1 = 0.0
	I2 = 0.0
	R = 0.0
	V1 = 0.0
	Vw = 0.0
	D = 0.0
	Dc = 0.0

	return u0 = [S,C,E1,E2,I1,I2,R,V1,Vw,D,Dc]
end

function condition_cut(u,t,integrator)
	param = integrator.p[1]	
	if all(u[param.ODEinfectIndex] .< 1/param.N₀)
		return true
	else 
		return false
	end
end

function affect_cut!(integrator)
	param = integrator.p[1]	
	for k in param.ODEinfectIndex	# = C,E1,E2,I1,I2
		if integrator.u[k]  < 1/param.N₀
			integrator.u[k] = 0.0
		end
	end
end




