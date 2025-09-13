function reset_condition(u, t, integrator)
    Θ, ϖ, dt = integrator.p.Θ, integrator.p.ϖ, integrator.dt
    return abs(u[1] - u[2]) >= Θ && rand() < ϖ * dt
end

function reset_affect!(integrator)
    integrator.u[1] = integrator.u[2]
end
