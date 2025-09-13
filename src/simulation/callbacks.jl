function reset_condition(u, t, integrator)
    u[1] - u[2]
end

function reset_affect!(integrator)
    Θ, ϖ = integrator.p.Θ, integrator.p.ϖ
    u = integrator.u[1] - integrator.u[2]
    if abs(u) >= Θ && rand() < ϖ * integrator.dt
        integrator.u[1] = integrator.u[2]
    end
end
