@inline function solve_mass_parameter(v::T, power::Int) where {T}
    power == 0 && return one(T)
    if v <= zero(T)
        return zero(T)
    elseif v >= one(T)
        return one(T)
    end

    lo = zero(T)
    hi = one(T)
    u = clamp(v, T(1e-12), one(T) - T(1e-12))
    for _ = 1:80
        up = u^power
        f = (power + one(T)) * up - power * up * u - v
        abs(f) < T(32) * eps(T) && return clamp(u, lo, hi)
        df = power * (power + one(T)) * u^(power - 1) * (one(T) - u)
        candidate = ifelse(iszero(df), (lo + hi) / 2, u - f / df)
        if f > zero(T)
            hi = u
        else
            lo = u
        end
        if !(lo < candidate < hi)
            candidate = (lo + hi) / 2
        end
        u = candidate
    end
    return (lo + hi) / 2
end
