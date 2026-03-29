const Vec3{T} = SVector{3, T}

struct PhaseSpacePoint{T}
    momenta::Vector{FourVector{T}}
    weight::T
end

Base.length(point::PhaseSpacePoint) = length(point.momenta)
Base.getindex(point::PhaseSpacePoint, i::Int) = point.momenta[i]
Base.iterate(point::PhaseSpacePoint, state...) = iterate(point.momenta, state...)

phase_space_weight(point::PhaseSpacePoint) = point.weight

@inline spatial(p::FourVector{T}) where {T} = Vec3{T}(
    LorentzVectorBase.px(p),
    LorentzVectorBase.py(p),
    LorentzVectorBase.pz(p),
)

@inline function fourvector(v::Vec3{T}, energy::T) where {T}
    return FourVector(v[1], v[2], v[3]; E = energy)
end

@inline function add_fourvectors(a::FourVector{T}, b::FourVector{T}) where {T}
    return FourVector(
        LorentzVectorBase.px(a) + LorentzVectorBase.px(b),
        LorentzVectorBase.py(a) + LorentzVectorBase.py(b),
        LorentzVectorBase.pz(a) + LorentzVectorBase.pz(b);
        E = LorentzVectorBase.E(a) + LorentzVectorBase.E(b),
    )
end

function total_momentum(momenta)
    T = eltype(eltype(momenta))
    total = FourVector(zero(T), zero(T), zero(T); E = zero(T))
    for p in momenta
        total = add_fourvectors(total, p)
    end
    return total
end

@inline kallen(x, y, z) = x^2 + y^2 + z^2 - 2 * (x * y + y * z + z * x)

@inline function breakup_momentum(M, m1, m2)
    threshold = (m1 + m2)^2
    upper = (m1 - m2)^2
    λ = (M^2 - threshold) * (M^2 - upper)
    λ = λ < zero(λ) && !isapprox(λ, zero(λ); atol = sqrt(eps(real(float(M))))) ? λ : max(zero(λ), λ)
    return sqrt(λ) / (2 * M)
end

@inline function two_body_density(M, m1, m2)
    return breakup_momentum(M, m1, m2) / (4 * M)
end

@inline function massless_phase_space_volume(s, n::Integer)
    return (π / 2)^(n - 1) * s^(n - 2) / (factorial(big(n - 1)) * factorial(big(n - 2)))
end

@inline function boost(p::FourVector{T}, parent::FourVector{T}) where {T}
    parent_mass = LorentzVectorBase.mass(parent)
    if iszero(parent_mass) || isapprox(parent_mass, zero(T); atol = sqrt(eps(T)))
        return p
    end
    β = spatial(parent) / LorentzVectorBase.E(parent)
    β2 = dot(β, β)
    if iszero(β2)
        return p
    end
    γ = LorentzVectorBase.E(parent) / parent_mass
    pvec = spatial(p)
    βdotp = dot(β, pvec)
    factor = ((γ - one(T)) / β2) * βdotp + γ * LorentzVectorBase.E(p)
    boosted_vec = pvec + factor * β
    boosted_energy = γ * (LorentzVectorBase.E(p) + βdotp)
    return fourvector(boosted_vec, boosted_energy)
end

@inline function unit_direction(cosθ::T, ϕ::T) where {T}
    sinθ = sqrt(max(zero(T), one(T) - cosθ^2))
    sinϕ, cosϕ = sincos(ϕ)
    return Vec3{T}(sinθ * cosϕ, sinθ * sinϕ, cosθ)
end

function mandelstam(momenta, i::Int, j::Int)
    p = add_fourvectors(momenta[i], momenta[j])
    return LorentzVectorBase.mass2(p)
end

function invariant_masses(momenta)
    @assert length(momenta) == 3 "Dalitz invariants are defined here for 3-body final states."
    return (
        s23 = mandelstam(momenta, 2, 3),
        s31 = mandelstam(momenta, 3, 1),
        s12 = mandelstam(momenta, 1, 2),
    )
end
