mutable struct GTS3C
    TLB::Float64
    S::Float64
    DB04::Float64
    DB16::Float64
    DB28::Float64
    DB32::Float64
    DB40::Float64
    DB48::Float64
    DB01::Float64
    ZA::Float64
    T0::Float64
    Z0::Float64
    G0::Float64
    RL::Float64
    DD::Float64
    DB14::Float64
    TR12::Float64
end

mutable struct MESO7
    TN1::MVector{5, Float64}
    TN2::MVector{4, Float64}
    TN3::MVector{5, Float64}
    TGN1::MVector{2, Float64}
    TGN2::MVector{2, Float64}
    TGN3::MVector{2, Float64}
end

mutable struct DMIX
    DM04::Float64
    DM16::Float64
    DM28::Float64
    DM32::Float64
    DM40::Float64
    DM01::Float64
    DM14::Float64
end

mutable struct LPOLY
    PLG::Matrix{Float64}
    CTLOC::Float64
    STLOC::Float64
    C2TLOC::Float64
    S2TLOC::Float64
    C3TLOC::Float64
    S3TLOC::Float64
    DAY::Float64
    DF::Float64
    DFA::Float64
    APD::Float64
    APDF::Float64
    APT::Vector{Float64}
    XLONG::Float64
end

mutable struct PARMB
    GSURF::Float64
    RE::Float64
end

#=
=#
