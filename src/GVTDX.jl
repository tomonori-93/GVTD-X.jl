"""
Author: Satoki Tsujino
Date: 2026/09/15
"""
module GVTDX

開発中
const lib = joinpath(@__DIR__, "..", "lib", "libGVTDX.so")

export Retrieval_velocity, Retrieve_velocity_GVTDX, Retrieve_velocity_GVTD, Retrieve_velocity_GBVTD

"""
"""
function Retrieval_velocity(
                    nrot::Int32,
                    ndiv::Int32,
                    nrdiv::Int32,
                    r_t::Vector{Float64},
                    theta_t::Vector{Float64}, 
                    rdiv_t::Vector{Float64},
                    lon_tc::Float64,
                    lat_tc::Float64,
                    usp::Float64,
                    vsp::Float64,
                    nthres_undef::Vector{Int32},
                    skip_min_t::Int32,
                    flag_GVTDX::Int32,
                    missing_value::Float64,
                    lon_rdr::Float64,
                    lat_rdr::Float64,
                    Vra_in::Array{Float64},
                    flag_datagap::Bool,
                    umd::Float64,
                    vmd::Float64)

    n, m, l = size(Vra_in)

    @assert stride(Vra_in,1) == 1   # column-major保証

    VTtot      = Array{Float64}(undef,n,m,l)
    VRtot      = Array{Float64}(undef,n,m,l)
    VRT0       = Array{Float64}(undef,n,m,l)
    VDR0       = Array{Float64}(undef,n,m,l)
    Vra        = Array{Float64}(undef,n,m,l)
    Vra_ret    = Array{Float64}(undef,n,m,l)
    VRTn       = Array{Float64}(undef,nrot,n,m,l)
    VRRn       = Array{Float64}(undef,nrot,n,m,l)
    phin       = Array{Float64}(undef,nrot,n,m,l)
    zetan      = Array{Float64}(undef,nrot,n,m,l)
    VDTm       = Array{Float64}(undef,ndiv,n,m,l)
    VDRm       = Array{Float64}(undef,ndiv,n,m,l)
    Vn_0       = Array{Float64}(undef,n,m,l)
    Vra_Er     = Array{Float64}(undef,n,m,l)
    Vra_Er_ret = Array{Float64}(undef,n,m,l)
    VTtot_Er   = Array{Float64}(undef,n,m,l)
    VRtot_Er   = Array{Float64}(undef,n,m,l)
    Uxtot_Er   = Array{Float64}(undef,n,m,l)
    Vytot_Er   = Array{Float64}(undef,n,m,l)
    lond       = Matrix{Float64}(undef,n,m)
    latd       = Matrix{Float64}(undef,n,m)
    zeta0      = Array{Float64}(undef,n,m,l)
    zetatot    = Array{Float64}(undef,n,m,l)

    ccall((:c_Retrieval_velocity, lib),
          Cvoid,
          (Cint, Cint, Cint, Cint, Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Float64, Float64, Float64, Float64, Ptr{Int32}, Int32, Int32, Float64, Float64, Float64, 
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Bool, Ptr{Float64}, Ptr{Float64}),
          n, m, l, nrot, ndiv, nrdiv, r_t, theta_ref_t, rdiv_t, lon_tc, lat_tc,
          usp, vsp, nthres_undef, skip_min_t, flag_GVTDX, missing_value,
          lon_rdr, lat_rdr, Vra_in,
          VTtot, VRtot, VRT0, VDR0, Vra, Vra_ret, VRTn, VRRn, phin, zetan,
          VDTm, VDRm, Vn_0, Vra_Er, Vra_Er_ret, VTtot_Er, VRtot_Er, Uxtot_Er, Vytot_Er,
          lond, latd, zeta0, zetatot, flag_datagap, umd, vmd )

    return VTtot, VRtot, VRT0, VDR0, Vra, Vra_ret, VRTn, VRRn, phin, zetan,
           Vra_Er, Vra_ret, VTtot_Er, VRtot_Er, Uxtot_Er, Vytot_Er, lond, latd, zeta0, zetatot

end

function Retrieval_velocity_GVTDX(
                    nrot::Int32,
                    ndiv::Int32,
                    nrdiv::Int32,
                    r::Vector{Float64},
                    t::Vector{Float64},
                    rh::Vector{Float64},
                    td::Array{Float64},
                    rdiv::Vector{Float64},
                    Vd::Array{Float64},
                    Vn::Float64,
                    RadTC::Float64,
                    missing_value::Float64)

    n, m = size(Vd)

    @assert stride(Vd,1) == 1   # column-major保証

    VT         = Array{Float64}(undef,n,m)
    VR         = Array{Float64}(undef,n,m)
    VRT0       = Array{Float64}(undef,n,m)
    VDR0       = Array{Float64}(undef,n,m)
    VRTn       = Array{Float64}(undef,nrot,n,m)
    VRRn       = Array{Float64}(undef,nrot,n,m)
    phin       = Array{Float64}(undef,nrot,n,m)
    zetan      = Array{Float64}(undef,nrot,n,m)
    VDTm       = Array{Float64}(undef,ndiv,n,m)
    VDRm       = Array{Float64}(undef,ndiv,n,m)

    ccall((:c_Retrieval_velocity_GVTDX, lib),
          Cvoid,
          (Cint, Cint, Cint, Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Ptr{Float64}, Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64, 
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
          n, m, nrot, ndiv, nrdiv, r, t, rh, td, rdiv, 
          Vd, Vn, RadTC, VT, VR, VRT0, VDR0, 
          VRTn, VRRn, VDTm, VDRm, missing_value, 
          phin, zetan, VRT0_GVTD, VDR0_GVTD, Vn_0,
          VRTns_r, VRTnc_r, VRRns_r, VRRnc_r, zetans_r, zetanc_r)

    return VT, VR, VRT0, VDR0, VRTn, VRRn, VDTm, VDRm, phin, zetan

end

function Retrieval_velocity_GVTD(
                    nasym::Int32,
                    r::Vector{Float64},
                    t::Vector{Float64},
                    td::Array{Float64},
                    Vd::Array{Float64},
                    RadTC::Float64,
                    missing_value::Float64,
                    flag_datagap::Bool)

    n, m = size(Vd)

    @assert stride(Vd,1) == 1   # column-major保証

    VT         = Array{Float64}(undef,n,m)
    VR         = Array{Float64}(undef,n,m)
    VT0        = Array{Float64}(undef,n,m)
    VR0        = Array{Float64}(undef,n,m)
    VTSn       = Array{Float64}(undef,nasym,n,m)
    VTCn       = Array{Float64}(undef,nasym,n,m)

    ccall((:c_Retrieval_velocity_GVTD, lib),
          Cvoid,
          (Cint, Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64, 
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Float64, Bool),
          n, m, nasym, r, t, td, Vd, RadTC,
          VT, VR, VT0, VR0, VTSn, VTCn, missing_value )

    return VT, VR, VT0, VR0, VTSn, VTCn

end

function Retrieval_velocity_GBVTD(
                    nasym::Int32,
                    r::Vector{Float64},
                    t::Vector{Float64},
                    td::Array{Float64},
                    Vd::Array{Float64},
                    RadTC::Float64,
                    missing_value::Float64,
                    flag_datagap::Bool)

    n, m = size(Vd)

    @assert stride(Vd,1) == 1   # column-major保証

    VT         = Array{Float64}(undef,n,m)
    VR         = Array{Float64}(undef,n,m)
    VT0        = Array{Float64}(undef,n,m)
    VR0        = Array{Float64}(undef,n,m)
    VTSn       = Array{Float64}(undef,nasym,n,m)
    VTCn       = Array{Float64}(undef,nasym,n,m)

    ccall((:c_Retrieval_velocity_GBVTD, lib),
          Cvoid,
          (Cint, Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64, 
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Float64, Bool),
          n, m, nasym, r, t, td, Vd, RadTC,
          VT, VR, VT0, VR0, VTSn, VTCn, missing_value )

    return VT, VR, VT0, VR0, VTSn, VTCn

end

end
