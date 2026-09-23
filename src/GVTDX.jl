"""
Author: Satoki Tsujino
Date: 2026/09/15
"""
module GVTDX

開発中
const lib = joinpath(@__DIR__, "..", "lib", "libGVTDX.so")

export GVTDX_main, Retrieve_GVTDX, Retrieve_GVTD, Retrieve_GBVTD

"""
"""
function Retrieval_control(
    integer(c_int), value :: nrot              !! the rotating maximum wavenumber used in the retrieval
    integer(c_int), value :: ndiv              !! the divergent maximum wavenumber used in the retrieval
    Float64        :: r_t(n)     !! radial coordinate on which Vra_in is defined [m]
    Float64        :: theta_ref_t(m)     !! azimuthal coordinate on which Vra_in is defined [rad]
    Float64, value :: lon_tc   !! longitude of the TC center [degree]
    Float64, value :: lat_tc   !! latitude of the TC center [degree]
    Float64, value :: usp      !! zonal component of the moving velocity of the TC [m/s]
    Float64, value :: vsp      !! meridional component of the moving velocity of the TC [m/s]
    integer(c_int)        :: nthres_undef(2)   !! thresholds of the azimuthal sampling number to determine the (1) innermost and (2) outermost radii, respectively
    integer(c_int), value :: skip_min_t        !! threshold of the azimuthal sampling number to determine unused radii
    integer(c_int), value :: flag_GVTDX        !! 1: GVTD-X, 2: GVTD, 3: GBVTD
    real,           value :: missing_value             !! undefined value for the original Doppler radar data
    Float64, value :: lon_rdr !! longitudinal position of the radar (degree)
    Float64, value :: lat_rdr !! latitudinal position of the radar (degree)
    real                  :: Vra_in(n,m,nz)  !! input Doppler velocity [m/s]
    logical(c_bool),       optional :: flag_datagap !! Flag for use of optimal wavenumber from data gap (Lee et al. 2000)
    Float64,        optional :: umd(nz)  !! zonal component of mean wind [m s-1]
    Float64,        optional :: vmd(nz)  !! meridional component of mean wind [m s-1]
    integer(c_int), value, optional :: nrdiv             !! radial grid number where the divergence is defined
    Float64,        optional :: rdiv_t(nrdiv)  !! radial grids where the divergence is defined
                    nrot::Int32,
                    ndiv::Int32,
                    r_t::Vector{Float64},
                    theta_t::Vector{Float64}, 
                    lon_tc::Float64,
                    lat_tc::Float64,
                    usp::Float64,
                    vsp::Float64,
                    nthres_undef::Vector{Int32},
                    h_pix::Matrix{Float64},
                    missing_value::Float64)

    n, m, l = size(Vra_in)

    @assert stride(Vra_in,1) == 1   # column-major保証

    VTtot      = Array{Float64}(undef,n,m,l)  !! retrieved total tangential wind [m s-1]
    VRtot      = Array{Float64}(undef,n,m,l)  !! retrieved total radial wind [m s-1]
    VRT0       = Array{Float64}(undef,n,m,l)  !! retrieved axisymmetric tangential wind [m s-1]
    VDR0       = Array{Float64}(undef,n,m,l)  !! retrieved axisymmetric radial wind [m s-1]
    Vra        = Array{Float64}(undef,n,m,l)      !! storm-relative input Doppler velocity [m s-1]
    Vra_ret    = Array{Float64}(undef,n,m,l)  !! storm-relative retrieved Doppler velocity [m s-1]
    VRTn       = Array{Float64}(undef,nrot,n,m,l)  !! retrieved wavenumber-N rotational-tangential wind [m s-1]
    VRRn       = Array{Float64}(undef,nrot,n,m,l)  !! retrieved wavenumber-N rotational-radial wind [m s-1]
    phin       = Array{Float64}(undef,nrot,n,m,l)  !! retrieved wavenumber-N streamfunction [m2 s-1]
    zetan      = Array{Float64}(undef,nrot,n,m,l)  !! retrieved wavenumber-N vorticity [s-1]
    VDTm       = Array{Float64}(undef,ndiv,n,m,l)  !! retrieved wavenumber-M divergent-tangential wind [m s-1]
    VDRm       = Array{Float64}(undef,ndiv,n,m,l)  !! retrieved wavenumber-M divergent-radial wind [m s-1]
    Vn_0       = Array{Float64}(undef,n,m,l)  !! storm-relative mean wind normal to line of sight [m s-1]
    Vra_Er     = Array{Float64}(undef,n,m,l)  !! earth-relative Doppler velocity [m s-1]
    Vra_Er_ret = Array{Float64}(undef,n,m,l)  !! earth-relative retrieved Doppler velocity [m s-1]
    VTtot_Er   = Array{Float64}(undef,n,m,l)  !! earth-relative retrieved total tangential wind [m s-1]
    VRtot_Er   = Array{Float64}(undef,n,m,l)  !! earth-relative retrieved total radial wind [m s-1]
    Uxtot_Er   = Array{Float64}(undef,n,m,l)  !! earth-relative retrieved total zonal wind [m s-1]
    Vytot_Er   = Array{Float64}(undef,n,m,l)  !! earth-relative retrieved total meridional wind [m s-1]
    lond       = Matrix{Float64}(undef,n,m)  !! Longitude [degree]
    latd       = Matrix{Float64}(undef,n,m)  !! Latitude [degree]
    zeta0      = Array{Float64}(undef,n,m,l)  !! retrieved axisymmetric vorticity [s-1]
    zetatot    = Array{Float64}(undef,n,m,l)  !! retrieved total vorticity [s-1]

    ccall((:c_Retrieval_control, lib),
          Cvoid,
          (Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Ptr{Float64}, Ptr{Float64}, Float64, Float64, Float64, 
           Float64, Float64, Float64),
          n, m, lon_pix, lat_pix, h_pix, lon_cor, lat_cor,
          re, rp, hsat, psat, lsat, missing_value)

    return VTtot, VRtot, VRT0, VDR0, Vra, Vra_ret, VRTn, VRRn, phin, zetan,
           Vra_Er, Vra_ret, VTtot_Er, VRtot_Er, Uxtot_Er, Vytot_Er, lond, latd, zeta0, zetatot

end

function Retrieve_GVTDX(x_in::Matrix{Float64},
                        y_in::Matrix{Float64}, 
                        iv::Matrix{Float64},
                        ivad::Matrix{Float64},
                        x_out::Vector{Float64},
                        y_out::Vector{Float64},
                        missing_value::Float64)

    n, m = size(x_in)
    l = size(x_out)[1]
    k = size(y_out)[1]

    @assert stride(x_in,1) == 1   # column-major保証

    ov = Matrix{Float64}(undef,l,k)
    ovad = Matrix{Float64}(undef,l,k)

    ccall((:c_tri_interpolation_2d, lib),
          Cvoid,
          (Cint, Cint, Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Float64),
          n, m, l, k, x_in, y_in, iv, ivad, x_out, y_out, ov, ovad, missing_value)

    return ov, ovad

end

function convert_Tbb2Zph(tval::Matrix{Float64},
                         t1d::Vector{Float64},
                         z1d::Vector{Float64},
                         missing_value::Float64)

    n, m = size(tval)
    l = size(t1d)[1]

    @assert stride(tval,1) == 1   # column-major保証

    zval = Matrix{Float64}(undef,n,m)

    ccall((:c_convert_Tbb2Zph, lib),
          Cvoid,
          (Cint, Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Float64),
          n, m, l, tval, zval, t1d, z1d, missing_value)

    return zval

end

end
