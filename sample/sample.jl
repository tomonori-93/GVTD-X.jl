include("../src/GVTDX.jl")
include("./wrap_netcdf_GVTDX.jl")

# For Himawari-8/9 satellite, a sample script of parallax correction
using NCDatasets  # Reading NetCDF
using DelimitedFiles
using DataStructures: OrderedDict  # Creating NetCDF attributes
using .GVTDX  # GVTDX retrieval module
using .wrap_netcdf_GVTDX

#-- setting section
infile_list = "input_files.txt"  # 1st column: NetCDF, 2nd column: Temperature vs height profile data (ASCII)
nrot = 3  # azimuthal wavenumber for retrieved rotating wind component
read_vname = "dealiased_velocity"
rad_vname = "radius"
azi_vname = "azimuth"
alt_vname = "sweep"
tim_vname = "sweep_time"
_fillval_name = "_FillValue"
lon_tc_name = "lon_tc"
lat_tc_name = "lat_tc"
lon_rdr_name = "lon_radar"
lat_rdr_name = "lat_radar"
usp_name = "mv_lon"
vsp_name = "mv_lat"
#-----------------------

d2r = π/180.0
r2d = 180.0/π

# Read infile_list
infiles = readlines(infile_list)
nl = size(infiles)[1]

# Start loop to read and perform parallax correction
for i in 1:nl
    local infile_nc = infiles[i]

    # Input NetCDF satellite data
    local al = NCDataset(infile_nc,"r")
    println("Read : $infile_nc")

    # Set arrays and parameters for reading Himawari data
    lon_tc = Float64(al[lon_tc_name].var[1])  # degrees
    lat_tc = Float64(al[lat_tc_name].var[1])  # degrees
    lon_rdr = Float64(al[lon_rdr_name].var[1])  # degrees
    lat_rdr = Float64(al[lat_rdr_name].var[1])  # degrees
    usp = Float64(al[usp_name].var[1])  # m/s
    vsp = Float64(al[vsp_name].var[1])  # m/s

    # In NCDatasets, _FillValue is forced to missing type in JuliaLang, so if set Float, you need it as follows.
    missing_value = Float64(al[read_vname].attrib[_fillval_name])
    nr, nt, nz = size(al[read_vname].var)
    Vra_in = fill(missing_value,nr,nt,nz)  # NOTE: x,y
    r_t = fill(convert(Float64,missing_value),nr)
    theta_t = fill(convert(Float64,missing_value),nt)
    alt = fill(convert(Float64,missing_value),nz)

    # Substitute NetCDF data into ds_{val,lon,lat}
    Vra_in  = map(x -> convert(Float64,x), al[read_vname].var)
    r_t     = map(x -> convert(Float64,x), al[rad_vname].var)
    theta_t = map(x -> convert(Float64,x), al[azi_vname].var)
    alt     = map(x -> convert(Float64,x), al[alt_vname].var)

    # Set the same values as moving velocity in mean wind velocity
    umd = fill(convert(Float64,usp),nz)
    vmd = fill(convert(Float64,vsp),nz)

    nthres_undef = [div(nt,2); div(nt,2)]
    skip_min_t = div(nt,2)
    flag_GVTDX = 1  # 1: GVTDX, 2: GVTD, 3: GBVTD
    flag_datagap = false  # Not used in GVTDX, used in GVTD and GBVTD

#println("check, ", nthres_undef, skip_min_t, Vra_in[div(nr,2),div(nt,2),2])
    # Retrieve wind components from Doppler velocity
    VTtot, VRtot, VRT0, VDR0, Vra, Vra_ret, VRTn, VRRn,
#    VTtot, VRtot, VRT0, VDR0, Vra, Vra_ret, VRTn, VRRn, phin, zetan,
           Vra_Er, Vra_ret, VTtot_Er, VRtot_Er, Uxtot_Er, Vytot_Er, 
           lond, latd = Retrieve_velocity(
#           lond, latd, zeta0, zetatot = Retrieve_velocity(
                    nrot, r_t, theta_t.*d2r, lon_tc, lat_tc,
                    usp, vsp, nthres_undef, skip_min_t, flag_GVTDX,
                    missing_value, lon_rdr, lat_rdr, Vra_in,
                    flag_datagap, umd, vmd
    )

    VTtot    = Array{Real}(VTtot)
    VRtot    = Array{Real}(VRtot)
    VRT0     = Array{Real}(VRT0)
    VDR0     = Array{Real}(VDR0)
    VRTn     = Array{Real}(VRTn)
    VRRn     = Array{Real}(VRRn)
    Vra      = Array{Real}(Vra)
    Vra_ret  = Array{Real}(Vra_ret)
    VTtot_Er = Array{Real}(VTtot_Er)
    VRtot_Er = Array{Real}(VRtot_Er)
    Uxtot_Er = Array{Real}(Uxtot_Er)
    Vytot_Er = Array{Real}(Vytot_Er)

    outfile = chopsuffix(infile_nc,".nc") * ".GVTDX.nc"
    ncdump_GVTDX(outfile,String(infile_nc),rad_vname,azi_vname,alt_vname,tim_vname,read_vname,nrot,
#                 r_t,theta_t,
                 VTtot,VRtot,VRT0,VDR0,VRTn,VRRn,
                 Vra,Vra_ret,VTtot_Er,VRtot_Er,Uxtot_Er,Vytot_Er)

    println("Output file: $outfile ...")

end

