include("../src/GVTDX.jl")
include("../src/wrap_netcdf_GVTDX.jl")

# For Himawari-8/9 satellite, a sample script of parallax correction
using NCDatasets  # Reading NetCDF
using DelimitedFiles
using DataStructures: OrderedDict  # Creating NetCDF attributes
using .GVTDX  # GVTDX retrieval module
using .wrap_netcdf_GVTDX

#-- setting section
infile_list = "input_files.txt"  # 1st column: NetCDF, 2nd column: Temperature vs height profile data (ASCII)
hskip_zprof = 4  # array number to skip reading in the temp/height profile data
nrot = 3  # azimuthal wavenumber for retrieved rotating wind component
read_vname = "dealiased_velocity"
rad_vname = "radius"
azi_vname = "azimuth"
alt_vname = "altitude"
_fillval_name = "_FillValue"
lon_tc_name = "lon_tc"
lat_tc_name = "lat_tc"
lon_rdr_name = "lon_rdr"
lat_rdr_name = "lat_rdr"
#-----------------------

d2r = π/180.0
r2d = 180.0/π

# Read infile_list
infiles = readlines(infile_list)
nl = size(infiles)[1]

# Start loop to read and perform parallax correction
for i in 1:nl
    local infile_nc = split(infiles[i])

    # Input NetCDF satellite data
    local al = NCDataset(infile_nc,"r")
    println("Read : $infile_nc")

    # Set arrays and parameters for reading Himawari data
    lon_tc = al[lon_tc_name].var  # degrees
    lat_tc = al[lat_tc_name].var  # degrees
    lon_rdr = al[lon_rdr_name].var  # degrees
    lat_rdr = al[lat_rdr_name].var  # degrees
    # In NCDatasets, _FillValue is forced to missing type in JuliaLang, so if set Float, you need it as follows.
    missing_value = al[read_vname].attrib[_fillval_name]
    nr, nt, nz = size(al[read_vname].var)
    ds_val = fill(missing_value,nr,nt)  # NOTE: x,y
    ds_rad = fill(convert(Float64,missing_value),nr)
    ds_azi = fill(convert(Float64,missing_value),nt)
    ds_alt = fill(convert(Float64,missing_value),nz)

    # Substitute NetCDF data into ds_{val,lon,lat}
    ds_val = map(x -> convert(Float64,x), al[read_vname].var)
    ds_rad = map(x -> convert(Float64,x), al[rad_vname].var)
    ds_azi = map(x -> convert(Float64,x), al[azi_vname].var)
    ds_alt = map(x -> convert(Float64,x), al[alt_vname].var)

    # Retrieve wind components from Doppler velocity
    VTtot, VRtot, VRT0, VDR0, VRTn, VRRn = Retrieval_control(
        aaa開発中
    )
    gph_grid, val_grid = ParallaxCorrect(ds_lon.*d2r, ds_lat.*d2r, ds_gph, ds_val,
                            lon_grid.*d2r, lat_grid.*d2r, R_e, R_p, h_sat, p_sat, l_sat, missing_value)

    lon_grid = Array{Real}(lon_grid)
    lat_grid = Array{Real}(lat_grid)
    gph_grid = Array{Real}(gph_grid)
    val_grid = Array{Real}(val_grid)

    outfile = chopsuffix(infile_nc,".nc") * ".GVTDX.nc"
    ncdump_gvtdx(outfile,String(infile_nc),ds_rad,ds_azi,VTtot,VRtot,VRT0,VDR0,VRTn,VRRn,
                 Vra,Vra_ret,Vra_Er,Vra_Er_ret,missing_value)

    println("Output file: $outfile ...")

end

