"""
Author: Satoki Tsujino
Date: 2026/09/14
"""
module wrap_netcdf_GVTDX  # sub module for handling NetCDF
# For Himawari-8/9 satellite, wrapper library
using NCDatasets  # Reading NetCDF
using DataStructures: OrderedDict  # Creating NetCDF attributes

export ncdump_GVTDX

#-- Define function to create NetCDF for VTTrac
"""
   ncdump_GVTDX(fname_nc,ref_fname_nc,n_grid,rad_grid,azi_grid,gph_grid,u_grid,v_grid,score_grid,_fillvalue)

Create a NetCDF file and store the VTTrac results

# Arguments
- `fname_nc::String` : Created NetCDF file name
- `ref_fname_nc::String` : Input NetCDF file name (used to copy global attributes and obs_time)
- `rad_grid::Array{Real,1}` : [nx_grid] longitude [degree]
- `azi_grid::Array{Real,1}` : [ny_grid] latitude [degree]
- `gph_grid::Array{Real,2}` : [nx_grid,ny_grid] height [m] at each rad_grid/azi_grid
- `val_grid::Array{Real,2}` : [nx_grid,ny_grid] variable [unit] at each rad_grid/azi_grid
- `read_vname::String` : varname to be output
"""
function ncdump_GVTDX(fname_nc::String,ref_fname_nc::String,
                      rad_name::String,azi_name::String,alt_name::String,tim_name::String,read_vname::String,nrot::Int,
#                      rad_grid::Array{Real,1},azi_grid::Array{Real,1},alt_grid::Array{Real,1},
                      VTtot::Array{Real,3},VRtot::Array{Real,3},
                      VRT0::Array{Real,3},VDR0::Array{Real,3},
                      VRTn::Array{Real,4},VRRn::Array{Real,4},
#                      VDTm::Array{Real,4},VDRm::Array{Real,4},
                      Vra::Array{Real,3},Vra_ret::Array{Real,3},
                      VTtot_Er::Array{Real,3},VRtot_Er::Array{Real,3}, 
                      Uxtot_Er::Array{Real,3},Vytot_Er::Array{Real,3})

    ds_old = NCDataset(ref_fname_nc,"r")
    #fname_nc_new = fname_nc[1:end-3] * ".GVTDX.nc"
    ds_new = NCDataset(fname_nc,"c")

    rad_grid = ds_old[rad_name].var
    azi_grid = ds_old[azi_name].var
    alt_grid = ds_old[alt_name].var
    #drad_grid = rad_grid[2] - rad_grid[1]
    #dazi_grid = azi_grid[2] - azi_grid[1]

    # Define the dimensions.
    defDim(ds_new,rad_name,size(rad_grid[:])[1])
    defDim(ds_new,azi_name,size(azi_grid[:])[1])
    defDim(ds_new,alt_name,size(alt_grid[:])[1])

    # Define global attributes
    ds_new.attrib["title"] = "GVTD-X-retrieved winds"
    ds_new.attrib["source"] = "Satoki Tsujino"
    ds_new.attrib["institution"] = "Meteorological Research Institute"
    _fillvalue = ds_old[read_vname].attrib["_FillValue"]

    # Define the variables with the attribute units
    #println(typeof(ds_old["obs_time"].attrib["units"]))
    #println(typeof(ds_old["obs_time"].attrib["long_name"]))
    rad = defVar(ds_new,rad_name,Float32,(rad_name,), attrib = OrderedDict(
            "units" => ds_old[rad_name].attrib["units"], 
            "long_name" => ds_old[rad_name].attrib["long_name"]
            ))
    azi = defVar(ds_new,azi_name,Float32,(azi_name,), attrib = OrderedDict(
            "units" => ds_old[azi_name].attrib["units"], 
            "long_name" => ds_old[azi_name].attrib["long_name"]
            ))
    alt = defVar(ds_new,alt_name,Float32,(alt_name,), attrib = OrderedDict(
            "units" => ds_old[alt_name].attrib["units"], 
            "long_name" => ds_old[alt_name].attrib["long_name"]
            ))
    rad[:] = rad_grid
    azi[:] = azi_grid
    alt[:] = alt_grid

    tim = defVar(ds_new,tim_name,Float32,(alt_name,), attrib = OrderedDict(
            "units" => ds_old[tim_name].attrib["units"], 
            "long_name" => ds_old[tim_name].attrib["long_name"]
            ))
    tim[:] = ds_old[tim_name].var

    VTtot_nc = defVar(ds_new,"VTtot",Float32,(rad_name,azi_name,alt_name), attrib = OrderedDict(
                 "units" => "m s-1", 
                 "long_name" => "retrieved total tangential wind",
                 "_FillValue" => Float32(_fillvalue),
                 ))
    VRtot_nc = defVar(ds_new,"VRtot",Float32,(rad_name,azi_name,alt_name), attrib = OrderedDict(
                 "units" => "m s-1", 
                 "long_name" => "retrieved total radial wind",
                 "_FillValue" => Float32(_fillvalue),
                 ))
    VRT0_nc = defVar(ds_new,"VRT0",Float32,(rad_name,azi_name,alt_name), attrib = OrderedDict(
                "units" => "m s-1", 
                "long_name" => "retrieved axisymmetric tangential wind",
                "_FillValue" => Float32(_fillvalue),
                ))
    VDR0_nc = defVar(ds_new,"VDR0",Float32,(rad_name,azi_name,alt_name), attrib = OrderedDict(
                "units" => "m s-1", 
                "long_name" => "retrieved axisymmetric radial wind",
                "_FillValue" => Float32(_fillvalue),
                ))
    Vra_org_nc = defVar(ds_new,"Vra",Float32,(rad_name,azi_name,alt_name), attrib = OrderedDict(
                   "units" => "m s-1", 
                   "long_name" => "Storm-relative Doppler velocity",
                   "_FillValue" => Float32(_fillvalue),
                   ))
    Vra_ret_nc = defVar(ds_new,"Vra_ret",Float32,(rad_name,azi_name,alt_name), attrib = OrderedDict(
                   "units" => "m s-1", 
                   "long_name" => "Storm-relative retrieved Doppler velocity",
                   "_FillValue" => Float32(_fillvalue),
                   ))
    VTtot_Er_nc = defVar(ds_new,"VTtot_Er",Float32,(rad_name,azi_name,alt_name), attrib = OrderedDict(
                    "units" => "m s-1", 
                    "long_name" => "Earth-relative total tangential wind",
                    "_FillValue" => Float32(_fillvalue),
                    ))
    VRtot_Er_nc = defVar(ds_new,"VRtot_Er",Float32,(rad_name,azi_name,alt_name), attrib = OrderedDict(
                    "units" => "m s-1", 
                    "long_name" => "Earth-relative total radial wind",
                    "_FillValue" => Float32(_fillvalue),
                    ))
    Uxtot_Er_nc = defVar(ds_new,"Uxtot_Er",Float32,(rad_name,azi_name,alt_name), attrib = OrderedDict(
                    "units" => "m s-1", 
                    "long_name" => "Earth-relative total zonal wind",
                    "_FillValue" => Float32(_fillvalue),
                    ))
    Vytot_Er_nc = defVar(ds_new,"Vytot_Er",Float32,(rad_name,azi_name,alt_name), attrib = OrderedDict(
                    "units" => "m s-1", 
                    "long_name" => "Earth-relative total meridional wind",
                    "_FillValue" => Float32(_fillvalue),
                    ))
    VTtot_nc[:,:,:] = VTtot
    VRtot_nc[:,:,:] = VRtot
    VRT0_nc[:,:,:] = VRT0
    VDR0_nc[:,:,:] = VDR0
    Vra_org_nc[:,:,:] = Vra
    Vra_ret_nc[:,:,:] = Vra_ret
    VTtot_Er_nc[:,:,:] = VTtot_Er
    VRtot_Er_nc[:,:,:] = VRtot_Er
    Uxtot_Er_nc[:,:,:] = Uxtot_Er
    Vytot_Er_nc[:,:,:] = Vytot_Er

    if nrot > 0
       #VRTn_nc = Vector{String}(undef,nrot)
       #VRRn_nc = Vector{String}(undef,nrot)
       for i=1:nrot
          VRTn_nc = defVar(ds_new,"VRT$i",Float32,(rad_name,azi_name,alt_name), attrib = OrderedDict(
                       "units" => "m s-1", 
                       "long_name" => "retrieved wavenumber-$i rotational-tangential wind",
                       "_FillValue" => Float32(_fillvalue),
                    ))
          VRRn_nc = defVar(ds_new,"VRR$i",Float32,(rad_name,azi_name,alt_name), attrib = OrderedDict(
                       "units" => "m s-1", 
                       "long_name" => "retrieved wavenumber-$i rotational-radial wind",
                       "_FillValue" => Float32(_fillvalue),
                    ))
          VRTn_nc[:,:,:] = VRTn[i,:,:,:]
          VRRn_nc[:,:,:] = VRRn[i,:,:,:]
       end
    end

#    if ndiv > 0
#       VRTn_nc = Vector{String}(undef,ndiv)
#       VRRn_nc = Vector{String}(undef,ndiv)
#       for i=1:ndiv
#          VDTm_nc[i] = defVar(ds_new,"VDT"*string(i),Float32,(rad_name,azi_name,alt_name), attrib = OrderedDict(
#                          "units" => "m s-1", 
#                          "long_name" => "retrieved wavenumber-"*string(i)*" divergent-tangential wind",
#                          "_FillValue" => Float32(_fillvalue),
#                       ))
#          VDRm_nc[i] = defVar(ds_new,"VDR"*string(i),Float32,(rad_name,azi_name,alt_name), attrib = OrderedDict(
#                          "units" => "m s-1", 
#                          "long_name" => "retrieved wavenumber-"*string(i)*" divergent-radial wind",
#                          "_FillValue" => Float32(_fillvalue),
#                       ))
#          VDTm_nc[i].var[:,:,:] = VDTm_nc[:,:,:,i]
#          VDRm_nc[i].var[:,:,:] = VDRm_nc[:,:,:,i]
#       end
#    end

    # add additional attributes
    #v.attrib["comments"] = "this is a string attribute with Unicode Ω ∈ ∑ ∫ f(x) dx"

    close(ds_new)
end

end
