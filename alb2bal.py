# Script to interpolate from ANHA4 grid onto nemo-nordic grid

import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import griddata
import netCDF4 as nc
from netCDF4 import Dataset

# files from U Alberta
alb = xr.open_dataset('ANHA4-EPM151_y2015m01d05_gridT.nc')
#print(alb)
alb_lat = alb.nav_lat_grid_T.values
#alb_lat = alb_lat[:,0]
alb_lon = alb.nav_lon_grid_T.values
#alb_lon = alb_lon[0,:]
alb_lev = alb.deptht.values
alb_dic = alb.vodic.values
alb_dic = alb_dic[0,:,:,:]
alb_alk = alb.voalk.values
alb_alk = alb_alk[0,:,:,:]
alb_oxy = alb.vooxy.values
alb_oxy = alb_oxy[0,:,:,:]
alb_fe  = alb.vofed.values
alb_fe  = alb_fe[0,:,:,:]
alb_po4 = alb.vopo4.values
alb_po4 = alb_po4[0,:,:,:]
alb_dop = alb.vodop.values
alb_dop = alb_dop[0,:,:,:]

# loop over sub-surface
#for k in range(1,30):
    # set zero values equal to value in cell above 
#    alb_dic_lyr = alb_dic[k,:,:]
#    alb_dic_lyr[alb_dic_lyr==0] = alb_dic[k-1,:,:]
#    alb_dic[k,:,:] = alb_dic_lyr

#for k in range(1,30):
    # interpolate to depth 
#    alb_dic_lev = np.interp

# read nordic domain
bal = xr.open_dataset('nordic_domain_cfg.nc')
#print(bal)
model_lat = bal.nav_lat.values
bal_lat = model_lat[:,0]
model_lon = bal.nav_lon.values
bal_lon = model_lon[0,:]
bal_lev = bal.nav_lev.values
bal_sea = bal.top_level.values
bal_sea = bal_sea[0,:,:]
bot_lev = bal.bottom_level.values
bot_lev = bot_lev[0,:,:]
bal_y, bal_x = np.meshgrid(bal_lat, bal_lon)

bal_oxy = np.zeros((56,1046,1238))
bal_oxy_tmp2 = np.zeros((50,1046,1238))

for k in range(0,35):
    print(k)
    alb_y = np.ravel(alb_lat)
    alb_x = np.ravel(alb_lon)
    alb_c = np.ravel(alb_dic[k,:,:])
    alb_a = np.ravel(alb_alk[k,:,:])
    alb_o = np.ravel(alb_oxy[k,:,:])
    alb_f = np.ravel(alb_fe[k,:,:])
    alb_p = np.ravel(alb_po4[k,:,:])
    alb_d = np.ravel(alb_dop[k,:,:])

    alb_y = alb_y[alb_c!=0]
    alb_x = alb_x[alb_c!=0]
    alb_c = alb_c[alb_c!=0]
    alb_a = alb_a[alb_a!=0]
    alb_o = alb_o[alb_o!=0]
    alb_f = alb_f[alb_f!=0]
    alb_p = alb_p[alb_p!=0]
    alb_d = alb_d[alb_d!=0]

    bal_oxy_tmp1 = griddata((alb_x, alb_y), alb_o, (bal_x, bal_y), method='nearest')
    bal_oxy_tmp1 = np.transpose(bal_oxy_tmp1)
    bal_oxy_tmp2[k,:,:] = bal_oxy_tmp1
    
for j in range(0,1046):
    for i in range(0,1238):
        bal_oxy[:,j,i] = griddata(alb_lev[0:35], bal_oxy_tmp2[0:35,j,i], bal_lev, method='linear')

for k in range(0,56):
    net_lev = bot_lev - k
    bal_oxy_lyr = bal_oxy[k,:,:]
    bal_oxy_lyr[net_lev<=0] = 0
    bal_oxy[k,:,:] = bal_oxy_lyr

try: ncfile.close()
except: pass
ncfile = Dataset('./o2_init.nc', mode='w')
print(ncfile)

depth_dim = ncfile.createDimension('depth', 56)     # depth axis
lat_dim = ncfile.createDimension('latD', 1046)     # latitude axis
lon_dim = ncfile.createDimension('lonD', 1238)    # longitude axis
for dim in ncfile.dimensions.items():
    print(dim)

depth = ncfile.createVariable('depth', np.float32, ('depth',))
depth.units = 'm'
depth.long_name = 'depth'
lat = ncfile.createVariable('lat', np.float32, ('latD','lonD'))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lon = ncfile.createVariable('lon', np.float32, ('latD','lonD'))
lon.units = 'degrees_east'
lon.long_name = 'longitude'

oxy = ncfile.createVariable('oxy',np.float64,('depth','latD','lonD'))
oxy.units = 'mol/m³'
oxy.long_name = 'oxygen_concentration'

depth[:]     = bal_lev
lat[:,:]     = model_lat[0:1046,0:1238]
lon[:,:]     = model_lon[0:1046,0:1238]
oxy[:,:,:]   = bal_oxy

print(ncfile)
ncfile.close(); print('Dataset is closed')

#bal_po4 = np.zeros(1046,1238)

#bal_dic = griddata((alb_x, alb_y), alb_c, (bal_x, bal_y), method='nearest')
#bal_dic = np.transpose(bal_dic)
#bal_dic[bal_sea==0]=np.nan
#fig=plt.figure(figsize=(20,20)); pcol = plt.pcolormesh(bal_dic); plt.clim(1,3); cbar=plt.colorbar(pcol); plt.savefig('interp_carbon')

#bal_alk = griddata((alb_x, alb_y), alb_a, (bal_x, bal_y), method='nearest')
#bal_alk = np.transpose(bal_alk)
#bal_alk[bal_sea==0]=np.nan
#fig=plt.figure(figsize=(20,20)); pcol = plt.pcolormesh(bal_alk); plt.clim(1,3); cbar=plt.colorbar(pcol); plt.savefig('interp_alkalinity')

#bal_oxy = griddata((alb_x, alb_y), alb_o, (bal_x, bal_y), method='nearest')
#bal_oxy = np.transpose(bal_oxy)
#bal_oxy[bal_sea==0]=np.nan
fig=plt.figure(figsize=(20,20)); pcol = plt.pcolormesh(bal_oxy[0,:,:]); plt.clim(0,0.4); cbar=plt.colorbar(pcol); plt.show()
fig=plt.figure(figsize=(20,20)); pcol = plt.pcolormesh(bal_oxy[23,:,:]); plt.clim(0,0.4); cbar=plt.colorbar(pcol); plt.show()

#bal_fe = griddata((alb_x, alb_y), alb_f, (bal_x, bal_y), method='nearest')
#bal_fe = np.transpose(bal_fe)
#bal_fe[bal_sea==0]=np.nan
#fig=plt.figure(figsize=(20,20)); pcol = plt.pcolormesh(bal_fe); plt.clim(0,2e-6); cbar=plt.colorbar(pcol); plt.savefig('interp_iron')

#bal_po4 = griddata((alb_x, alb_y), alb_p, (bal_x, bal_y), method='nearest')
#bal_po4 = np.transpose(bal_po4)
#bal_po4[bal_sea==0]=np.nan
#fig=plt.figure(figsize=(20,20)); pcol = plt.pcolormesh(bal_po4); plt.clim(0,0.001); cbar=plt.colorbar(pcol); plt.savefig('interp_phosphate')

#bal_dop = griddata((alb_x, alb_y), alb_d, (bal_x, bal_y), method='nearest')
#bal_dop = np.transpose(bal_dop)
#bal_dop[bal_sea==0]=np.nan
#fig=plt.figure(figsize=(20,20)); pcol = plt.pcolormesh(bal_dop); plt.clim(0,0.001); cbar=plt.colorbar(pcol); plt.savefig('interp_dop')

