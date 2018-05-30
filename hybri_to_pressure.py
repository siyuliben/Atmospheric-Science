####################### hybrid to pressure 1 ##################
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset,date2num,num2date
from mpl_toolkits.basemap import Basemap
import sys
sys.setrecursionlimit(1000000)

import csv
##################### reading filenames #####################
columns = [[] for _ in range(15)]

with open("Tracefilename.csv",'r') as f:
    c = csv.reader(f, delimiter=',', skipinitialspace=True)
    for row in c:
        for i, col in enumerate(row):
            columns[i].append(col)
PSS=np.asarray(columns[0][:])
GPH=np.asarray(columns[1][:])
GPH500=np.asarray(columns[2][:])
zonalwind=np.asarray(columns[3][:])
zonalwind300=np.asarray(columns[4][:])
### This CSV file has the list of file name that we want to input or output.
### Create this file first then run the script
### Each column have is different input variable name or output file name.
###################### reading filenames #####################

filenameindex=10 ############ number that you have to change to get the input file.
                 ############ TRACE data have 36 files. therefore the number is from 1 to 36.


##################### input file with file number.
dpath1='TRACE21/SOURCE/'
ncfile1=PSS[filenameindex-1]

###################### set pressure level ####
n=np.array([100000.0,92500.0,85000.0,70000.0,60000.0,50000.0,40000.0,\
            30000.0,25000.0,20000.0,15000.0,10000.0,7000.0])
pindex=3 #this index is to set the pressure level index in array 'n'. Here we set it to 3 which is 700hPa.

######################## this is input filename #######
dpath2='TRACE21/SOURCE/'
ncfile2=zonalwind[filenameindex-1]

######################## this is output filename #######
dpath3='elisontimmlab_rit/siyu/TRACE/U300/'
filename=zonalwind300[filenameindex-1]

####################### variable name
var='U'

########### one more thing that has to be changed if variable is not zonal wind.
########### At the every end. We have to add some description for each variable.
########### please go to the last few lines to write something about the variable you are going to select.
########### the line starts with 'w_nc_var.setncatts'.


################# now reading the hybrid pressure level.
print ncfile1
nc_f1 = dpath1+ncfile1 # Your filename
nc_fid1 = Dataset(nc_f1, )

lats = nc_fid1.variables['lat'][:]
lons = nc_fid1.variables['lon'][:]
time1 = nc_fid1.variables['time']
time = nc_fid1.variables['time'][:]
field = nc_fid1.variables['PS'][:]
time_units = time1.units

testilev=nc_fid1.variables['ilev'][:]
testhybi=nc_fid1.variables['hybi'][:]
testhyai=nc_fid1.variables['hyai'][:]
testP0=nc_fid1.variables['P0'][:]

bigarray = np.zeros([testilev.shape[0],time.shape[0],lats.shape[0],lons.shape[0]])

i=0
j=0
k=0
l=0


while i < testilev.shape[0]:
    while j < time.shape[0]:
        while k < lats.shape[0]:
            while l < lons.shape[0]:
                PP = testhyai[i]*testP0+testhybi[i]*field[j,k,l]
                bigarray[i,j,k,l]=PP
                l=l+1
            k=k+1;l=0
        j=j+1;k=0
    i=i+1;j=0

############################ hybrid to pressure 2 ############################
### conversion from hybrid-pressure to stander-pressure

nc_f1 = dpath2+ncfile2 # Your filename
nc_fid1 = Dataset(nc_f1, )

nc_f1 = dpath1+ncfile1 # Your filename
nc_fid1 = Dataset(nc_f1, )

lats = nc_fid1.variables['lat'][:]
lons = nc_fid1.variables['lon'][:]
time1 = nc_fid1.variables['time']
time = nc_fid1.variables['time'][:]
field = nc_fid1.variables[var][:]
time_units = time1.units

testilev=nc_fid1.variables['ilev'][:]
testhybi=nc_fid1.variables['hybi'][:]
testhyai=nc_fid1.variables['hyai'][:] 
testP0=nc_fid1.variables['P0'][:]



final = np.zeros([time.shape[0],lats.shape[0],lons.shape[0]])

i=0
j=0
k=0
l=0


while j < time.shape[0]:
    while k < lats.shape[0]:
        while l < lons.shape[0]:
            a,b=bigarray[:,j,k,l][np.argsort(np.abs(bigarray[:,j,k,l]-n[pindex]))[0:2]]
            Ia=np.where(bigarray[:,j,k,l]==a)[0][0]
            Ib=np.where(bigarray[:,j,k,l]==b)[0][0]
            if n[pindex]<bigarray[Ia,j,k,l] and n[pindex]<bigarray[Ib,j,k,l]:
                    final[j,k,l]=np.nan
            if n[pindex]>bigarray[Ia,j,k,l] and n[pindex]>bigarray[Ib,j,k,l]:
                    final[j,k,l]=np.nan
            if Ia>=field[0,:,0,0].shape or Ib>=field[0,:,0,0].shape:
                    final[j,k,l]=np.nan
            else:
                if a > b:
                    slope = (field[j,Ia-1,k,l]-field[j,Ib-1,k,l])/(np.log(bigarray[Ia,j,k,l])-np.log(bigarray[Ib,j,k,l]))
                    final[j,k,l]=(np.log(n[pindex])-np.log(bigarray[Ib,j,k,l]))*slope + field[j,Ib-1,k,l]
                else:
                    slope = (field[j,Ib-1,k,l]-field[j,Ia-1,k,l])/(np.log(bigarray[Ib,j,k,l])-np.log(bigarray[Ia,j,k,l]))
                    final[j,k,l]=(np.log(n[pindex])-np.log(bigarray[Ia,j,k,l]))*slope + field[j,Ia-1,k,l]
        
            l=l+1
        k=k+1;l=0
    j=j+1;k=0


#############################################################################
####################### Writing netcdf4 (3-dim data fields)###################


w_nc_fid = Dataset(dpath3+filename, "w", format="NETCDF4")


########### lon ##################
w_nc_fid.createDimension('lat',lats.size)
w_nc_var = w_nc_fid.createVariable("lat","f4",("lat"))
w_nc_var.setncatts({'long_name': "latitude",\
                    'units': "degrees_north"})
########### lat #####################
w_nc_fid.createDimension('lon',lons.size)
w_nc_var = w_nc_fid.createVariable("lon","f4",("lon"))
w_nc_var.setncatts({'long_name': "longitude",\
                    'units': "degrees_east"})
########### time ##########################
w_nc_fid.createDimension('time',None)
w_nc_var = w_nc_fid.createVariable("time","f8",("time"))
w_nc_var.setncatts({'long_name': "time",\
                    'units': "ka BP", 'calendar': "noleap",\
                    'history':"ka BP"})


w_nc_fid.variables['time'][:] = time
w_nc_fid.variables['lat'][:] = lats
w_nc_fid.variables['lon'][:] = lons



w_nc_var = w_nc_fid.createVariable(var,"f8",("time","lat","lon"))
w_nc_var.setncatts({'long_name': 'U at 300hPa',\
                    'units': "m/s", 'level_desc': '300hPa',\
                    'var_desc': "U"})
w_nc_fid.variables[var][:,:,:] = final[:,:,:]
w_nc_fid.close()

