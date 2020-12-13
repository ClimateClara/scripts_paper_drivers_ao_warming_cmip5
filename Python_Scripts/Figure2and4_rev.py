#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 18:28:01 2016

Compare energy budget in CMIP5 models

@author: Clara Burgard
"""

from pylab import *
from mpl_toolkits.basemap import Basemap, addcyclic
import sys,os,glob
import scipy.io as sio
import scipy.stats as stats
import scipy.ndimage as ndi
import numpy.ma as ma
import matplotlib.dates as mdates
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

#Model(s) used
model_list=['ACCESS1-0','ACCESS1-3','CanESM2','CCSM4','CESM1-BGC','CESM1-CAM5','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','FGOALS-g2',\
'GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','GISS-E2-R','HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR',\
'MIROC5','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M','NorESM1-ME']#,'Ensmean','GFDL-ESM2G','GISS-E2-R-CC','MRI-CGCM3','HadGEM2-ES']

allmodel_list=append(model_list,'Ensemble mean')

#open the dictionary that contains all variables
d={}

#Amount of seconds in the year for the conversion from W to J
yearinsec=3600*365*24.

#path where to store output
outputpath='/work/mh0033/m300411/DataEB/RESULTS/PAPER/Prelim2/'

################## READ IN ARCTIC OCEAN AREA (computed in AOmask.sh) ###################################

d["AOarea_ocean_mod"]=zeros((len(model_list)+1))

for i,mod in enumerate(model_list):
  print mod,'AO area'

  inputpath_time='/work/mh0033/m300411/DataEB/WORK_DATA/'
    
  file0=glob.glob(inputpath_time+'AOarea2/Arctic_ocean_totalarea_%s.nc' %mod)
  fid0=sio.netcdf_file(file0[0])
  d["AOarea_ocean_mod"][i]=fid0.variables['tos'][0]
  
# Ensemble mean AO area
d["AOarea_ocean_mod"][len(model_list)]=nanmean(d["AOarea_ocean_mod"][0:len(model_list)])

#############################################################################################################

#Define axis labels and titles for each variable
tit={}
lab={}
tit['rsds']='ISW'
lab['rsds']='Incoming SW [W/m$^2$]'
tit['rsus']='OSW'
lab['rsus']='Outgoing SW [W/m$^2$]'
tit['rlds']='ILW'
lab['rlds']='Incoming LW [W/m$^2$]'
tit['rlus']='OLW'
lab['rlus']='Outgoing LW [W/m$^2$]'
tit['hfss']='SHF'
lab['hfss']='Sensible Heat Flux [W/m$^2$]'
tit['hfls']='Latent Heat Flux'
lab['hfls']='Latent Heat Flux [W/m$^2$]'
tit['sit']='SIVPA'
lab['sit']='Sea-ice Volume Per Area [m$^3$/m$^2$]'
tit['sic']='Sea-ice Concentration'
lab['sic']='Sea-ice Concentration [%]'
tit['hfds']='OHF'
lab['hfds']='Oceanic Heat Flux [W/m$^2$]'    
tit['netlw']='Net LW'
lab['netlw']='Net LW [W/m$^2$]'
tit['netsw']='Net SW'
lab['netsw']='Net SW [W/m$^2$]'
tit['turb']='Turbulent Heat Fluxes'
lab['turb']='Turbulent Heat Fluxes [W/m$^2$]'
tit['sum_atmos']='Sum Atmo'
lab['sum_atmos']='Sum Atmospheric Fluxes [W/m$^2$]'
tit['totsiv']='Total Sea-ice Volume'
lab['totsiv']='Total Sea-ice Volume [x 10$^{12}$ m$^3$]'
tit['totsia']='Total Sea-ice Area'
lab['totsia']='Total Sea-ice Area [x 10$^{12}$ m$^2$]'
tit['totsnv']='Snow Depth on Sea Ice'
lab['totsnv']='Snow Depth on Sea Ice [m]'
tit['heat_content']='Oceanic Heat Content'
lab['heat_content']='Oceanic Heat Content [x 10$^{25}$ J]'

##############################################################################################################

############################# READ IN ATMOSPHERIC VARIABLES (prepared in work_data_formatting.sh) ##############

#Variables read in from a file
variables=['rlds','rsds','rlus','rsus','hfss','hfls','totsiv','totsia']
#Variables read in from a file + additional variables computed here
variables_all=['rlds','rsds','rlus','rsus','hfss','hfls','netlw','netsw','turb','sum_atmos','totsiv','totsia']

#Create the arrays
for x in variables:
        d["{0}_mod".format(x)]=zeros((2868,len(model_list)))
        d["{0}_mod2".format(x)]=zeros((2868,len(model_list)))


#Read in the variables
for i,mod in enumerate(model_list):
  print mod
 
  
  inputpath_time='/work/mh0033/m300411/DataEB/WORK_DATA/'
  p=0
  for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
    for var in variables:
      if var=='totsiv':
        var1='sit'
      elif var=='totsia':
        var1='sic'
      elif var=='totsnv':
        var1='snd'
      else:
        var1=var
      
      if var=='totsia' or var=='totsiv':
        file0=glob.glob(inputpath_time+mod+'/Method170516/%s_Arctic_fldmean_%s_ensmean_186101-209912.nc' %(var,mod))
      elif var=='totsnv':
        file0=glob.glob(inputpath_time+mod+'/hist_rcp/totsnv_%s_ensmean_186101-209912.nc' %(mod))
      else:
        file0=glob.glob(inputpath_time+mod+'/Method170516/%s_Arctic_monfldmean_remap_%s_ensmean_186101-209912.nc' %(var,mod))
      if not file0:
        print mod+' prob!'
      else:
        fid0=sio.netcdf_file(file0[0])
        time0=fid0.variables['time'][:]
        year0=(time0/365.)+1850
        mon0=(year0-np.floor(year0))*12
        day0=(mon0-np.floor(mon0))*30
        date0=transpose(array([np.floor(year0),np.ceil(mon0),np.floor(day0)]))
        date_new=date0.astype(int)
	if var=='hfss' or var=='hfls':
	  d["{0}_mod".format(var)][:,i]=-fid0.variables[var1][:,0,0]
	else:
	  d["{0}_mod".format(var)][:,i]=fid0.variables[var1][:,0,0]
	d["{0}_mod2".format(var)][:,i]=d["{0}_mod".format(var)][:,i]
	#d["{0}_mod".format(var)][:,i]=d["{0}_mod2".format(var)][:,i]-mean(d["{0}_mod2".format(var)][tref:tref+30,i])

            
#    ##SUMS
#    d["netlw_mod"]=d["rlds_mod"]-d["rlus_mod"]
#    d["netsw_mod"]=d["rsds_mod"]-d["rsus_mod"]
#    d["turb_mod"]=d["hfss_mod"]+d["hfls_mod"]
#    d["sum_atmos_mod"]=d["netlw_mod"]+d["netsw_mod"]+d["turb_mod"]
    d["netlw_mod2"]=d["rlds_mod2"]-d["rlus_mod2"]
    d["netsw_mod2"]=d["rsds_mod2"]-d["rsus_mod2"]
    d["turb_mod2"]=d["hfss_mod2"]+d["hfls_mod2"]
    d["sum_atmos_mod2"]=d["netlw_mod2"]+d["netsw_mod2"]+d["turb_mod2"]

#build the ensemble mean over the model(s)
for var in variables_all:
#  d["{0}_ensmean".format(var)]=nanmean(d["{0}_mod".format(var)],axis=1)
#  d["{0}_ensstd".format(var)]=nanstd(d["{0}_mod".format(var)],axis=1)
  d["{0}_ensmean2".format(var)]=nanmean(d["{0}_mod2".format(var)],axis=1)
  d["{0}_ensstd2".format(var)]=nanstd(d["{0}_mod2".format(var)],axis=1)

  
#append ensmean
for var in variables_all:
	d['{0}_orig'.format(var)]=zeros((2868,len(allmodel_list)))
	for n in range(2868) :
		d['{0}_orig'.format(var)][n,:] = append(d["{0}_mod2".format(var)][n,:], nanmean(d["{0}_mod2".format(var)][n,:]))

#decompose in month and year
for var in variables_all:
	d['{0}_orig_mon'.format(var)] = zeros((13,237,len(allmodel_list)))
	for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
		d['{0}_orig_mon'.format(var)][k,:,:]=d['{0}_orig'.format(var)][1+k:237*12+1:12,:]
	k=12
	d['{0}_orig_mon'.format(var)][k,:]=nanmean(d['{0}_orig_mon'.format(var)][0:12,:],axis=0)

#build the anomaly
for var in variables_all:
  d['{0}_orig_mon_ano'.format(var)] = zeros((13,237,len(allmodel_list)))
  for n in range(237):
    d['{0}_orig_mon_ano'.format(var)][:,n,:] = d['{0}_orig_mon'.format(var)][:,n,:] - nanmean( d['{0}_orig_mon'.format(var)][:,0:100,:] ,axis=1 )

    
## CHECK IF IT HAS BEEN READ IN WELL
#for var in variables_all:
#    plt.figure()
#    plt.plot(d["{0}_mod2".format(var)])
#    plt.plot(d["{0}_ensmean2".format(var)][:],'k-',linewidth=3)
#    plt.title(tit[var]+str(k))

#figure()
#plot(nanmean(d["totsia_orig_mon"][0:12,0:100,:],axis=1),'b-')
#plot(nanmean(d["totsia_orig_mon"][0:12,207:237,:],axis=1),'r-')
#
#figure()
#plot(nanmean(d["totsiv_orig_mon"][0:12,0:100,:],axis=1),'b-')
#plot(nanmean(d["totsiv_orig_mon"][0:12,207:237,:],axis=1),'r-')

########### COMPUTE OCEANIC HEAT CONTENT ###################################

#Create the array
d["txv_mod"]=np.zeros((2868,len(model_list))) #sum of (temperature x cell volume) over depth
d["txv_mod2"]=np.zeros((2868,len(model_list))) #sum of (temperature x cell volume) over depth
heat_content=zeros((2868,len(model_list))) # J 
heat_content2=zeros((2868,len(model_list))) # J

#density of ocean water
rho=1025. #kg m-3

#heat capacity of ocean water in the different models (see http://es-doc.org/#Documentation), if not given, default = 4000.
cp = np.zeros((len(model_list)))
for i,mod in enumerate(model_list):
  print mod
  if mod=='CSIRO-Mk3-6-0':
    cp[i]=4186.
  elif mod=='GFDL-CM3' or mod=='GFDL-ESM2M':
    cp[i]=3992.1
  elif mod=='HadGEM2-CC' or mod=='HadGEM2-ES':
    cp[i]=3988.
  elif mod=='MPI-ESM-LR' or mod=='MPI-ESM-MR':
    cp[i]=3902.
  else:
    cp[i]=4000.

#read in temperature x cell volumes
for i,mod in enumerate(model_list):
  print mod

  inputpath_time='/work/mh0033/m300411/DataEB/WORK_DATA/'
  for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
    file0=glob.glob(inputpath_time+mod+'/Method170516/heatcolumn_Arctic_monfldsum_%s_ensmean_186101-209912.nc' %mod)
    fid0=sio.netcdf_file(file0[0])
    if mod=='MPI-ESM-LR':
      time0=fid0.variables['time'][:]
      year0=(time0/365)+1850
      mon0=(year0-floor(year0))*12
      day0=(mon0-floor(mon0))*30
      date0=transpose(array([floor(year0),ceil(mon0),floor(day0)]))
      date_new=date0.astype(int)
    d["txv_mod"][:,i]=fid0.variables['thetao'][:,0,0]
    #d["txv_mod"][:,i][abs(d["txv_mod"][:,i])>999999]=nan
    d["txv_mod2"][:,i]=d["txv_mod"][:,i]
    #d["txv_mod"][:,i]=d["txv_mod"][:,i]-nanmean(d["txv_mod"][tref:tref+30,i])
   
    
#compute the heat content
for i,mod in enumerate(model_list):
  heat_content[:,i]=rho*cp[i]*d["txv_mod"][:,i]
  heat_content2[:,i]=rho*cp[i]*d["txv_mod2"][:,i]

#Ensemble mean
heat_content_ensmean = np.nanmean(heat_content[:,:],axis=1)
heat_content2_ensmean = np.nanmean(heat_content2[:,:],axis=1)

### CHECK IF IT HAS BEEN READ IN WELL
#plt.figure()
#plt.plot(heat_content2[:,0],label=k)
#plt.legend(loc='best')

#######################################################################################

########### ADD MULTI-MODEL ENSEMBLE MEAN TO THE ARRAYS

#append multi-model ensemble mean values to the other models
for i,var in enumerate(variables_all):
  d['{0}'.format(var)] = zeros((2868,len(allmodel_list)))
  for n in range(2868) :
    d['{0}'.format(var)][n,:] = append( d['{0}_mod2'.format(var)][n,:], nanmean( d['{0}_mod2'.format(var)][n,:] ))

#append multi-model ensemble mean values to the other models
var = 'heat_content'
d['{0}'.format(var)] = zeros((2868,len(allmodel_list)))
for n in range(2868) :
  d['{0}'.format(var)][n,:] = append( heat_content2[n,:], nanmean( heat_content2[n,:] ))

variables_all.append('heat_content')


########### CONVERT EVERYTHING TO ENERGY

##### Constants
rho_i=910.    #ice density kg/m3
rho_sn=250.   #snow density kg/m3
lat_fus=334774.   #latent heat of fusion J/kg
yearinsec=3600*365*24.  #amount of seconds in a year
moninsec=3600*24*30.    #amount of seconds in a month
#####


seaice=-d["totsiv_mod2"]*rho_i*lat_fus #convert sea-ice volume to latent energy

######
#Absolute sources
# Energy exchanged through atmospheric fluxes
E_a = zeros((2868,len(allmodel_list)))

for n in range(2868) :
  E_a[n,:] = append(d["sum_atmos_mod2"][n,:]*moninsec*d['AOarea_ocean_mod'][:-1:], nanmean(d["sum_atmos_mod2"][n,:]*moninsec*d['AOarea_ocean_mod'][-1]))


######
# Absolute sinks
# Heat content : Ocean, Ice volume

H_oc = zeros((2867,len(allmodel_list)))
H_ic = zeros((2867,len(allmodel_list)))

#compute difference between one month and the month before
for n in range(2867) :
  H_oc[n,:] = append(heat_content2[n+1,:] - heat_content2[n,:], nanmean(heat_content2[n+1,:] - heat_content2[n,:]))
  H_ic[n,:] = append(seaice[n+1,:] - seaice[n,:], nanmean(seaice[n+1,:] - seaice[n,:]))


# to compare directly to the fluxes : take the value in the middle of the month
H_oc_comp = H_oc[0:-1:,:] + (H_oc[1::,:] - H_oc[0:-1:,:]) * 0.5 
H_ic_comp = H_ic[0:-1:,:] + (H_ic[1::,:] - H_ic[0:-1:,:]) * 0.5 

######

#### test how it looks like
#z=50*12
#for m,mod in enumerate(allmodel_list):
#  figure()
#  plot(E_a[1:z+1,m] , 'g-')
#  plot(H_oc_comp[0:z,m] + H_ic_comp[0:z,m], 'b-')
#  plot(H_oc_comp[0:z,m] + H_ic_comp[0:z,m] - E_a[1:z+1,m], 'r-')
#  grid()
#  title(mod)
###################################################

##### DECOMPOSE IN MONTHLY AND ANNUAL

var_name = ['H_oc', 'H_ic']
for i,var in enumerate([H_oc, H_ic]):
  d['{0}_mon'.format(var_name[i])] = zeros((13,237,len(allmodel_list)))
  for m,mod in enumerate(allmodel_list):
    print mod
    for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
      if k==11:
        d['{0}_mon'.format(var_name[i])][k,:,:]=var[k:237*12:12,:]+(var[12:238*12:12,:]-var[k:237*12:12,:]) * 0.5
      else:
        d['{0}_mon'.format(var_name[i])][k,:,:]=var[k:237*12:12,:]+(var[k+1:237*12:12,:]-var[k:237*12:12,:]) * 0.5

var_name = ['E_a']
for i,var in enumerate([E_a]):
  d['{0}_mon'.format(var_name[i])] = zeros((13,237,len(allmodel_list)))
  for m,mod in enumerate(model_list):
    print mod
    for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
      d['{0}_mon'.format(var_name[i])][k,:,:]=var[1+k:237*12+1:12,:]

#annual value (in J)
k=12
var_name = ['E_a', 'H_oc', 'H_ic']
for i,var in enumerate([E_a, H_oc, H_ic]):
  d['{0}_mon'.format(var_name[i])][k,:]=nansum(d['{0}_mon'.format(var_name[i])][0:12,:],axis=0)
#for i,var in enumerate(flux_Wm2):
#  d['{0}_mon'.format(var)][k,:]=nansum(d['{0}_mon'.format(var)][0:12,:],axis=0)

d['H_tot_mon']=d['H_oc_mon']+d['H_ic_mon']

##### test if it works
##general formula
#for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','year']):  
#  figure()
#  plot((d['E_a_mon'][k,:] ) / 10**21, 'g-')
#  plot((d['H_oc_mon'][k,:] + d['H_ic_mon'][k,:] ) / 10**21, 'b-')
#  ylim(-4,4)
#  title(month)

#######################################################

#Anomalies
var_name = ['E_a', 'H_oc', 'H_ic', 'H_tot']
for i,var in enumerate(['E_a', 'H_oc', 'H_ic', 'H_tot']):
  d['{0}_mon_ano'.format(var_name[i])] = zeros((13,237,len(allmodel_list)))
  for n in range(237):
    d['{0}_mon_ano'.format(var_name[i])][:,n,:] = d['{0}_mon'.format(var_name[i])][:,n,:] - nanmean( d['{0}_mon'.format(var_name[i])][:,0:100,:] ,axis=1 )

d['E_o_mon_ano'] = d['H_tot_mon_ano']-d['E_a_mon_ano']

########################################

############## DIFFERENT MODEL GROUPS

group_atm = ['ACCESS1-0','ACCESS1-3','CCSM4','CESM1-BGC','CESM1-CAM5','CSIRO-Mk3-6-0','GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-LR','IPSL-CM5B-LR','MIROC5']
group_lat = ['CanESM2','CMCC-CM','CMCC-CMS','CNRM-CM5','GISS-E2-R','HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-MR','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3']
group_both = ['FGOALS-g2','GFDL-CM3','NorESM1-M','NorESM1-ME']



#### Figure 2

matplotlib.rcParams.update({'font.size': 28})
#annual values
k=12

#number of rows
rows=3
#number of models
n=len(allmodel_list)
#x-axis
x=date_new[100*12:-24:12,0]
if n % rows == 0:
     f, axs = plt.subplots(rows, n/rows, sharex=True, sharey=True)
else:
     f, axs = plt.subplots(rows, n/rows+1, sharex=True, sharey=True)

#define format     
f.set_size_inches(33,16)
f.subplots_adjust(bottom=0.15)
f.subplots_adjust(left=0.07)
f.subplots_adjust(right=0.95)
f.subplots_adjust(top=0.92)

for m in range(n):
     y1 = (d['E_a_mon_ano'][k,100::,m] ) / 10**23
     y2 = (d['H_oc_mon_ano'][k,100::,m] + d['H_ic_mon_ano'][k,100::,m] ) / 10**23 - (d['E_a_mon_ano'][k,100::,m] ) / 10**23
     y3 = (d['H_oc_mon_ano'][k,100::,m] + d['H_ic_mon_ano'][k,100::,m] ) / 10**23

     if m % rows == 0:
	  axs[0, m/rows].plot(x,cumsum(y1), color='red', linewidth=5)
	  axs[0, m/rows].plot(x,cumsum(y2), color='turquoise', linewidth=5)
	  axs[0, m/rows].plot(x,cumsum(y3), color='b', linewidth=5)
	  if allmodel_list[m] in group_atm:
	  	axs[0,m/rows].set_title(allmodel_list[m], color='red')
	  elif allmodel_list[m] in group_lat:
	  	axs[0,m/rows].set_title(allmodel_list[m], color='turquoise')
	  elif allmodel_list[m] in group_both:
	  	axs[0,m/rows].set_title(allmodel_list[m], color='grey')
	  elif allmodel_list[m] == 'Ensemble mean':
	  	axs[0,m/rows].set_title(allmodel_list[m], color='k')
	  axs[0,m/rows].plot((x[0],x[-1]),(0,0),'k-',linewidth=1.5)
	  axs[0,m/rows].set_xlim(x[0],x[-1])
	  #axs[0,m/rows].grid()
     elif m % rows == 1:
	  axs[1, m/rows].plot(x,cumsum(y1), color='red', linewidth=5)
	  axs[1, m/rows].plot(x,cumsum(y2), color='turquoise', linewidth=5)
	  axs[1, m/rows].plot(x,cumsum(y3), color='b', linewidth=5)
	  if allmodel_list[m] in group_atm:
	  	axs[1,m/rows].set_title(allmodel_list[m], color='red')
	  elif allmodel_list[m] in group_lat:
	  	axs[1,m/rows].set_title(allmodel_list[m], color='turquoise')
	  elif allmodel_list[m] in group_both:
	  	axs[1,m/rows].set_title(allmodel_list[m], color='grey')
	  elif allmodel_list[m] == 'Ensemble mean':
	  	axs[1,m/rows].set_title(allmodel_list[m], color='k')
	  axs[1,m/rows].plot((x[0],x[-1]),(0,0),'k-',linewidth=1.5)
	  axs[1,m/rows].set_xlim(x[0],x[-1])
	  #axs[1,m/rows].grid()
     elif m % rows == 2: 
	  axs[2, m/rows].plot(x,cumsum(y1), color='red', linewidth=5)
	  axs[2, m/rows].plot(x,cumsum(y2), color='turquoise', linewidth=5)
	  axs[2, m/rows].plot(x,cumsum(y3), color='b', linewidth=5)
	  if allmodel_list[m] in group_atm:
	  	axs[2,m/rows].set_title(allmodel_list[m], color='r')
	  elif allmodel_list[m] in group_lat:
	  	axs[2,m/rows].set_title(allmodel_list[m], color='turquoise')
	  elif allmodel_list[m] in group_both:
	  	axs[2,m/rows].set_title(allmodel_list[m], color='grey')
	  elif allmodel_list[m] == 'Ensemble mean':
	  	axs[2,m/rows].set_title(allmodel_list[m], color='k')
	  axs[2,m/rows].plot((x[0],x[-1]),(0,0),'k-',linewidth=1.5)
	  axs[2,m/rows].set_xlim(x[0],x[-1])
	  #axs[2,m/rows].grid()
    

     
plt.setp([a.get_xticklabels() for a in axs[0, :]], visible=False)
plt.setp([a.get_yticklabels() for a in axs[:, 1]], visible=False)
plt.setp([a.get_yticklabels() for a in axs[:, 2]], visible=False)
plt.setp([a.get_yticklabels() for a in axs[:, 3]], visible=False)
plt.setp([a.get_xticklabels() for a in axs[:, 0]], rotation=90)
plt.setp([a.get_xticklabels() for a in axs[:, 1]], rotation=90)
plt.setp([a.get_xticklabels() for a in axs[:, 2]], rotation=90)
plt.setp([a.get_xticklabels() for a in axs[:, 3]], rotation=90)
plt.setp([a.get_xticklabels() for a in axs[:, 4]], rotation=90)
plt.setp([a.get_xticklabels() for a in axs[:, 5]], rotation=90)
plt.setp([a.get_xticklabels() for a in axs[:, 6]], rotation=90)
plt.setp([a.get_xticklabels() for a in axs[:, 7]], rotation=90)
plt.setp([a.get_xticklabels() for a in axs[:, 8]], rotation=90)
#plt.setp(axs[1,8].get_xticklabels(), rotation=90, visible=True)

f.text(0.5,0.01,'Time [Years]' ,ha='center')
f.text(0.02,0.7, 'Cumulated Energy Anomalies [x 10$^{23}$ J]' ,ha='center',rotation='vertical')
#f.savefig(outputpath+'all_fluxes_cumsum.pdf',bbox_inches='tight',orientation='landscape')
f.show()



#####################################################################################

#### Figure 4

blax=arange(-3,3)
blay=arange(3,-3,-1)

fig, ax = plt.subplots()

fig1 = matplotlib.pyplot.gcf()
fig1.set_size_inches(17,15)
fig1.subplots_adjust(bottom=0.08)
fig1.subplots_adjust(left=0.05)
fig1.subplots_adjust(right=0.95)
fig1.subplots_adjust(top=0.92)

for m,mod in enumerate(allmodel_list):
	if mod in group_atm:
		scatter(sum(d['E_o_mon_ano'.format(var)][k,100::,m]) / 10**23, sum(d['E_a_mon_ano'][12,100::,m]) / 10**23, c='r', marker='o', edgecolors='none',s=150)
	elif mod in group_lat:
		scatter(sum(d['E_o_mon_ano'.format(var)][k,100::,m]) / 10**23, sum(d['E_a_mon_ano'][12,100::,m]) / 10**23, c='turquoise', marker='o', edgecolors='none',s=150)
	elif mod in group_both:
		scatter(sum(d['E_o_mon_ano'.format(var)][k,100::,m]) / 10**23, sum(d['E_a_mon_ano'][12,100::,m]) / 10**23, c='grey', marker='o', edgecolors='none',s=150)
	elif mod=='Ensemble mean':
		scatter(sum(d['E_o_mon_ano'.format(var)][k,100::,m]) / 10**23, sum(d['E_a_mon_ano'][12,100::,m]) / 10**23, c='k', marker='o', edgecolors='none',s=150)
#plot(blax,blay,'k--')
#for m,mod in enumerate(allmodel_list):
#     scatter(sum(d['E_o_mon_ano'.format(var)][k,100::,m]) / 10**23, sum(d['E_a_mon_ano'][12,100::,m]) / 10**23, c=colm[m], marker=markerstyle[m], label=mod, edgecolors='none',s=70)
#      plot(blax,blay,'k--')
plot(blax,blay,'k--')
xlim(-2,2.5)
ylim(-2,2.5)
ylabel('Cumulated energy due to $\Delta H_{\t{sfc}}$ end of 21st century [x 10$^{23}$ J]')
xlabel('Cumulated energy due to $\Delta H_{\t{lat}}$ end of 21st century [x 10$^{23}$ J]')
grid()
fit = np.polyfit(sum(d['E_o_mon_ano'.format(var)][12,100::,:],axis=0), sum(d['E_a_mon_ano'][12,100::,:],axis=0), deg=1)
plot(arange(-1.2,2.5,0.2), fit[0] * arange(-1.2,2.5,0.2) + fit[1]/10**23, color='black')
#plot(sum(d['E_o_mon_ano'.format(var)][k,100::,:]/10**23,axis=0), fit[0] * sum(d['E_o_mon_ano'.format(var)][k,100::,:]/10**23,axis=0) + fit[1]/10**23, color='black')
slope, intercept, r_value, p_value, std_err = stats.linregress(sum(d['E_o_mon_ano'.format(var)][k,100::,:],axis=0), sum(d['E_a_mon_ano'][12,100::,:],axis=0))
legend(loc=1, scatterpoints = 1,ncol=2,fontsize=18)
text(-1.4,-1.5,'R = '+str(round(corrcoef(sum(d['E_o_mon_ano'.format(var)][k,100::,:],axis=0), sum(d['E_a_mon_ano'][12,100::,:],axis=0))[0,1],2)))
#text(min(sum(d['E_o_mon_ano'.format(var)][k,100::,:],axis=0))/10**23 + 0.2,-1.0,'R$^2$ ='+str(round(r_value**2,2)))

#savefig(outputpath+'EaVSEo_rev.pdf', bbox_inches='tight', orientation='landscape')




