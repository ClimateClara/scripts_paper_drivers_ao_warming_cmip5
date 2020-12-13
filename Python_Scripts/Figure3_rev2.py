#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 18:14:18 2016

Comparison of data from the 26 CMIP5 models
Figure 3

#modified : inserted some more variables, like meridional atmospheric transport, humidity

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
import pandas as pd

#Model(s) used
model_list=['ACCESS1-0','ACCESS1-3','CanESM2','CCSM4','CESM1-BGC','CESM1-CAM5','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','FGOALS-g2',\
'GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','GISS-E2-R','HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR',\
'MIROC5','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3','NorESM1-M','NorESM1-ME']#,'Ensmean','GFDL-ESM2G','GISS-E2-R-CC','MRI-CGCM3','HadGEM2-ES']



allmodel_list=append(model_list,'Ensemble mean')


#################

  
group_atm = ['ACCESS1-0','ACCESS1-3','CCSM4','CESM1-BGC','CESM1-CAM5','CSIRO-Mk3-6-0','GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-LR','IPSL-CM5B-LR','MIROC5']
group_lat = ['CanESM2','CMCC-CM','CMCC-CMS','CNRM-CM5','GISS-E2-R','HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-MR','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3']
group_both = ['FGOALS-g2','GFDL-CM3','NorESM1-M','NorESM1-ME']



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

for m,mod in enumerate(model_list):
	if mod in group_atm:
		scatter(m,d["AOarea_ocean_mod"][m],c='r')
	elif mod in group_lat:
		scatter(m,d["AOarea_ocean_mod"][m],c='turquoise')
	elif mod in group_both:
		scatter(m,d["AOarea_ocean_mod"][m],c='grey')

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
variables=['rlds','rsds','rlus','rsus','hfss','hfls','totsiv','totsia','tas','psxprw','prw','rlut','rsut','rsdt']
#Variables read in from a file + additional variables computed here
variables_all=['rlds','rsds','rlus','rsus','hfss','hfls','netlw','netsw','turb','sum_atmos','totsiv','totsia','tas','psxprw','prw','rlut','rsut','rsdt','rad_toa','atmmer']

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
      elif var=='psxprw':
        var1='ps'
      else:
        var1=var
      
      if var=='totsia' or var=='totsiv':
        file0=glob.glob(inputpath_time+mod+'/Method170516/%s_Arctic_fldmean_%s_ensmean_186101-209912.nc' %(var,mod))
      elif var=='totsnv':
        file0=glob.glob(inputpath_time+mod+'/hist_rcp/totsnv_%s_ensmean_186101-209912.nc' %(mod))
      else:
        file0=glob.glob(inputpath_time+mod+'/Method170516/%s_Arctic_monfldmean_remap_%s_ensmean_186101-209912.nc' %(var,mod))
      if not file0 or (mod=='FGOALS-g2' and var=='rlut'):
        print mod+' prob!'
        d["{0}_mod".format(var)][:,i] = nan
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
    d["rad_toa_mod2"]=d["rsdt_mod2"]-d["rsut_mod2"]-d["rlut_mod2"]
    d["atmmer_mod2"]=-d["rad_toa_mod2"]+d["sum_atmos_mod2"]
    

#build the ensemble mean over the model(s)
for var in variables_all:
#  d["{0}_ensmean".format(var)]=nanmean(d["{0}_mod".format(var)],axis=1)
#  d["{0}_ensstd".format(var)]=nanstd(d["{0}_mod".format(var)],axis=1)
  d["{0}_ensmean2".format(var)]=nanmean(d["{0}_mod2".format(var)],axis=1)
  d["{0}_ensstd2".format(var)]=nanstd(d["{0}_mod2".format(var)],axis=1)
  
## CHECK IF IT HAS BEEN READ IN WELL
#for var in variables_all:
#    plt.figure()
#    plt.plot(d["{0}_mod2".format(var)])
#    plt.plot(d["{0}_ensmean2".format(var)][:],'k-',linewidth=3)
#    plt.title(tit[var]+str(k))

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

#figure()
#for m,mod in enumerate(allmodel_list):
#  #figure()
#  if mod in group_atm:
#		plot(d['psxprw_orig_mon_ano'][12,100::,m],'r-')
#  elif mod in group_lat:
#		plot(d['psxprw_orig_mon_ano'][12,100::,m],'c-')
#  elif mod in group_both:
#		plot(d['psxprw_orig_mon_ano'][12,100::,m],'k-')
#  title(mod)
#
#figure()
#for m,mod in enumerate(allmodel_list):
#  #figure()
#  if mod in group_atm:
#		plot(d['prw_orig_mon_ano'][12,100::,m],'r-')
#  elif mod in group_lat:
#		plot(d['prw_orig_mon_ano'][12,100::,m],'c-')
#  elif mod in group_both:
#		plot(d['prw_orig_mon_ano'][12,100::,m],'k-')
#  title(mod)
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

#append ensmean
for var in ['heat_content']:
	d['{0}_orig'.format(var)]=zeros((2868,len(allmodel_list)))
	for n in range(2868) :
		d['{0}_orig'.format(var)][n,:] = append(heat_content2[n,:], nanmean(heat_content2[n,:]))

#decompose in month and year
for var in ['heat_content']:
	d['{0}_orig_mon'.format(var)] = zeros((13,237,len(allmodel_list)))
	for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
		d['{0}_orig_mon'.format(var)][k,:,:]=d['{0}_orig'.format(var)][1+k:237*12+1:12,:]
	k=12
	d['{0}_orig_mon'.format(var)][k,:]=nanmean(d['{0}_orig_mon'.format(var)][0:12,:],axis=0)

#build the anomaly
for var in ['heat_content']:
  d['{0}_orig_mon_ano'.format(var)] = zeros((13,237,len(allmodel_list)))
  for n in range(237):
    d['{0}_orig_mon_ano'.format(var)][:,n,:] = d['{0}_orig_mon'.format(var)][:,n,:] - nanmean( d['{0}_orig_mon'.format(var)][:,0:100,:] ,axis=1 )


#######################################################################################

########### ICE PRODUCTION (computed in understand_Hlat.sh) ###################################

#declare variables
d["icep_mod"]=zeros((238,len(model_list)))
d["icep_mod2"]=zeros((238,len(model_list)))


for i,mod in enumerate(model_list):
  print mod

  inputpath_time='/work/mh0033/m300411/DataEB/WORK_DATA/'
  for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
    file0=glob.glob(inputpath_time+mod+'/Method170516/iceprod0903_Arctic_fldmean_%s_ensmean_186101-209912.nc' %mod)
    fid0=sio.netcdf_file(file0[0])
    if mod=='MPI-ESM-LR':
      time0=fid0.variables['time'][:]
      year0=(time0/365)+1850
      mon0=(year0-floor(year0))*12
      day0=(mon0-floor(mon0))*30
      date0=transpose(array([floor(year0),ceil(mon0),floor(day0)]))
      date_new=date0.astype(int)
    d["icep_mod"][:,i]=fid0.variables['sit'][:,0,0]
    d["icep_mod2"][:,i]=d["icep_mod"][:,i]
#    d["icep_mod"][:,i]=d["icep_mod"][:,i]-nanmean(d["icep_mod"][tref:tref+30,i])
   
d['icep_ensmean']=nanmean(d['icep_mod'][:,:],axis=1)
d['icep_ensmean2']=nanmean(d['icep_mod2'][:,:],axis=1)

##CHECK
#figure()
#plot(d['icep_mod'][:,:])
#legend(loc='best')

#append ensmean
for var in ['icep']:
	d['{0}_orig'.format(var)]=zeros((238,len(allmodel_list)))
	for n in range(238) :
		d['{0}_orig'.format(var)][n,:] = append(d["{0}_mod2".format(var)][n,:], nanmean(d["{0}_mod2".format(var)][n,:]))

#build the anomaly
for var in ['icep']:
  d['{0}_orig_mon_ano'.format(var)] = zeros((237,len(allmodel_list)))
  for n,j in enumerate(range(1,238)):
    d['{0}_orig_mon_ano'.format(var)][n,:] = d['{0}_orig'.format(var)][j,:] - nanmean( d['{0}_orig'.format(var)][0:100,:] ,axis=0 )


#######################################################################################

################## SEA SURFACE TEMPERATURE #################

#declare variables
d["sst_mod"]=zeros((2868,len(model_list)))
d["sst_mod2"]=zeros((2868,len(model_list)))


for i,mod in enumerate(model_list):
  print mod

  inputpath_time='/work/mh0033/m300411/DataEB/WORK_DATA/'
  for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
	if mod == 'FGOALS-g2' or mod=='HadGEM2-ES':
		d["sst_mod"][:,i]=nan
		d["sst_mod2"][:,i]=nan
	else:
		file0=glob.glob(inputpath_time+mod+'/Method170516/tos_Arctic_monfldmean_%s_ensmean_186101-209912.nc' %mod)
		fid0=sio.netcdf_file(file0[0])
		if mod=='MPI-ESM-LR':
			time0=fid0.variables['time'][:]
			year0=(time0/365)+1850
			mon0=(year0-floor(year0))*12
			day0=(mon0-floor(mon0))*30
			date0=transpose(array([floor(year0),ceil(mon0),floor(day0)]))
			date_new=date0.astype(int)
		d["sst_mod"][:,i]=fid0.variables['tos'][:,0,0]
		d["sst_mod2"][:,i]=d["sst_mod"][:,i]
		#d["sst_mod"][:,i]=d["sst_mod"][:,i]-nanmean(d["sst_mod"][tref:tref+30,i])
   
d['sst_ensmean']=nanmean(d['sst_mod'][:,:],axis=1)
d['sst_ensmean2']=nanmean(d['sst_mod2'][:,:],axis=1)

d['vert_grad_mod2'] = d["sst_mod2"] - d["tas_mod2"]
##CHECK
#figure()
#plot(d['sst_mod'][:,:])
#legend(loc='best')


#append ensmean
for var in ['sst','vert_grad']:
	d['{0}_orig'.format(var)]=zeros((2868,len(allmodel_list)))
	for n in range(2868) :
		d['{0}_orig'.format(var)][n,:] = append(d["{0}_mod2".format(var)][n,:], nanmean(d["{0}_mod2".format(var)][n,:]))

#decompose in month and year
for var in ['sst','vert_grad']:
	d['{0}_orig_mon'.format(var)] = zeros((13,237,len(allmodel_list)))
	for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
		d['{0}_orig_mon'.format(var)][k,:,:]=d['{0}_orig'.format(var)][1+k:237*12+1:12,:]
	k=12
	d['{0}_orig_mon'.format(var)][k,:]=nanmean(d['{0}_orig_mon'.format(var)][0:12,:],axis=0)

#build the anomaly
for var in ['sst','vert_grad']:
  d['{0}_orig_mon_ano'.format(var)] = zeros((13,237,len(allmodel_list)))
  for n in range(237):
    d['{0}_orig_mon_ano'.format(var)][:,n,:] = d['{0}_orig_mon'.format(var)][:,n,:] - nanmean( d['{0}_orig_mon'.format(var)][:,0:100,:] ,axis=1 )

#figure()
#for m,mod in enumerate(allmodel_list):
#	if mod in group_atm:
#		plot(d['tas_orig_mon_ano'][12,100::,m],'r-')
#	elif mod in group_lat:
#		plot(d['tas_orig_mon_ano'][12,100::,m],'c-')
#	elif mod in group_both:
#		plot(d['tas_orig_mon_ano'][12,100::,m],'k-')
#
#figure()
#for m,mod in enumerate(allmodel_list):
#	if mod in group_atm:
#		plot(d['vert_grad_orig_mon_ano'][12,100::,m],'r-')
#	elif mod in group_lat:
#		plot(d['vert_grad_orig_mon_ano'][12,100::,m],'c-')
#	elif mod in group_both:
#		plot(d['vert_grad_orig_mon_ano'][12,100::,m],'k-')

##################################################################################

################## SURFACE OCEAN TEMPERATURE GLOBAL #################

#declare variables
d["sstglob_mod"]=zeros((2868,len(model_list)))
d["sstglob_mod2"]=zeros((2868,len(model_list)))


for i,mod in enumerate(model_list):
  	print mod

  	inputpath_time='/work/mh0033/m300411/DataEB/WORK_DATA/'
	if mod == 'FGOALS-g2' or mod=='HadGEM2-ES':
		d["sstglob_mod"][:,i]=nan
		d["sstglob_mod2"][:,i]=nan
	else:
    		file0=glob.glob(inputpath_time+mod+'/Method170516/AA_ocean_%s_ensmean_186101-209912.nc' %mod)
    		fid0=sio.netcdf_file(file0[0])
    		if mod=='MPI-ESM-LR':
			time0=fid0.variables['time'][:]
			year0=(time0/365)+1850
			mon0=(year0-floor(year0))*12
			day0=(mon0-floor(mon0))*30
			date0=transpose(array([floor(year0),ceil(mon0),floor(day0)]))
			date_new=date0.astype(int)
    		d["sstglob_mod"][:,i]=fid0.variables['tos'][:,0,0]
    		d["sstglob_mod2"][:,i]=d["sstglob_mod"][:,i]
#    		d["sstglob_mod"][:,i]=d["sstglob_mod"][:,i]-nanmean(d["sstglob_mod"][tref:tref+30,i])
   
d['sstglob_ensmean']=nanmean(d['sstglob_mod'][:,:],axis=1)
d['sstglob_ensmean2']=nanmean(d['sstglob_mod2'][:,:],axis=1)

#append ensmean
for var in ['sstglob']:
	d['{0}_orig'.format(var)]=zeros((2868,len(allmodel_list)))
	for n in range(2868) :
		d['{0}_orig'.format(var)][n,:] = append(d["{0}_mod2".format(var)][n,:], nanmean(d["{0}_mod2".format(var)][n,:]))

#decompose in month and year
for var in ['sstglob']:
	d['{0}_orig_mon'.format(var)] = zeros((13,237,len(allmodel_list)))
	for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
		d['{0}_orig_mon'.format(var)][k,:,:]=d['{0}_orig'.format(var)][1+k:237*12+1:12,:]
	k=12
	d['{0}_orig_mon'.format(var)][k,:]=nanmean(d['{0}_orig_mon'.format(var)][0:12,:],axis=0)

#build the anomaly
for var in ['sstglob']:
  d['{0}_orig_mon_ano'.format(var)] = zeros((13,237,len(allmodel_list)))
  for n in range(237):
    d['{0}_orig_mon_ano'.format(var)][:,n,:] = d['{0}_orig_mon'.format(var)][:,n,:] - nanmean( d['{0}_orig_mon'.format(var)][:,0:100,:] ,axis=1 )

##CHECK
#figure()
#for m,mod in enumerate(model_list):
#	figure()
#	plot(d["sstglob_orig_mon"][12,:,m])
#	title(mod)


################## SURFACE TEMPERATURE #################

#global
d["tasglob_mod"]=zeros((2868,len(model_list)))
d["tasglob_mod2"]=zeros((2868,len(model_list)))
#tropical
d["tas_trop_mod"]=zeros((2868,len(model_list)))
d["tas_trop_mod2"]=zeros((2868,len(model_list)))
#lower latitudes
d["tas_low_mod"]=zeros((2868,len(model_list)))
d["tas_low_mod2"]=zeros((2868,len(model_list)))
#Arctic (land+ocean)
d["tas_arc_mod"]=zeros((2868,len(model_list)))
d["tas_arc_mod2"]=zeros((2868,len(model_list)))

for var in ['tasglob','tas_low','tas_trop','tas_arc']:
	for i,mod in enumerate(model_list):
  		print mod

  		inputpath_time='/work/mh0033/m300411/DataEB/WORK_DATA/'
  		if var=='tasglob':
			file0=glob.glob(inputpath_time+mod+'/Method170516/tas_glob_ocean_%s_ensmean_186101-209912.nc' %mod)
		elif var=='tas_trop':
			file0=glob.glob(inputpath_time+mod+'/Method170516/AAtrop_%s_ensmean_186101-209912.nc' %mod)
		elif var=='tas_arc':
			file0=glob.glob(inputpath_time+mod+'/Method170516/AAarc_%s_ensmean_186101-209912.nc' %mod)
   		elif var=='tas_low':
			file0=glob.glob(inputpath_time+mod+'/Method170516/tas_lowlat_ocean_%s_ensmean_186101-209912.nc' %mod)
		fid0=sio.netcdf_file(file0[0])
		if mod=='MPI-ESM-LR':
			time0=fid0.variables['time'][:]
			year0=(time0/365)+1850
			mon0=(year0-floor(year0))*12
			day0=(mon0-floor(mon0))*30
			date0=transpose(array([floor(year0),ceil(mon0),floor(day0)]))
			date_new=date0.astype(int)
		d["{0}_mod".format(var)][:,i]=fid0.variables['tas'][:,0,0]
		d["{0}_mod2".format(var)][:,i]=d["{0}_mod".format(var)][:,i]
#		d["{0}_mod".format(var)][:,i]=d["{0}_mod".format(var)][:,i]-nanmean(d["{0}_mod".format(var)][tref:tref+30,i])
   
	d['{0}_ensmean'.format(var)]=nanmean(d['{0}_mod'.format(var)][:,:],axis=1)
	d['{0}_ensmean2'.format(var)]=nanmean(d['{0}_mod2'.format(var)][:,:],axis=1)

#air temperature gradient between lower latitudes and Arctic
d['gradient_low_arc_atm_mod2'] = d['tas_low_mod2'] - d['tas_arc_mod2']

#append ensmean
for var in ['tasglob','tas_low','tas_trop','tas_arc','gradient_low_arc_atm']:
	d['{0}_orig'.format(var)]=zeros((2868,len(allmodel_list)))
	for n in range(2868) :
		d['{0}_orig'.format(var)][n,:] = append(d["{0}_mod2".format(var)][n,:], nanmean(d["{0}_mod2".format(var)][n,:]))

#decompose in month and year
for var in ['tasglob','tas_low','tas_trop','tas_arc','gradient_low_arc_atm']:
	d['{0}_orig_mon'.format(var)] = zeros((13,237,len(allmodel_list)))
	for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
		d['{0}_orig_mon'.format(var)][k,:,:]=d['{0}_orig'.format(var)][1+k:237*12+1:12,:]
	k=12
	d['{0}_orig_mon'.format(var)][k,:]=nanmean(d['{0}_orig_mon'.format(var)][0:12,:],axis=0)

#build the anomaly
for var in ['tasglob','tas_low','tas_trop','tas_arc','gradient_low_arc_atm']:
  d['{0}_orig_mon_ano'.format(var)] = zeros((13,237,len(allmodel_list)))
  for n in range(237):
    d['{0}_orig_mon_ano'.format(var)][:,n,:] = d['{0}_orig_mon'.format(var)][:,n,:] - nanmean( d['{0}_orig_mon'.format(var)][:,0:100,:] ,axis=1 )

#Arctic Amplification
AAmp0 = (nanmean(d['tas_orig_mon_ano'][12,209:239,:],axis=0)) / (nanmean(d['tasglob_orig_mon_ano'][12,209:239,:],axis=0))
AAmp = (nanmean(d['tas_arc_orig_mon_ano'][12,209:239,:],axis=0)) / (nanmean(d['tas_trop_orig_mon_ano'][12,209:239,:],axis=0))




##################################################################################


    
################## MEAN TEMPERATURE OVER ALL DEPTHS #################

#ocean temperature Arctic
d["thetao_mod"]=zeros((2868,len(model_list)))
d["thetao_mod2"]=zeros((2868,len(model_list)))
#ocean temperature lower latitudes
d["thetao_low_mod"]=zeros((2868,len(model_list)))
d["thetao_low_mod2"]=zeros((2868,len(model_list)))

for var in ['thetao','thetao_low']:
	for i,mod in enumerate(model_list):
	  print mod

	  inputpath_time='/work/mh0033/m300411/DataEB/WORK_DATA/'
	  if var=='thetao':
	    	file0=glob.glob(inputpath_time+mod+'/Method170516/thetao_fldmean_Arctic_%s_ensmean_186101-209912.nc' %mod)
	  else : 
		file0=glob.glob(inputpath_time+mod+'/Method170516/thetao_fldmean_lowlat_%s_ensmean_186101-209912.nc' %mod)
	  if not file0 : 
		  d["{0}_mod".format(var)][:,i] = nan
		  d["{0}_mod2".format(var)][:,i] = nan
	  else:
		  fid0=sio.netcdf_file(file0[0])
		  if mod=='MPI-ESM-LR':
		      time0=fid0.variables['time'][:]
		      year0=(time0/365)+1850
		      mon0=(year0-floor(year0))*12
		      day0=(mon0-floor(mon0))*30
		      date0=transpose(array([floor(year0),ceil(mon0),floor(day0)]))
		      date_new=date0.astype(int)
		  d["{0}_mod".format(var)][:,i]=fid0.variables['thetao'][:,0,0]
		  d["{0}_mod".format(var)][:,i][abs(d["{0}_mod".format(var)][:,i])>999999]=nan
		  d["{0}_mod2".format(var)][:,i]=d["{0}_mod".format(var)][:,i]
#		  d["{0}_mod".format(var)][:,i]=d["{0}_mod".format(var)][:,i]-nanmean(d["{0}_mod".format(var)][tref:tref+30,i])
	   
	d['{0}_ensmean'.format(var)]=nanmean(d['{0}_mod'.format(var)][:,:],axis=1)
	d['{0}_ensmean2'.format(var)]=nanmean(d['{0}_mod2'.format(var)][:,:],axis=1)

 #ocean temperature gradient between lower latitudes and Arctic
d['gradient_low_arc_oc_mod2'] = d['thetao_low_mod2'] - d['thetao_mod2']

#append ensmean
for var in ['thetao','thetao_low','gradient_low_arc_oc']:
	d['{0}_orig'.format(var)]=zeros((2868,len(allmodel_list)))
	for n in range(2868) :
		d['{0}_orig'.format(var)][n,:] = append(d["{0}_mod2".format(var)][n,:], nanmean(d["{0}_mod2".format(var)][n,:]))

#decompose in month and year
for var in ['thetao','thetao_low','gradient_low_arc_oc']:
	d['{0}_orig_mon'.format(var)] = zeros((13,237,len(allmodel_list)))
	for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
		d['{0}_orig_mon'.format(var)][k,:,:]=d['{0}_orig'.format(var)][1+k:237*12+1:12,:]
	k=12
	d['{0}_orig_mon'.format(var)][k,:]=nanmean(d['{0}_orig_mon'.format(var)][0:12,:],axis=0)

#build the anomaly
for var in ['thetao','thetao_low','gradient_low_arc_oc']:
  d['{0}_orig_mon_ano'.format(var)] = zeros((13,237,len(allmodel_list)))
  for n in range(237):
    d['{0}_orig_mon_ano'.format(var)][:,n,:] = d['{0}_orig_mon'.format(var)][:,n,:] - nanmean( d['{0}_orig_mon'.format(var)][:,0:100,:] ,axis=1 )


#################################################################################

##################  MASS TRANSPORT

#look at all the straits
for line in ['barents','bering','canada','fram','denmark','icefar','farsco']:
	d['mfo_{0}'.format(line)] = zeros((2868,len(model_list)))

for m,mod in enumerate(model_list):
	file0=glob.glob(inputpath_time+mod+'/hist_rcp2/mfo_%s_ensmean_186101-209912.nc' %mod)
	if not file0 or mod=='MPI-ESM-LR' or mod=='MRI-CGCM3' or mod=='GFDL-ESM2M':# or mod=='MRI-CGCM3' or mod=='GFDL-ESM2M':
		d["mfo_barents"][:,m] = nan
		d["mfo_bering"][:,m] = nan
		d["mfo_canada"][:,m] = nan
		d["mfo_fram"][:,m] = nan
		d["mfo_denmark"][:,m] = nan
		d["mfo_icefar"][:,m] = nan
		d["mfo_farsco"][:,m] = nan
	else:
		fid0=sio.netcdf_file(file0[0])
		time0=fid0.variables['time'][:]
		year0=(time0/365)+1850
		mon0=(year0-floor(year0))*12
		day0=(mon0-floor(mon0))*30
		date0=transpose(array([floor(year0),ceil(mon0),floor(day0)]))
		date_new=date0.astype(int)
		if mod=='MPI-ESM-LR':
			d["mfo_barents"][:,m]=-fid0.variables['mfo'][:,0]
		else:
			d["mfo_barents"][:,m]=fid0.variables['mfo'][:,0]
		d["mfo_bering"][:,m]=fid0.variables['mfo'][:,1]
		if mod=='MPI-ESM-LR':
			d["mfo_canada"][:,m]=-fid0.variables['mfo'][:,2]
		else:
			d["mfo_canada"][:,m]=fid0.variables['mfo'][:,2]
		if mod=='MRI-CGCM3' or mod=='GFDL-ESM2M':
			d["mfo_fram"][:,m]=-fid0.variables['mfo'][:,9]
		else:
			d["mfo_fram"][:,m]=fid0.variables['mfo'][:,9]
		d["mfo_denmark"][:,m]=fid0.variables['mfo'][:,3]
		d["mfo_icefar"][:,m]=fid0.variables['mfo'][:,10]
		d["mfo_farsco"][:,m]=fid0.variables['mfo'][:,7]


d['mfo_total1'] = d["mfo_barents"] + d["mfo_bering"] + d["mfo_canada"] + d["mfo_fram"]
d['mfo_total2'] = d["mfo_denmark"] + d["mfo_icefar"] + d["mfo_farsco"] + d["mfo_canada"] + d["mfo_bering"]

#append ensmean
for var in ['barents','bering','canada','fram','denmark','icefar','farsco','total1','total2']:
	d['mfo_{0}_orig'.format(var)]=zeros((2868,len(allmodel_list)))
	for n in range(2868) :
		d['mfo_{0}_orig'.format(var)][n,:] = append(d["mfo_{0}".format(var)][n,:], nanmean(d["mfo_{0}".format(var)][n,:]))

#decompose in month and year
for var in ['barents','bering','canada','fram','denmark','icefar','farsco','total1','total2']:
	d['mfo_{0}_orig_mon'.format(var)] = zeros((13,237,len(allmodel_list)))
	for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
		d['mfo_{0}_orig_mon'.format(var)][k,:,:]=d['mfo_{0}_orig'.format(var)][1+k:237*12+1:12,:]
	k=12
	d['mfo_{0}_orig_mon'.format(var)][k,:]=nanmean(d['mfo_{0}_orig_mon'.format(var)][0:12,:],axis=0)

#build the anomaly
for var in ['barents','bering','canada','fram','denmark','icefar','farsco','total1','total2']:
  d['mfo_{0}_orig_mon_ano'.format(var)] = zeros((13,237,len(allmodel_list)))
  for n in range(237):
    d['mfo_{0}_orig_mon_ano'.format(var)][:,n,:] = d['mfo_{0}_orig_mon'.format(var)][:,n,:] - nanmean( d['mfo_{0}_orig_mon'.format(var)][:,0:100,:] ,axis=1 )


# Denmark Strait looks weird, we do not use it anymore    
    
#############################################################################################

################## MEAN TEMPERATURE ACROSS THE STRAITS #################

for line in ['barents','bering','canada','fram','denmark','icefar','farsco','fram_l','fram_r']:
	d['thetafo_{0}'.format(line)] = zeros((2868,len(model_list)))

short={}
short['barents'] = 'BO'
short['bering'] = 'BS'
short['canada'] = 'CA'
short['fram'] = 'FS'
short['denmark'] = 'DS'
short['icefar'] = 'IF'
short['farsco'] = 'FaSc'
short['fram_l'] = 'FS_l'
short['fram_r'] = 'FS_r'

for m,mod in enumerate(model_list):
	print mod
	for line in ['barents','bering','canada','fram','denmark','icefar','farsco','fram_l','fram_r']:
		file0=glob.glob(inputpath_time+mod+'/Method170516/thetao_fldmean_%s_%s_ensmean_186101-209912.nc' %(short['{0}'.format(line)],mod))
		if not file0:#or mod=='GFDL-ESM2G' : #or mod=='MPI-ESM-LR':# or mod=='MRI-CGCM3' or mod=='GFDL-ESM2M':
			d["thetafo_{0}".format(line)][:,m] = nan
			print 'prob '+mod
		else:
			fid0=sio.netcdf_file(file0[0])
			time0=fid0.variables['time'][:]
			year0=(time0/365)+1850
			mon0=(year0-floor(year0))*12
			day0=(mon0-floor(mon0))*30
			date0=transpose(array([floor(year0),ceil(mon0),floor(day0)]))
			date_new=date0.astype(int)
			d["thetafo_{0}".format(line)][:,m]=fid0.variables['thetao'][:,0,0]
			d["thetafo_{0}".format(line)][:,m][abs(d["thetafo_{0}".format(line)][:,m])>999999]=nan


#append ensmean
for var in ['barents','bering','canada','fram','denmark','icefar','farsco','fram_l','fram_r']:
	d['thetafo_{0}_orig'.format(var)]=zeros((2868,len(allmodel_list)))
	for n in range(2868) :
		d['thetafo_{0}_orig'.format(var)][n,:] = append(d["thetafo_{0}".format(var)][n,:], nanmean(d["thetafo_{0}".format(var)][n,:]))

#decompose in month and year
for var in ['barents','bering','canada','fram','denmark','icefar','farsco','fram_l','fram_r']:
	d['thetafo_{0}_orig_mon'.format(var)] = zeros((13,237,len(allmodel_list)))
	for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
		d['thetafo_{0}_orig_mon'.format(var)][k,:,:]=d['thetafo_{0}_orig'.format(var)][1+k:237*12+1:12,:]
	k=12
	d['thetafo_{0}_orig_mon'.format(var)][k,:]=nanmean(d['thetafo_{0}_orig_mon'.format(var)][0:12,:],axis=0)

#build the anomaly
for var in ['barents','bering','canada','fram','denmark','icefar','farsco','fram_l','fram_r']:
  d['thetafo_{0}_orig_mon_ano'.format(var)] = zeros((13,237,len(allmodel_list)))
  for n in range(237):
    d['thetafo_{0}_orig_mon_ano'.format(var)][:,n,:] = d['thetafo_{0}_orig_mon'.format(var)][:,n,:] - nanmean( d['thetafo_{0}_orig_mon'.format(var)][:,0:100,:] ,axis=1 )

################### SEA ICE EXPORT THROUGH FRAM STRAIT

#sea ice export through Fram Strait
d["exp_mod"]=zeros((2868,len(model_list)))
d["exp_mod2"]=zeros((2868,len(model_list)))

for i,mod in enumerate(model_list):
  print mod

  inputpath_time='/work/mh0033/m300411/DataEB/WORK_DATA/'
  for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
	if mod=='GISS-E2-R' or mod=='FGOALS-g2':
		file0=glob.glob(inputpath_time+mod+'/Method170516/transiyfs_%s_ensmean_186101-209912.nc' %mod)
		fid0=sio.netcdf_file(file0[0])
		d["exp_mod"][:,i]=fid0.variables['transiy'][:,0,0]
	else:
		file0=glob.glob(inputpath_time+mod+'/hist_rcp2/transifs_%s_ensmean_186101-209912.nc' %mod)
		if not file0:#or mod=='GFDL-ESM2G' : #or mod=='MPI-ESM-LR':# or mod=='MRI-CGCM3' or mod=='GFDL-ESM2M':
			d["exp_mod"][:,i] = nan
			print 'prob '+mod
		else:
			fid0=sio.netcdf_file(file0[0])
			d["exp_mod"][:,i]=fid0.variables['transifs'][:]
			if mod=='CNRM-CM5' or mod=='NorESM1-M' or mod=='NorESM1-ME':
				d["exp_mod"][:,i] = -d["exp_mod"][:,i]
	if mod=='MPI-ESM-LR':
		time0=fid0.variables['time'][:]
		year0=(time0/365)+1850
		mon0=(year0-floor(year0))*12
		day0=(mon0-floor(mon0))*30
		date0=transpose(array([floor(year0),ceil(mon0),floor(day0)]))
		date_new=date0.astype(int)    
	d["exp_mod"][:,i][abs(d["exp_mod"][:,i])>10**30]=nan
	d["exp_mod2"][:,i]=d["exp_mod"][:,i]*334774.
#	d["exp_mod"][:,i]=d["exp_mod"][:,i]-nanmean(d["exp_mod"][tref:tref+30,i])
   


d['exp_ensmean']=nanmean(d['exp_mod'][:,:],axis=1)
d['exp2_ensmean']=nanmean(d['exp_mod2'][:,:],axis=1)

##CHECK
#figure()
#for m,mod in enumerate(model_list):
#	#figure()
#	plot(d['exp_mod2'][:,m])
#	legend(loc='best')
#	title(mod)

#append ensmean
for var in ['exp']:
	d['{0}_orig'.format(var)]=zeros((2868,len(allmodel_list)))
	for n in range(2868) :
		d['{0}_orig'.format(var)][n,:] = append(d['{0}_mod2'.format(var)][n,:], nanmean(d['{0}_mod2'.format(var)][n,:]))

#decompose in month and year
for var in ['exp']:
	d['{0}_orig_mon'.format(var)] = zeros((13,237,len(allmodel_list)))
	for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
		d['{0}_orig_mon'.format(var)][k,:,:]=d['{0}_orig'.format(var)][1+k:237*12+1:12,:]
	k=12
	d['{0}_orig_mon'.format(var)][k,:]=nanmean(d['{0}_orig_mon'.format(var)][0:12,:],axis=0)

#build the anomaly
for var in ['exp']:
  d['{0}_orig_mon_ano'.format(var)] = zeros((13,237,len(allmodel_list)))
  for n in range(237):
    d['{0}_orig_mon_ano'.format(var)][:,n,:] = d['{0}_orig_mon'.format(var)][:,n,:] - nanmean( d['{0}_orig_mon'.format(var)][:,0:100,:] ,axis=1 )

#######################################################

################################ PRESSURE BARENTS SEA #################################################

d["ps_mod"]=zeros((2868,len(model_list)))
d["ps_mod2"]=zeros((2868,len(model_list)))


for i,mod in enumerate(model_list):
	print mod

	inputpath_time='/work/mh0033/m300411/DataEB/WORK_DATA/'
	file0=glob.glob(inputpath_time+mod+'/Method170516/ps_BarentsSea_remap_fldmean_%s_ensmean_186101-209912.nc' %mod)
	if not file0:
		d["ps_mod"][:,i] = nan
		print 'prob '+mod
	else:
		fid0=sio.netcdf_file(file0[0])
		if mod=='MPI-ESM-LR':
			time0=fid0.variables['time'][:]
			year0=(time0/365)+1850
			mon0=(year0-floor(year0))*12
			day0=(mon0-floor(mon0))*30
			date0=transpose(array([floor(year0),ceil(mon0),floor(day0)]))
			date_new=date0.astype(int)
	d["ps_mod"][:,i]=fid0.variables['ps'][:,0,0]/100
	d["ps_mod2"][:,i]=d["ps_mod"][:,i]
#	d["ps_mod"][:,i]=d["ps_mod"][:,i]-nanmean(d["ps_mod"][tref:tref+30,i])
   
d['ps_ensmean']=nanmean(d['ps_mod'][:,:],axis=1)
d['ps_ensmean2']=nanmean(d['ps_mod2'][:,:],axis=1)

##CHECK
#figure()
#plot(d['ps_mod'][:,:])
#legend(loc='best')

#append ensmean
for var in ['ps']:
	d['{0}_orig'.format(var)]=zeros((2868,len(allmodel_list)))
	for n in range(2868) :
		d['{0}_orig'.format(var)][n,:] = append(d["{0}_mod2".format(var)][n,:], nanmean(d["{0}_mod2".format(var)][n,:]))

#decompose in month and year
for var in ['ps']:
	d['{0}_orig_mon'.format(var)] = zeros((13,237,len(allmodel_list)))
	for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
		d['{0}_orig_mon'.format(var)][k,:,:]=d['{0}_orig'.format(var)][1+k:237*12+1:12,:]
	k=12
	d['{0}_orig_mon'.format(var)][k,:]=nanmean(d['{0}_orig_mon'.format(var)][0:12,:],axis=0)

#build the anomaly
for var in ['ps']:
  d['{0}_orig_mon_ano'.format(var)] = zeros((13,237,len(allmodel_list)))
  for n in range(237):
    d['{0}_orig_mon_ano'.format(var)][:,n,:] = d['{0}_orig_mon'.format(var)][:,n,:] - nanmean( d['{0}_orig_mon'.format(var)][:,0:100,:] ,axis=1 )

###################################

################### CONVERT EVERYTHING TO ENERGY

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
E_am = zeros((2868,len(allmodel_list)))

for n in range(2868) :
  E_a[n,:] = append(d["sum_atmos_mod2"][n,:]*moninsec*d['AOarea_ocean_mod'][:-1:], nanmean(d["sum_atmos_mod2"][n,:]*moninsec*d['AOarea_ocean_mod'][-1]))
  E_am[n,:] = append(d["atmmer_mod2"][n,:]*moninsec*d['AOarea_ocean_mod'][:-1:], nanmean(d["atmmer_mod2"][n,:]*moninsec*d['AOarea_ocean_mod'][-1]))


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

var_name = ['E_a','E_am']
for i,var in enumerate([E_a, E_am]):
  d['{0}_mon'.format(var_name[i])] = zeros((13,237,len(allmodel_list)))
  for m,mod in enumerate(model_list):
    print mod
    for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
      d['{0}_mon'.format(var_name[i])][k,:,:]=var[1+k:237*12+1:12,:]

#annual value (in J)
k=12
var_name = ['E_a', 'E_am', 'H_oc', 'H_ic']
for i,var in enumerate([E_a, E_am, H_oc, H_ic]):
  d['{0}_mon'.format(var_name[i])][k,:]=nansum(d['{0}_mon'.format(var_name[i])][0:12,:],axis=0)


d['H_tot_mon']=d['H_oc_mon']+d['H_ic_mon']

################### Anomalies ################################

var_name = ['E_a', 'E_am', 'H_oc', 'H_ic', 'H_tot']
for i,var in enumerate(var_name):
  d['{0}_mon_ano'.format(var_name[i])] = zeros((13,237,len(allmodel_list)))
  for n in range(237):
    d['{0}_mon_ano'.format(var_name[i])][:,n,:] = d['{0}_mon'.format(var_name[i])][:,n,:] - nanmean( d['{0}_mon'.format(var_name[i])][:,0:100,:] ,axis=1 )

d['E_o_mon_ano'] = d['H_tot_mon_ano']-d['E_a_mon_ano']

#################################################

group_atm = ['ACCESS1-0','ACCESS1-3','CCSM4','CESM1-BGC','CESM1-CAM5','CSIRO-Mk3-6-0','GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-LR','IPSL-CM5B-LR','MIROC5']
group_lat = ['CanESM2','CMCC-CM','CMCC-CMS','CNRM-CM5','GISS-E2-R','HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-MR','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3']
group_both = ['FGOALS-g2','GFDL-CM3','NorESM1-M','NorESM1-ME']

#################### SCATTER CHANGES

#append multi-model ensemble mean
for i,var in enumerate(variables_all):
  d['{0}'.format(var)] = zeros((2868,len(allmodel_list)))
  for n in range(2868) :
    d['{0}'.format(var)][n,:] = append( d['{0}_mod2'.format(var)][n,:], nanmean( d['{0}_mod2'.format(var)][n,:] ))

var = 'heat_content'
d['{0}'.format(var)] = zeros((2868,len(allmodel_list)))
for n in range(2868) :
  d['{0}'.format(var)][n,:] = append( heat_content2[n,:], nanmean( heat_content2[n,:] ))

variables_all.append('heat_content')

#create variable = sum of two radiative fluxes
d['radtot_mod2'] = d['netsw_mod2'] + d['netlw_mod2']
var = 'radtot'
d['{0}'.format(var)] = zeros((2868,len(allmodel_list)))
for n in range(2868) :
	d['{0}'.format(var)][n,:] = append( d['{0}_mod2'.format(var)][n,:], nanmean( d['{0}_mod2'.format(var)][n,:] ))

variables_all.append('radtot')

variables_all.append('vert_grad')
var = 'vert_grad'
d['{0}'.format(var)] = zeros((2868,len(allmodel_list)))
for n in range(2868) :
	d['{0}'.format(var)][n,:] = append( d['{0}_mod2'.format(var)][n,:], nanmean( d['{0}_mod2'.format(var)][n,:] ))


for i,var in enumerate(variables_all):
  d['{0}_mon'.format(var)] = zeros((13,237,len(allmodel_list)))
  d['diff_{0}'.format(var)] = zeros((13,len(allmodel_list)))
  for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
      d['{0}_mon'.format(var)][k,:,:] = d['{0}'.format(var)][1+k:237*12+1:12,:]
      d['diff_{0}'.format(var)][k,:] = nanmean( d['{0}_mon'.format(var)][k,207:237,:], axis=0) - nanmean(d['{0}_mon'.format(var)][k,0:100,:], axis=0 )
  d['{0}_mon'.format(var)][12,:,:]  = nanmean( d['{0}_mon'.format(var)][0:12,:,:], axis=0 )
  d['diff_{0}'.format(var)][12,:] = nanmean( d['{0}_mon'.format(var)][12,207:237,:], axis=0) - nanmean(d['{0}_mon'.format(var)][12,0:100,:], axis=0 )

#var_simple = variables_all.pop(-3)
  
#correlation coefficient
variables_ordered = ['rsds','rsus','netsw','rlds','rlus','netlw','radtot','hfss','hfls','turb','atmmer','sum_atmos','totsia','totsiv','heat_content']
for var in variables_ordered:
#  figure() 
#  for m,mod in enumerate(allmodel_list):
#    if mod in group_atm:
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_sum_atmos'][12,m],edgecolors='None',c='r')
#    elif mod in group_lat:
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_sum_atmos'][12,m],edgecolors='None',c='c')
#    elif mod in group_both:
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_sum_atmos'][12,m],edgecolors='None',c='grey')
#    elif mod=='Ensemble mean':
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_sum_atmos'][12,m],edgecolors='None',c='k')
#  title(var)
#  xlabel(var)
#  ylabel('Sum Atmos')
  
  fit = np.polyfit(d['diff_{0}'.format(var)][12,:],d['diff_sum_atmos'][12,:], deg=1)
#  plot(arange(min(d['diff_{0}'.format(var)][12,:]),max(d['diff_{0}'.format(var)][12,:])), fit[0] * arange(min(d['diff_{0}'.format(var)][12,:]),max(d['diff_{0}'.format(var)][12,:])) + fit[1], color='black')
  #plot(sum(d['E_o_mon_ano'.format(var)][k,100::,:]/10**23,axis=0), fit[0] * sum(d['E_o_mon_ano'.format(var)][k,100::,:]/10**23,axis=0) + fit[1]/10**23, color='black')
  slope, intercept, r_value, p_value, std_err = stats.linregress(d['diff_{0}'.format(var)][12,:],d['diff_sum_atmos'][12,:])
#  legend(loc=1, scatterpoints = 1,ncol=2,fontsize=18)
#  if var=='vert_grad':
#    ccoeff=corrcoef(concatenate((d['diff_{0}'.format(var)][12,0:10],d['diff_{0}'.format(var)][12,11:16],d['diff_{0}'.format(var)][12,17::]),axis=0),concatenate((d['diff_sum_atmos'][12,0:10],d['diff_sum_atmos'][12,11:16],d['diff_sum_atmos'][12,17::]),axis=0))[0,1]  
#  else:
#    ccoeff=corrcoef(d['diff_{0}'.format(var)][12,:],d['diff_sum_atmos'][12,:])[0,1]
  #print 'corrcoeff '+var+' ='+str(ccoeff)
  print 'r_value '+var+' ='+str(round(r_value,2))
  print 'p_value '+var+' ='+str(round(p_value,3))
  
##CHECK  
#for var in variables_all:
#  #for m,mod in enumerate(allmodel_list):
#    figure()
#    #plot(d['{0}'.format(var)][:,:])
#    plot(d['{0}_mon'.format(var)][12,:,:])
#    plot(d['diff_{0}'.format(var)][12,:])

#for var in ['totsia','totsiv','heat_content']:
#	d['diff_{0}'.format(var)][12,:] = d['diff_{0}'.format(var)][12,:]/d['start_{0}'.format(var)][12,:]
var_list = ['rsds','rsus','netsw','rlds','rlus','netlw','radtot','hfss','hfls','turb','atmmer','sum_atmos','totsia','totsiv','heat_content']#,'vert_grad','tas','psxprw','prw','rad_toa']
for m,mod in enumerate(allmodel_list):
  d['diff_{0}'.format(mod)] = array([d['diff_rsds'][:,m],d['diff_rsus'][:,m],d['diff_netsw'][:,m],d['diff_rlds'][:,m],d['diff_rlus'][:,m],d['diff_netlw'][:,m],d['diff_radtot'][:,m],d['diff_hfss'][:,m],d['diff_hfls'][:,m],d['diff_turb'][:,m],d['diff_atmmer'][:,m],d['diff_sum_atmos'][:,m],d['diff_totsia'][:,m],d['diff_totsiv'][:,m],d['diff_heat_content'][:,m]])#,d['diff_vert_grad'][:,m],d['diff_tas'][:,m],d['diff_psxprw'][:,m]/10**5,d['diff_prw'][:,m]])


 
###### Figure 2b is divided in several parts
matplotlib.rcParams.update({'font.size': 24})  
figure()

fig1 = matplotlib.pyplot.gcf()
fig1.set_size_inches(18,12)
fig1.subplots_adjust(bottom=0.08)
fig1.subplots_adjust(left=0.18)
fig1.subplots_adjust(right=0.95)
fig1.subplots_adjust(top=0.92)
ax = fig1.add_subplot(111) 

#xlabel("Component")
ylabel("$\Delta$ Flux [W/m$^2$]")

comp_list1=['SW Down','SW Up','Net SW','LW Down','LW Up','Net LW','Net SW+Net LW','Sensible HF','Latent HF','Turbulent HF','Atm mer','Net Atmosph. Flux']#,'TOA','Vertical T grad','Surface Air Temp','Water vap press']

width=0.2
size=150

for m,mod in enumerate(allmodel_list):
    if mod in group_atm :
       scatter(range(len(comp_list1[0:12])), d['diff_{0}'.format(mod)][0:12,12], c='r', edgecolors='none', s=size)
       #scatter(range(13,16), d['diff_{0}'.format(mod)][16:19,12], c='r', edgecolors='none', s=size)
    elif mod in group_lat :
       scatter(arange(len(comp_list1[0:12]))+width, d['diff_{0}'.format(mod)][0:12,12], color='turquoise', edgecolors='none', s=size)	
       #scatter(arange(13,16)+width, d['diff_{0}'.format(mod)][16:19,12], color='turquoise', edgecolors='none', s=size)	
    elif mod in group_both :
       scatter(arange(len(comp_list1[0:12]))+2*width, d['diff_{0}'.format(mod)][0:12,12], c='grey', edgecolors='none', s=size)
       #scatter(arange(13,16)+2*width, d['diff_{0}'.format(mod)][16:19,12], c='grey', edgecolors='none', s=size)
    elif mod == 'Ensemble mean' :
       scatter(arange(len(comp_list1[0:12]))+3*width, d['diff_{0}'.format(mod)][0:12,12], c='k', edgecolors='none', s=size)
       #scatter(arange(13,16)+3*width, d['diff_{0}'.format(mod)][16:19,12], c='k', edgecolors='none', s=size)
axhline(xmin=0,xmax=len(comp_list1),y=0, c='k',linewidth=2)

ylim(-30, 60)
xlim(-1,len(comp_list1))

grid()
ax.set_xticks( arange(len(comp_list1))+1.5*width)
ax.set_xticklabels( comp_list1, rotation=90)
#ax.axis["bottom"].major_ticklabels.set_axis_direction("left")

plt.draw()
plt.show()
#savefig(outputpath+'changes_scatter_atmofluxes_rev60.pdf', bbox_inches='tight', orientation='landscape')

#######
matplotlib.rcParams.update({'font.size': 24})  
figure()

fig1 = matplotlib.pyplot.gcf()
fig1.set_size_inches(5,12)
fig1.subplots_adjust(bottom=0.08)
fig1.subplots_adjust(left=0.18)
fig1.subplots_adjust(right=0.95)
fig1.subplots_adjust(top=0.92)
ax = fig1.add_subplot(111) 

#xlabel("Component")
ylabel("$\Delta$ SIA [x 10$^{12}$ m$^2$] or $\Delta$ SIV [x 10$^{13}$ m$^3$] or $\Delta$ OHC [x 10$^{22}$ J]")


comp_list2=['Sea Ice Area','Sea Ice Volume','Oceanic Sensible Heat Content']

width=0.2
size=150

for m,mod in enumerate(allmodel_list):
	if mod in group_atm :
	  scatter(range(0,1), d['diff_{0}'.format(mod)][12:13,12]/10**12, c='r', edgecolors='none', s=size)
	  scatter(range(1,2), d['diff_{0}'.format(mod)][13:14,12]/10**13, c='r', edgecolors='none', s=size)
	  scatter(range(2,3), d['diff_{0}'.format(mod)][14:15,12]/10**22, c='r', edgecolors='none', s=size)
	elif mod in group_lat :
	  scatter(arange(0,1)+width, d['diff_{0}'.format(mod)][12:13,12]/10**12, color='turquoise', edgecolors='none', s=size)
	  scatter(arange(1,2)+width, d['diff_{0}'.format(mod)][13:14,12]/10**13, color='turquoise', edgecolors='none', s=size)
	  scatter(arange(2,3)+width, d['diff_{0}'.format(mod)][14:15,12]/10**22, color='turquoise', edgecolors='none', s=size)	

	elif mod in group_both :
	  scatter(arange(0,1)+2*width, d['diff_{0}'.format(mod)][12:13,12]/10**12, c='grey', edgecolors='none', s=size)
	  scatter(arange(1,2)+2*width, d['diff_{0}'.format(mod)][13:14,12]/10**13, c='grey', edgecolors='none', s=size)
	  scatter(arange(2,3)+2*width, d['diff_{0}'.format(mod)][14:15,12]/10**22, c='grey', edgecolors='none', s=size)		
	elif mod == 'Ensemble mean' :
	  scatter(arange(0,1)+3*width, d['diff_{0}'.format(mod)][12:13,12]/10**12, c='k', edgecolors='none', s=size)
	  scatter(arange(1,2)+3*width, d['diff_{0}'.format(mod)][13:14,12]/10**13, c='k', edgecolors='none', s=size)
	  scatter(arange(2,3)+3*width, d['diff_{0}'.format(mod)][14:15,12]/10**22, c='k', edgecolors='none', s=size)	

axhline(xmin=0,xmax=len(comp_list2),y=0, c='k',linewidth=2)


xlim(-0.5,len(comp_list2))
ylim(-6, 12)

grid()
ax.set_xticks( arange(len(comp_list2))+1.5*width)
ax.set_xticklabels( comp_list2, rotation=90)
#host.axis["bottom"].major_ticklabels.set_axis_direction("left")

plt.draw()
plt.show()
#savefig(outputpath+'changes_scatter_atmoother.pdf', bbox_inches='tight', orientation='landscape')

##### prepare Figure 2a

d['E_o_orig_mon'] = d['H_tot_mon'] - d['E_a_mon']
d['E_o_W_orig_mon'] = d['E_o_orig_mon']/moninsec
d['E_o_W_orig_mon'][12,:] = d['E_o_orig_mon'][12,:]/yearinsec


for var in ['E_o_W']:
  d['{0}_orig_mon_ano'.format(var)] = zeros((13,237,len(allmodel_list)))
  for n in range(237):
    d['{0}_orig_mon_ano'.format(var)][:,n,:] = d['{0}_orig_mon'.format(var)][:,n,:] - nanmean( d['{0}_orig_mon'.format(var)][:,0:100,:] ,axis=1 )

    
######################## ONE PLOT WITH ALL CAUSES
AAmp_ocean = (nanmean(d['sst_orig_mon_ano'][12,209:239,:],axis=0)) / (nanmean(d['sstglob_orig_mon_ano'][12,209:239,:],axis=0))


variables_tot=['mfo_barents','mfo_bering','mfo_fram','mfo_canada','mfo_denmark','mfo_icefar','mfo_farsco','thetafo_barents','thetafo_bering','thetafo_fram','thetafo_canada','thetafo_denmark','thetafo_icefar','thetafo_farsco','thetafo_fram_l','thetafo_fram_r','gradient_low_arc_atm','gradient_low_arc_oc','sst','AAmp','AAmp_ocean','icep','exp','E_o_W','ps']

for i,var in enumerate(variables_tot):
  d['diff_{0}'.format(var)] = zeros((13,len(allmodel_list)))
  for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','year']):
  	if var == 'icep' :
	      d['diff_{0}'.format(var)][k,:] = nanmean( d['{0}_orig_mon_ano'.format(var)][207:237,:], axis=0)
	elif var=='AAmp':
	      d['diff_{0}'.format(var)][k,:] = AAmp
	elif var=='AAmp_ocean':
	      d['diff_{0}'.format(var)][k,:] = AAmp_ocean
  	else : 
	      d['diff_{0}'.format(var)][k,:] = nanmean( d['{0}_orig_mon_ano'.format(var)][k,207:237,:], axis=0)
  
#correlation coefficient between ocean variables and Flat
var_order=['thetafo_barents','thetafo_bering','thetafo_fram','thetafo_fram_l','thetafo_fram_r','thetafo_canada','thetafo_denmark','thetafo_icefar','thetafo_farsco','gradient_low_arc_atm','gradient_low_arc_oc','sst','AAmp','icep','mfo_barents','mfo_bering','mfo_fram','mfo_canada','ps','exp']
for var in var_order:
#  figure() 
#  for m,mod in enumerate(allmodel_list):
#    if mod in group_atm:
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_E_o_W'][12,m],edgecolors='None',c='r')
#    elif mod in group_lat:
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_E_o_W'][12,m],edgecolors='None',c='c')
#    elif mod in group_both:
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_E_o_W'][12,m],edgecolors='None',c='grey')
#    elif mod=='Ensemble mean':
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_E_o_W'][12,m],edgecolors='None',c='k')
#  title(var)
#  xlabel(var)
#  ylabel('E_o_W')
  
  x=d['diff_{0}'.format(var)][12,:]
  y=d['diff_E_o_W'][12,:]
  idx = np.isfinite(x) & np.isfinite(y)
  fit = np.polyfit(x[idx], y[idx], 1)
  #fit = np.polyfit(d['diff_{0}'.format(var)][12,:],d['diff_E_o_W'][12,:], deg=1)
  #plot(arange(min(d['diff_{0}'.format(var)][12,:]),max(d['diff_{0}'.format(var)][12,:])), fit[0] * arange(min(d['diff_{0}'.format(var)][12,:]),max(d['diff_{0}'.format(var)][12,:])) + fit[1], color='black')
  #plot(sum(d['E_o_mon_ano'.format(var)][k,100::,:]/10**23,axis=0), fit[0] * sum(d['E_o_mon_ano'.format(var)][k,100::,:]/10**23,axis=0) + fit[1]/10**23, color='black')
  slope, intercept, r_value, p_value, std_err = stats.linregress(x[idx],y[idx])
#  legend(loc=1, scatterpoints = 1,ncol=2,fontsize=18)
  
  ccoeff=corrcoef(d['diff_{0}'.format(var)][12,:],d['diff_E_o_W'][12,:])[0,1]
  #print 'corrcoeff '+var+' ='+str(ccoeff)
  print 'r_value '+var+' ='+str(round(r_value,2))
  #print 'R**2 '+var+' ='+str(r_value**2)
  print 'p_value '+var+' ='+str(round(p_value,3))
  
#correlation coefficient between ocean variables and Fatm  
var_order=['thetafo_barents','thetafo_bering','thetafo_fram','thetafo_fram_l','thetafo_fram_r','thetafo_canada','thetafo_denmark','thetafo_icefar','thetafo_farsco','gradient_low_arc_atm','gradient_low_arc_oc','sst','AAmp','icep','mfo_barents','mfo_bering','mfo_fram','mfo_canada','ps','exp']
for var in var_order:
#  figure() 
#  for m,mod in enumerate(allmodel_list):
#    if mod in group_atm:
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_sum_atmos'][12,m],edgecolors='None',c='r')
#    elif mod in group_lat:
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_sum_atmos'][12,m],edgecolors='None',c='c')
#    elif mod in group_both:
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_sum_atmos'][12,m],edgecolors='None',c='grey')
#    elif mod=='Ensemble mean':
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_sum_atmos'][12,m],edgecolors='None',c='k')
#  title(var)
#  xlabel(var)
#  ylabel('Sum Atmos')
  
  x=d['diff_{0}'.format(var)][12,:]
  y=d['diff_sum_atmos'][12,:]
  idx = np.isfinite(x) & np.isfinite(y)
  fit = np.polyfit(x[idx], y[idx], 1)
  slope, intercept, r_value, p_value, std_err = stats.linregress(x[idx],y[idx])
#  legend(loc=1, scatterpoints = 1,ncol=2,fontsize=18)
  
  ccoeff=corrcoef(d['diff_{0}'.format(var)][12,:],d['diff_sum_atmos'][12,:])[0,1]
  #print 'corrcoeff '+var+' ='+str(ccoeff)
  print 'r_value '+var+' ='+str(round(r_value,2))
  #print 'R**2 '+var+' ='+str(r_value**2)
  print 'p_value '+var+' ='+str(round(p_value,3))

#correlation coefficient between atmos variables and Flat
for var in variables_ordered:
#  figure() 
#  for m,mod in enumerate(allmodel_list):
#    if mod in group_atm:
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_E_o_W'][12,m],edgecolors='None',c='r')
#    elif mod in group_lat:
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_E_o_W'][12,m],edgecolors='None',c='c')
#    elif mod in group_both:
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_E_o_W'][12,m],edgecolors='None',c='grey')
#    elif mod=='Ensemble mean':
#      scatter(d['diff_{0}'.format(var)][12,m],d['diff_E_o_W'][12,m],edgecolors='None',c='k')
#  title(var)
#  xlabel(var)
#  ylabel('E_o_W')
  
  x=d['diff_{0}'.format(var)][12,:]
  y=d['diff_E_o_W'][12,:]
  idx = np.isfinite(x) & np.isfinite(y)
  fit = np.polyfit(x[idx], y[idx], 1)
  slope, intercept, r_value, p_value, std_err = stats.linregress(x[idx],y[idx])
#  legend(loc=1, scatterpoints = 1,ncol=2,fontsize=18)
  
  ccoeff=corrcoef(d['diff_{0}'.format(var)][12,:],d['diff_E_o_W'][12,:])[0,1]
  #print 'corrcoeff '+var+' ='+str(ccoeff)
  print 'r_value '+var+' ='+str(round(r_value,2))
  #print 'R**2 '+var+' ='+str(r_value**2)
  print 'p_value '+var+' ='+str(round(p_value,3))
  
#for var in variables_tot:	
#    figure()
#	#title(var)
#    plot(d['diff_{0}'.format(var)][12,:],label=var)
#    legend()

#var_list2 = ['mfo_barents','mfo_bering','mfo_fram','mfo_canada','mfo_denmark','mfo_icefar','mfo_farsco','thetafo_barents','thetafo_bering','thetafo_fram','thetafo_canada','thetafo_denmark','thetafo_icefar','thetafo_farsco','thetafo_fram_l','thetafo_fram_r','gradient_low_arc_atm','gradient_low_arc_oc','sst','AAmp','AAmp_ocean','icep','exp','E_o_W','ps']
var_list2 = ['mfo_barents','mfo_bering','mfo_fram','mfo_canada','mfo_denmark','mfo_icefar','mfo_farsco','thetafo_barents','thetafo_bering','thetafo_fram','thetafo_fram_l','thetafo_fram_r','thetafo_canada','thetafo_denmark','thetafo_icefar','thetafo_farsco','gradient_low_arc_atm','gradient_low_arc_oc','sst','AAmp','icep','exp','E_o_W','ps']

for m,mod in enumerate(allmodel_list):
  d['diff_{0}'.format(mod)] = array([d['diff_mfo_barents'][:,m],d['diff_mfo_bering'][:,m],d['diff_mfo_fram'][:,m],d['diff_mfo_canada'][:,m],d['diff_mfo_denmark'][:,m],d['diff_mfo_icefar'][:,m],d['diff_mfo_farsco'][:,m],d['diff_thetafo_barents'][:,m],d['diff_thetafo_bering'][:,m],d['diff_thetafo_fram'][:,m],d['diff_thetafo_fram_l'][:,m],d['diff_thetafo_fram_r'][:,m],d['diff_thetafo_canada'][:,m],d['diff_thetafo_denmark'][:,m],d['diff_thetafo_icefar'][:,m],d['diff_thetafo_farsco'][:,m],d['diff_gradient_low_arc_atm'][:,m],d['diff_gradient_low_arc_oc'][:,m],d['diff_sst'][:,m],d['diff_AAmp'][:,m],d['diff_icep'][:,m],-d['diff_exp'][:,m],d['diff_E_o_W'][:,m],d['diff_ps'][:,m]])

#figure()
#scatter(d['diff_E_o_W'][12,:],d['diff_sum_atmos'][12,:])
#
#figure()
#scatter(d['diff_E_o_W'][12,:],d['diff_atmmer'][12,:])


#################

  
group_atm = ['ACCESS1-0','ACCESS1-3','CCSM4','CESM1-BGC','CESM1-CAM5','CSIRO-Mk3-6-0','GFDL-ESM2G','GFDL-ESM2M','IPSL-CM5A-LR','IPSL-CM5B-LR','MIROC5']
group_lat = ['CanESM2','CMCC-CM','CMCC-CMS','CNRM-CM5','GISS-E2-R','HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-MR','MPI-ESM-LR','MPI-ESM-MR','MRI-CGCM3']
group_both = ['FGOALS-g2','GFDL-CM3','NorESM1-M','NorESM1-ME']


########## divide into three different plots

matplotlib.rcParams.update({'font.size': 26})  
figure()

fig1 = matplotlib.pyplot.gcf()
fig1.set_size_inches(9,15)
fig1.subplots_adjust(bottom=0.08)
fig1.subplots_adjust(left=0.18)
fig1.subplots_adjust(right=0.95)
fig1.subplots_adjust(top=0.92)

ax = fig1.add_subplot(111) 

#xlabel("Component")
ylabel("$\Delta$ Mass transport [x 10$^9$ kg/s] or Surface pressure Barents Sea [hPa]")


comp_list1=['Mass Barents','Mass Bering','Mass Fram','Mass Canada', 'Surface Pressure Barents Sea']

width=0.2
size=150
for m,mod in enumerate(allmodel_list):
#	if mod in model_mfo:
		if mod in group_atm :
		  scatter(arange(0,4), d['diff_{0}'.format(mod)][0:4,12]/10**9, c='r', edgecolors='none',s=size)	
		  scatter(arange(4,5), d['diff_{0}'.format(mod)][23,12], c='r', edgecolors='none', s=size)
		elif mod in group_lat :
		  scatter(arange(0,4)+width, d['diff_{0}'.format(mod)][0:4,12]/10**9, c='turquoise', edgecolors='none', s=size)
		  scatter(arange(4,5)+width, d['diff_{0}'.format(mod)][23,12], c='turquoise', edgecolors='none', s=size)
		elif mod in group_both :
		  scatter(arange(0,4)+2*width, d['diff_{0}'.format(mod)][0:4,12]/10**9, c='grey', edgecolors='none', s=size)
		  scatter(arange(4,5)+2*width, d['diff_{0}'.format(mod)][23,12], c='grey', edgecolors='none', s=size)
		elif mod == 'Ensemble mean' :
		  scatter(arange(0,4)+3*width, d['diff_{0}'.format(mod)][0:4,12]/10**9, c='k', edgecolors='none', s=size)
		  scatter(arange(4,5)+3*width, d['diff_{0}'.format(mod)][23,12], c='k', edgecolors='none', s=size)

axhline(xmin=0,xmax=len(comp_list1),y=0, c='k',linewidth=2)

ylim(-3, 3)
xlim(-1,len(comp_list1))

grid()
ax.set_xticks(arange(len(comp_list1))+1.5*width)
ax.set_xticklabels(comp_list1,rotation=90)

plt.draw()
plt.show()
#savefig(outputpath+'changes_oc_scatter_mass_ps.pdf', bbox_inches='tight', orientation='landscape')

####


matplotlib.rcParams.update({'font.size': 26})  
figure()

fig1 = matplotlib.pyplot.gcf()
fig1.set_size_inches(22,15)
fig1.subplots_adjust(bottom=0.08)
fig1.subplots_adjust(left=0.18)
fig1.subplots_adjust(right=0.95)
fig1.subplots_adjust(top=0.92)

ax = fig1.add_subplot(111) 

#xlabel("Component")
ylabel("$\Delta$ Temperature [K] or Arctic Amplification or $\Delta$ Ice production [x 10$^{12}$ m$^3$/year]")

comp_list2=['Temp Barents','Temp Bering','Temp Fram','Temp Fram W','Temp Fram E','Temp Canada','Temp Denmark','Temp IceFar','Temp FarSco','Gradient Atmo','Gradient Ocean','Sea Surface Temperature','Arctic Amplification','Ice production']

width=0.2
size=150
for m,mod in enumerate(allmodel_list):
#	if mod in model_mfo:
        if mod in group_atm :
          scatter(arange(0,13), d['diff_{0}'.format(mod)][7:20,12], c='r', edgecolors='none', s=size)
          scatter(arange(13,14), d['diff_{0}'.format(mod)][20,12]/10**12, c='r', edgecolors='none', s=size)
#          scatter(arange(11,12), d['diff_{0}'.format(mod)][21:22,12]/10**12, c='r', edgecolors='none', s=size)
#          scatter(arange(12,15), d['diff_{0}'.format(mod)][18:21,12], c='r', edgecolors='none', s=size)
        elif mod in group_lat :
          scatter(arange(0,13)+width, d['diff_{0}'.format(mod)][7:20,12], c='turquoise', edgecolors='none', s=size)
          scatter(arange(13,14)+width, d['diff_{0}'.format(mod)][20,12]/10**12, c='turquoise', edgecolors='none', s=size)
#          scatter(arange(11,12)+width, d['diff_{0}'.format(mod)][21:22,12]/10**12, c='turquoise', edgecolors='none', s=size)
#          scatter(arange(12,15)+width, d['diff_{0}'.format(mod)][18:21,12], c='turquoise', edgecolors='none', s=size)
        elif mod in group_both :
		 scatter(arange(0,13)+2*width, d['diff_{0}'.format(mod)][7:20,12], c='grey', edgecolors='none', s=size)
		 scatter(arange(13,14)+2*width, d['diff_{0}'.format(mod)][20,12]/10**12, c='grey', edgecolors='none', s=size)
#		 scatter(arange(11,12)+2*width, d['diff_{0}'.format(mod)][21:22,12]/10**12, c='grey', edgecolors='none', s=size)			
#		 scatter(arange(12,15)+2*width, d['diff_{0}'.format(mod)][18:21,12], c='grey', edgecolors='none', s=size)			
        elif mod == 'Ensemble mean' :
		 scatter(arange(0,13)+3*width, d['diff_{0}'.format(mod)][7:20,12], c='k', edgecolors='none', s=size)
		 scatter(arange(13,14)+3*width, d['diff_{0}'.format(mod)][20,12]/10**12, c='k', edgecolors='none', s=size)
#		 scatter(arange(11,12)+3*width, d['diff_{0}'.format(mod)][21:22,12]/10**12, c='k', edgecolors='none', s=size)		
#		 scatter(arange(12,15)+3*width, d['diff_{0}'.format(mod)][18:21,12], c='k', edgecolors='none', s=size)			

    
axhline(xmin=0,xmax=len(comp_list2),y=0, c='k',linewidth=2)

xlim(-1,len(comp_list2))
ylim(-6, 6)

grid()
ax.set_xticks(arange(len(comp_list2))+1.5*width)
ax.set_xticklabels(comp_list2, rotation=90)

plt.draw()
plt.show()
#savefig(outputpath+'changes_oc_scatter_temp_rev.pdf', bbox_inches='tight', orientation='landscape')

#####

matplotlib.rcParams.update({'font.size': 26})  
figure()

fig1 = matplotlib.pyplot.gcf()
fig1.set_size_inches(4,15)
fig1.subplots_adjust(bottom=0.08)
fig1.subplots_adjust(left=0.18)
fig1.subplots_adjust(right=0.95)
fig1.subplots_adjust(top=0.92)

ax = fig1.add_subplot(111) 

#xlabel("Component")
ylabel("$\Delta$ Meridional total and latent heat flux [x 10$^{14}$ W]")

comp_list3=['Sea Ice Export Flux','Meridional oceanic heat flux']

width=0.2
size=150
for m,mod in enumerate(allmodel_list):
#	if mod in model_mfo:
		if mod in group_atm :
		  scatter(arange(0,2), d['diff_{0}'.format(mod)][21:23,12]/10**14, c='r', edgecolors='none', s=size)	
		elif mod in group_lat :
		  scatter(arange(0,2)+width, d['diff_{0}'.format(mod)][21:23,12]/10**14, c='turquoise', edgecolors='none', s=size)		
		elif mod in group_both :	
		  scatter(arange(0,2)+2*width, d['diff_{0}'.format(mod)][21:23,12]/10**14, c='grey', edgecolors='none', s=size)		
		elif mod == 'Ensemble mean' :
		  scatter(arange(0,2)+3*width, d['diff_{0}'.format(mod)][21:23,12]/10**14, c='k', edgecolors='none', s=size)		

axhline(xmin=0,xmax=len(comp_list3),y=0, c='k',linewidth=2)

xlim(-1,len(comp_list3))
#par2.set_ylim(-0.4, 1.2)


grid()
ax.set_xticks(arange(len(comp_list3))+1.5*width)
ax.set_xticklabels(comp_list3,rotation=90)

plt.draw()
plt.show()
#savefig(outputpath+'changes_oc_scatter_merhf.pdf', bbox_inches='tight', orientation='landscape')

###############################################



