#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Close the Arctic Ocean energy budget for MPI-ESM-LR

Last modification : 01.09.2016
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
model_list=['MPI-ESM-LR']

#open the dictionary that contains all variables
d={}

#Amount of seconds in the year for the conversion from W to J
yearinsec=3600*365*24.

#path where to store output
outputpath='/work/mh0033/m300411/DataEB/RESULTS/PAPER/Prelim2/'

################## READ IN ARCTIC OCEAN AREA (computed in AOmask.sh) ###################################

d["AOarea_ocean_mod"]=np.zeros((len(model_list)+1))

for i,mod in enumerate(model_list):
  print mod,'AO area'

  inputpath_time='/work/mh0033/m300411/DataEB/WORK_DATA/'
  
  file0=glob.glob(inputpath_time+'AOarea2/Arctic_ocean_totalarea_%s.nc' %mod)
  fid0=sio.netcdf_file(file0[0])
  d["AOarea_ocean_mod"][i]=fid0.variables['tos'][0]

d["AOarea_ocean_mod"][len(model_list)]=np.nanmean(d["AOarea_ocean_mod"][0:len(model_list)])

###############################################################################################################

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
variables=['rlds','rsds','rlus','rsus','hfss','hfls','totsiv','totsia','totsnv']
#Variables read in from a file + additional variables computed here
variables_all=['rlds','rsds','rlus','rsus','hfss','hfls','netlw','netsw','turb','sum_atmos','totsiv','totsia','totsnv']

#Create the arrays
for x in variables:
        d["{0}_mod".format(x)]=np.zeros((2868,len(model_list)))
        d["{0}_mod2".format(x)]=np.zeros((2868,len(model_list)))

#Read in the variables
for i,mod in enumerate(model_list):
  print mod
   
  inputpath_time='/work/mh0033/m300411/DataEB/WORK_DATA/'
  p=0
  for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
    for var in variables:
      print var
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
        #file0=glob.glob(inputpath_time+mod+'/hist_rcp/totsnv_%s_ensmean_186101-209912.nc' %(mod))
        file0=glob.glob(inputpath_time+mod+'/Method170516/totsnv_Arctic_fldmean_%s_ensmean_186101-209912.nc' %(mod))
      else:
        file0=glob.glob(inputpath_time+mod+'/Method170516/%s_Arctic_monfldmean_remap_%s_ensmean_186101-209912.nc' %(var,mod))
      if not file0:
        print mod+' '+var+' prob!'
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

      
    ##SUMS
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
  d["{0}_ensmean2".format(var)] = np.nanmean(d["{0}_mod2".format(var)],axis=1)
  d["{0}_ensstd2".format(var)] = np.nanstd(d["{0}_mod2".format(var)],axis=1)

## CHECK IF IT HAS BEEN READ IN WELL
#for var in variables_all:
#    plt.figure()
#    plt.plot(d["{0}_mod2".format(var)])
#    plt.plot(d["{0}_ensmean2".format(var)][:],'k-',linewidth=3)
#    plt.title(tit[var]+str(k))
##################################################################################################


########### COMPUTE OCEANIC HEAT CONTENT ###################################

#Create the array
d["txv_mod"] = np.zeros((2868,len(model_list))) #sum of (temperature x cell volume) over depth
d["txv_mod2"] = np.zeros((2868,len(model_list))) #sum of (temperature x cell volume) over depth
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



############# OCEAN HEAT TRANSPORT IN Y DIRECTION ############################################

d["hfx_mod"]=zeros((2868,len(model_list)))
d["hfx_mod2"]=zeros((2868,len(model_list)))

for i,mod in enumerate(model_list):
  print mod

  #Pfad ändern
  inputpath_time='/work/mh0033/m300411/DataEB/WORK_DATA/'
  for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
    file0=glob.glob(inputpath_time+'/hf_ocean/hf_tot_fldsum_mon_%s_ensmean_186101-209912.nc' %mod)
    fid0=sio.netcdf_file(file0[0])
    d["hfx_mod"][:,i]=fid0.variables['hfx'][:,0,0]
    d["hfx_mod2"][:,i]=d["hfx_mod"][:,i]

### CHECK IF IT HAS BEEN READ IN WELL
#plt.figure()
#plt.plot(d["hfx_mod2"][:,0])
#plt.legend(loc='best')

############################################################################

########### CONVERT EVERYTHING TO ENERGY

##### Constants
rho_i=910.    #ice density kg/m3
rho_sn=250.   #snow density kg/m3
lat_fus=334774.   #latent heat of fusion J/kg
yearinsec=3600*365*24.  #amount of seconds in a year
moninsec=3600*24*30.    #amount of seconds in a month
#####

seaice=-d["totsiv_mod2"]*rho_i*lat_fus #convert sea-ice volume to latent energy
snow=-d["totsnv_mod2"]*rho_sn*lat_fus  #convert snow volume to latent energy

#lateral (meridional) oceanic heat flux
lathf = d["hfx_mod2"]

######
#Absolute sources
# Energy exchanged through atmospheric and oceanic fluxes
E_a = (d["sum_atmos_mod2"][:,0]*moninsec*d['AOarea_ocean_mod'][0]) # in J now
E_o = (lathf[:,0]*moninsec) #in J now

######

######
# Absolute sinks
# Heat content : Ocean, Ice volume, Snow volume

#compute difference between one month and the month before
H_oc = heat_content2[1::,0] - heat_content2[0:-1:,0]
H_ic = seaice[1::,0] - seaice[0:-1:,0]
H_sn = snow[1::,0] - snow[0:-1:,0]


# to compare directly to the fluxes : take the value in the middle of the month
H_oc_comp = H_oc[0:-1:] + (H_oc[1::] - H_oc[0:-1:]) * 0.5 
H_ic_comp = H_ic[0:-1:] + (H_ic[1::] - H_ic[0:-1:]) * 0.5 
H_sn_comp = H_sn[0:-1:] + (H_sn[1::] - H_sn[0:-1:]) * 0.5 

######

#### test if it closes : YEP ! :)
#z=30*12
#plt.figure()
#plt.plot(E_a[1:z+1] + E_o[1:z+1], 'r-')
#plt.plot(H_oc_comp[0:z] + H_ic_comp[0:z] + H_sn_comp[0:z], 'b-')
#plt.grid()

###################################################

##### DECOMPOSE IN MONTHLY AND ANNUAL

var_name = ['H_oc', 'H_ic', 'H_sn'] #, 'H_ex']
for i,var in enumerate([H_oc, H_ic, H_sn]): #, H_ex]):
  d['{0}_mon'.format(var_name[i])] = zeros((13,237))
  for m,mod in enumerate(model_list):
    print mod
    for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
      if k==11:
	d['{0}_mon'.format(var_name[i])][k,:]=var[k:237*12:12]+(var[12:238*12:12]-var[k:237*12:12]) * 0.5
      else:
	d['{0}_mon'.format(var_name[i])][k,:]=var[k:237*12:12]+(var[k+1:237*12:12]-var[k:237*12:12]) * 0.5

    
var_name = ['E_a', 'E_o']
for i,var in enumerate([E_a, E_o]):
  d['{0}_mon'.format(var_name[i])] = zeros((13,237))
  for m,mod in enumerate(model_list):
    print mod
    for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
      d['{0}_mon'.format(var_name[i])][k,:]=var[1+k:237*12+1:12]

#monthly energy (in J)
flux_Wm2=['rsds','rlds','rsus','rlus','hfss','hfls','netsw','netlw','turb','sum_atmos']
for var in flux_Wm2:
  d['{0}_energy'.format(var)]=zeros((2868,len(model_list)))
  for n in range(2868) :
    d['{0}_energy'.format(var)][n,:] = d["{0}_mod2".format(var)][n,:]*moninsec*d['AOarea_ocean_mod'][:-1:]

#monthly energy fluxes    (in W/m2)   
for i,var in enumerate(flux_Wm2):
  d['{0}_mon'.format(var)] = zeros((13,237,len(model_list)))
  for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):
      d['{0}_mon'.format(var)][k,:,:]=d['{0}_energy'.format(var)][1+k:237*12+1:12,:]

#annual value (in J)
k=12
var_name = ['E_a', 'E_o', 'H_oc', 'H_ic', 'H_sn']
for i,var in enumerate([E_a, E_o, H_oc, H_ic, H_sn]):
  d['{0}_mon'.format(var_name[i])][k,:]=nansum(d['{0}_mon'.format(var_name[i])][0:12,:],axis=0)
for i,var in enumerate(flux_Wm2):
  d['{0}_mon'.format(var)][k,:]=nansum(d['{0}_mon'.format(var)][0:12,:],axis=0)

d['H_tot_mon']=d['H_oc_mon']+d['H_ic_mon']+d['H_sn_mon']

##### test if it works
##general formula
##for each month
#for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','year']):  
#  figure()
#  plot((d['E_a_mon'][k,:] + d['E_o_mon'][k,:]) / 10**21, 'g-')
#  plot((d['H_oc_mon'][k,:] + d['H_ic_mon'][k,:] + d['H_sn_mon'][k,:]) / 10**21, 'b-')
#  ylim(-4,4)
#  title(month)
#
##only part of it
#z=15
#for k,month in enumerate(['Jan','Feb', 'Mar', 'Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','year']):  
#  figure()
#  plot((d['E_a_mon'][k,0:z] + d['E_o_mon'][k,0:z]) / 10**21, 'g-')
#  plot((d['H_oc_mon'][k,0:z] + d['H_ic_mon'][k,0:z] + d['H_sn_mon'][k,0:z]) / 10**21, 'b-')
#  ylim(-4,4)
#  title(month)


#######################################################

###### COMPARE COMPUTED OCEANIC HEAT FLUX WITH REAL ONE
###### FIGURE 1 ######################################

k = 12
blax=range(-12,30)
blay=range(-12,30)

matplotlib.rcParams.update({'font.size': 28})  

fig, ax = plt.subplots()

fig1 = matplotlib.pyplot.gcf()
fig1.set_size_inches(22,15)
fig1.subplots_adjust(bottom=0.08)
fig1.subplots_adjust(left=0.05)
fig1.subplots_adjust(right=0.95)
fig1.subplots_adjust(top=0.92)

scatter((d['H_tot_mon'][k,:]-d['H_sn_mon'][k,:]-d['E_a_mon'][k,:])/(yearinsec*d['AOarea_ocean_mod'][0]),(d['E_o_mon'][k,:])/(yearinsec*d['AOarea_ocean_mod'][0]),c='r',edgecolors='None',s=50)
plot(blax,blay,'k-')
#title(month)
grid()
xlim(16,28)
ylim(16,28)  
xlabel('Computed F$_\t{lat}$ [W/m$^2$]')
ylabel('Model output F$_\t{lat}$ [W/m$^2$]')
#fig.savefig(outputpath+'compare_EoWm2_MPI-ESM-LR.pdf',bbox_inches='tight',orientation='landscape')
#####

#figure()
#plot(d['H_tot_mon'][k,:]-d['H_sn_mon'][k,:]-d['E_a_mon'][k,:])
#plot(d['E_o_mon'][k,:])




