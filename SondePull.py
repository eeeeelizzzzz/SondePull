#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 13:10:56 2021

@author: elizabeth.smith

specifics in here are tailored to pull from wyo archive
and find a pbl top height. the dates/times/locs/loop are
written for a specific field campaign so those would
need modification for others
"""
import metpy 
from datetime import datetime, timezone
from siphon.simplewebservice.igra2 import IGRAUpperAir
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import metpy.calc as mpcalc
import scipy.stats as stat
import scipy as sp
from scipy import signal
from metpy.plots import add_metpy_logo, SkewT
from metpy.units import units
import netCDF4 as nc
import glob


from siphon.simplewebservice.wyoming import WyomingUpperAir
#from siphon.simplewebservice.igra2 import IGRAUpperAir
#station = 'USM00072357' #koun for IGRA

# big fonts are the best fonts
plt.rcParams['font.size'] = 20


# times/locs for PBLTops campaign  or CHEESEHEAD
#years = [2019]
years = [2020]
months = [8,9]
#months = [9,10]
days = np.arange(1,32,1)
hours = np.arange(0,24,3) # checking hourly for any special sondes, could change if you don't care
stations = ['OUN','SHV']
#stations = ['CHS']
grade = 0 # parameter defining sounding type, 0 = public quality NWS sonde, 1 = research quality NCAR sonde

def NCAR_sounding(file_name):
    f = nc.Dataset(file_name,'r')
    z = f.variables['alt'][:] #alt above MSL, m
    z = z-z[0]
    p = f.variables['pres'][:] # pressure, hPa
    t = f.variables['tdry'][:] + 273.15 # dry bulb T, K
    td = f.variables['dp'][:] + 273.15 # dew point T, K
    ws = f.variables['wspd'][:] # windspeed, m/s
    wd = f.variables['wdir'][:] # winddirections, deg
    u = f.variables['u_wind'][:] # u component, m/s
    v = f.variables['v_wind'][:] # v component, m/s
    rh = f.variables['rh'][:] # rel hum., %
    pt = f.variables['theta'][:] + 273.15 # pot temp, K
    mr = f.variables['mr'][:] # mixing ratio, #g/kg
    vpt = f.variables['theta_v'][:] + 273.15 #virt. pot. temp, K
    
    # Note: theta_v in the file IS NOT in potential temperature units!!!
    # #this is important! 
    # Fix needed:
    vpt = vpt * (1000. / p)**(0.286)

    
    return [z,p,t,td,pt,vpt,mr,rh,u,v,ws,wd]

def fill_finder(field, value=-999):
    field[np.where(field==value)[0]]=np.nan
    return field

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

for S in range(len(stations)):
    BLHGT = []
    BL25 = []
    BL75 = []
    SND_TIME = []
    
    Parcel_bl  = []
    Parcel_a_bl = []
    Pt_bl = []
    Q_bl = []
    Rh_bl = []
    Sfc_inv_bl  = []
    Elv_inv_bl = []
    
    
    
    if S == 1: # these start points clip the dates we dont need for SHV
        s_d = 4
        s_m = 1
    else:
        s_d = 0
        s_m = 0
    for M in range(s_m,len(months)):
        if M == 1: # if sept(30 days)
            e_d = 30
        else:
            e_d = len(days)
        for D in range(s_d,e_d-1):
            for H in range(len(hours)):
                station = stations[S]
                month = months[M]
                day = days[D]
                hour = hours[H]
                year = years[0]
                # set the date you want to pull a sonde
                date = datetime(year, month, day, hour).replace(tzinfo=timezone.utc)
                
                # this creates a string for printing and for file writeout
                file_name = station+'_'+str(year)+'_0'+str(month)+'_'+str(day)+'_'+str(hour)+'Z'
                
                # if using NWS public sonde data
                if grade == 0:
                    # see if a data file exisits on the wyoming archive for your given date
                    try:
                        df_orig = WyomingUpperAir.request_data(date, station)
                    except:
                        # if not, tell us and moce on
                        print('No data located for '+file_name)
                        continue
                    # get rid of fields we don't need (strings make issues in interpolation)
                    print('Data located for '+file_name+'.......')
                    df_orig1 = df_orig
                    df_orig = df_orig.drop(['station','station_number','time','latitude','longitude','elevation'],axis=1)
                    # interpolating z to get reg spacing for differencing and get rid of nans
                    # this may be problematic for some application, use caution
                    intz = np.arange(df_orig['height'].min(),4000,20) 
                    # in one we step pull data out of dataframe, interpolate to intz, put back in dataframe
                    df = pd.DataFrame(np.array([np.interp(intz,df_orig['height'].values,df_orig[d].values) for d in df_orig]).T)
                    # replace column headers on new datafram with the original ones
                    df.columns = df_orig.columns
                    
                    # get variables from sounding
                    z = df['height'].to_numpy() 
                    z = z-z[0] #adjusting for surface elevation 
                    z = z * units.meter
                    p = df['pressure'].to_numpy() 
                    p = p * units.hPa
                    t = df['temperature'].to_numpy()+273.15 
                    t = t * units.K
                    td = df['dewpoint'].to_numpy()+273.15
                    td = td * units.K
                    wd = df['direction'].to_numpy()
                    wd = wd * units.degrees
                    ws = df['speed'].to_numpy()
                    ws = ws * units.knot
                    u = df['u_wind'].to_numpy()
                    u = u * units.knot
                    v = df['v_wind'].to_numpy()
                    v = v * units.knot
                    rh = mpcalc.relative_humidity_from_dewpoint(t,td)
                    rh = rh * units.percent
                    rh = rh * 100.
                    
                    #compute mixing ratio, pot. temp, spec. humidity, virtual pot. temp. 
                    pt = mpcalc.potential_temperature(p,t) 
                    # this is written for metpy 1.0, older versions of metpy want the arguments reversed in the below q to dp function
                    # check your version with metpy.__version__
                    q = mpcalc.specific_humidity_from_dewpoint(p,td) 
                    mr = mpcalc.mixing_ratio_from_specific_humidity(q)
                    vpt = mpcalc.virtual_potential_temperature(p, t, mr)
                    mr =  mr * 1000 * (units.gram/units.kilogram) #adjust units for viz
                
                # if using NCAR research sondes
                if grade == 1:
                    file_path = '/Users/elizabeth.smith/Documents/CHEESEHEAD/sondes/'
                    snd_name = 'CHEESEHEAD_ISS1_RS41_v1_'+f"{date:%Y%m%d_%H}"
                    snd_file = glob.glob(file_path+snd_name+'*.nc')
                    if len(snd_file)<1: #glob found none
                        #print("Missing obs on "+f"{date:%Y%m%d_%H}")
                        continue
                    else:
                        snd_file = np.sort(glob.glob(file_path+snd_name+'*.nc'))[0]
                        print("Found obs on "+f"{date:%Y%m%d_%H}")
                        ### JTC 4/8/21  Adding full date out to seconds and 
                        # fixing timezone to be UTC
                        full_date = datetime(int(snd_file[-18:-14]),
                                             int(snd_file[-14:-12]),
                                             int(snd_file[-12:-10]),
                                             int(snd_file[-9:-7]),
                                             int(snd_file[-7:-5]),
                                             int(snd_file[-5:-3])).replace(tzinfo=timezone.utc)
                        snd_vars1 = NCAR_sounding(snd_file)# funtion returns [z0,p1,t2,td3,pt4,vpt5,mr6,rh7,u8,v9,ws10,wd11]
                        snd_vars=[]
                        for VAR in range(len(snd_vars1)):
                            Field = fill_finder(snd_vars1[VAR].data) #turning -999 into nans?
                            snd_vars.append(Field)
                    t = snd_vars[2] * units.K
                    z = snd_vars[0] * units.meter
                    p = snd_vars[1] * units.hPa
                    td = snd_vars[3] * units.K
                    pt = snd_vars[4] * units.K
                    vpt = snd_vars[5] * units.K
                    mr = snd_vars[6] * (units.gram/units.kilogram)
                    rh = snd_vars[7] * units.percent
                    u = snd_vars[8] * units.meter_per_second
                    v = snd_vars[9] * units.meter_per_second
                    ws = snd_vars[10] * units.meter_per_second
                    wd = snd_vars[11] * units.meter_per_second
                    q = mpcalc.specific_humidity_from_dewpoint(p,td)
                    
                    # Interpolate through NaNs prior to interpolation
    
                    nans, x = nan_helper(t)
                    t[nans] = np.interp(x(nans), x(~nans), t[~nans])
                    
                    nans, x = nan_helper(z)
                    z[nans] = np.interp(x(nans), x(~nans), z[~nans])
                    
                    nans, x = nan_helper(p)
                    p[nans] = np.interp(x(nans), x(~nans), p[~nans])
                    
                    nans, x = nan_helper(td)
                    td[nans] = np.interp(x(nans), x(~nans), td[~nans])
                    
                    nans, x = nan_helper(pt)
                    pt[nans] = np.interp(x(nans), x(~nans), pt[~nans])
                    
                    nans, x = nan_helper(vpt)
                    vpt[nans] = np.interp(x(nans), x(~nans), vpt[~nans])
                    
                    nans, x = nan_helper(mr)
                    mr[nans] = np.interp(x(nans), x(~nans), mr[~nans])
                    
                    nans, x = nan_helper(rh)
                    rh[nans] = np.interp(x(nans), x(~nans), rh[~nans])
                    
                    nans, x = nan_helper(u)
                    u[nans] = np.interp(x(nans), x(~nans), u[~nans])
                    
                    nans, x = nan_helper(v)
                    v[nans] = np.interp(x(nans), x(~nans), v[~nans])
                    
                    nans, x = nan_helper(ws)
                    ws[nans] = np.interp(x(nans), x(~nans), ws[~nans])
                    
                    nans, x = nan_helper(wd)
                    wd[nans] = np.interp(x(nans), x(~nans), wd[~nans])
                    
                    nans, x = nan_helper(q)
                    q[nans] = np.interp(x(nans), x(~nans), q[~nans])
                    
                    # Interpolate to 20-m grid
                    z_int =  np.arange(np.nanmin(z.magnitude), 4000, 20)
                    t_int =  np.interp(z_int, z.magnitude, t)
                    p_int = np.interp(z_int, z.magnitude, p)
                    td_int = np.interp(z_int, z.magnitude, td)
                    pt_int = np.interp(z_int, z.magnitude, pt)
                    vpt_int = np.interp(z_int, z.magnitude, vpt)
                    mr_int = np.interp(z_int, z.magnitude, mr)
                    rh_int = np.interp(z_int, z.magnitude, rh)
                    u_int = np.interp(z_int, z.magnitude, u)
                    v_int = np.interp(z_int, z.magnitude, v)
                    wd_int = np.interp(z_int, z.magnitude, wd)
                    ws_int = np.interp(z_int, z.magnitude, ws)
                    q_int = np.interp(z_int, z.magnitude, q)
                    
                    # Assign units for MetPy calculations and rewrite original variables
                    z = z_int * units.meter
                    t = t_int 
                    p = p_int 
                    td = t_int 
                    pt = pt_int 
                    vpt = vpt_int 
                    mr = mr_int 
                    u = u_int 
                    v = v_int 
                    wd = wd_int  
                    ws = ws_int 
                    rh = rh_int 
                    q = q_int 
                    
                    # Smooth data over 500-m window
                    t_smooth = signal.savgol_filter(t, 15, 5)
                    pt_smooth = signal.savgol_filter(pt, 15, 5)
                    q_smooth = signal.savgol_filter(q, 15, 5)
                    rh_smooth = signal.savgol_filter(rh, 15, 5)
                    vpt_smooth = signal.savgol_filter(vpt, 15, 5)
                    mr_smooth = signal.savgol_filter(mr, 15, 5)
                
                
                    # Compute differences over 200-m intervals because of minor noise in hi-res data
                    
                    d_mr = np.gradient(mr_smooth)
                    d_pt = np.gradient(pt_smooth)
                    d_t = np.gradient(t_smooth)
                    d_q = np.gradient(q_smooth)
                    d_vpt = np.gradient(vpt_smooth)
                    d_rh = np.gradient(rh_smooth)
                    
                    
                    
                
                # vertical gradients (actually differences, but since interpolated to constant spacing /dz would cancel anyway)
                # we care about value of gradient
                if grade == 0: #not using research grade sondes
                    d_mr = np.diff(mr)
                    d_pt = np.diff(pt)
                    d_q = np.diff(q)
                    d_vpt = np.diff(vpt)
                    d_rh = np.diff(rh)
                    d_t = np.diff(t)
                
                
                # range over which to search for the boundary layer top
                top_search = np.where(z.magnitude>2500)[0][0] # we don't care above here
                bot_search = np.where(z.magnitude>50)[0][0] # we don't care below here
                
                # we will be applying several methods (M) as described in Seidel et al (2010) JGR-Atmos.
                # M1 - parcel - height theta_v = theta_v(sfc)
                # M1_a - parcel_conig - height theta_v = theta_v(sfc) + .6K from Coniglio et al (2013) Wea. Forecasting
                # M2 - max pt gradient
                # M3 - min spec. hum gradient
                # M4 - min RH gradient
                # M5 - min refractivity gradient -- do not have this varible so skipping
                # M6 - base of elevated temp inversion
                # M7 - top of sfc based temp inversion
                
                # M1
                # where in the search range does vpt become greater than sfc value(mean from sfc to 100m)
                sfc_vpt = (vpt[0]+vpt[bot_search])/2. # mean of sfc value and the value at top of 100 m layer 
                # to account for low-level data anomalies, missing data points, pre-launch data, etc.
                parcel_idx = np.where(vpt[bot_search:top_search] > sfc_vpt)[0]
                if parcel_idx.shape[0] < 1: #no points were found where vpt>sfc value
                    parcel_bl = np.nan # so put a nan in
                elif parcel_idx[0]==1: # this means vpt is immediately increasing so lets ignore this
                    parcel_bl = np.nan # so put a nan in
                else:
                    parcel_bl = z[parcel_idx[0]+bot_search].magnitude #height of the first index where vpt is great than the sfc value
                    
                # M1_a
                # where in the search range does vpt become greater than sfc value +.6K
                parcel_a_idx = np.where(vpt[bot_search:top_search] > (sfc_vpt+.6*units.K))[0]
                if parcel_a_idx.shape[0] < 1: #no points were found where vpt>sfc value +.6K
                    parcel_a_bl = np.nan # so put a nan in
                elif parcel_a_idx[0]==1: # this means vpt is immediately increasing so lets ignore this
                    parcel_a_bl = np.nan # so put a nan in
                else:
                    parcel_a_bl = z[parcel_a_idx[0]+bot_search].magnitude #height of the first index where vpt is great than the sfc value +.6K
            
                
                # M2
                # max pt grad
                pt_bl = z[np.argmax(d_pt[bot_search:top_search])+bot_search].magnitude
                
                # M3
                # min spec hum gradient
                q_bl = z[np.argmin(d_q[bot_search:top_search])+bot_search].magnitude
                
                # M4
                # min RHgradient
                rh_bl = z[np.argmin(d_rh[bot_search:top_search])+bot_search].magnitude
                
                # for M6 and M7 we need temperature inversions, so lets get a profile of sign changes of t
                tsign = np.sign(d_t) #find sign of gradient vals
                signchange = ((np.roll(tsign, 1) - tsign) != 0).astype(int) #return 1 for sign change, 0 for not
                signchange[0]=0 #setting since numpy.roll does a circular shift, so if the last element has different 
                # sign than the first, the first element in the signchange array will be 1
                #note -1, +1, and 0 are all considered unique signs here
                
                # M6 and M7 - doing 7 and 6 together searching from bottom up
                #  top of sfc based temp inversion, bottom of elevated temp inversion
                sfc_inv_idx = np.where(signchange==1)[0] # first time the sign changes! 
                if sfc_inv_idx.shape[0] < 1: # no points found, no inversions exist
                    sfc_inv_bl = np.nan
                    elv_inv_bl = np.nan # by extension this is also a nan and we can move on
                else:
                    #we actually want the top of this inversion so we need to keep looking up
                    top_inv_idx = np.where(signchange[sfc_inv_idx[0]:]==1)[0]
                    if top_inv_idx.shape[0] < 1: # no points found, no inversions exist
                        sfc_inv_bl = z[sfc_inv_idx[0]].magnitude #inversion was 1 point deep
                        # there are no inversions aloft so 
                        elv_inv_bl = np.nan
                    elif z[top_inv_idx[0]+sfc_inv_idx[0]] - z[sfc_inv_idx[0]] > 300*units.meter: 
                    #if the next sign change is too far away to be the top of an inversion layer
                        sfc_inv_bl = z[sfc_inv_idx[0]].magnitude
                        elv_inv_bl = z[top_inv_idx[0]+sfc_inv_idx[0]].magnitude
                    else:
                        sfc_inv_bl = z[top_inv_idx[0]+sfc_inv_idx[0]].magnitude # top of inversion layer
                        next_inv_idx = np.where(signchange[(top_inv_idx[0]+sfc_inv_idx[0]+1):]==1)[0] #start at the next level and keep looking up
                        if next_inv_idx.shape[0] < 1: # no points found, no inversions exist
                            elv_inv_bl = np.nan
                        else:
                            elv_inv_bl = z[next_inv_idx[0]+top_inv_idx[0]+sfc_inv_idx[0]].magnitude #bottom of the layer
                
                
                
                # put em all together
                list_bl = [parcel_bl, parcel_a_bl, pt_bl, q_bl, rh_bl, sfc_inv_bl, elv_inv_bl]
                bls = np.asarray(list_bl)
                      
                #bl height is the median of all bl heights 
                bl_hgt = np.nanmedian(bls)*units.meter 
                
                #append to station lists for writeout later
                BLHGT.append(bl_hgt.magnitude) # bl height in meters
                BL25.append(np.percentile(bls,25))
                BL75.append(np.percentile(bls,75))
                if grade == 1: #research grade sonde have higher res time info
                    SND_TIME.append(int(full_date.timestamp())) #epoch seconds for the launch time
                if grade == 0:
                    SND_TIME.append(int(date.timestamp()))
                Parcel_bl.append(parcel_bl)
                Parcel_a_bl.append(parcel_a_bl) #.6K
                Pt_bl.append(pt_bl)
                Q_bl.append(q_bl)
                Rh_bl.append(rh_bl)
                Sfc_inv_bl.append(sfc_inv_bl)
                Elv_inv_bl.append(elv_inv_bl)
                
                
                
                # plot it up
                if grade==1:
                    file_name = 'CHS_'+snd_file[-18:-3]
                
                # profiles
                fig = plt.figure(figsize=(12,12))
                fig.suptitle(date.strftime("%d %b %Y %H%M")+' '+station)
                ax1 = fig.add_subplot(111)
                ax2 = ax1.twiny()
                ax1.plot(pt,z,'r',lw=3,label='Pot. Temp.')
                ax2.plot(rh,z,'b',lw=3,label='RH')
                ax1.axhline(parcel_bl,color='c',label='M1')
                ax1.axhline(parcel_a_bl,color='b',label='M1_a')
                ax1.axhline(pt_bl,color='r',label='M2')
                ax1.axhline(q_bl,color='m',label='M3')
                ax2.axhline(rh_bl,color='g',label='M4')
                ax2.axhline(elv_inv_bl,color='darkred',label='M6')
                ax2.axhline(sfc_inv_bl,color='orange',label='M7')
                ax2.axhline(bl_hgt,color='k',lw=5,ls='--')
                ax1.set_xlim(275,320)
                ax2.set_xlim(0,100)
                plt.ylim(0,3000)
                ax1.legend(loc='upper right')
                ax2.legend(loc='upper left')
                #plt.grid()
                plt.savefig('/Users/elizabeth.smith/Documents/PBLTops/sonde_plots/'+file_name+'.png')
                #plt.show()
                plt.close()
                
                
                
                ## gradients
                # fig = plt.figure(figsize=(12,12))
                # ax1 = fig.add_subplot(111)
                # ax2 = ax1.twiny()
                # ax1.plot(d_pt,z[1:],'r',lw=3,label='Pot. Temp.')
                # ax1.plot(d_vpt,z[1:],'g',lw=3,label='V Pot. Temp')
                # ax2.plot(d_mr,z[1:],'b',lw=3,label='Mix. Ratio')
                # plt.axhline(pt_bl,color='r')
                # plt.axhline(vpt_bl,color='g')
                # plt.axhline(mr_bl,color='b')
                # # plt.axhline(bl_hgt,color='k',lw=5,ls='--')
                # #ax1.set_xlim(-5,5)
                # #ax2.set_xlim(-5,5)
                # plt.ylim(0,3000)
                # ax1.legend(loc='lower right')
                # ax2.legend(loc='lower left')
                # #plt.grid()
                # plt.show()


                #skewT plotting
                # skew = SkewT()
                # # Plot the data using normal plotting functions, in this case using
                # # log scaling in Y, as dictated by the typical meteorological plot
                # skew.plot(p, t, 'r')
                # skew.plot(p, td, 'g')
                # skew.plot_barbs(p, u, v)
                
                # # Add the relevant special lines
                # skew.plot_dry_adiabats()
                # skew.plot_moist_adiabats()
                # skew.plot_mixing_lines()
                # skew.ax.set_ylim(1000, 100)
                
                # # Add the MetPy logo!
                # fig = plt.gcf()
                # add_metpy_logo(fig, 115, 100)
                # plt.show()
                print('End')
                
                
                
    # write out
    snd_time = np.asarray(SND_TIME)
    bl_height = np.asarray(BLHGT)
    bl_25 = np.asarray(BL25)
    bl_75 = np.asarray(BL75)
    parcelb = np.asarray(Parcel_bl)
    parcelb_a = np.asarray(Parcel_a_bl) #.6K
    ptb = np.asarray(Pt_bl)
    qb = np.asarray(Q_bl)
    rhb = np.asarray(Rh_bl)
    sfcb = np.asarray(Sfc_inv_bl)
    elvb = np.asarray(Elv_inv_bl)
    
    
    # create output file nc4.Dataset(name, write mode, clear if it exists, file format)
    output_file = nc.Dataset('/Users/elizabeth.smith/Documents/PBLTops/output_files/'+str(station)+'.nc', 'w', clobber=True, format='NETCDF3_64BIT')
    
    # global attributes
    output_file.title = 'PBL Heights estimated from '+str(station)
    output_file.author = 'Elizabeth Smith, NSSL'
    output_file.contact = 'Elizabeth.smith@noaa.gov'
    output_file.reference = 'https://github.com/eeeeelizzzzz/SondePull/'
    
    # define dimensions       (name,value)
    output_file.createDimension('t', len(snd_time)) #time dimension
    
    # create a variable file.createVariable(name, precision, dimensions) = values (usually some array)    
    output_file.createVariable('t','f8',('t'))[:] = snd_time
    setattr(output_file.variables['t'],'units','seconds epoch time (since 00UTC on 1/1/1970)')
    setattr(output_file.variables['t'],'description','valid times for sondes (+1hr since launch) epoch seconds')

    # create a variable file.createVariable(name, precision, dimensions) = values (usually some array)    
    output_file.createVariable('pbl_h','f8',('t'))[:] = bl_height
    setattr(output_file.variables['pbl_h'],'units','m AGL')
    setattr(output_file.variables['pbl_h'],'description','PBL height estimated from radiosonde profiles. See github reference for details.')   
    
    output_file.createVariable('pbl_25perc','f8',('t'))[:] = bl_25
    setattr(output_file.variables['pbl_25perc'],'units','m AGL')
    setattr(output_file.variables['pbl_25perc'],'description','25th percentile of all considered PBL height estimates')
    
    output_file.createVariable('pbl_75perc','f8',('t'))[:] = bl_75
    setattr(output_file.variables['pbl_75perc'],'units','m AGL')
    setattr(output_file.variables['pbl_75perc'],'description','75th percentile of all considered PBL height estimates')
    
    output_file.createVariable('pbl_h_parcel','f8',('t'))[:] = parcelb
    setattr(output_file.variables['pbl_h_parcel'],'units','m AGL')
    setattr(output_file.variables['pbl_h_parcel'],'description','PBL height estimated from the parcel method')
    
    output_file.createVariable('pbl_h_extparcel','f8',('t'))[:] = parcelb_a
    setattr(output_file.variables['pbl_h_extparcel'],'units','m AGL')
    setattr(output_file.variables['pbl_h_extparcel'],'description','PBL height estimated from the Coniglio extension of the parcel method')
    
    output_file.createVariable('pbl_h_pt','f8',('t'))[:] = ptb
    setattr(output_file.variables['pbl_h_pt'],'units','m AGL')
    setattr(output_file.variables['pbl_h_pt'],'description','PBL height estimated from max pt gradient')
    
    output_file.createVariable('pbl_h_qb','f8',('t'))[:] = qb
    setattr(output_file.variables['pbl_h_qb'],'units','m AGL')
    setattr(output_file.variables['pbl_h_qb'],'description','PBL height estimated from min q gradient')
    
    output_file.createVariable('pbl_h_rh','f8',('t'))[:] = rhb
    setattr(output_file.variables['pbl_h_rh'],'units','m AGL')
    setattr(output_file.variables['pbl_h_rh'],'description','PBL height estimated from min rh gradient')
    
    output_file.createVariable('pbl_h_sfcinv','f8',('t'))[:] = sfcb
    setattr(output_file.variables['pbl_h_pt'],'units','m AGL')
    setattr(output_file.variables['pbl_h_pt'],'description','PBL height estimated from top of sfc based inversion')
    
    output_file.createVariable('pbl_h_elvinv','f8',('t'))[:] = elvb
    setattr(output_file.variables['pbl_h_pt'],'units','m AGL')
    setattr(output_file.variables['pbl_h_pt'],'description','PBL height estimated from base of elevated inversion')
    output_file.close()
    print("File written")
    
    
    
                
                
                
                
                