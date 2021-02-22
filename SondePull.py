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
from metpy.plots import add_metpy_logo, SkewT
from metpy.units import units


from siphon.simplewebservice.wyoming import WyomingUpperAir
#from siphon.simplewebservice.igra2 import IGRAUpperAir
#station = 'USM00072357' #koun for IGRA

# big fonts are the best fonts
plt.rcParams['font.size'] = 20

# times/locs for PBLTops campaign 
# months = [8,9]
# days = np.arange(1,32,1)
# hours = np.arange(0,24,1) # checking hourly for any special sondes, could change if you don't care

# for testing
months = [8,9]
days = np.arange(1,3,1)
hours = np.arange(0,24,6) # checking hourly for any special sondes, could change if you don't care
stations = ['OUN','SHV']
for S in range(len(stations)):
    if S == 1: # these start points clip the dates we dont need for SHV
        s_d = 4
        s_m = 1
    else:
        s_d = 0
        s_m = 0
    for M in range(s_m,len(months)):
        if M == 1: # if sept(30 days)
            e_d = len(days)-1
        else:
            e_d = len(days)
        for D in range(s_d,e_d):
            for H in range(len(hours)):
                station = stations[S]
                month = months[M]
                day = days[D]
                hour = hours[H]
                # set the date you want to pull a sonde
                date = datetime(2020, month, day, hour)
                
                # this creates a string for printing and for file writeout
                file_name = station+'_2020_0'+str(month)+'_'+str(day)+'_'+str(hour)+'Z'
                
                # see if a data file exisits on the wyoming archive for your given date
                try:
                    df_orig = WyomingUpperAir.request_data(date, station)
                except:
                    # if not, tell us and moce on
                    print('No data located for '+file_name)
                    continue
                # get rid of fields we don't need (strings make issues in interpolation)
                print('Data located for '+file_name+'.......')
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
                
                #compute mixing ratio, pot. temp, spec. humidity, virtual pot. temp. 
                pt = mpcalc.potential_temperature(p,t) 
                q = mpcalc.specific_humidity_from_dewpoint(p,td)
                mr = mpcalc.mixing_ratio_from_specific_humidity(q)
                vpt = mpcalc.virtual_potential_temperature(p, t, mr)
                mr =  mr * 1000 * (units.gram/units.kilogram) #adjust units for viz
                
                # vertical gradients (actually differences, but since interpolated to constant spacing /dz would cancel anyway)
                # we care about value of gradient
                d_mr = np.diff(mr)
                d_pt = np.diff(pt)
                d_q = np.diff(q)
                d_vpt = np.diff(vpt)
                d_rh = np.diff(rh)
                d_t = np.diff(t)
                
                
                # range over which to search for the boundary layer top
                top_search = np.where(z.magnitude>3000)[0][0] # we don't care above here
                bot_search = np.where(z.magnitude>150)[0][0] # we don't care below here
                
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
                # where in the search range does vpt become greater than sfc value
                parcel_idx = np.where(vpt[bot_search:top_search] > vpt[0])[0]
                if parcel_idx.shape[0] < 1: #no points were found where vpt>sfc value
                    parcel_bl = np.nan # so put a nan in
                if parcel_idx[0]==1: # this means vpt is immediately increasing so lets ignore this
                    parcel_bl = np.nan # so put a nan in
                else:
                    parcel_bl = z[parcel_idx[0]+bot_search].magnitude #height of the first index where vpt is great than the sfc value
                    
                # M1_a
                # where in the search range does vpt become greater than sfc value +.6K
                parcel_a_idx = np.where(vpt[bot_search:top_search] > (vpt[0]+.6*units.K))[0]
                if parcel_a_idx.shape[0] < 1: #no points were found where vpt>sfc value +.6K
                    parcel_a_bl = np.nan # so put a nan in
                if parcel_a_idx[0]==1: # this means vpt is immediately increasing so lets ignore this
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
                tsign = np.sign(d_t.magnitude) #find sign of gradient vals
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
                        next_inv_idx = np.where(signchange[(top_inv_idx[0]+sfc_inv_idx[0]):]==1)[0] #keep looking up
                        if next_inv_idx.shape[0] < 1: # no points found, no inversions exist
                            elv_inv_bl = np.nan
                        else:
                            elv_inv_bl = z[next_inv_idx[0]+top_inv_idx[0]+sfc_inv_idx[0]].magnitude #bottom of the layer
                # now we want to do some screening for temps increasing from the sfc up -- problem spot.
                if sfc_inv_bl < 10: # if it is less than 10 m
                    print('hit')
                    sfc_inv_bl = np.nan # this 0 would pull down the stats.
                
                
                
                # put em all together
                list_bl = [parcel_bl, parcel_a_bl, pt_bl, q_bl, rh_bl, sfc_inv_bl, elv_inv_bl]
                bls = np.asarray(list_bl)
                #bl height is the median of all bl heights 
                bl_hgt = np.nanmedian(bls)*units.meter 
                
                # plot it up
                
                
                # profiles
                fig = plt.figure(figsize=(12,12))
                fig.suptitle(date.strftime("%d %b %Y %H%M")+' '+station)
                ax1 = fig.add_subplot(111)
                ax2 = ax1.twiny()
                ax1.plot(pt,z,'r',lw=3,label='Pot. Temp.')
                ax2.plot(rh*(100),z,'b',lw=3,label='RH')
                ax1.axhline(parcel_bl,color='c',label='M1')
                ax1.axhline(parcel_a_bl,color='b',label='M1_a')
                ax1.axhline(pt_bl,color='r',label='M2')
                ax1.axhline(q_bl,color='m',label='M3')
                ax2.axhline(rh_bl,color='g',label='M4')
                ax2.axhline(elv_inv_bl,color='darkred',label='M6')
                ax2.axhline(sfc_inv_bl,color='orange',label='M7')
                ax2.axhline(bl_hgt,color='k',lw=5,ls='--')
                ax1.set_xlim(280,320)
                ax2.set_xlim(0,100)
                plt.ylim(0,3000)
                ax1.legend(loc='upper right')
                ax2.legend(loc='upper left')
                #plt.grid()
                #plt.savefig('/Users/elizabeth.smith/Documents/PBLTops/sonde_plots/'+file_name+'.png')
                plt.show()
                #plt.close()
                
                
                
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