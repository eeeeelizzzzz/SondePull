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
months = [8,9]
days = np.arange(1,32,1)
hours = np.arange(0,24,1) # checking hourly for any special sondes, could change if you don't care
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
                
                # vertical gradients
                # we care about value of gradient
                d_mr = np.abs(np.diff(mr))
                d_pt = np.abs(np.diff(pt))
                d_q = np.abs(np.diff(q))
                d_vpt = np.abs(np.diff(vpt))
                
                
                # range over which to search for the boundary layer top
                top_search = np.where(z.magnitude>3000)[0][0] # we don't care above here
                bot_search = np.where(z.magnitude>150)[0][0] # we don't care below here
                
                # do the searchin
                pt_bl = z[np.argmax(d_pt[bot_search:top_search])+bot_search]
                vpt_bl = z[np.argmax(d_vpt[bot_search:top_search])+bot_search]
                
                # confine moisture variable to within some range of vtemp variable ID'd level
                new_llim = np.where(z>vpt_bl-vpt_bl*.25)[0]
                new_ulim = np.where(z>vpt_bl+vpt_bl*.25)[0]
                if len(new_llim)>0: 
                    new_bot = max(bot_search, new_llim[0])
                else: #did not find one, beyond bounds
                    new_bot = bot_search
                if len(new_ulim)>0:
                    new_top = min(top_search, new_ulim[0])
                else: #did not find one, beyond bounds
                    new_top = top_search
            
                q_bl = z[np.argmax(d_q[new_bot:new_top])+new_bot]
                mr_bl = z[np.argmax(d_mr[new_bot:new_top])+new_bot]
                
                # put em all together
                list_bl = [mr_bl.magnitude, pt_bl.magnitude, q_bl.magnitude, vpt_bl.magnitude]
                bls = np.asarray(list_bl)
                bl_hgt = np.median(bls)*units.meter 
                
                # plot it up
                
                
                # profiles
                fig = plt.figure(figsize=(12,12))
                fig.suptitle(date.strftime("%d %b %Y %H%M")+' '+station)
                ax1 = fig.add_subplot(111)
                ax2 = ax1.twiny()
                ax1.plot(pt,z,'r',lw=3,label='Pot. Temp.')
                plt.axhline(pt_bl,color='r')
                ax1.plot(vpt,z,'g',lw=3,label='V Pot. Temp.')
                plt.axhline(vpt_bl,color='g')
                ax2.plot(mr,z,'b',lw=3,label='Mix. Ratio')
                plt.axhline(mr_bl,color='b')
                plt.axhline(bl_hgt,color='k',lw=5,ls='--')
                ax1.set_xlim(300,320)
                ax2.set_xlim(0,20)
                plt.ylim(0,3000)
                ax1.legend(loc='upper right')
                ax2.legend(loc='upper left')
                #plt.grid()
                plt.savefig('/Users/elizabeth.smith/Documents/PBLTops/sonde_plots/'+file_name+'.png')
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