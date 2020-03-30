#Convert counts to SV 

from netCDF4 import Dataset
import numpy as np
from numpy import transpose
import pandas as pd
import datetime
from datetime import timedelta
import matplotlib.pyplot as plt
from cmocean import cm
from dtki import d_n_calc, matlab2datetime, datetime2matlabdn
from dateutil import parser
from geopy.distance import vincenty as distance #pip install vincenty
import os
import glob
import json
from scipy import stats
import scipy.io as sio
import pytz
import math
import statistics

#WORKHORSE 300kHz Sentinel
# import .mat file of CHP_02_03 data
# WH sentinel CHP_02 S/N 2128
# Kc; beam1: 0.3926 / beam2: 0.3933 / beam3: 0.4311 / beam4: 0.3922
# Er; beam1: 72 / beam2: 71 / beam3: 69 / beam4: 71

filelist = glob.glob('/Users/dong/Downloads/CHP_02_03/*.mat')

#1. calculate the avg Sali and Temp of each Timestamp
Sali_list = []
Temp_list = []
for mat in filelist:
    if 'CHP02_SBE' in mat:
        adcp =sio.loadmat(mat)
        time = adcp['ctd']['date'][0][0].tolist()
        depth = adcp['ctd']['dpth'][0][0].tolist()
        adcp_lat = adcp['ctd']['LatN'][0][0].tolist()
        adcp_lon = adcp['ctd']['LonE'][0][0].tolist()
        sali = adcp['ctd']['sali'][0][0].tolist()
        temp = adcp['ctd']['temp'][0][0].tolist()
        adcp_time = [(str(int(i[0]))+'-'+str(int(i[1]))+'-'+str(int(i[2]))+' '+str(int(i[3]))+'-'+str(int(i[4]))+'-'+str(int(i[5]))) for i in time]
        adcp_datetime = [datetime.datetime.strptime(i, '%Y-%m-%d %H-%M-%S') for i in adcp_time]

        Sali_list.append(sali)
        Temp_list.append(temp)
        continue

    else:
        continue
    
#caculate average of Sali, Temp of each timestamp from CTD data
Sali_avg = []
Temp_avg = []
iterator = range(0, len(Sali_list[0]))
for k in iterator:
    each_sali_element = [Sali_list[0][k], Sali_list[1][k], Sali_list[2][k], Sali_list[3][k], Sali_list[4][k]]
    each_temp_element = [Temp_list[0][k], Temp_list[1][k], Temp_list[2][k], Temp_list[3][k], Temp_list[4][k]]
    each_sali = np.average(each_sali_element)
    each_temp = np.average(each_temp_element)

    Sali_avg.append(each_sali)
    Temp_avg.append(each_temp)

#2. Calculate each Sv value
for mat in filelist:
    if 'CHP02_WH' in mat:
        adcp =sio.loadmat(mat)
        time = adcp['adcp']['date'][0][0].tolist()
        adcp_lat = adcp['adcp']['LatN'][0][0].tolist()
        adcp_lon = adcp['adcp']['LonE'][0][0].tolist()
        temp = adcp['adcp']['ET'][0][0].tolist()
        depth = adcp['adcp']['Zgrd'][0][0].tolist() #Bin depth, based on the depth of instrument. upward positive, downward negative
        instrument_Depth = adcp['adcp']['ED'][0][0].tolist()

        echo1 =  adcp['adcp']['CorMag1'][0][0].tolist()
        echo2 =  adcp['adcp']['CorMag2'][0][0].tolist()
        echo3 =  adcp['adcp']['CorMag3'][0][0].tolist()
        echo4 =  adcp['adcp']['CorMag4'][0][0].tolist()

        adcp_time = [(str(int(i[0]))+'-'+str(int(i[1]))+'-'+str(int(i[2]))+' '+str(int(i[3]))+'-'+str(int(i[4]))+'-'+str(int(i[5]))) for i in time]
        adcp_datetime = [datetime.datetime.strptime(i, '%Y-%m-%d %H-%M-%S') for i in adcp_time]

        adcp_depth = adcp['adcp']['Zgrd'][0][0][-10].tolist()

        echo_List = [echo1, echo2, echo3, echo4]
        Kc_list = [0.3926, 0.3933, 0.4311, 0.3922]
        Er_list = [72, 71, 69, 71]


        beam_iterator = range(0, len(echo_List)) #set iterator for each Beam (total 4 beams)
        for b in beam_iterator: 
            echo = echo_List[b]
            iterator= range(0, len(echo))
            for e in iterator:
                
                C = -143.5
                Tx = temp[e][0] #real time temp of trasducer(℃)
                tpl = 8.13 #transmit pulse length(m)
                Ldbm = 10*math.log10(tpl)
                #tp =  #transmit power(Watts)
                Pdbw = 14.0 #10*math.log10(tp) (dbWatt)
                Kc = Kc_list[b]#beam specific sensitivity coefficient (serial number needed/ ask RDI with S/N)
                Er = Er_list[b]#need to ask

                alpha = 0.068 #absorption coefficient of water(dB/m)

                iterator= range(0,len(echo[e]))
                for i in iterator:

                    E = echo[e][i] #echo intensity(counts)
                    B = 1.76 #blank after trasmit(m)
                    L = 8.13 #transmit pulse length(m)
                    D = 8 #depth cell length(m)
                    N = 24 #depth cell number of the scattering layer being measured
                    theta = 20 #beam angle from the system vertical
                    c1 = 1475.1 #speed of sound used by the instrument (preset) (m/s)
                    
                    
                    Sali = Sali_avg[e]
                    Temp = Temp_avg[e]

                    if (np.isnan(instrument_Depth[e][0]) == True) or (instrument_Depth[e][0] < 140):
                        echo[e][i] = 0
                    else:
                        Depth = instrument_Depth[e][0] - depth[e][i]
                        c_ =  1448.96 + 4.591*Temp - 5.304*(10**(-2))*Temp**2 + 2.374*(10**(-4))*Temp**3 + 1.340*(Sali-35) + 1.630*(10**(-2))*Depth + 1.675*(10**(-7))*Depth**2 - 1.025*(10**(-2))*Temp*(Sali-35) - 7.139*(10**(-13))*Temp*(Depth**3)
                            #sound speed at each depth cell for each ensemble (m/s), Eq by Mackenzie(1981), Adaptable range (temp )
                        

                        R = ((B+(L+D)/2 +((N-1)*D)+(D/4))/np.cos(theta))*(c_/c1) # range along the beam to the scatter(m)
                        SV = C + 10*np.log10((Tx + 273.16)*(R**2)) - Ldbm - Pdbw + 2*alpha*R + Kc*(E - Er)
                        echo[e][i] = SV

            echo_List[b] = echo
        
        echo = []
        iterator = range(0, len(echo_List[0]))
        for i in iterator:
            iterator = range(0, len(echo_List[0][i]))
            one_col_echo = []
            for e in iterator:
                all_echo = [int(echo_List[0][i][e]), int(echo_List[1][i][e]), int(echo_List[2][i][e]), int(echo_List[3][i][e])]
                if np.sum(all_echo) == 0:
                    one_col_echo.append(np.nan)
                else:
                    cell_echo = statistics.mean(all_echo)
                    one_col_echo.append(cell_echo)
            echo.append(one_col_echo)
        print(len(echo))
        
        print(echo)
        #plotting
        plt.figure(figsize= (40, 10))
        echo = np.transpose(echo)
        plt.pcolormesh(adcp_datetime[4000:4300], adcp_depth, echo[:, 4000:4300], cmap=cm.thermal)
        cbar = plt.colorbar(pad=0.005, aspect=20)
        cbar.set_label('SV', fontsize=25)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
        cbar.ax.tick_params(labelsize=25)

        os.chdir('/Users/dong/Downloads')
        plt.savefig('echo' + 'sample' + '.pdf')
        plt.close()
            