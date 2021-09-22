#!/usr/bin/python3
#DATA T2, RHv1, WSP10, WDIR10, TSK
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from datetime import datetime, timedelta
import os
import psycopg2
import re

#local time translation in data files calculation
def to_loc_time_calc():
  cols = [0, 1]
  date_str = []
  with open(file_list[0]) as file:
    hh = np.loadtxt(file, dtype='str', skiprows=1, usecols=cols)
    file.close()
  for i in range(len(hh)):
    date_str.append(hh[i][0][1:] + ' ' + hh[i][1][:-1])
  date_const = datetime.strptime(date_str[0], "%Y-%m-%d %H:%M:%S")
  hh = []
  for i in range(len(date_str)):
    date_temp = datetime.strptime(date_str[i], "%Y-%m-%d %H:%M:%S")
    date_zero = datetime(date_const.year, date_const.month, date_const.day, 17, 0, 0)
    hh.append(int((date_temp - date_zero).total_seconds() / 3600))
  skip_rows = hh.index(0)+1
  max_rows = hh.index(24) - hh.index(0)+1
  #max_rows = hh.index(6) - hh.index(0)+1
  hh_calc = []
  for i in range(len(hh)):
    if hh[i] >= 0 and hh[i] < 25:
      hh_calc.append(hh[i])
  return hh_calc, skip_rows, max_rows

#local time translation in data  files observation
def to_loc_time_obs(date_list):
    date_const = date_list[0]
    hh = []
    for i in range(len(date_list)):
        date_temp = date_list[i]
        date_zero = datetime(date_const.year, date_const.month, date_const.day, 0, 0, 0)
        hh.append(float((date_temp - date_zero).total_seconds() / 3600))
    hh_obs = []
    for i in range(len(hh)):
        if hh[i] >= 0 and hh[i] <= 24:
           hh_obs.append(hh[i])
    return hh_obs

#selection of observation data from the database
def export_data_obs_from_db(t1, t2, datetime, table_name, temp=None, rh=None, wsp=None, wdir=None):
    # list name column in data base
    par_list_input = [temp, rh, wsp, wdir]
    # dict name columns that is not none
    dict_par_output = {}
    flag = 1
    for par in par_list_input:
        if par is not None:
            t1 = str(t1)
            t2 = str(t2)
            conn = psycopg2.connect(dbname='test', user='test', password='abcdef',host='127.0.0.1')
            cur = conn.cursor()
            #one time date conversion
            if flag:
                cur.execute("SELECT " + datetime + " FROM " + table_name + " WHERE " + datetime + " >= '" + t1 + "' and " + datetime + " <= '" + t2 + "';")
                hh = []
                for row in cur:
                    hh.append(row[0])
                hh_loc_obs = to_loc_time_obs(hh)
                dict_par_output.update({'hh': hh_loc_obs})
            flag = 0
            par_obs = []
            cur.execute("SELECT " + par + " FROM " + table_name + " WHERE " + datetime + " >= '" + t1 + "' and " + datetime + " <= '" + t2 + "';")
            for row in cur:
                par_obs.append(row[0])
            dict_par_output.update({par: par_obs})
            cur.close()
            conn.close()
    return dict_par_output


#minimum and maximum temperature per month for graphs scale bounds - values taken from the internet
def temp_bounds(month_number):
  bounds = {1: [-30,-10],2: [-30,-5],3: [-20,0],4: [-5,10],5: [0,20],6: [7,30],7: [10,30],8: [8,25],9: [4,20], 10: [-5,5],11: [-15,-5],12: [-30,-10]}
  return min(bounds.get(month_number)), max(bounds.get(month_number))

#reading meteorological parameters
def read_values(path_to_file, cols):
  with open(path_to_file) as file:
    values = np.loadtxt(file, skiprows=skip_rows, usecols=cols, max_rows=max_rows)
    file.close()
  return values

#calculation standard deviation
def sigma(hh_calc,hh_obs,par_calc,par_obs):
    hh_par_calc_dict = dict(zip(hh_calc, par_calc))
    hh_par_obs_dict = dict(zip(hh_obs, par_obs))
    hh_par_obs_dict = {k: v for k, v in hh_par_obs_dict.items() if (k).is_integer()}
    list_obs = []
    list_calc = []
    for k, v in hh_par_obs_dict.items():
        if k in hh_par_calc_dict.keys():
            list_obs.append(hh_par_obs_dict[k])
            list_calc.append(hh_par_calc_dict[k])
    standard_deviation = 0
    for obs, calc in zip(list_obs, list_calc):
        standard_deviation += (obs - calc) ** 2
    standard_deviation = round(np.sqrt(standard_deviation / len(list_calc)),5)
    return standard_deviation

#a graphs of calculations and observations are drawn for an individual location
def draw_graph_single(date,graph_csl_key,hh_calc, par_calc,ylabel,y_min,y_max,y_step,par_name, hh_obs=None, par_obs=None,standard_deviation=None):
    fig, ax = plt.subplots(figsize=(15, 10))
    if standard_deviation is not None:
        fig.suptitle(date.strftime("%d-%m-%Y")+'  Standard deviation = '+str(standard_deviation), fontsize=20)
    else:
        fig.suptitle(date.strftime("%d-%m-%Y"), fontsize=20)
    ax.plot(hh_calc, par_calc, color=graph_csl[graph_csl_key][0], label=graph_csl[graph_csl_key][1], linewidth=0.7)
    if par_obs is not None:
      ax.scatter(hh_obs, par_obs, color=graph_csl[graph_csl_key][2], marker=graph_csl[graph_csl_key][3], label=graph_csl[graph_csl_key][4])
    ax.set_xlabel('Local Time [h]', fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which='major', length=10)
    ax.tick_params(which='minor', length=4)
    ax.set_xlim(min(hh_calc), max(hh_calc))
    ax.set_xticks(np.arange(min(hh_calc), max(hh_calc)+1, step=4))
    #definition y scale bounds
    if y_min > min(par_calc):
       y_min = min(par_calc)-1
    if y_max < max(par_calc):
       y_max = max(par_calc)
    try:
        if y_min > min(par_obs):
           y_min = min(par_obs)-1
    except:
           pass
           #No observations were found to compare or no minimum
    try:
        if y_max < max(par_obs):
           y_max = max(par_obs)
    except:
          pass
          #No observations were found to compare or no maximum
    ax.set_ylim(int(y_min), int(y_max))
    ax.set_yticks(np.arange(int(y_min), int(y_max) + 1, y_step))
    fig.legend(loc='center right')
    fig.savefig('.../plots_2018/ts/single/' + par_name + '_' + file_name + '_single.png')
    #plt.show()
    return 0

#a graphs of calculations and observations are drawn for an all locations
def draw_graph_all(date, graph_csl,data_all,y_min,y_max,y_step,y_label,par_name,key_hh_calc,key_hh_obs,key_par_calc,key_par_obs):
    fig, ax = plt.subplots(figsize=(15, 10))
    fig.suptitle(date.strftime("%d-%m-%Y"), fontsize=20)
    locations = ['airport','bec','center','gmc','tor']
    for data_loc, loc in zip(data_all,locations):
        ax.plot(data_loc[key_hh_calc], data_loc[key_par_calc], color=graph_csl[loc][0], label=graph_csl[loc][1], linewidth=0.7)
        if data_loc[key_par_obs] is not None:
           ax.scatter(data_loc[key_hh_obs], data_loc[key_par_obs], color=graph_csl[loc][2], marker=graph_csl[loc][3], label=graph_csl[loc][4])
    ax.set_xlabel('Local Time [h]', fontsize=14)
    ax.set_ylabel(y_label, fontsize=14)
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which='major', length=10)
    ax.tick_params(which='minor', length=4)
    ax.set_xlim(min(hh_calc), max(hh_calc))
    ax.set_xticks(np.arange(min(hh_calc), max(hh_calc)+1, step=4))
    #definition y scale bounds
    for data_loc in data_all:
        if y_min > min(data_loc[key_par_calc]):
           y_min = min(data_loc[key_par_calc])-1
        if y_max < max(data_loc[key_par_calc]):
           y_max = max(data_loc[key_par_calc])
    try:
        if y_min > min(data_loc[key_par_obs]):
           y_min = min(data_loc[key_par_obs])-1
    except:
           pass
           #No observations were found to compare or no minimum
    try:
        if y_max < max(data_loc[key_par_obs]):
           y_max = max(data_loc[key_par_obs])
    except:
          pass
          #No observations were found to compare or no maximum

    ax.set_ylim(int(y_min), int(y_max))
    ax.set_yticks(np.arange(int(y_min), int(y_max) + 1, y_step))
    fig.legend(loc='center right')
    fig.savefig('.../plots_2018/ts/all/' + par_name + '_' +date.strftime("%d-%m-%Y")+'_all.png')
    #plt.show()
    return 0

#read the data and draw your own graph for each area
def read_data_and_draw_graph(file,column_rh, date,location):
    wsp10_calc = read_values(file, 4)
    wdir10_calc = read_values(file, 5)
    tmp2_calc = read_values(file, 6)
    rhv1_calc = read_values(file, column_rh)
    tsk_calc = read_values(file, 10)
    temp_min, temp_max = temp_bounds(date.month)
    standard_deviation_tmp = None
    standard_deviation_rhv1 = None
    standard_deviation_wsp10 = None
    standard_deviation_wdir10 = None
    standard_deviation_tsk = None
    # check for files with observations
    if location == 'airport' or location == 'gmc':
       if location == 'airport':
          data_obs = export_data_obs_from_db(date, date + timedelta(days=1), datetime='"Datetime"',table_name='"rp5_Bogashevo"', temp='"T"', rh='"U"')
       else:
           data_obs = export_data_obs_from_db(date, date + timedelta(days=1), datetime='"Datetime"', table_name='"rp5_Tomsk"', temp='"T"', rh='"U"')
       hh_obs = data_obs['hh']
       tmp_obs = data_obs['"T"']
       rhv1_obs = data_obs['"U"']
       standard_deviation_tmp = sigma(hh_calc,hh_obs,tmp2_calc,tmp_obs)
       standard_deviation_rhv1 = sigma(hh_calc,hh_obs,rhv1_calc,rhv1_obs)
       wsp10_obs = None
       wdir10_obs = None
       tsk_obs = None
    if location == 'bec' or location == 'tor':
       try:
            if location == 'bec':
                data_obs = export_data_obs_from_db(date, date + timedelta(days=1), datetime='"Datetime"', table_name='"iao_Wind_UZM_dd_mm_year"', wsp='"Vh_UZM_BEC_10"',wdir='"Fi_UZM_BEC_10"')
                hh_obs = data_obs['hh']
                wsp10_obs = data_obs['"Vh_UZM_BEC_10"']
                wdir10_obs = data_obs['"Fi_UZM_BEC_10"']
            else:
                data_obs = export_data_obs_from_db(date, date + timedelta(days=1), datetime='"Datetime"', table_name='"iao_Wind_UZM_dd_mm_year"', wsp='"Vh_UZM_Academ"', wdir='"Fi_UZM_Academ"')
                hh_obs = data_obs['hh']
                wsp10_obs = data_obs['"Vh_UZM_Academ"']
                wdir10_obs = data_obs['"Fi_UZM_Academ"']
            standard_deviation_wsp10 = sigma(hh_calc, hh_obs, wsp10_calc, wsp10_obs)
            standard_deviation_wdir10 = sigma(hh_calc, hh_obs, wdir10_calc, wdir10_obs)
            if location == 'bec':
                data_obs = export_data_obs_from_db(date, date + timedelta(days=1), datetime='"Datetime"', table_name='"iao_temperature_UZM_dd_mm_year"', temp='"T_UZM_BEC_10"')
                tmp_obs = data_obs['"T_UZM_BEC_10"']
            else:
                data_obs = export_data_obs_from_db(date, date + timedelta(days=1), datetime='"Datetime"', table_name='"iao_temperature_UZM_dd_mm_year"', temp='"T_UZM_Academ"')
                tmp_obs = data_obs['"T_UZM_Academ"']
            standard_deviation_tmp = sigma(hh_calc, hh_obs, tmp2_calc, tmp_obs)
        except:
            hh_obs = None
            tmp_obs = None
            rhv1_obs = None
            wsp10_obs = None
            wdir10_obs = None
            tsk_obs = None
            #No observations were found to compare
    if location == 'center':
       hh_obs = None
       tmp_obs = None
       rhv1_obs = None
       wsp10_obs = None
       wdir10_obs = None
       tsk_obs = None

    draw_graph_single(date, location, hh_calc, tmp2_calc, ylabel[0], temp_min, temp_max, 1, par_name[0],hh_obs=hh_obs, par_obs=tmp_obs,standard_deviation=standard_deviation_tmp)
    draw_graph_single(date, location, hh_calc, rhv1_calc, ylabel[1], 0, 101, 20, par_name[1],  hh_obs= hh_obs, par_obs=rhv1_obs,standard_deviation=standard_deviation_rhv1)
    draw_graph_single(date, location, hh_calc, wsp10_calc, ylabel[2], 0, 10, 1, par_name[2],hh_obs= hh_obs, par_obs=wsp10_obs,standard_deviation=standard_deviation_wsp10)
    draw_graph_single(date, location, hh_calc, wdir10_calc, ylabel[3], 0, 360, 60, par_name[3],hh_obs= hh_obs, par_obs=wdir10_obs,standard_deviation=standard_deviation_wdir10)
    draw_graph_single(date, location, hh_calc, tsk_calc, ylabel[4], int(min(tsk_calc))-1,int(max(tsk_calc))+1 , 1, par_name[4],hh_obs= None, par_obs=None,standard_deviation=None)
    #save read data to the dictionary for future use
    data_calc_obs = {'hh_calc': hh_calc, 'hh_obs': hh_obs, 'tmp2_calc': tmp2_calc, 'rhv1_calc': rhv1_calc,
                     'wsp10_calc': wsp10_calc, 'wdir10_calc': wdir10_calc, 'tsk_calc': tsk_calc, 'tmp_obs': tmp_obs, 'rhv1_obs': rhv1_obs, 'wsp10_obs': wsp10_obs, 'wdir10_obs': wdir10_obs, 'tsk_obs':tsk_obs}

    return data_calc_obs

data_path = '.../Data/4plots/wrfout_d03_2018-01-02/'
file_list = glob.glob(data_path+'*ts.dat')

#local time and range limits with positive values in data files calculation
hh_calc, skip_rows, max_rows = to_loc_time_calc()

#dictionary with parameters for displaying the values of calculations and observations on graphs
#color, symbol,label
graph_csl = {'airport':['blue','airport_calc','red','^','airport_obs'],'bec':['red','bec_calc','black','o','bec_obs'],
             'center':['orange','center_calc','white','*','center_obs'],'gmc':['green','gmc_calc','brown','D','gmc_obs'], 'tor':['black','tor_calc','green','o','tor_obs']}

ylabel = ['Temperature $[^\circ$C]','Relative humidity v1 [%]','Wind Speed '+r'$[\frac{m}{s}]$','Wind Direction $[^\circ$]','Temperature skin $[^\circ$C]']

#for output file name
par_name = ['t2','rhv1','wsp10','wdir10','tsk']


#date data definition
date = re.findall(r'\d{4}-\d{2}-\d{2}', file_list[0])
date = datetime.strptime(date[0], "%Y-%m-%d")+timedelta(days=1)

#draw graphs individually
for file in file_list:

    file_name = file.split('/')[-1]
    if re.findall('airport', file):
       data_airport = read_data_and_draw_graph(file,21,date,'airport')
    if re.findall('bec', file):
       data_bec = read_data_and_draw_graph(file,20,date,'bec')
    if re.findall('center', file):
       data_center = read_data_and_draw_graph(file, 21, date, 'center')
    if re.findall('gmc', file):
       data_gmc = read_data_and_draw_graph(file, 21, date, 'gmc')
    if re.findall('tor', file):
       data_tor = read_data_and_draw_graph(file, 20, date, 'tor')

#draw all the graphs at once
data_all = [data_airport,data_bec,data_center,data_gmc,data_tor]
temp_min, temp_max = temp_bounds(date.month)
draw_graph_all(date, graph_csl,data_all,temp_min,temp_max,1,ylabel[0],par_name[0],'hh_calc','hh_obs','tmp2_calc','tmp_obs')
draw_graph_all(date, graph_csl,data_all,0,101,20,ylabel[1],par_name[1],'hh_calc','hh_obs','rhv1_calc','rhv1_obs')
draw_graph_all(date, graph_csl,data_all,0,10,1,ylabel[2],par_name[2],'hh_calc','hh_obs','wsp10_calc','wsp10_obs')
draw_graph_all(date, graph_csl,data_all,0,360,60,ylabel[3],par_name[3],'hh_calc','hh_obs','wdir10_calc','wdir10_obs')
draw_graph_all(date, graph_csl,data_all,temp_min,temp_max, 1, ylabel[4],par_name[4],'hh_calc', 'hh_obs', 'tsk_calc', 'tsk_obs')



