#!/usr/bin/python3
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from datetime import datetime, timedelta
import os
import re
import psycopg2

#minimum and maximum temperature per month for graphs scale bounds
def temp_bounds(month_number):
  bounds = {1: [-30,-10],2: [-30,-5],3: [-20,0],4: [-5,10],5: [0,20],6: [7,30],7: [10,30],8: [8,25],9: [4,20], 10: [-5,5],11: [-15,-5],12: [-30,-10]}
  return min(bounds.get(month_number)), max(bounds.get(month_number))

#reading meteorological parameters
def read_values(path_to_file, cols):
  with open(path_to_file) as file:
    values = np.loadtxt(file, skiprows=1, usecols=cols, max_rows=6)
    file.close()
  return values

def export_data_obs_from_db(t, table_name):
    # dict name columns that is not none
    #dict_par_output = {}
    t = repr(str(t))
    conn = psycopg2.connect(dbname='meteo_db', user='postgres', password='postgres', host='localhost')
    cur = conn.cursor()
    if re.findall('temperature_mtr5', table_name):
       temp_obs = []
       cur.execute('SELECT "T_0","T_50","T_100","T_150","T_200","T_250","T_300","T_350","T_400","T_450","T_500" FROM '+table_name+' WHERE "Datetime" = %s;' % t)
       for row in cur:
           for r in row:
               temp_obs.append(r)
       print("temp_obs")
       print(temp_obs)
       cur.close()
       conn.close()
       return temp_obs
    if re.findall('sodar_wind_speed', table_name):
        wsp10_obs = []
        cur.execute('SELECT "Vh_75","Vh_100","Vh_125","Vh_150","Vh_175","Vh_200","Vh_225","Vh_250","Vh_275","Vh_300","Vh_325", "Vh_350","Vh_375","Vh_400","Vh_425","Vh_450","Vh_475","Vh_500" FROM ' + table_name + ' WHERE "Datetime" = %s;' % t)
        for row in cur:
            for r in row:
                wsp10_obs.append(r)
        print("wsp10_obs")
        print(wsp10_obs)
        cur.close()
        conn.close()
        return wsp10_obs
    if re.findall('sodar_wind_dir', table_name):
        wdir10_obs = []
        cur.execute('SELECT "Fi_75","Fi_100","Fi_125","Fi_150","Fi_175","Fi_200","Fi_225","Fi_250","Fi_275","Fi_300","Fi_325", "Fi_350","Fi_375","Fi_400","Fi_425","Fi_450","Fi_475","Fi_500" FROM ' + table_name + ' WHERE "Datetime" = %s;' % t)
        for row in cur:
            for r in row:
                wdir10_obs.append(r)
        print("wdir10_obs")
        print(wdir10_obs)
        cur.close()
        conn.close()
        return wdir10_obs


#nearest value to current for hgt
def nearest(lst, target):
  return min(lst, key=lambda x: abs(target-x))

def sigma(hgt_calc,hgt_obs,par_calc,par_obs):
    hh_par_calc_dict = dict(zip(hgt_calc, par_calc))
    hh_par_obs_dict = dict(zip(hgt_obs, par_obs))
    #print(hh_par_calc_dict)
    #print(hh_par_obs_dict)
    list_obs = []
    list_calc = []
    list_k_temp = []
    for k in hh_par_calc_dict.keys():
        k_temp = nearest(hh_par_obs_dict.keys(), k)#+4 so that the nearest value is 50
        #print(k,' ',k_temp)
        if k_temp not in list_k_temp:
            list_k_temp.append(k_temp)
        #list_calc.append(hh_par_calc_dict[k])
        #list_obs.append(hh_par_obs_dict[k_temp])
    #print(list_k_temp)
    for k in list_k_temp:
        k_temp = nearest(hh_par_calc_dict.keys(), k)#+4 so that the nearest value is 50
        #print(k,' ',k_temp)
        list_calc.append(hh_par_calc_dict[k_temp])
        list_obs.append(hh_par_obs_dict[k])
    #print("Observations")
    #print(list_obs)
    #print("Calculations")
    #print(list_calc)
    num_none = 0;
    standard_deviation = 0
    for obs, calc in zip(list_obs, list_calc):
        if (obs != None) and (calc != None):
            standard_deviation += (obs - calc) ** 2
        else:
            num_none+=1
    if num_none == 0:
        standard_deviation = round(np.sqrt(standard_deviation / len(list_calc)),5)
    else:
        standard_deviation = round(np.sqrt(standard_deviation / (len(list_calc)-num_none)), 5)
    #print('std=',standard_deviation)
    return standard_deviation


#a graphs of calculations and observations are drawn for an individual location
def draw_graph_single(date, time, graph_csl_key, x_label, par_calc, hgt_calc,par_obs, hgt_obs, x_min, x_max, y_min, y_max, y_step, obs_found, par_name, standard_deviation=None):
  fig, ax = plt.subplots(figsize=(15, 10))
  if standard_deviation is not None:
    fig.suptitle(date.strftime("%d-%m-%Y") +' '+time+ '  Standard deviation = ' + str(standard_deviation), fontsize=20)
  else:
    fig.suptitle(date.strftime("%d-%m-%Y") +' '+ time, fontsize=20)
  ax.plot(par_calc, hgt_calc, color=graph_csl[graph_csl_key][0], label=graph_csl[graph_csl_key][1], linewidth=0.7)
  #obs exist
  if obs_found and par_obs != []:
    ax.scatter(par_obs, hgt_obs, color=graph_csl[graph_csl_key][2], marker=graph_csl[graph_csl_key][3],
               label=graph_csl[graph_csl_key][4])
  ax.set_xlabel(x_label, fontsize=15)
  ax.set_ylabel('Height above surface [m]', fontsize=15)
  ax.xaxis.set_minor_locator(AutoMinorLocator(3))
  ax.yaxis.set_minor_locator(AutoMinorLocator(2))
  ax.tick_params(which='major', length=10)
  ax.tick_params(which='minor', length=4)
  # definition x scale bounds
  if x_min > min(par_calc):
      x_min = min(par_calc)-1
  if x_max < max(par_calc):
      x_max = max(par_calc)
  try:
    if x_min > min(par_obs):
       x_min = min(par_obs) - 1
  except:
    pass
    # No observations were found to compare or no minimum
  try:
    if x_max < max(par_obs):
       x_max = max(par_obs)
  except:
    pass
    # No observations were found to compare or no maximum
  ax.set_xlim(int(x_min), int(x_max))
  if re.findall('Temperature', x_label):
    ax.set_xticks(np.arange(int(x_min), int(x_max) + 1, step=3))
  elif re.findall('Wind Speed',x_label):
      ax.set_xticks(np.arange(int(x_min), int(x_max) + 1, step=1))
  else:
      ax.set_xticks(np.arange(int(x_min), int(x_max) + 1, step=60))
  ax.set_ylim(y_min, y_max)
  ax.set_yticks(np.arange(y_min, y_max + 1, y_step))
  fig.legend(loc='center right')
  fig.savefig('.../plots_2018/point/single/'+par_name+'_'+file_name+'_single.png')
  #plt.show()
  return 0

def draw_graph_all(date,time, graph_csl,data_all, x_label,x_min,x_max, y_min,y_max, y_step,key_hgt_calc,key_hgt_obs,key_par_calc,key_par_obs,par_name):
    fig, ax = plt.subplots(figsize=(15, 10))
    fig.suptitle(date.strftime("%d-%m-%Y")+' '+ time, fontsize=20)
    locations = ['airport','bec','center','gmc','tor']
    for data_loc, loc in zip(data_all,locations):
        ax.plot(data_loc[key_par_calc], data_loc[key_hgt_calc], color=graph_csl[loc][0], label=graph_csl[loc][1], linewidth=0.7)
        if data_loc[key_par_obs] is not None and data_loc[key_par_obs] != []:
           ax.scatter(data_loc[key_par_obs], data_loc[key_hgt_obs], color=graph_csl[loc][2], marker=graph_csl[loc][3], label=graph_csl[loc][4])
    ax.set_xlabel(x_label, fontsize=15)
    ax.set_ylabel('Height above surface [m]', fontsize=15)
    ax.xaxis.set_minor_locator(AutoMinorLocator(3))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which='major', length=10)
    ax.tick_params(which='minor', length=4)
    #definition y scale bounds
    for data_loc in data_all:
        if x_min > min(data_loc[key_par_calc]):
           x_min = min(data_loc[key_par_calc])-1
        if x_max < max(data_loc[key_par_calc]):
           x_max = max(data_loc[key_par_calc])
        try:
            if x_min > min(data_loc[key_par_obs]):
               x_min = min(data_loc[key_par_obs])-1
        except:
               pass
               #No observations were found to compare or no minimum
        try:
            if x_max < max(data_loc[key_par_obs]):
               x_max = max(data_loc[key_par_obs])
        except:
              pass
              #No observations were found to compare or no maximum
    ax.set_xlim(int(x_min), int(x_max))
    if re.findall('Temperature', x_label):
        ax.set_xticks(np.arange(int(x_min), int(x_max) + 1, step=3))
    elif re.findall('Wind Speed', x_label):
        ax.set_xticks(np.arange(int(x_min), int(x_max) + 1, step=1))
    else:
        ax.set_xticks(np.arange(int(x_min), int(x_max) + 1, step=60))
    ax.set_ylim(y_min, y_max)
    ax.set_yticks(np.arange(y_min, y_max + 1, y_step))
    fig.legend(loc='center right')
    fig.savefig('.../plots_2018/point/all/'+par_name+'_'+file_name+'_all.png')
    #plt.show()
    return 0

#observation tor mtp5
def define_temp_obs(f,skip_rows_temp_obs):
  with open(data_path + f) as fin:
    values = np.loadtxt(fin, dtype='str', skiprows=skip_rows_temp_obs, max_rows=1)[4:]
  fin.close()
  val = np.zeros(values.shape)
  fin.close()
  for i in range(values.shape[0]):
    val[i] = float(values[i])
  return val

#i = 0

def define_area_and_draw_graph(file,graph_csl_key):
  '''
  #if data in files
  if not re.findall('tor', file):
    tmp_obs = None
    hgt_obs = None
  else:
    try:
        data_obs = export_data_obs_from_db(date, date + timedelta(days=1), datetime='"Datetime"',
                                           table_name='"iao_temperature_UZM_dd_mm_year"', temp='"T_UZM_BAK_10"')
        tmp_obs = data_obs['"T_UZM_BAK_10"']
    for f in os.listdir(data_path):
      if f.startswith('mtp5'):
        obs_found = True
        f_obs_found = f
  if not re.findall('gmc', file) and (re.split('_',file))[-1][:-4] not in gmc_list:
    tmp_calc = read_values(file, 5)
    hgt_calc = read_values(file, 11)
  else:
      with open(file) as f:
          values = np.loadtxt(f, skiprows=3, usecols=[5, 11], max_rows=6)
      f.close()
      tmp_calc = values[:,0]
      hgt_calc = values[:,1]
  ######################
  #conversion to local time
  date = (file[97:120]).split('_')[0]
  time = ''
  for i in range(len((file[97:120]).split('_')[1].split('%3A'))):
      time = time + (file[97:120]).split('_')[1].split('%3A')[i] + ':'
  date = date + ' ' + time[:-1]
  date_temp = datetime.strptime(date, "%Y-%m-%d %H:%M:%S")
  loctime = int((date_temp - date_const).total_seconds() / 3600)
  #if ((loctime >= 0) and (date_temp.day < date_const.day + 2)):
  if (loctime >= 0):
      date_temp = date_const + timedelta(days=1)
  #########################################
  date = re.findall(r'\d{4}-\d{2}-\d{2}', file)
  time = re.findall(r'\d{2}%3A\d{2}%3A\d{2}', file)
  time = time[-1].replace('%3A', ':')
  date = datetime.strptime(date[-1], "%Y-%m-%d")
  temp_min, temp_max = temp_bounds(date.month)
  if obs_found:
      tmp_obs = define_temp_obs(f_obs_found, 0)
      hgt_obs = np.arange(0, 501, 50)
  else:
      tmp_obs = None
      hgt_obs = np.arange(0, 501, 50)'''
  #if data in data base
  #obs_found = False
  standard_deviation_temp = None
  standard_deviation_wsp10 = None
  standard_deviation_wdir10 = None
  temp_min, temp_max = temp_bounds(date.month)
  tmp_obs = None
  wsp10_obs = None
  wdir10_obs = None
  hgt_obs_temp = None
  hgt_obs_wsp10 = None
  hgt_obs_wdir10 = None
  obs_found_temp = False
  obs_found_wsp10 = False
  obs_found_wdir10 = False
  if re.findall('tor', file):
      #obs_found = True
      hgt_calc = read_values(file, 11)
      tmp_calc = read_values(file, 5)
      wsp10_calc = read_values(file, 3)
      wdir10_calc = read_values(file, 4)
      try:
          tmp_obs = export_data_obs_from_db(date, table_name='"iao_temperature_mtr5_05_10_2019"')
          obs_found_temp = True
          hgt_obs_temp = [hgt for hgt in range(0, 501, 50)]
          standard_deviation_temp = sigma(hgt_calc, hgt_obs_temp, tmp_calc, tmp_obs)
      except:
          #print('bad temp')
          tmp_obs = None
          hgt_obs_temp = None
      try:
          #obs_found = True
          wsp10_obs = export_data_obs_from_db(date, table_name='"iao_sodar_wind_speed_06_09_2019"')
          obs_found_wsp10 = True
          hgt_obs_wsp10 = [hgt for hgt in range(75, 501, 25)]
          if wsp10_obs != []:
            standard_deviation_wsp10 = sigma(hgt_calc, hgt_obs_wsp10, wsp10_calc, wsp10_obs)
      except:
          #print('bad wsp10')
          wsp10_obs = None
          hgt_obs_wsp10 = None
      try:
          #obs_found = True
          wdir10_obs = export_data_obs_from_db(date, table_name='"iao_sodar_wind_dir_06_09_2019"')
          obs_found_wdir10 = True
          hgt_obs_wdir10 = [hgt for hgt in range(75, 501, 25)]
          if wdir10_obs != []:
            standard_deviation_wdir10 = sigma(hgt_calc, hgt_obs_wdir10, wdir10_calc, wdir10_obs)
      except:
          #print('bad wdir10')
          wdir10_obs = None
          hgt_obs_wdir10 = None
      #wsp10_obs = None
      #wsp10_obs = [hgt for hgt in range(0, 6, 1)]
      #hgt_obs_wsp10 = None

  if not re.findall('gmc', file):
      tmp_calc = read_values(file, 5)
      hgt_calc = read_values(file, 11)
      wsp10_calc = read_values(file, 3)
      wdir10_calc = read_values(file, 4)
  else:
      with open(file) as f:
          values = np.loadtxt(f, skiprows=3, usecols=[3, 4, 5, 11], max_rows=6)
      f.close()
      wsp10_calc = values[:, 0]
      wdir10_calc = values[:, 1]
      tmp_calc = values[:, 2]
      hgt_calc = values[:, 3]

  #if there are no observations
  draw_graph_single(date, time, graph_csl_key, 'Temperature $[^\circ$C]', tmp_calc, hgt_calc, tmp_obs, hgt_obs_temp, temp_min, temp_max, 0, 500, 100, obs_found_temp,'temp',standard_deviation_temp)
  draw_graph_single(date, time, graph_csl_key, 'Wind Speed '+r'$[\frac{m}{s}]$',  wsp10_calc, hgt_calc, wsp10_obs, hgt_obs_wsp10, 0, 10, 0, 500, 100, obs_found_wsp10,'wsp10', standard_deviation_wsp10)
  draw_graph_single(date, time, graph_csl_key, 'Wind Direction $[^\circ$]',  wdir10_calc, hgt_calc, wdir10_obs, hgt_obs_wdir10, 0, 360, 0, 500, 100, obs_found_wdir10,'wdir10', standard_deviation_wdir10)
  #if there are observations
  if obs_found_temp:
      draw_graph_single(date, time, graph_csl_key, 'Temperature $[^\circ$C]', tmp_calc, hgt_calc, tmp_obs, hgt_obs_temp, temp_min, temp_max, 0, 500, 100, obs_found_temp, 'temp', standard_deviation_temp)
  if obs_found_wsp10:
      draw_graph_single(date, time, graph_csl_key, 'Wind Speed ' + r'$[\frac{m}{s}]$', wsp10_calc, hgt_calc, wsp10_obs, hgt_obs_wsp10, 0, 10, 0, 500, 100, obs_found_wsp10, 'wsp10', standard_deviation_wsp10)
  if obs_found_wdir10:
      draw_graph_single(date, time, graph_csl_key, 'Wind Direction $[^\circ$]', wdir10_calc, hgt_calc, wdir10_obs, hgt_obs_wdir10, 0, 360, 0, 500, 100, obs_found_wdir10, 'wdir10', standard_deviation_wdir10)
  data_calc_obs_temp = {'hgt_calc': hgt_calc, 'hgt_obs': hgt_obs_temp, 'tmp_calc': tmp_calc, 'tmp_obs': tmp_obs}
  data_calc_obs_wsp10 = {'hgt_calc': hgt_calc, 'hgt_obs': hgt_obs_wsp10, 'wsp10_calc': wsp10_calc, 'wsp10_obs': wsp10_obs}
  data_calc_obs_wdir10 = {'hgt_calc': hgt_calc, 'hgt_obs': hgt_obs_wdir10, 'wdir10_calc': wdir10_calc, 'wdir10_obs': wdir10_obs}
  return data_calc_obs_temp,data_calc_obs_wsp10,data_calc_obs_wdir10


folder_name = 'wrfout_d03_2018-01-02'
data_path = '.../Data/4plots/'+folder_name+'/'

#data point
data_list = glob.glob(data_path+'*dat')
file_list = []
for f in data_list:
  if re.findall('point', f):
     file_list.append(f)

#dictionary with parameters for displaying the values of calculations and observations on graphs
#color, symbol,label
graph_csl = {'airport':['blue','airport_calc','red','^','airport_obs'],'bec':['red','bec_calc','black','o','bec_obs'],
             'center':['orange','center_calc','white','*','center_obs'],'gmc':['green','gmc_calc','brown','D','gmc_obs'], 'tor':['black','tor_calc','green','o','tor_obs']}

#range for files with the same time interval
i = 0
j = 5
while j <= len(file_list):
    for file in file_list[i:j]:
        file_name = file.split('/')[-1]
        date = re.findall(r'\d{4}-\d{2}-\d{2}', file)[-1]
        time = re.findall(r'\d{2}:\d{2}:\d{2}', file)
        time = time[-1]
        date_time = date+' '+time
        date = datetime.strptime(date_time, "%Y-%m-%d %H:%M:%S")
        if re.findall('airport', file):
           data_airport = define_area_and_draw_graph(file, 'airport')
        if re.findall('bec', file):
           data_bec = define_area_and_draw_graph(file, 'bec')
        if re.findall('center', file):
           data_center = define_area_and_draw_graph(file, 'center')
        if re.findall('gmc', file):
           data_gmc = define_area_and_draw_graph(file, 'gmc')
        if re.findall('tor', file):
           data_tor = define_area_and_draw_graph(file, 'tor')
    data_all_temp = [data_airport[0], data_bec[0], data_center[0], data_gmc[0], data_tor[0]]
    data_all_wsp10 = [data_airport[1], data_bec[1], data_center[1], data_gmc[1], data_tor[1]]
    data_all_wdir10 = [data_airport[2], data_bec[2], data_center[2], data_gmc[2], data_tor[2]]
    temp_min, temp_max = temp_bounds(date.month)
    draw_graph_all(date, time, graph_csl, data_all_temp, 'Temperature $[^\circ$C]', temp_min, temp_max, 0, 500, 100,
                   'hgt_calc', 'hgt_obs', 'tmp_calc', 'tmp_obs', 'temp')
    draw_graph_all(date, time, graph_csl, data_all_wsp10, 'Wind Speed ' + r'$[\frac{m}{s}]$', 0, 10, 0, 500, 100,
                   'hgt_calc', 'hgt_obs', 'wsp10_calc', 'wsp10_obs', 'wsp10')
    draw_graph_all(date, time, graph_csl, data_all_wdir10, 'Wind Direction $[^\circ$]', 0, 360, 0, 500, 100, 'hgt_calc',
                   'hgt_obs', 'wdir10_calc', 'wdir10_obs', 'wdir10')
    i = j
    j += 5

