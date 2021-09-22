#!/usr/bin/python3

#DATA T2, TSK
import glob
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.ticker import MaxNLocator
from datetime import datetime, timedelta
import numpy as np
import os
import re
from mpl_toolkits.basemap import Basemap

def parName(file_list,i):
    with open(file_list[0]) as file:
        par_name = np.loadtxt(file, dtype='str', usecols=i)[0]
    file.close()
    return par_name

#read coordinates lon,lat
def read_bounds_lonlat(file_list,cols):
    with open(file_list[0]) as file:
         coord = np.loadtxt(file,skiprows=1,usecols=cols)
    file.close()
    xlong_min = coord[0, 0]
    xlong_max = coord[-1, 0]
    xlat_min = coord[0, 1]
    xlat_max = coord[-1, 1]
    return xlong_min, xlong_max, xlat_min, xlat_max
#read coordinates x,y
def read_coord_lonlat(file_list,cols):
    with open(file_list[0]) as file:
         coord = np.loadtxt(file,skiprows=1,usecols=cols)
    file.close()
    x = coord[:, 0].reshape(51, 51)
    y = coord[:, 1].reshape(51, 51)
    return x, y

def scale_bounds(file_list,i):
    with open(file_list[0]) as file:
        z = np.loadtxt(file,skiprows=1,usecols=i)
    file.close()
    minv = min(z)
    maxv = max(z)

    for f in file_list[1:]:
        with open(f) as file:
            z = np.loadtxt(file,skiprows=1,usecols=i)
        file.close()
        if minv>min(z):
            minv=min(z)
        if maxv<max(z):
            maxv = max(z)
    return minv, maxv

def draw_picture(file_list,x,y,par_name,scale_name,path):

    for f in file_list:
        file_name = f.split('/')[-1]
        print(file_name)
        #read values meteo parameter
        with open(f) as file:
            data = np.loadtxt(file,skiprows=1,usecols=[6,7,8])
        file.close()
        u10 = data[:,0].reshape(51,51)
        v10 = data[:,1].reshape(51,51)
        t2 = data[:, 2].reshape(51, 51)
        date = re.findall(r'\d{4}-\d{2}-\d{2}', f)
        time = re.findall(r'\d{2}:\d{2}:\d{2}', f)
        time = time[-1]
        
        date = datetime.strptime(date[-1], "%Y-%m-%d")

        #levels = MaxNLocator(nbins='auto').tick_values(minv, maxv)

        levels = np.arange(round(minv)-1, round(maxv)+2,1)
        fig = plt.figure(figsize=(22, 12))
        fig.suptitle(par_name + ' ' + date.strftime("%d-%m-%Y") + ' ' + time, fontsize=20)
        m = Basemap(llcrnrlon=nllcrnrlon, llcrnrlat=nllcrnrlat, urcrnrlon=nurcrnrlon, urcrnrlat=nurcrnrlat, projection='lcc', lat_1=50., lat_2=60., lat_0=56.5, lon_0=85)
        #River
        m.readshapefile('.../Data/shapefiles/River/RUS_water_areas_dcw', name='RUS_water_areas_dcw', drawbounds=True)
        #Buildings
        #m.readshapefile('.../Data/shapefiles/Buildings/gis_osm_buildings_a_free_1', name='gis_osm_buildings_a_free_1', drawbounds=True)
        x, y = m(lons,lats)
        # Academ
        xpt, ypt = m(85.04886, 56.47301)
        #lonpt, latpt = m(xpt, ypt)
        m.plot(xpt, ypt, 'bo')
        plt.text(xpt,ypt,'ACADEM', fontsize = 15, color = 'w',bbox=dict(facecolor='k', alpha=0.3))
        # Aeroport
        xpt, ypt = m(85.21121, 56.38289)
        #lonpt, latpt = m(xpt, ypt)
        m.plot(xpt, ypt, 'bo')
        plt.text(xpt,ypt,'Bogashevo', fontsize = 15, color = 'w',bbox=dict(facecolor='k', alpha=0.3))
        # GMC
        xpt, ypt = m(84.96747, 56.44603)
        #lonpt, latpt = m(xpt, ypt)
        m.plot(xpt, ypt, 'bo')
        plt.text(xpt,ypt,'GMC', fontsize = 15, color = 'w',bbox=dict(facecolor='k', alpha=0.3))
        # BEC
        xpt, ypt = m(85.09772, 56.48197)
        #lonpt, latpt = m(xpt, ypt)
        m.plot(xpt, ypt, 'bo')
        plt.text(xpt,ypt,'BEC', fontsize = 15, color = 'w', bbox=dict(facecolor='k', alpha=0.3))
        cs = m.contourf(x, y, t2, levels, cmap='jet')
        m.drawparallels(np.arange(nllcrnrlat, nurcrnrlat, 0.2), labels=[True, False, False, False])
        m.drawmeridians(np.arange(nllcrnrlon, nurcrnrlon, 0.1), labels=[False, False, False, True])
        m.colorbar(cs, label=scale_name)
        Q = m.quiver(x, y, u10, v10, scale=110)
        plt.quiverkey(Q, 0.5, 0.9, 1.0, r'$1.0 \frac{m}{s}$', labelpos='E',coordinates='figure',fontproperties={'weight': 'bold'})
        fig.savefig(path+par_name+'_'+file_name+'.png')
        #plt.show()
    return 0


file_list = glob.glob('.../Data/4plots/wrfout_d03_2018-01-02/*surf.dat')
path1 = '.../plots_2018/surf/t2/'
path2 = '.../plots_2018/surf/tsk/'

#T2
i = 8
scale_name = 'Temperature, $[^\circ$C]'
par_name = parName(file_list,i)
#column numbers XLONG, XLAT
cols=[2,3]
nllcrnrlon, nurcrnrlon, nllcrnrlat, nurcrnrlat  = read_bounds_lonlat(file_list,cols)
#print(nllcrnrlon)
#print(nurcrnrlon)
#print(nllcrnrlat)
#print(nurcrnrlat)
lons, lats = read_coord_lonlat(file_list,cols)
minv, maxv = scale_bounds(file_list,i)
#print(minv,maxv)
draw_picture(file_list,lons,lats,par_name,scale_name,path1)

#TSK
i = 12
scale_name = 'Temperature skin $[^\circ$C]'
par_name = parName(file_list,i)
cols=[2,3]
nllcrnrlon, nurcrnrlon, nllcrnrlat, nurcrnrlat  = read_bounds_lonlat(file_list,cols)
lons, lats = read_coord_lonlat(file_list,cols)
minv, maxv = scale_bounds(file_list,i)
#print(minv,maxv)
draw_picture(file_list,lons,lats,par_name,scale_name,path2)

