#!/usr/bin/python

#import ALL THE THINGS! (most are probably not needed. I'm too lazy to purge.)
import matplotlib
matplotlib.use('Agg')
from matplotlib.collections import PatchCollection
import urllib2
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import colors, ticker, cm
import numpy as np
import scipy
import math
from scipy.stats import gaussian_kde
from scipy.stats import kde
from scipy.interpolate import griddata
import scipy.ndimage as ndi
from scipy.spatial import KDTree
from math import radians, cos, sin, asin, sqrt
import datetime
import urllib
from matplotlib.colors import LogNorm
from itertools import groupby
import pyart
import sys
import warnings
warnings.filterwarnings("ignore")

#User Inputish area *****
home_directory = '/Users/blumberg/'
projection_type='lcc'
map_width=600000
map_height=400000
UND_names = ['Kevin Mahoney']
#************************

#Time stuff for the current time.
start_time = datetime.datetime.now()
c = datetime.datetime.now() + datetime.timedelta(hours=5)
year = c.year
month = c.month
day = c.day
hour = c.hour
minute = c.minute
if year<10:
    year_str="0"+str(year)
else:
    year_str=str(year)
if month<10:
    month_str="0"+str(month)
else:
    month_str=str(month)
if day<10:
    day_str="0"+str(day)
else:
    day_str=str(day)
if hour<10:
    hour_str="0"+str(hour)
else:
    hour_str=str(hour)
if minute<10:
    minute_str="0"+str(minute)
else:
        minute_str=str(minute)
chaser_time_string =year_str+month_str+day_str+"_"+hour_str+minute_str
chaser_time = c.strftime(" %B %d, %Y\n%H:%M UTC")


#Central Plains Radar Names and Locations
#Could not find a good dataset for this...so I made my own. haha
radars = ['BIS','MVX','MBX','ABR','UDX','FSD','OAX','LNX','UEX','TWX','DDC','GLD','ICT','FDR','INX','TLX','VNX',
          'AMA','BRO','CRP','FWS','HGX','GRK','DFX','LBB','MAF','EWX','SJT','DYX','SRX','LZK','FTG','GJX','PUX',
          'CYS','RIW','EAX','SGF','LSX','DVN','DMX','MPX','DLH']
radars_lats = [46.64,47.32,48.24,45.27,44.08,43.35,41.19,41.57,40.19,39.00,37.46,39.22,37.39,34.21,36.11,35.20,36.44,
               35.13,25.55,27.46,32.34,29.28,30.43,29.16,33.40,31.57,29.42,31.22,32.32,35.17,34.50,39.47,39.04,38.28,
               41.09,43.04,38.49,37.14,38.42,41.37,41.44,44.51,46.51]
radars_lons = [-100.45,-97.20,-100.52,-98.25,-102.50,-96.45,-96.22,-100.35,-98.27,-96.14,-99.58,-101.42,-97.26,-98.59,-95.34,-97.17,-98.08,
               -101.43,-97.25,-97.30,-97.18,-95.05,-97.24,-96.47,-101.49,-102.12,-98.02,-100.30,-99.15,-94.22,-92.16,-104.33,-108.13,-104.11,
               -104.48,-108.28,-94.16,-93.23,-90.41,-90.35,-93.43,-93.34,-92.12]
               
#Define functions
def second_smallest(numbers):
     m1, m2 = float('inf'), float('inf')
     for x in numbers:
         if x <= m1:
             m1, m2 = x, m1
         elif x < m2:
             m2 = x
     return m2
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367. * c
    return km
def nanargmax(a):
    idx = np.argmax(a, axis=None)
    multi_idx = np.unravel_index(idx, a.shape)
    if np.isnan(a[multi_idx]):
        nan_count = np.sum(np.isnan(a))
        # In numpy < 1.8 use idx = np.argsort(a, axis=None)[-nan_count-1]
        idx = np.argpartition(a, -nan_count-1, axis=None)[-nan_count-1]
        multi_idx = np.unravel_index(idx, a.shape)
    return multi_idx

#Pull in SN locations
lats=[]
lons=[]
names=[]
url = 'http://www.spotternetwork.org/feeds/gr-no.txt'
data = urllib2.urlopen(url)
for line in data:
    split_line = line.split( )
    if len(split_line) > 0:
        if split_line[0]=="Object:":
            latlons = split_line[1].split(',')
            lats.extend([float(latlons[0])])
            lons.extend([float(latlons[1])])
	if (split_line[0]=="Icon:") and (split_line[1].split(',')[2]=='000'):
            temporary = (split_line[1].rsplit(',')[5]+' '+split_line[2].split(',')[0])
            temporary = temporary.split('2')[0]
            temporary = temporary.rsplit(' ')
            name1 = temporary[0]
            name2 = temporary[1]
            name2 = name2.rsplit("\\")[0]
            names.extend([name1[1::]+' '+name2])

#UND Chaser locations ***********
UND_lats=[37.37545]
UND_lons=[-100.1950]
#i=0
#for UND_name in UND_names:
#    try:
#        index = names.index(UND_name)
#        UND_lats.extend([lats[index]])
#        UND_lons.extend([lons[index]])
#        if i==0:
#            with open("Mahoney_location.txt", "w") as text_file:
#                text_file.write("%s %s" % (UND_lats[0], UND_lons[0]))
#        i=i+1
#    except ValueError:
#        with open('Mahoney_location.txt') as f:
#            lines = f.readlines()
#            lines = lines[0].split()
#            UND_lats.extend([float(lines[0])])
#            UND_lons.extend([float(lines[1])])
#        continue
#**********************************

#Find nearest radar to first UND chaser***********
radar_dist = []
for i in range(len(radars)):
    radar_dist.extend([haversine(UND_lons[0], UND_lats[0], radars_lons[i], radars_lats[i])])
nearest_radar_idx = radar_dist.index(min(radar_dist))
nearest_radar = 'K'+radars[nearest_radar_idx]
second_nearest_radar_idx = radar_dist.index(second_smallest(radar_dist))
second_nearest_radar = 'K'+radars[second_nearest_radar_idx]
#*************************************************

#Obnoxious datetime string stuff to download recent radar data***********
d = datetime.datetime.now() + datetime.timedelta(hours=5)
date_time = d.strftime(" %B %d, %Y\n%I:%M%p CDT")
date_time_text = date_time
year = d.year
month = d.month
day = d.day
hour = d.hour
minute = d.minute
restart = True
i=0
while restart:
    if minute<0:
        minute = minute+60
        hour = hour-1
        print hour
    if hour<0:
        hour=hour+23
        day = day-1
        print hour
    if day<1:
        if month in [5,7,8,10,12]:
            day=day+30
            month=month-1
        if month in [1,2,4,6,9,11]:
            day=day+31
            month=month-1
        if month == 2:
            day=day+28
            month=month-1
    if month<1:
        month=month+12
        year = year-1
        
    if year<10:
        year_str="0"+str(year)
    else:
        year_str=str(year)
    month = d.month
    if month<10:
        month_str="0"+str(month)
    else:
        month_str=str(month)
    day = d.day
    if day<10:
        day_str="0"+str(day)
    else:
        day_str=str(day)
    if hour<10:
        hour_str="0"+str(hour)
    else:
        hour_str=str(hour)
    if minute<10:
        minute_str="0"+str(minute)
    else:
        minute_str=str(minute)
        
    date_string = year_str+month_str+day_str
    file_string = nearest_radar+'_'+date_string+'_'+hour_str+minute_str
    
    #where to get radar data and where to save it
    file_to_get = "http://mesonet-nexrad.agron.iastate.edu/level2/raw/"+nearest_radar+"/"+file_string
    file_save_name = "/Users/blumberg/Desktop/location_sa/"+file_string
    
    if i<1:
	i=i+1
	minute=minute-1
        continue
    try:
        testfile = urllib.URLopener()
        testfile.retrieve(file_to_get,file_save_name)
        print file_string+" EXISTS!!"
        break
    except IOError:
        i=i+1
        minute=minute-1
        if i==15:
            nearest_radar==second_nearest_radar
        continue
        
radar_time = hour_str+":"+minute_str+" UTC"
#*******************************************************************************************

#Have to place this here b/c of the basemap function that reads the following shapefiles
m = Basemap(width=map_width,height=map_height,projection=projection_type,resolution='i',
            lat_0=UND_lats[0],lon_0=UND_lons[0])

#Grab state shapefiles within nearest 200km*************************************
"""
cities_lons=[]
cities_lats=[]
cities=[]
states=[]
dist_to_cities=[]
cities_info = m.readshapefile(home_directory+'/spotter_network_density/shapefiles/cities','cities',drawbounds=False) #This is why I had to set basemap projection earlier
for info, att in zip(m.cities_info, m.cities): 
    dist = haversine(info['LON'], info['LAT'], UND_lons[0], UND_lats[0])
    if dist < 200 and len(info['ST'])==2:
        dist_to_cities.extend([dist])
        cities_lons.extend([info['LON']])
        cities_lats.extend([info['LAT']])
        cities.extend([info['NAME']])
        states.extend([info['ST']])
UND_states = list(set(states))
print UND_states
"""
#************************************************************************************

#Pull in Watch Data from SPC*********************************************************
sock = urllib.urlopen("http://www.spc.noaa.gov/products/watch/")
htmlSource = sock.readlines()
sock.close()
tstorm_websites = []
tornado_websites = []
pds_tornado_websites = []
pds_tstorm_websites = []
for line in htmlSource:
    if '.html">Severe Thunderstorm Watch' in line:
        new_line = line.rsplit('"', 1)[0]
        new_line = new_line.rsplit('"', 1)[1]
        new_line = new_line[18::]
        tstorm_website = 'http://www.spc.noaa.gov/products/watch/wou'+new_line
        tstorm_websites.extend([tstorm_website])
    if '.html">Tornado Watch' in line:
        new_line = line.rsplit('"', 1)[0]
        new_line = new_line.rsplit('"', 1)[1]
        new_line = new_line[18::]
        tornado_website = 'http://www.spc.noaa.gov/products/watch/wou'+new_line
        tornado_websites.extend([tornado_website])
    if '.html">PDS Tornado Watch' in line:
        new_line = line.rsplit('"', 1)[0]
        new_line = new_line.rsplit('"', 1)[1]
        new_line = new_line[18::]
        pds_tornado_website = 'http://www.spc.noaa.gov/products/watch/wou'+new_line
        pds_tornado_websites.extend([pds_tornado_website])
    if '.html">PDS Severe Thunderstorm Watch' in line:
        new_line = line.rsplit('"', 1)[0]
        new_line = new_line.rsplit('"', 1)[1]
        new_line = new_line[18::]
        pds_tstorm_website = 'http://www.spc.noaa.gov/products/watch/wou'+new_line
        pds_tstorm_websites.extend([pds_tstorm_website])
#**************************************************************************************
        
#use watch websites to pull in current watches*****************************************
tstorm_lat_list = []
tstorm_lon_list = []
tornado_lat_list = []
tornado_lon_list = []
for website in tstorm_websites:
    sock = urllib.urlopen(website)
    htmlSource = sock.readlines()
    sock.close()
    for line in htmlSource:
        if 'LAT...LON' in line:
            tstorm_lat_list.extend([[float(line.rsplit()[1][0:2]+'.'+line.rsplit()[1][2:4]),float(line.rsplit()[2][0:2]+'.'+line.rsplit()[2][2:4]),float(line.rsplit()[3][0:2]+'.'+line.rsplit()[3][2:4]),float(line.rsplit()[4][0:2]+'.'+line.rsplit()[4][2:4])]])
            tstorm_lon_list.extend([[float('-'+line.rsplit()[1][4:6]+'.'+line.rsplit()[1][6::]),float('-'+line.rsplit()[2][4:6]+'.'+line.rsplit()[2][6::]),float('-'+line.rsplit()[3][4:6]+'.'+line.rsplit()[3][6::]),float('-'+line.rsplit()[4][4:6]+'.'+line.rsplit()[4][6::])]])
for website in tornado_websites:
    sock = urllib.urlopen(website)
    htmlSource = sock.readlines()
    sock.close()
    for line in htmlSource:
        if 'LAT...LON' in line:
            tornado_lat_list.extend([[float(line.rsplit()[1][0:2]+'.'+line.rsplit()[1][2:4]),float(line.rsplit()[2][0:2]+'.'+line.rsplit()[2][2:4]),float(line.rsplit()[3][0:2]+'.'+line.rsplit()[3][2:4]),float(line.rsplit()[4][0:2]+'.'+line.rsplit()[4][2:4])]])
            tornado_lon_list.extend([[float('-'+line.rsplit()[1][4:6]+'.'+line.rsplit()[1][6::]),float('-'+line.rsplit()[2][4:6]+'.'+line.rsplit()[2][6::]),float('-'+line.rsplit()[3][4:6]+'.'+line.rsplit()[3][6::]),float('-'+line.rsplit()[4][4:6]+'.'+line.rsplit()[4][6::])]])
#****************************************************************************************

#Plot text
created_text = 'Created by David Goines\nUniversity of North Dakota\naero.und.edu/~dgoines'
title_text = 'OU/NSSL CLAMPS Trailer Location\nNearest Radar: '+nearest_radar+'\nLatest 0.5$^\circ$ scan: '+radar_time

#CREATE PLOT****************************************************************************************************************************************
print "Creating plot"
new_UND_lons, new_UND_lats = m(UND_lons,UND_lats)
new_lons, new_lats = m(lons,lats)
fig = plt.figure(figsize=(18, 11))
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8])
#*******************Map boundaries
m.drawcoastlines() # Draws the coastlines
m.drawcountries(linewidth=2.5) # Draw the country borders
m.drawmapboundary(fill_color='aqua') # This will plot the oceans the given color
m.drawstates(linewidth=1.5, zorder=3) # Draw the state borders
m.fillcontinents(color='#FFFFFF', lake_color='aqua', zorder=0) # Fill the continents with a color
#*******************Highways
#print "pull in road shapefiles..."
#for UND_state in UND_states:
#    highways_info = m.readshapefile(home_directory+'/spotter_network_density/shapefiles/roads/tl_2013_'+UND_state+'_prisecroads','highways',drawbounds=False)
#    for info, highway in zip(m.highways_info, m.highways):
#        if info['RTTYP'] == 'S':
#            x,y = zip(*highway)
#            plt.plot(x,y, marker=None,color='#A59A2A',linewidth=0.35, zorder=2)
#        if info['RTTYP'] == 'U':
#            x,y = zip(*highway)
#            plt.plot(x,y, marker=None,color='#A54F2A',linewidth=0.55, zorder=2)
#        if info['RTTYP'] == 'I':
#            x,y = zip(*highway)
#            plt.plot(x,y, marker=None,color='b',linewidth=1.0, zorder=2)       
#print "road shapefiles pulled in."
#********************Watches
for i in xrange(len(tstorm_lon_list)):
    tstorm_lat_list_new = []
    tstorm_lon_list_new = []
    for j in xrange(len(tstorm_lon_list[i])):
	if tstorm_lon_list[i][j] > -50.:
            tstorm_lon_list[i][j] = tstorm_lon_list[i][j]-100
	tstorm_lon_temp,tstorm_lat_temp = m(tstorm_lon_list[i][j],tstorm_lat_list[i][j])
        tstorm_lat_list_new.extend([tstorm_lat_temp])
        tstorm_lon_list_new.extend([tstorm_lon_temp])
    tstorm_lon_list_new.extend([tstorm_lon_list_new[0]])
    tstorm_lat_list_new.extend([tstorm_lat_list_new[0]])
    plt.plot(tstorm_lon_list_new,tstorm_lat_list_new,'b--',linewidth=1.5, zorder=6)
for i in xrange(len(tornado_lon_list)):
    tornado_lat_list_new = []
    tornado_lon_list_new = []
    for j in xrange(len(tornado_lon_list[i])):
        if tornado_lon_list[i][j] > -50.:
            tornado_lon_list[i][j] = tornado_lon_list[i][j]-100
	tornado_lon_temp,tornado_lat_temp = m(tornado_lon_list[i][j],tornado_lat_list[i][j])
        tornado_lat_list_new.extend([tornado_lat_temp])
        tornado_lon_list_new.extend([tornado_lon_temp])
    tornado_lon_list_new.extend([tornado_lon_list_new[0]])
    tornado_lat_list_new.extend([tornado_lat_list_new[0]])
    plt.plot(tornado_lon_list_new,tornado_lat_list_new,'r--',linewidth=2.0, zorder=6)
#******************Warnings
warnings_info = m.readshapefile('./current_ww','warnings',drawbounds=False)
for info, warning in zip(m.warnings_info, m.warnings):
    if (info['TYPE'] == 'TO' and info['GTYPE'] == 'P') and (info['STATUS'] != 'CAN' and info['STATUS'] != 'EXP' and info['STATUS'] != 'EXA' and info['STATUS'] != 'EXB' and info['ISSUED'] != 'None'):
        x,y = zip(*warning)
        plt.fill(x,y, edgecolor='k', facecolor='#FF0000',linewidth=0.75, alpha=0.2, zorder=5)
	plt.plot(x,y, color='#FF0000',linewidth=1.2, zorder=5)
    if (info['TYPE'] == 'SV' and info['GTYPE'] == 'P') and (info['STATUS'] != 'CAN' and info['STATUS'] != 'EXP' and info['STATUS'] != 'EXA' and info['STATUS'] != 'EXB' and info['ISSUED'] != 'None'):
        x,y = zip(*warning)
        plt.fill(x,y, edgecolor='k', facecolor='y',linewidth=0.65, zorder=4, alpha=0.15 )
	plt.plot(x,y, color='y',linewidth=1.1, zorder=4)
#*****************The Plot
ax1.annotate(title_text, xy=(0.5, 1.0), xytext=(-15, -15), fontsize=14,
    xycoords='axes fraction', textcoords='offset points',
    bbox=dict(facecolor='#F1F0F2', edgecolor='black', alpha=1.0, boxstyle='round,pad=0.5'),
    horizontalalignment='center', verticalalignment='top',zorder=9)
ax1.annotate(chaser_time, xy=(0.95, 1.0), xytext=(-15, -15), fontsize=11,
    xycoords='axes fraction', textcoords='offset points',
    bbox=dict(facecolor='#F1F0F2', edgecolor='black', alpha=1.0, boxstyle='round,pad=0.5'),
    horizontalalignment='center', verticalalignment='top',zorder=9)
ax1.annotate(created_text, xy=(0.11, 0.10), xytext=(-15, -15), fontsize=10,
    xycoords='axes fraction', textcoords='offset points',
    bbox=dict(facecolor='#F1F0F2', edgecolor='black', alpha=0.8, boxstyle='round,pad=0.5'),
    horizontalalignment='center', verticalalignment='top',zorder=9)
#*****************Legend
#plt.plot([1000,1000],[1000,1000], marker=None,color='#A59A2A',linewidth=0.35, zorder=10,label='State Highway')
#plt.plot([1000,1000],[1000,1000], marker=None,color='#A54F2A',linewidth=0.55, zorder=10,label='U.S. Highway')
#plt.plot([1000,1000],[1000,1000], marker=None,color='b',linewidth=1.0, zorder=10,label='Interstate Highway') 
plt.plot([1000,1000],[1000,1000],'b--',linewidth=1.5, zorder=10,label="Svr. T-Storm Watch")
plt.plot([1000,1000],[1000,1000],'r--',linewidth=2.0, zorder=10,label="Tornado Watch")
plt.fill([1000,1000],[1000,1000], edgecolor='k', facecolor='y',linewidth=0.65, zorder=10, alpha=0.15, label='Svr. T-Storm Warning')
plt.fill([1000,1000],[1000,1000], edgecolor='k', facecolor='#FF0000',linewidth=0.75, alpha=0.2, zorder=10, label='Tornado Warning')
#plt.scatter(new_lons,new_lats,s=8., marker='o',c='k',zorder=7, alpha=0.75, label='Chasers/Spotters')
if len(new_UND_lons)>0:
    plt.scatter(new_UND_lons,new_UND_lats,s=280., marker='*',c='#FF33FF',edgecolor='k',linewidth=0.125,zorder=7, alpha=1.0, label='CLAMPS')
#atts_info = m.readshapefile(home_directory+'/spotter_network_density/shapefiles/nexrad_attributes/current_nexattr','atts',drawbounds=False)
#TVS_count=0
#for info, att in zip(m.atts_info, m.atts):
#    if info['TVS'] == 'TVS' and info['MESO'] != 'NONE' and int(info['RANGE']) <= 100:
#        TVS_count=TVS_count+1
#        if TVS_count==0:
#            plt.scatter(att[0],att[1],s=90., marker='v',c='r',edgecolor='r',zorder=8, alpha=1.0, label='Tor. Vortex Signature')
#        else:
#            plt.scatter(att[0],att[1],s=90., marker='v',c='r',edgecolor='r',zorder=8, alpha=1.0)
#plt.legend( loc=2, borderaxespad=0.,prop={'size':10})
leg = plt.legend( loc=2, borderaxespad=0.,prop={'size':10})
leg = leg.get_frame()
leg.set_facecolor('#F1F0F2')
m.drawcounties(zorder=3)

#PLOT RADAR DATA*******************************************************************
filename = "/Users/blumberg/Desktop/location_sa/"+file_string
radar = pyart.io.read(filename)
display = pyart.graph.RadarMapDisplay(radar)
#set projection and map for PyArt. It should be the same as the first projection
CS = display.plot_ppi_map('reflectivity', 0, vmin=10, vmax=75,
                     lat_0=UND_lats[0],lon_0=UND_lons[0],
                     width=map_width,height=map_height, projection=projection_type,
                     resolution='i',colorbar_flag=False,title_flag=False,
                     lat_lines=None, lon_lines=None, embelish=False,
                     cmap = pyart.graph.cm.NWSRef)

#Did this to get my own colorbar******
cint_refl = 5.0
cbar_refl_min = 10
cbar_refl_max = 75 + (cint_refl)
cflevs_refl = np.arange(cbar_refl_min, cbar_refl_max, cint_refl)
refl_ticks = np.arange(cbar_refl_min, cbar_refl_max, cint_refl)
blah1 = np.random.rand(3,2)
blah2 = np.random.rand(3,2)
blah3 = np.random.rand(3,2)
CS1 = plt.contourf(blah1,blah2,blah3,levels=cflevs_refl, cmap = pyart.graph.cm.NWSRef)   
cbar_ax = fig.add_axes([0.85, 0.05, 0.05, 0.85],frameon=False)
cbar = plt.colorbar(CS1, pad=0.05, fraction=0.40, orientation='vertical')
cbar_ax.set_xticks([]) 
cbar_ax.set_yticks([]) 
cbar.set_ticks(refl_ticks)
cbar.set_label('Radar Reflectivity (dBZ)',size=16)
clabs = ['%i' % f for f in refl_ticks]
#***************************************

plt.show()
#Save file stuff
save_name = "/Users/blumberg/Desktop/location_sa/UND_zoom_"+chaser_time_string+".png"
plt.savefig(save_name, bbox_inches='tight')
print "Plot created!!!"
#****************************************************************************************************************************************

end_time = datetime.datetime.now()

total_time = end_time-start_time
print "Total time: "+str(total_time)


sys.exit()
