
import os 
import glob
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from herbie.archive import Herbie

os.chdir("R:/landscape_ecology/projects/hurrisk/validation/")
tracks = gpd.read_file("track_segments.shp")
tracksOld = gpd.read_file("old_tracks/track_segments.shp")

# Find the new tracks that were not in the old tracks shapefile.
tracks = tracks[~tracks.track_id.isin(tracksOld.track_id)]

def decrementDate(date):
    if int(date[11:13]) == 00 and (int(date[8:10]) == 1):
        dateNew = date[:11] + '23' + date[13:]
        dateNew = dateNew[:5] + str(int(dateNew[5:7])-1).zfill(2) + dateNew[7:]
        if int(date[5:7]) in [2,4,6,8,9,11]:
            dateNew = dateNew[:8] + '31' + dateNew[10:]
        if int(date[5:7]) in [1,5,7,8,10,12]:
            dateNew = dateNew[:8] + '30' + dateNew[10:]
        if int(date[5:7]) == 3:
            dateNew = dateNew[:8] + '28' + dateNew[10:]
    if int(date[11:13]) > 0:
        dateNew = date[:11] + str(int(date[11:13])-1).zfill(2) + date[13:]
    if int(date[11:13]) == 0 and int(date[8:10]) > 1:
        dateNew = date[:11] + '23' + date[13:]
        dateNew = dateNew[:8] + str(int(dateNew[8:10])-1).zfill(2) + dateNew[10:]
    return(dateNew)

        
# Format the start and end dates for Herbie for each hurricane
starts = [str(tracks['int_s'][i]) for i in tracks.index]#[1:] # The first date is too early
startDates = [starts[i][0:4]+'-'+starts[i][4:6]+'-'+starts[i][6:8]+' '+starts[i][8:10]+':00' for i in range(len(starts))]
ends = [str(tracks['int_e'][i]) for i in tracks.index]#[1:] # The first date is too early
endDates = [ends[i][0:4]+'-'+ends[i][4:6]+'-'+ends[i][6:8]+' '+ends[i][8:10]+':00' for i in range(len(ends))]

del startDates[15] # 2016-11-18 was problematic
del endDates[15]

wsMaps = {}
dirMaps = {}
for i in range(len(startDates)):
    date = startDates[i]
    wsMaps[startDates[i]] = '' # place holder. The raster will go here
    dirMaps[startDates[i]] = ''
    
    dar = [] # empty list for wind speed rasters
    wdar = [] # empth list for wind direction rasters
    forecast = 0
    while date != endDates[i]:
        H = Herbie(date, model='hrrr', product='sfc', fxx=forecast)
        H.download()
        
        """ Error handling. Catches missing data, and finds a previous forecast as an alternative. """
        try: 
            dsx = H.xarray('UGRD:10 m') # I use u wind because it's convenient. Could use v too.
        except ValueError: # if the grib2 file is missing, decrement time by 1 hour and get forecast
            forecast += 1 #
            date = decrementDate(date)
            print("fart")
            continue
        except AttributeError: # if the grib2 file is missing, decrement time by 1 hour and get forecast
            forecast += 1 # 
            date = decrementDate(date)
            print("ding")
            continue
        except IndexError:
            forecast += 1
            date = decrementDate(date)
            print("poop")
            continue
        finally:
            print("Oops, I sharted again.")
            
        date = date[:11] + str(int(date[11:13])+forecast).zfill(2) + date[13:] # put date back
        if int(date[11:13]) > 23:
            diff = int(date[11:13]) - 23 -1
            date = date[:8] + str(int(date[8:10])+1).zfill(2) + " " + str(diff).zfill(2) + date[13:]
        
        forecast = 0 # reset forecast

        """ Extract the of the data. """
        dsx = dsx.data_vars['u10'].to_numpy()
        dsy = H.xarray('VGRD:10 m')
        dsy = dsy.data_vars['v10'].to_numpy()
        
        
        """ Calculate wind speed, make a raster, and append it to the list. """
        dar.append(np.sqrt(dsx**2+dsy**2))
        
        """ Calculate wind direction """
        wd = np.arctan2(dsx, dsy) # angle in radians
        wd = np.degrees(wd) # convert to degrees
        wd[wd<0] = wd[wd<0]+360 # convert negative angles to positive
        wdar.append(wd)

        """ Delete these files. They take up too much space. """
        pathTmp = H.get_localFilePath().as_posix()[:40]
        os.chdir(pathTmp)
        [os.remove(file) for file in os.listdir()]

        """ Increment the date """
        if int(date[11:13]) <= 23:
            date = date[:11] + str(int(date[11:13])+1).zfill(2) + date[13:]
        if int(date[11:13]) == 23 and (int(date[5:7]) in [1,3,5,7,8,10,12]) and int(date[8:10]) == 31:
            date = date[:11] + '00' + date[13:]
            date = date[:5] + str(int(date[5:7])+1).zfill(2) + date[7:]
            date = date[:8] + "01" + date[10:]
        if int(date[11:13]) == 23 and (int(date[5:7]) in [2,4,6,9,11]) and int(date[8:10]) == 30:
            date = date[:11] + '00' + date[13:]
            date = date[:5] + str(int(date[5:7])+1).zfill(2) + date[7:]
            date = date[:8] + "01" + date[10:]
        if int(date[11:13]) == 24:
            date = date[:11] + '00' + date[13:]
            date = date[:8] + str(int(date[8:10])+1).zfill(2) + date[10:]
        
    dar = np.dstack(dar)
    maxws = dar.max(axis=2)
    maxws = np.flipud(maxws)
    wsMaps[startDates[i]] = maxws
    
    wdar = np.dstack(wdar)
    m,n = dar.shape[:2]
    maxdir = wdar[np.arange(m)[:,None],np.arange(n),dar.argmax(2)]
    maxdir = np.flipud(maxdir)
    dirMaps[startDates[i]] = maxdir
    
    print(wsMaps.keys())
 
"""Save the maps as arrays"""
os.chdir("R:/landscape_ecology/projects/hurrisk/validation/HRRR_outputs/")
for date in wsMaps.keys():
    np.save(file = "WS_" + date[:10] + ".npy", arr=wsMaps[date])
    
for date in dirMaps.keys():
    np.save(file = "DIR_" + date[:10] + ".npy", arr=dirMaps[date])


""" Should have saved maps in another format. Here I save to text. The lat, lon,
WS, and dir files should match sequentially, and we'll be able to get WS and/or direction
from the nearest lat and lon combination. """
H = Herbie(startDates[0], model='hrrr', product='sfc', fxx=0)
H.download()
t = H.xarray("UGRD:10 m")
t("latitude")
tLat = t.coords['latitude'].to_numpy()
tLon = t.coords['longitude'].to_numpy()
np.savetxt("latitude.txt", tLat, fmt="%f") # Latitude text file
np.savetxt("longitude.txt", tLon, fmt="%f") # longitude text file

# convert saved arrays to .txt files
files = glob.glob('*.npy', recursive=False)
for file in files:
    tmp = np.load(file)
    np.savetxt(file + '.txt', tmp, fmt='%f')
    
    
