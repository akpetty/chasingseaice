#!/usr/bin/python
#;--------------------------------------------------------------------------------
#; Program: icedrift.pro
#; Author: Linette Boisvert
#; Date created: August 29, 2018
#;
#; Python-translated by Rob Russell, 10/26/2018
#' Modified to read in a .ini file with the parameters needed, and fixed a couple of things
#; Note: If you want to disable the plot (old and new sequences) you can comment out the 
#; plotting code near the bottom.
#; Requires: numpy (python-numpy) and optionally matplotlib if you want to plot stuff.
#;
#;This program takes in John Sonntag's sequence and flight plan files, extracts the longitude and latitude (from 
#; sequence file) and elapsed time (from plan file) and calculates the new longitude and latiude waypoints with
#; drift correction. New waypoints are output into the same file type and formatting as John's sequence file for 
#; him to injest and change the flight path of the plane. This is done so that roughly the same sea ice that 
#; ICESat-2 flew over is also sampled by IceBridge at a slightly later time.
#; 
#; Code is made specifically for IceBridge Mid Weddell and West Weddell Sea Ice flights
#; 
#; *note currently the .plan file needs to be modified before it can be injested. this is done by just deleting the lines until after the 
#; headers for the 8 column section
#; ex: Y04211 Y04210    6.1     6.0   250      30.5  0.12  2609.4  8.50
#; * also note: the 180 turn lines must also be deleted from the plane file for this to work.
#; 
#; * note: zulu time is 3 hours ahead  of punta & ushuaia time, 
#;
#; inputs: latitude (deg), longitude (deg), take off time (hr), take off time minute (min), icesat2 cross over time hour (hr), icesat2 cross over time minutes (min), wind speed (m/s) and wind direction (degrees), 
#; pressure (mb) measurements taken from the dc8.
#; 
#; latitude and longitude from the plane: in degrees, -90,90 and -180 to 180 E
#; 
#; outputs: latitude and longitude adjusted for drift of sea ice.
#;
#;function icedrift, take_off_time, IS2_CROSSover_TIME, WIND_SPEED_PLANE, WIND_ANGLE_PLANE, PLANE_PRESSURE, PLANE_LATITUDE, PLANE_LONGITUDE
#;
#;---------------------------------------------------------------------------------------------
#function icedrift, take_off_hour,take_off_minutes, IS2_CROSSover_hour, IS2_crossover_minutes, WIND_SPEED_PLANE, WIND_ANGLE_PLANE, PLANE_PRESSURE, PLANE_LATITUDE, PLANE_LONGITUDE, sequence_file, plan_file, output_file
import numpy as np
import math
import os
import sys 

def icedrift(take_off_hour, take_off_minutes,IS2_CROSSOVER_HOUR,IS2_CROSSOVER_MINUTES,WIND_SPEED_PLANE,WIND_ANGLE_PLANE,PLANE_PRESSURE,PLANE_LATITUDE, PLANE_LONGITUDE, sequence_file, plan_file, output_file):
    #drift_corrected=np.array((5,100))
    #;put in variables (if you want) 
    #take_off_hour=13
    #take_off_minutes=38
    #IS2_CROSSOVER_HOUR=18
    #IS2_CROSSOVER_MINUTES=00
    #;19:36 in zulu time 
    #WIND_SPEED_PLANE=26
    #WIND_ANGLE_PLANE=315
    #PLANE_PRESSURE=964.9
    #PLANE_LATITUDE=-70.25 
    #PLANE_LONGITUDE=-046.06
    #sequence_file='siseelyeloop_10_19_18.sequence'
    #plan_file='siseelyeloop_10_19_18.plan'
    #output_file='siseelyeloop_drift_10_19_18_a.sequence'

    #;-----------------------------------------------------------------------
    #;convert knots to m/s for windspeed
    WIND_SPEED_PLANE=WIND_SPEED_PLANE*0.51444444444

    #; set pi constant
    pi=(4.0*math.atan(1.0))
    #;---------------------------------------------------------
    #; read in sequence file

    file_latlon= open(sequence_file,'r')# ;'siwestweddell.sequence'
    #;file_latlon= 'simidweddell.sequence'
    #openr,1,file_latlon
    temp1 = file_latlon.readlines()

    #; read in flight plan file
    file_time= open(plan_file,'r')  #siwestweddell.plan'
    #;file_time= 'simidweddell1.plan'
    # openr,2,file_time
    temp = file_time.readlines() 

    #;----------------------------------------------------------
    #;read in lat, longitude from ERA fields and the all of the slopes, intercepts and angles for converting 
    #; wind measurements taken from the plane to surface winds and directions

    # 3d float array 
    all_slope_xy = np.fromfile('all_slope_xy_OIB', '<f4').reshape(111, 51,6)
    #openr,10,'all_slope_xy_OIB'
    #readu,10,all_slope_xy
    #close,10
    
    # 3d float array 
    #all_intercept_xy=fltarr(6,51,111)
    #openr,11,'all_intercept_xy_OIB'
    all_intercept_xy = np.fromfile('all_intercept_xy_OIB', '<f4').reshape(111, 51,6)
    #readu,11,all_intercept_xy
    #lose,11
    
    # 3d float array 
    #all_angle_xy=dblarr(6,51,111) is this a double array? Yes, it is.
    all_angle_xy = np.fromfile('all_angle_xy_OIB', '<f8').reshape(111, 51,6)
    #openr,12,'all_angle_xy_OIB'
    #readu,12,all_angle_xy
    #close,12

    # Net CDF format, unless otherwise noted - use the downloaded ERA files (32 bit floats)
    #filen= 'october_2010_weddell_surface.nc'
    #ncid = NCDF_OPEN(filen)
    era_lon = np.fromfile("ERA_lon",'f4')
    #era_lon = np.fromfile("era_lon1")
    #f1.close()
    #result1 = ncdf_varid(ncid,'longitude') ;read in lons
    #ncdf_varget,ncid,result1,ERA_lon
    
    era_lat = np.fromfile("ERA_lat",'f4')
    #result2 = ncdf_varid(ncid,'latitude') ; read in lats
    #ncdf_varget,ncid,result2, ERA_lat

    #;----------------------
    #;STEP1: Calculate the windspeed and direction at the surface from height
    #;
    #;equation for wind speed at surface will be provided by Ed Blanchard-Wrigglesworth
    #;equation for direction of ice motion will be provided by Ron Kwok
    #;
    #;inputs:   windspeed, direction and pressure
    #;output: windspeed and direction at surface (aka. velocity and direction of ice)
    #;---------------------
    ilevs=np.array([850,875,900,925,950,975,1000])    #;%pressure levels

    #;compute the closest pressure level from the pressure level of the plane
    pres_diff=PLANE_PRESSURE-ilevs
    min_pres = min(abs(pres_diff))
    #lev_min=where(abs(pres_diff) eq min_pres)

    #lev_min=abs(pres_diff)==min_pres 
    # This returns the INDEX of the closest value 
    lev_min = (np.abs( abs(pres_diff) - abs(min_pres))).argmin()

    #;compute the nearest latitude and longitude from the plane to those given by ERA and the equations  
    lat_n = PLANE_LATITUDE - era_lat 
    # Closest lat
    #near_lat=(where(abs(lat_n) eq min(abs(lat_n))))
    near_latn = min(abs(lat_n))
    near_lat = (np.abs( abs(lat_n) - abs(near_latn))).argmin()
    
    lon_n = PLANE_LONGITUDE - era_lon 
    # Closest lon
    #near_lon=where(abs(lon_n) eq min(abs(lon_n)))
    near_lonn = min(abs(lon_n))
    near_lon = (np.abs( abs(lon_n) - abs(near_lonn))).argmin()
    
    #; surface winds and direction calculations
    # 100 -1.07985 is Z
    # 010 -0.8167  is Y
    # 001 -0.07738 is X
    
    #surface_wind=(WIND_SPEED_PLANE-all_intercept_xy(lev_min,near_lat, near_lon))/all_slope_xy(lev_min, near_lat,near_lon) #;calculate winds at surface
    surface_wind=(WIND_SPEED_PLANE-all_intercept_xy[near_lon,near_lat, lev_min])/all_slope_xy[near_lon, near_lat, lev_min] #;calculate winds at surface
    surface_direction=WIND_ANGLE_PLANE-all_angle_xy[near_lon ,near_lat,lev_min] #; calculate angle of winds at surface

    #; find the closest longitude point to the ones given by ron
    #ron_lon=360.0-[0,15,30,45,60]-360.0 #; to get it into -180 to 180 E
    ron_lon1 = np.array([0,15,30,45,60])
    ron_lon = 360 - ron_lon1 - 360
    lon_ron=PLANE_LONGITUDE-ron_lon
    min_diff=min(abs(lon_ron))

    #close_lon=where(abs(lon_ron) eq min_diff)
    close_lon = (np.abs( abs(lon_ron) - abs(min_diff))).argmin()
    
    #;velocity scaling factors (turn winds into ice velocity):
    #scaling_factor=[0.012,0.011,0.010,0.009,0.008]*2.0
    scaling_factor=np.array([0.012,0.011,0.010,0.009,0.008]) * 2
    scaling=scaling_factor[close_lon] 

    #;this will be computed via an equation from ron kwok
    ice_velocity=  (surface_wind/1000.0)*scaling #;km/s

    #;wind direction from wind angle taken from the plane into radians
    ice_dir=((surface_direction-25.0))*(pi/180)
    print "Velocity : " + str(ice_velocity) + " " + "Direction " + str(ice_dir/(pi/180))
    #;stop
    
    #;convert to u and v components of the wind
    #;U_component=wind_speed*cos(wind_angle)
    #;V_component=wind_speed*sin(wind_angle)

    #;-----------------------------------------------------------------
    #; create/open the sequence file that will be output to John Sonntag

    #openw,3, output_file
    outf = open(output_file,"w")

    #; read in the first 3 lines from the sequence file that have PUQ (starting point) and record them in the 
    #; output sequence file
    #readf,1,temp1

    # temp1 is the sequence file
    # temp is the plan file
    outf.write(temp1[0])
    outf.write(temp1[1])
    outf.write(temp1[2])

    tsplit_1st_WP = temp1[0]
    # Make sure that the below code is doing the above code
    #tsplit_1st_WP = strsplit(temp1,/extract)
    #printf,3, tsplit_1st_WP
    #readf,1,temp1
    #tsplit = strsplit(temp1,/extract)
    #printf,3, tsplit
    #readf,1,temp1
    #tsplit = strsplit(temp1,/extract)
    #printf,3, temp1
    allNewLats=[]
    allNewLons=[]
    allOldLats=[]
    allOldLons=[]

    count=0
    #;-------------------------------------------------------------------
    #; go through each line of the sequence file (starting with line 4 of the sequence file and line one 
    #; of the plan file) and adjust the waypoints according to drift, until the end of the file is reached
    #WHILE ~EOF(1) DO BEGIN
    for i in range(3,len(temp1)): #temp1[3:]:
        #readf,1,temp1
        #readf,2,temp
        #print, temp1
        #print, temp
        #tsplit = strsplit(temp,/extract)
        #tsplit1=strsplit(temp1,/extract)
        tsplit = temp[i-3].split()
        tsplit1 = temp1[i].split()
        # If the line has 180 at the beginning then account for that by selecting the right fields
        #; when 180 degree turns are encountered in the files, elasped time must still be accounted for
        #; latitude and longitude must remain the same for the turns? I believe this is correct...check with John  
        #if tsplit(0) eq '180' then elapsed_time=float(tsplit(5)) else  elapsed_time=float(tsplit(8))
        if tsplit[0]=="180":
            #elapsed_time=float(tsplit[5]) 
            #lat=lat_prev
            #lon=lon_prev ???? THESE DONT EXIST YET
            print "180 degree turn, skipping..." 
        elif tsplit1[0]=="PUQ" or tsplit1[0]=="SAWH":
            outf.write(temp1[i])
        else:
            elapsed_time = float(tsplit[8]) 
            lat=float(tsplit1[1])
            allOldLats.append(lat)
            lon=float(tsplit1[2]) 
            allOldLons.append(lon)
            #if tsplit(0) eq '180' then lat=lat_prev else lat=float(tsplit1(1))
            #if tsplit(0) eq '180' then lon=lon_prev else lon=float(tsplit1(2))

            #;---------------------
            #; STEP2: Calculate time difference between IceSat-2 crossing and when OIB reaches that point
            #; 
            #; this will come from John Sonntag
            #; 
            #; input: IS2 time, oib take off time from punta
            #; output: total time difference
            #;----------------------

            #;icebridge time in minutes:
            take_off_time=take_off_hour*60.0*60.0+take_off_minutes*60.0

            #;is2 time in minutes:
            is2_crossover_time=IS2_CROSSOVER_HOUR*60.0*60.0+IS2_CROSSOVER_MINUTES*60.0

            #;compute time difference between takeoff and is2 crossover

            IS2_time=is2_crossover_time-take_off_time

            #;convert time to seconds
            time_diff_to_IS2=(elapsed_time*60.0*60-IS2_time)
            #;time_diff_to_IS2=(IS2_time+elapsed_time*60.*60.)
            print "Time diff to IS2 " + str(time_diff_to_IS2)


            #;---------------------
            #; STEP3: Using time difference from step2 and velocty and direction of sea ice drift from step1 calculate the total distance moved
            #; 
            #; input: time, drift speed
            #; output: sea ice travel distance
            #;---------------------

            #;calculate the distance traveled based on ice drift and time: time must be seconds, and velocity must be in km/s
            dis_traveled= time_diff_to_IS2*ice_velocity

            #;--------------------
            #; STEP4: Use John's code with lat/lon and angle to get new location of way point
            #; 
            #; input: original lat long, distance, angle
            #; output: new lat lon
            #;--------------------
            #;
            #;PURPOSE:  Computes ending waypoint lat2/lon2 given starting waypoint
            #;lat1/lon1, starting true course tc1 and distance dist. Lat/lon
            #;are in radians N latitude/E longitude, tc1 is in radians,
            #;and dist is in kilometers.
            #;
            #;#define PI (4.0*atan((double)(1.0)))
            #;#define DEG2RAD (PI/180.0)
            #;#define RAD2NM  (180.0*60.0/PI)
            #;#define RAD2KM  (RAD2NM*6076.1*12.0*2.54/100.0/1000.0)
            #;
            #;AUTHOR:   John GARY Sonntag
            #;
            #;DATE:     July 1999
            #;*------------------------------------------------------------------------*/

            #; convert original latitude and longitude to radians
            original_latitude=lat*(pi/180.0)
            original_longitude=lon*(pi/180.0)

            original_longitude=-original_longitude

            #;produce constants
            pi=(4.0*math.atan(float(1.0)))
            rad2nm=(180.0*60.0/pi)
            rad2km=(rad2nm*6076.1*12.0*2.54/100.0/1000.0)

            #;***** subtract 360 degrees from the ice_dir to get which way the ice is moving towards... I think
            #;
            #;calculate new latitude waypoint
            #;changed it to negative cos...
            new_lat=math.asin(math.sin(original_latitude)*math.cos(dis_traveled/rad2km)-math.cos(original_latitude)*math.sin(dis_traveled/rad2km)*math.cos(ice_dir))

            #;how i calculated the new longitude
            templon=original_longitude+math.atan2(math.sin(ice_dir)*math.sin(dis_traveled/rad2km)*math.cos(original_latitude), math.cos(dis_traveled/rad2km)-math.sin(original_latitude)*math.sin(new_lat))

            #;change new lat and lon into degrees
            new_lat=new_lat/(pi/180)
            new_lon=templon/(pi/180)*(-1.0)

            #;these values are created in case the 180 turn is met
            lat_prev=lat
            lon_prev=lon

            #;print, lon, lat, new_lon, new_lat
            
            # This is not being used
            # dist_temp = MAP_2POINTS( lon, lat, new_lon, new_lat,radius=rad2km);*0.000539957
            
            #;print, lon, lat, new_lon, new_Lat, dis_traveled, dist_temp, surface_wind, surface_direction
            #;stop
            #;writing in the new latitude and longitudes into the file, once the PUQ waypoint is reached the original longitude and latitude are input 
            
            #if tsplit(0) eq '180' then goto, jump180
            #if tsplit(1) eq 'PUQ' then goto, jump1
            if tsplit[0]=="180":  # jump180
                print "180 turn"
                count=count+1.0
            elif tsplit[1]=="PUQ" or tsplit[1]=="SAWH": #jump1
                outf.write(tsplit_1st_WP)
            else:
                #if tsplit1(4) eq 6 then printf,3, tsplit(1),' ', new_lat,' ', new_lon, ' ',tsplit1(3), ' ',tsplit1(4),' ', tsplit1(5),' ', tsplit1(6) 
                #else printf,3,tsplit(1),' ', new_lat,' ', new_lon, ' ',tsplit1(3), ' ',tsplit1(4)
                if tsplit1[4] != "0":
                    outf.write(tsplit[1] + " " + str(new_lat) + " " + str(new_lon) + " " + tsplit1[3] + " " + tsplit1[4] + " " + tsplit1[5] + " " + tsplit1[6]+ "\n")
                    allNewLats.append(new_lat)
                    allNewLons.append(new_lon)
                else:
                    outf.write(tsplit[1] + " " + str(new_lat) + " " + str(new_lon) + " " + tsplit1[3] + " " + tsplit1[4] + "\n")  
                    allNewLats.append(new_lat)
                    allNewLons.append(new_lon) 
                count=count+1.0
    
    #;for my own personal testing
    #drift_corrected(0,count)=ice_velocity
    #drift_corrected(1,count)=lat
    #drift_corrected(2,count)=lon
    #drift_corrected(3,count)=new_lat
    #drift_corrected(4,count)=new_lon
    #count=count+1.0
    #Essentially skip writing out 
    #jump180: print, '180 turn'
    #ENDWHILE
    #;this is here to put in the first waypoint PUQ into the file
    # Rewrite the PUQ line to the file
    #jump1: printf,3,tsplit_1st_wp

    outf.close()
    file_latlon.close()
    file_time.close()
    #close,1
    #close,2
    #close,3

    from matplotlib import pyplot as plt
    from mpl_toolkits.basemap import Basemap

    m = Basemap(projection='spstere',boundinglat=-60,lon_0=90,resolution='l')
    #m.drawcoastlines()
    #m.fillcontinents(color='coral',lake_color='aqua')
    # draw parallels and meridians.
    
    x1,y1 = m(allNewLons,allNewLats)
    x2,y2 = m(allOldLons,allOldLats)
    
    m.plot(x1, y1, marker='o',color='m')
    m.plot(x2, y2, marker='o',color='g')
    m.etopo()
    #m.shadedrelief()
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))
    m.drawmapboundary(fill_color='aqua')
    #print str(len(allNewLons))
    #print str(len(allNewLats))
    plt.title("Original and Drift Sequences")
    plt.show()
    
    #  Dunno what this is, not needed
    # ?????
    #drift_corrected(3,61)=-53.0028
    #drift_corrected(4,61)=-70.8547
    # ?????

    # IDL Plotting code, not used
    # map=plot_pretty_v2(drift_corrected(3,*), drift_corrected(4,*),drift_corrected(0,*), minval=0, maxval=1, cbartitle='drift speed', colbar=0)
    #;-----------------------------------------------------------------------------------
    #!P.multi=0
    #!y.thick=40
    #!x.thick=40
    #!p.thick=50
    #!p.charthick=40
    #!p.charsize=25
    #loadct,13
    #;x=[1,2,3,4,5,6,7,8,9,10,11,12]
    #;x=[2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016]
    # The code below is not necessary 
    #set_plot, 'ps'
    #device, file = 'Sensible180.eps', /encap, /color, bits_per_pixel=8, xsize=350, ysize=350
    #plot, drift_corrected(2,0:60),drift_corrected(1,0:60), xrange=[300,350], yrange=[-50,-75];,ave_m,psym=2,min_value=[-500], xrange=[0,20], xtickinterval=[5],yrange=[0,20] ,xtitle='Buoy', ytitle='MERRRA-2'
    #oplot, drift_corrected(4,0:30),drift_corrected(3,0:30), color=250
    #oplot, drift_corrected(4,30:60), drift_corrected(3,30:60), color=150
    #device,/close

    #return,0
    #END

if __name__ == "__main__":
    # Read in the settings file
    try:
        opt1=sys.argv[1]
    except:
        "Please enter a config filename to use.  Exiting.."
        sys.exit()
    print "Reading settings from " + str(opt1)
    pini=open(str(opt1),'r')
    pdat=pini.readlines()
    pini.close()
    psettings={}
    for p in pdat:
      # Add to dictionary
      p=p.strip()
      if (len(p))>0:
        if (p[0] != "#" and p[0]!=" " and p[0]!=""):
          # Not a comment line
          var=p.split(" ")[0]
          val=p.split(" ")[1]
          psettings[var]=val 
  
    take_off_hour = int(psettings['take_off_hour'])
    print "Setting take_off_hour to " + str(take_off_hour)
    take_off_minutes = int(psettings['take_off_minutes'])
    print "Setting take_off_minutes to " + str(take_off_hour)
    IS2_CROSSover_hour = int(psettings['IS2_CROSSover_hour'])
    print "Setting IS2_CROSSover_hour to " + str(IS2_CROSSover_hour)
    IS2_crossover_minutes = int(psettings['IS2_crossover_minutes'])
    print "Setting IS2_crossover_minutes to " + str(IS2_crossover_minutes)
    WIND_SPEED_PLANE = float(psettings['WIND_SPEED_PLANE'])
    print "Setting WIND_SPEED_PLANE to " + str(WIND_SPEED_PLANE)
    WIND_ANGLE_PLANE = float(psettings['WIND_ANGLE_PLANE'])
    print "Setting WIND_ANGLE_PLANE to " + str(WIND_ANGLE_PLANE)
    PLANE_PRESSURE = float(psettings['PLANE_PRESSURE'])
    print "Setting PLANE_PRESSURE to " + str(PLANE_PRESSURE)
    PLANE_LATITUDE = float(psettings['PLANE_LATITUDE'])
    print "Setting PLANE_LATITUDE to " + str(PLANE_LATITUDE)
    PLANE_LONGITUDE = float(psettings['PLANE_LONGITUDE'])
    print "Setting PLANE_LONGITUDE to " + str(PLANE_LONGITUDE)
    sequence_file = str(psettings['sequence_file'])
    print "Setting sequence_file to " + str(sequence_file)
    plan_file = str(psettings['plan_file'])
    print "Setting plan_file to " + str(plan_file)
    output_file = str(psettings['output_file'])
    print "Setting output_file to " + str(output_file)

    icedrift(take_off_hour, take_off_minutes,IS2_CROSSover_hour,IS2_crossover_minutes,WIND_SPEED_PLANE,WIND_ANGLE_PLANE,PLANE_PRESSURE,PLANE_LATITUDE, PLANE_LONGITUDE, sequence_file, plan_file, output_file)
  