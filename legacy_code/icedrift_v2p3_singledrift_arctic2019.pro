PRO icedrift_v2p3_singledrift_arctic2019
  ;version 2p3, fixed error in sel_ronlat

  sdate = '20190401'
  
  sl = path_sep()
  ;dir = 'Y:\jharbeck\Working\drift\'
  ;dir = '/data/users/lnboisvert/data/Default/'
  dir = 'D:'+sl+'users'+sl+'jharbeck'+sl+'Working'+sl+'Drift_code'+sl
  flight_dir = dir+sdate+sl
  verbose = 1

  ;enter values here for this run
  take_off_hour = 13 ;hour (Z) ;valid time for correction (time the winds are measured)
  take_off_minutes = 10 ;minutes
  IS2_CROSSover_hour = 13 ;hour (Z) = Time at center waypoint on IS-2 track
  IS2_crossover_minutes = 03 ;minutes
  WIND_SPEED_PLANE = 4 ;knots (converted to m/s in code below)
  WIND_ANGLE_PLANE = 298 ;degrees (0 -> 360)
  ;PLANE_PRESSURE = 940.0 ;mb
  PLANE_ALTITUDE = 3500.0 ;feet, use to calc plane pressure below
  PLANE_LATITUDE = 85.813 ;degrees (-90 -> +90)
  PLANE_LONGITUDE = -119.84 ;degrees
  ;sequence_file = dir+'siis2arcticocean.sequence'
  ;plan_file = dir+'siis2arcticocean.plan'
  ;output_file = dir+file_basename(sequence_file,'.sequence')+'_drift_applied.sequence'
  ;output_image = dir+file_basename(sequence_file,'.sequence')+'_drift_applied.sequence.png'

  ;test for files
  ;if file_test(sequence_file) eq 0 then STOP,'cannot find sequence file'
  ;if file_test(plan_file) eq 0 then STOP,'cannot find plan file'

  drift_corrected=fltarr(5,100)
  ;=======================================================================

  ;-----------------------------------------------------------------------
  ;              adjust plane altitude to plane pressure
  ;-----------------------------------------------------------------------
  plane_altitude = plane_altitude/3.208 ;convert feet to meters 
  
  ;*assumption of adiabatic lapse rate and standard atmosphere profile
  aA = -6.5e-3
  aR = 287.0
  aG0 = 9.81
  aT0 = 288.16
  aP0 = 1013.25
  aH0 = 0.0
  aH2 = PLANE_ALTITUDE
  
  ;given altitude H2, compute pressure P2
  aT2 = aT0 + aA*(aH2-aH0)
  aX = -aG0/(aA*aR)
  aP2 = aP0*(aT2/aT0)^aX 
  print,'Calculated plane pressure is (mb) :',aP2
  PLANE_PRESSURE = aP2
  ;=======================================================================

  ;-----------------------------------------------------------------------
  ;                        Read in Ed's wind data
  ;-----------------------------------------------------------------------
  ;
  
  ;read in lat, lon from ERA fields and the all of the slopes, intercepts and angles for converting
  ; wind measurements taken from the plane to surface winds and directions

 

  ;              ****CONVERTED DIRECTLY FROM ED'S DATA****
  ilevs = [925,950,975,1000] ;pressure levels
  sz = [50,107] ;size of final array at each pressure level

  nlev = n_elements(ilevs)
  all_slope_xy = fltarr(nlev,sz[0],sz[1])
  all_intercept_xy = all_slope_xy
  all_angle_xy = all_slope_xy

  for f=0,nlev-1 do begin
    file = 'F_funcs_arctic_'+strcompress(ilevs[f],/remove_all)+'mb.h5'

    file_id = H5F_OPEN(file)

    ID = H5D_OPEN(file_id,'lat')
    lat = H5D_READ(ID)
    H5D_CLOSE,ID

    ID = H5D_OPEN(file_id,'lon')
    lon = H5D_READ(ID)
    H5D_CLOSE,ID

    ID = H5D_OPEN(file_id,'F_angle')
    F_angle = H5D_READ(ID)
    H5D_CLOSE,ID

    ID = H5D_OPEN(file_id,'F_slope')
    F_slope = H5D_READ(ID)
    H5D_CLOSE,ID

    ID = H5D_OPEN(file_id,'F_interc')
    F_interc = H5D_READ(ID)
    H5D_CLOSE,ID

    H5F_CLOSE,FILE_ID

    if f eq 0 then begin
      ERA_lat = (reform(lat,sz[0],sz[1]))[*,0]
      ERA_lon = (reform(lon,sz[0],sz[1]))[0,*] ;from Ed is 0->360
    endif
    all_slope_xy[f,*,*] = reform(F_slope,sz[0],sz[1])
    all_intercept_xy[f,*,*] = reform(F_interc,sz[0],sz[1])
    all_angle_xy[f,*,*] = reform(F_angle,sz[0],sz[1])

  endfor
  ;=======================================================================

  ;-----------------------------------------------------------------------
  ;STEP1: Calculate the windspeed and direction at the surface from height
  ;
  ;equation for wind speed at surface provided by Ed Blanchard-Wrigglesworth
  ;equation for direction of ice motion provided by Ron Kwok
  ;
  ;inputs: windspeed, direction and pressure
  ;output: windspeed and direction at surface (aka. velocity and direction of ice)
  ;-----------------------------------------------------------------------

  ;find the closest ERA pressure level, latitude and longitude to those of the aircraft
  pres_diff = abs(plane_pressure - ilevs)
  ;lev_min = where(pres_diff eq min(pres_diff))
  temp = min(pres_diff,lev_min)

  lat_n = abs(plane_latitude - era_lat)
  ;near_lat = where(lat_n eq min(lat_n))
  temp = min(lat_n,near_lat)

  if plane_longitude lt 0 then plane_longitude = plane_longitude+360.0 ;ensure the entry longitude is 0->360
  lon_n = abs(plane_longitude - era_lon)
  ;near_lon = where(lon_n eq min(lon_n))
  temp = min(lon_n,near_lon)
  
  ;convert knots to m/s for windspeed
  wind_speed_plane = wind_speed_plane*0.51444444444

  ;surface winds and direction calculations
  surface_wind = (wind_speed_plane - all_intercept_xy[lev_min,near_lat,near_lon]) / all_slope_xy[lev_min, near_lat,near_lon] ;calculate winds at surface
  surface_direction = wind_angle_plane - all_angle_xy[lev_min,near_lat,near_lon] ; calculate angle of winds at surface

  ;find the closest longitude point to the ones given by ron
  ronlon = [-45.0,-75.0,-105.0,-135.0,-165.0,-195.0]
  ronlon360 = ronlon+360.0 ;move it to 0->360, same as adjusted longitude above
  ronlat = [90.0,85.0,80.0,75.0]

  sel_londiff = min(abs(plane_longitude - ronlon360),sel_ronlon)
  sel_latdiff = min(abs(plane_longitude - ronlat),sel_ronlat)

  ;velocity scaling factors (turn winds into ice velocity):
  geo2surf_conv = 1.8 ;2.0 for Weddell Sea, 1.8 for Arctic
  scaling_factor=[[0.008,0.008,0.008,0.008,0.008,0.008],$  ;90N
    [0.008,0.008,0.007,0.008,0.009,0.010],$  ;85N
    [-99.0,-99.0,0.006,0.008,0.009,0.012],$  ;80N
    [-99.0,-99.0,-99.0,0.009,0.012,0.012]]   ;75N
  scaling = scaling_factor[sel_ronlon,sel_ronlat]
  if scaling lt -10.0 then STOP,'scaling factor selection issue'
  scaling = scaling*geo2surf_conv

  ;this will be computed via an equation from ron kwok
  ice_velocity =  (surface_wind[0]/1000.0)*scaling[0];km/s

  ;wind direction from wind angle taken from the plane into radians
  turning_factor = 2.0
  ;***subtracted 180 from surface_direction as we are calculating where it's moving to, not from***
  ice_dir = ((surface_direction[0]-180.0 + turning_factor))*(!pi/180)
  ice_dir_deg = ice_dir*180.0/ !pi
  if ice_dir_deg lt 0 then ice_dir_deg = ice_dir_deg+360.0
  print,'Ice velocity (km/h): ',strcompress(ice_velocity*3600.0),', Ice direction flowing to (deg):', strcompress(ice_dir_deg)

  ;=======================================================================

  ;-----------------------------------------------------------------------
  ;              open sequence, plan files and output files
  ;-----------------------------------------------------------------------
  ;open sequence file
;  openr,seqlun,sequence_file,/get_lun
;  seqtemp = ''
;
;  ;open plan file
;  openr,planlun,plan_file,/get_lun
;  plantemp = ''
;
;  ;create/open the sequence file that will be output to John Sonntag
;  openw,outlun,output_file,/get_lun
  ;=======================================================================

  ;-----------------------------------------------------------------------
  ;                    Work with input/output headers
  ;-----------------------------------------------------------------------
  ;read in the first 3 lines from the sequence file that have PUQ (starting point) and
  ;record them in the output sequence file
;  for n=0,2 do begin & readf,seqlun,seqtemp & printf,outlun,seqtemp & endfor
;  seqsplit_1st_WP = strsplit(seqtemp,' ',/extract)
;
;  ;move ahead in the plan file to the data we are looking for
;  while ~EOF(planlun) do begin
;    readf,planlun,plantemp
;    if strcmp(strmid(plantemp,0,7),'mission',/fold_case) eq 1 and strlen(plantemp) lt 20 then begin
;      for hdr=0,4 do readf,planlun,plantemp ;pass through Mission Analysis header and starting point entry line
;      break
;    endif
;  endwhile
  ;=======================================================================

  ;-----------------------------------------------------------------------
  ;  go through each line of the sequence file (starting with line 4 of
  ;   the sequence file and line one of the plan file) and adjust the
  ;  waypoints according to drift, until the end of the file is reached
  ;-----------------------------------------------------------------------
;
;  count=0L
;  WHILE ~EOF(seqlun) DO BEGIN
;
;    readf,seqlun,seqtemp
;    readf,planlun,plantemp
;
;    if verbose gt 0 then begin
;      print,'Sequence line: ',seqtemp
;      print,'Plan line: ',plantemp
;    endif
;
;    plansplit = strsplit(plantemp,/extract)
;    seqsplit = strsplit(seqtemp,/extract)
;
;    ;skip any line with a turn in it and move on to the next line
;    if plansplit[0] eq '180' then begin
;      readf,planlun,plantemp
;      plansplit = strsplit(plantemp,/extract)
;    endif
;
;    ;ensure we are comparing the same two lines
;    if strcmp(plansplit[1],seqsplit[0],/fold_case) eq 0 then STOP,'sequence and plan input lines do not agree'
;
;    ;assign data to variables from data string
;    elapsed_time = double(plansplit[8])
;    lat = double(seqsplit(1))
;    lon = double(seqsplit(2))


    ;---------------------
    ; STEP2: Calculate time difference between IceSat-2 crossing and when OIB reaches that point
    ;
    ; this will come from John Sonntag
    ;
    ; input: IS2 time, oib take off time from punta
    ; output: total time difference
    ;----------------------

    ;icebridge time in elapsed seconds:
    take_off_time = take_off_hour*3600.0+take_off_minutes*60.0

    ;is2 time in elapsed seconds:
    is2_crossover_time = is2_crossover_hour*3600.0+is2_crossover_minutes*60.0

    ;compute time difference between takeoff and is2 crossover
    IS2_time = IS2_crossover_time - take_off_time

    ;convert time to seconds
    ;time_diff_to_IS2 = (elapsed_time*60.0*60.0)-IS2_time
    time_diff_to_IS2= IS2_time ; this is used for single point adjustments
    ;time_diff_to_IS2=(IS2_time+elapsed_time*60.*60.)

    ;if verbose gt 0 then begin
      print,'Time difference (hours) between racetrack reference point and IS2 flyover: ',time_diff_to_is2/3600.0
    ;  print,'----------'
    ;endif

    ;---------------------
    ; STEP3: Using time difference from step2 and velocty and direction of sea ice drift from step1 calculate the total distance moved
    ;
    ; input: time, drift speed
    ; output: sea ice travel distance
    ;---------------------

    ;calculate the distance traveled based on ice drift and time: time must be seconds, and velocity must be in km/s
    dis_traveled = time_diff_to_IS2*ice_velocity
    print,'Distance reference point traveled (km): ',dis_traveled
    print,'----------'
    STOP
    ;--------------------
    ; STEP4: Use John's code with lat/lon and angle to get new location of way point
    ;
    ; input: original lat long, distance, angle
    ; output: new lat lon
    ;--------------------
    ;
    ;PURPOSE:  Computes ending waypoint lat2/lon2 given starting waypoint
    ;lat1/lon1, starting true course tc1 and distance dist. Lat/lon
    ;are in radians N latitude/E longitude, tc1 is in radians,
    ;and dist is in kilometers.
    ;
    ;#define PI (4.0*atan((double)(1.0)))
    ;#define DEG2RAD (PI/180.0)
    ;#define RAD2NM  (180.0*60.0/PI)
    ;#define RAD2KM  (RAD2NM*6076.1*12.0*2.54/100.0/1000.0)
    ;
    ;AUTHOR:   John Gary Sonntag
    ;
    ;DATE:     July 1999
    ;*------------------------------------------------------------------------*/
;
;    ; convert original latitude and longitude to radians
;    original_latitude = lat*(!pi/180.0)
;    original_longitude = lon*(!pi/180.0)
;
;    ;original_longitude = -1.0*original_longitude
;
;    ;produce constants
;    rad2nm = (180.0 * 60.0 / !pi)
;    rad2km = (rad2nm * 6076.1 * 12.0 * 2.54 / 100.0 / 1000.0)
;
;    ;***** subtract 360 degrees from the ice_dir to get which way the ice is moving towards... I think
;    ;
;    ;calculate new latitude waypoint
;    ;changed it to negative cos...
;    new_lat = asin(sin(original_latitude)*cos(dis_traveled/Rad2km) - cos(original_latitude)*sin(dis_traveled/Rad2km)*cos(ice_dir))
;
;    ;how i calculated the new longitude
;    templon = original_longitude + atan(sin(ice_dir)*sin(dis_traveled/Rad2km)*cos(original_latitude),cos(dis_traveled/Rad2km)-sin(original_latitude)*sin(new_lat))
;
;    ;change new lat and lon into degrees
;    new_lat = new_lat/(!pi/180.0)
;    new_lon = templon/(!pi/180.0);*(-1.0)
;
;    ;these values are created in case the 180 turn is met
;    lat_prev = lat
;    lon_prev = lon
;
;
;    ;dist_temp = MAP_2POINTS(lon, lat, new_lon, new_lat,radius=rad2km);*0.000539957
;    dist_temp = MAP_2POINTS(lon, lat, new_lon, new_lat, radius = rad2km*1000.0)/1000.0 ;track adjustment in kilometers
;    az_temp = (MAP_2POINTS(lon, lat, new_lon, new_lat))[1]
;    if az_temp lt 0 then az_temp = az_temp+360.0
;    az_temp = long(az_temp)
;
;    if verbose ge 1 then begin
;      print,'                 Lon:',strcompress(lon),',     Lat:',strcompress(lat)
;      print,'             Adj Lon:',strcompress(new_lon),', Adj Lat:',strcompress(new_lat)
;      print,' Calc Drift distance:',strcompress(dis_traveled)
;      print,'Dist old->new points:',strcompress(dist_temp)
;      print,'Point mvmt direction:',strcompress(az_temp)
;      print,'  Surface wind speed:',strcompress(surface_wind)
;      print,' Srfc Wind Direction:',strcompress(surface_direction)
;      bl=0
;    endif
;
;    ;writing in the new latitude and longitudes into the file, once the PUQ waypoint is reached the original longitude and latitude are input
;    if plansplit[1] eq 'BGTL' then goto, jump1
;    if seqsplit[4] eq 6 then begin
;      fmt = '(a6,a3,f12.8,a1,f12.8,a1,a3,a1,a1,a1,a3,a1,a4)'
;      printf,outlun,plansplit[1],' ', new_lat,' ',new_lon,' ',seqsplit[3], ' ',seqsplit[4],' ', seqsplit[5],' ', seqsplit[6],format=fmt
;    endif else begin
;      fmt = '(a6,a3,f12.8,a1,f12.8,a1,a3,a1,a1,a1)'
;      printf,outlun,plansplit[1],' ', new_lat,' ',new_lon,' ',seqsplit[3], ' ',seqsplit[4],format=fmt
;    endelse
;
;    ;for my own personal testing
;    drift_corrected[0,count] = ice_velocity
;    drift_corrected[1,count] = lat
;    drift_corrected[2,count] = lon
;    drift_corrected[3,count] = new_lat
;    drift_corrected[4,count] = new_lon
;
;    count++
;
;  ENDWHILE
;
;  ;this is here to put in the first waypoint BGTL into the file
;  jump1: printf,outlun,seqsplit_1st_wp
;
;  free_lun,seqlun
;  free_lun,planlun
;  free_lun,outlun
;  
;  
;  
;  ;=======================================================================
;
;  ;-----------------------------------------------------------------------
;  ;              check for file existence and length match
;  ;-----------------------------------------------------------------------
;  if file_test(output_file) eq 0 then STOP,'No output file'
;  result = query_ascii(output_file,out_info)
;  if result eq 0 then STOP,'Cannot query new sequence file'
;  result = query_ascii(sequence_file,orig_info)
;  if result eq 0 then STOP,'Cannot query original sequence file'
;  if out_info.lines eq orig_info.lines then begin
;    print,'New sequence file matches original length'
;  endif else STOP,'Issue with new sequence file length'
;
;  ;=======================================================================
;  
;  
;  ;-----------------------------------------------------------------------
;  ;Print out vector values for John Sonntag to use in his code for
;  ;the race track IS-2 flights
;  print, '  '
;  print, '------------------------------------------------------------'
;  print, 'Ice Direction    ', Ice_dir
;  print, 'Distance Traveled    ', dis_traveled
;  ;-----------------------------------------------------------------------
;  ;
;  ;
;
;  ;-----------------------------------------------------------------------
;  ;                        sanity plotting section
;  ;-----------------------------------------------------------------------
;  drift_corrected = drift_corrected[*,0:count-1]
;
;  STOP,'continue to plot graphs
;
;  ;new course
;  map1 = plot_pretty_v2p2(drift_corrected(3,*), drift_corrected(4,*),indgen(count),minval=0,maxval=count-1,title='Track Changes from Drift | Time:Purple->Red | O:New,Square:Orig)',colbar=0,palette=34,ssize=0.5);,wdim=[1200,800])
;  map1.plot.linestyle='-'
;  ;orig course
;  map2 = plot_pretty_v2p2(drift_corrected(1,*), drift_corrected(2,*),indgen(count),minval=0,maxval=count-1,overplot=1,ssize=0.8,palette=34)
;  map2.plot.symbol='square'
;  map2.plot.sym_filled=0
;  map2.plot.linestyle='--'
;  map2.plot.sym_thick=2.0
;
;  STOP
;  map1.plot.save,output_image
;  STOP
  ;=======================================================================
end