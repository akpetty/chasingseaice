""" driftcorrect.py
	
	Code to drift correct OIB underflight of ICESat-2 to account for sea ice drift.
	Takes in the wind measurement time, the plane position (at this time) and the IS-2 cross-over time
	Outputs a single drift correction which can optionally be applied for all rows of a given sequence file
	Still working on doing this iteratively for a long-line

	Current Python version written by Alek Petty (alek.a.petty@nasa.gov)
	Original code written by Linette Boisvert adapted by Jeremy Harbeck. Previous Python adaptation by Robbie Russell.
	
	Input:
		Wind data and correction factors
		Plane information
		
	Output:
		Sea ice drift corrections, including direction and distance travelled
		Updated sequence file
		Optional plot of a map showing both uncorrected/corrected waypoints and their lat/lon differences
	
	Code written and tested in Python 3.6, with python dependencies:
		numpy
		h5py
		matplotlib (for plotting)
		basemap (for plotting)

	Update history:
		30th July 2019: Version 1

	To do:
		Add more elegant way of grabbing nearest Ed and Ron data..
		Add more hemisphere/line-type options 
		Add options for just constant scalings, turning angles
		Need to add a time elapsed fucntion to do this correction on each row of sequence file for the long-line. 
    
    Questions:
		Shouldnt the ice moving to/from depend on if we are arriving before or after IS-2?
	
"""

import numpy as np
import math
import os
import sys 
import h5py

# Uncomment if you can't/don't want to plot
#from mpl_toolkits.basemap import Basemap
#from matplotlib import pyplot as plt

def PlotMapCheck(original_lons, original_lats, corrected_lons, corrected_lats, savefigure=0):
	""" Map original and corrected drifts as a sanity check

	"""
	fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(10, 3))
	ax=axs[0]
	plt.sca(ax)

	m = Basemap(projection='npstere',boundinglat=64,lon_0=-90,resolution='l')
	#m.drawcoastlines()
	#m.fillcontinents(color='coral',lake_color='aqua')
	# draw parallels and meridians.

	x1,y1 = m(corrected_lons, corrected_lats)
	x2,y2 = m(original_lons, original_lats)

	m.plot(x1, y1, marker='x',color='m')
	m.plot(x2, y2, marker='x',color='g')
	#m.etopo()
	m.drawcoastlines(linewidth=0.25, zorder=5)
	m.drawparallels(np.arange(90,-90,-5), linewidth = 0.25, zorder=10)
	m.drawmeridians(np.arange(-180.,180.,30.), latmax=85, linewidth = 0.25, zorder=10)
	m.fillcontinents(color='0.9',lake_color='grey', zorder=3)

	#m.shadedrelief()
	#m.drawparallels(np.arange(-80.,81.,20.))
	#m.drawmeridians(np.arange(-180.,181.,20.))
	#m.drawmapboundary(fill_color='aqua')
	#print str(len(allNewLons))
	#print str(len(allNewLats))
	#plt.title("Original and Drift Sequences")

	ax=axs[1]
	plt.sca(ax)
	plt.plot(np.array(original_lons)-np.array(corrected_lons), label='lon difference')
	#plt.plot(corrected_lons, label='corrected')
	plt.xlabel('waypoint')
	plt.ylabel('longitude difference (deg)')
	plt.tight_layout()

	ax=axs[2]
	plt.sca(ax)
	plt.plot(np.array(original_lats)-np.array(corrected_lats), label='lat difference')
	#plt.plot(corrected_lons, label='corrected')
	plt.xlabel('waypoint')
	plt.ylabel('latitude difference (deg)')
	plt.tight_layout()

	if (savefigure==1):
		plt.savefig('./mapcheck.pdf', dpi=300)

	plt.show()

def PlanePressureCalc(plane_altitudeT):
	"""adjust plane altitude to plane pressure

	returns: plane pressure in mb

	"""
	#convert altitiude in feet to meters 
	aH2 = plane_altitudeT/3.208 

	# assumption of adiabatic lapse rate and standard atmosphere profile
	# State what these constants are and the units!
	aA = -6.5e-3
	aR = 287.0
	aG0 = 9.81
	aT0 = 288.16
	aP0 = 1013.25
	aH0 = 0.0

	# for a given altitude H2, compute pressure P2
	aT2 = aT0 + aA*(aH2-aH0)
	aX = -aG0/(aA*aR)
	plane_pressure = aP0*(aT2/aT0)**aX 

	print('Calculated plane pressure is (mb) :',plane_pressure)

	return plane_pressure
 
def GrabWindData(filePath):
	"""read in values already generated from Ed wind data
	
	returns: slope intercept and angle of winds?

	"""
	# 3d float arrays
	all_slope_xy = np.fromfile(filePath+'all_slope_xy_OIB', '<f4').reshape(111, 51,6)

	all_intercept_xy = np.fromfile(filePath+'all_intercept_xy_OIB', '<f4').reshape(111, 51,6)

	all_angle_xy = np.fromfile(filePath+'all_angle_xy_OIB', '<f8').reshape(111, 51,6)

	return all_slope_xy, all_intercept_xy, all_angle_xy

def CalcWindData(filePath, pressure_levs, hemStr):
	"""read in Ed wind data and calculate slopes etc..
	
	returns: slope intercept and angle of winds?

	"""
	 #pressure levels
	
	
	#for f=0,nlev-1 do begin
	#	file = filePath+'F_funcs_arctic_'+strcompress(ilevs[f],/remove_all)+'mb.h5'

	sz = [107,50] #size of final array at each pressure level	
		

	nlev = len(pressure_levs)
	all_slope_xy = np.zeros((nlev,sz[0], sz[1]))
	all_intercept_xy = np.zeros((nlev,sz[0], sz[1]))
	all_angle_xy = np.zeros((nlev,sz[0], sz[1]))

	for x in range(len(pressure_levs)):
		#print(x)
		file = filePath+'F_funcs_'+hemStr+'_'+str(pressure_levs[x])+'mb.h5'
		h5file = h5py.File(file, 'r')
		h5file
		#print(h5file['F_angle'][x].shape)
		#print(h5file['lat'][0].reshape(sz[0], sz[1]))
		# DON"T SEE WHY YOU RESHAPE ED'S DATA IF JUST FINDING NEAREST INDEX??
		all_angle_xy[x] = h5file['F_angle'][0].reshape(sz[0], sz[1])
		all_slope_xy[x] = h5file['F_slope'][0].reshape(sz[0], sz[1])
		all_intercept_xy[x] = h5file['F_interc'][0].reshape(sz[0], sz[1])

		if (x==0):
			lat = h5file['lat'][0].reshape(sz[0], sz[1])[0, :]
			#print(lat)
			lon = h5file['lon'][0].reshape(sz[0], sz[1])[:, 0]
			#print(lon)


	return all_slope_xy, all_intercept_xy, all_angle_xy, lat, lon

def select_mindiff_level(value, array):
	""" Select the index of the array closest to the input value"""

	array_diff=np.abs(value-array)
	#print(array_diff)
	min_arraydiff = min(array_diff)
	# This returns the INDEX of the closest pressure level to plane 
	sel_index = np.abs(array_diff - min_arraydiff).argmin()	
	return sel_index

def CalcSurfaceWindsDrift(hemStr, plane_altitude, plane_latitude, plane_longitude,
		wind_angle_plane, wind_speed_plane, all_intercept_xy, all_angle_xy, 
		all_slope_xy, lon_reanalysis, lat_reanalysis, pressure_levs, time_diff):
	""" Calculate the windspeed and direction at the surface from height

	Ed data goes from plane to surface wind/direction. 
	We then need to convert this to geostrophic winds
	Then we apply geostrophic to ice drift..

	Notes: equation for wind speed at surface provided by Ed Blanchard-Wrigglesworth
	equation for direction of ice motion provided by Ron Kwok

	Inputs: windspeed, direction and pressure at plane altitude
	Returns: windspeed and direction at surface 
				(aka. velocity and direction of ice)
	"""

	# Convert plane altitude to air pressure
	plane_pressure = PlanePressureCalc(PLANE_ALTITUDE)

	# Convert knots to m/s for windspeed
	wind_speed_plane_ms = wind_speed_plane*0.51444444444

	# ensure the entry longitude is 0->360
	if (plane_longitude<0):
		plane_longitude = plane_longitude+360.0 


	#============== Surface Wind Calculation ======================	

	# Find the closest ERA pressure level
	lev_pressure = select_mindiff_level(plane_pressure, pressure_levs)
	#print('Plane air pressure:', plane_pressure)
	#print('Pressure levels:', pressure_levs)
	#print('Lev_pressure index:', lev_pressure)
	
	# Find the closest ERA latitude index
	lev_lat = select_mindiff_level(plane_latitude, lat_reanalysis)
	#print('Plane latitude:', plane_latitude)
	#print('Reanalysis latitudes:', lat_reanalysis)
	#print('Latitude index:', lev_lat)
	
	# Find the closest ERA longitude index
	lev_lon = select_mindiff_level(plane_longitude, lon_reanalysis)
	#print('Plane longitude:', plane_longitude)
	#print('Reanalysis longitudes:', lon_reanalysis)
	#print('Longitude index:', lev_lon)

	# Surface winds and direction calculations
	# Calculate winds at surface
	#print(all_intercept_xy.shape)
	surface_wind = (wind_speed_plane_ms - all_intercept_xy[lev_pressure,lev_lon, lev_lat]) / all_slope_xy[lev_pressure, lev_lon, lev_lat] 

	#calculate angle of winds at surface
	wind_direction_plane = wind_angle_plane - all_angle_xy[lev_pressure,lev_lon, lev_lat] 


	#============== Kwok Drift turning angle/scaling Calculation ======================	

	# Find the closest longitude point to the ones given by ron

	if (hemStr=='ARCTIC_SPRING'):

		# Geo winds from surface wind conversion
		surface_to_geo_conv = 2.0 # 2.0 for Weddell Sea/Arctic summer, 1.8 for Arctic spring

		lats_kwok = np.array([90.0, 85.0, 80.0, 75.0])
		lons_kwok = np.array([-45.0,-75.0,-105.0,-135.0,-165.0,-195.0])
		lons_kwok360=np.copy(lons_kwok)
		lons_kwok360[np.where(lons_kwok360<0)] = lons_kwok360[np.where(lons_kwok360<0)]+360.

		# Velocity scaling factors (geostrophic winds to ice drift):

		geo_to_ice_drift_scalings=([0.008,0.008,0.008,0.008,0.008,0.008],
			[0.008,0.008,0.007,0.008,0.009,0.010],
			[-99.0,-99.0,0.006,0.008,0.009,0.012],
			[-99.0,-99.0,-99.0,0.009,0.012,0.012]) # 90N, 85N, 80N, 75N 

		# Find the closest Kwok longitude index
		lev_lon_kwok = select_mindiff_level(plane_longitude, lons_kwok360)
		#print('Kwok lon index:', lev_lon_kwok)
		
		# Find the closest Ron latitude index
		lev_lat_kwok = select_mindiff_level(plane_latitude, lats_kwok)
		#print('Kwok lat index:', lev_lat_kwok)

		geo_to_ice_drift_scaling = geo_to_ice_drift_scalings[lev_lat_kwok][lev_lon_kwok]
		

	elif (hemStr=='ARCTIC_SUMMER'):
		surface_to_geo_conv = 2.0 # 2.0 for Weddell Sea/Arctic summer, 1.8 for Arctic spring
		geo_to_ice_drift_scaling=0.01

	print('Surface to geostrophic wind speed scaling factor:', surface_to_geo_conv)
	print('Geostrophic to ice drift scaling factor:', geo_to_ice_drift_scaling)
	if (geo_to_ice_drift_scaling < -10.0):
		# ADD LINE TO ALSO THROW AN EXCEPTION
		print('scaling factor selection issue')

	# Turning angle between surface winds and ice drift (+ve is to the right)
	# Earlier versions used 2 (I believe there was confusion about the surface and geostrophic winds)
	turning_angle = 27.

	ice_drift_scaling = surface_to_geo_conv * geo_to_ice_drift_scaling

	#============== Ice velocity calculation ======================	

	# ice_velocity in km/s
	ice_velocity =  (surface_wind*ice_drift_scaling)/1000.

	if (time_diff>0):
		# Measurement time is before IS-2 X-Over
		# Subtracted 180 from surface_direction as we are calculating where we need to move to correct for this***
		ice_dir_deg = (wind_direction_plane + turning_angle - 180.)
	else:
		ice_dir_deg = wind_direction_plane + turning_angle

	if (ice_dir_deg < 0):
		ice_dir_deg = ice_dir_deg+360.

	print('Ice velocity (cm/s): ',str(ice_velocity*100000.))
	print('Ice velocity (km/h): ',str(ice_velocity*3600.))
	#print('Ice direction - flowing to (deg):', str(ice_dir_deg))

	return ice_velocity, ice_dir_deg


def DriftCorrectlatLonPair(original_lat, original_lon, ice_direction, dist_traveled):
	"""
	Computes ending waypoint lat2/lon2 given starting waypoint
    lat1/lon1, starting true course tc1 and distance dist. Lat/lon
    are in radians N latitude/E longitude, tc1 is in radians,
    and dist is in kilometers.
     
	AUTHOR:   John Sonntag
	DATE:     July 1999

	"""

	# Define constant

	# Radius of the Earth in km
	rad2nm=180.0*60.0/np.pi
	radiusEarthkm=rad2nm*6076.1*12.0*2.54/100.0/1000.0
	# Or just use 6371

	# Calculate new latitude waypoint
	
	# Play around with setting these as constants to assess sensitivity
	#ice_direction=360.
	#dist_traveled=0.1

	#print('Radius of Earth:', radiusEarthkm)
	#print('Distance travelled:', dist_traveled)
	#print('Ice drift direction:', ice_direction)

	lat_rads=np.radians(original_lat)
	lon_rads=np.radians(original_lon)
   
    # Changed it to negative cos...?
	temp_lat=math.asin(math.sin(lat_rads)*math.cos(dist_traveled/radiusEarthkm)-math.cos(lat_rads)*math.sin(dist_traveled/radiusEarthkm)*math.cos(ice_direction))

	temp_lon=lon_rads + math.atan2(math.sin(ice_direction)*math.sin(dist_traveled/radiusEarthkm)*math.cos(original_lat), math.cos(dist_traveled/radiusEarthkm)-math.sin(original_lat)*math.sin(temp_lat))

	corrected_lat=np.degrees(temp_lat)
	corrected_lon=np.degrees(temp_lon)

	print('original_lat (deg):', original_lat)
	print('corrected_lat (deg):', corrected_lat)

	print('original_lon (deg):', original_lon)
	print('corrected_lon (deg):', corrected_lon)
	#print('lat_rads:', lat_rads)
	#print('lon_rads:', lon_rads)

	
	return corrected_lat, corrected_lon


def CalcSequenceDriftCorrection(distance_traveled, ice_direction, FILE_PATH):
	""" Output the drifts onto the sequence file
	
	Inputs: windspeed, direction and pressure at plane altitude
	Returns: windspeed and direction at surface 
				(aka. velocity and direction of ice)
	"""

	# Create file to write corrected data to
	output_file = open(FILE_PATH+'./siis2arcticocean.sequence.driftcorrected',"w")

	# Read in sequence file
	sequence_file= open(FILE_PATH+'./siis2arcticocean.sequence','r').readlines()
	#np.loadtxt('./siis2arcticocean.sequence', dtype=['str', 'f2', 'f2', 'f2', 'f2'])

    # Read in flight plan file
	plan_file= open(FILE_PATH+'./siis2arcticocean.plan','r').readlines()   #siwestweddell.plan'

	# Find index of Mission analysis row
	for pidx in range(len(plan_file)):
		if (plan_file[pidx]=='Mission Analysis\n'):
			break

	# Add those first three dead lines to the output file.
	output_file.write(sequence_file[0])
	output_file.write(sequence_file[1])
	output_file.write(sequence_file[2])
    
	# Declare empty arrays to store data in for later sanity checks
	corrected_lats=[]
	corrected_lons=[]
	original_lats=[]
	original_lons=[]
	
	for w_idx in range(3, len(sequence_file)):
		
		#waypoint_row_info = plan_file[w_idx].split()
		#elapsed_time = plan_file[pidx+6].split()[8]

		original_lat = float(sequence_file[w_idx].split()[1])
		original_lon = float(sequence_file[w_idx].split()[2])
		#print('original_lat:', original_lat)
		#print('original_lon:', original_lon)
		#print('Distance reference point traveled (km): ',distance_traveled)

		#FOR THE LONG LINE WE WANT TO ADD A CONDITION HERE THAT CALCULATES EXTRA TIME AND DISTANCE TRAVELLED
		
		corrected_lat, corrected_lon = DriftCorrectlatLonPair(original_lat, original_lon, ice_direction, distance_traveled)

		if (plan_file[pidx+w_idx+3].split()[0]=="180"):  
			print("180 turn")

		elif (sequence_file[w_idx].split()[0]=="BGTL"):
			output_file.write(sequence_file[w_idx])

		elif(sequence_file[w_idx].split()[4]!="0"):
			print("0 FLAG??")
			# WHAT IS THIS IF STATEMENT ABOUT  
			output_file.write(sequence_file[w_idx].split()[0] + "\t" + str(corrected_lat) + " " + str(corrected_lon) + " " + \
				sequence_file[w_idx].split()[3] + " " + sequence_file[w_idx].split()[4] + " " + sequence_file[w_idx].split()[5] \
				+ " " + sequence_file[w_idx].split()[6]+ "\n")
			corrected_lats.append(corrected_lat)
			corrected_lons.append(corrected_lon)
			original_lats.append(original_lat)
			original_lons.append(original_lon)
		
		else:
		    output_file.write(sequence_file[w_idx].split()[0] + "\t" + str(corrected_lat) + " " + str(corrected_lon) + " " + sequence_file[w_idx].split()[3] + " " + sequence_file[w_idx].split()[4] + "\n")  
		    corrected_lats.append(corrected_lat)
		    corrected_lons.append(corrected_lon) 
		    original_lats.append(original_lat)
		    original_lons.append(original_lon)

	return original_lons, original_lats, corrected_lats, corrected_lons

def IS2timeDiff(MEASUREMENT_HOUR, MEASUREMENT_MIN, IS2_CROSSOVER_MIN, IS2_CROSSOVER_HOUR):
	
	# Wind measurement time in seconds
    measurement_time = (MEASUREMENT_HOUR*60.0*60.0)+(MEASUREMENT_MIN*60.0)
    #print('OIB take-off time (seconds):', measurement_time)
    
    # IS2 x-over time in seconds
    is2_crossover_time = (IS2_CROSSOVER_HOUR*60.0*60.0) + (IS2_CROSSOVER_MIN*60.0)
    #print('IS2 X-over time (seconds):', is2_crossover_time)

    # Compute time difference between takeoff and IS2 crossover in seconds
    time_diff_to_IS2=is2_crossover_time-measurement_time
    print('OIB takeoff and IS-2 X-Over difference (hours):', time_diff_to_IS2/(60.*60.))
    
    #ADD ELAPSED TIME HERE IF WE DO LONG_LINE CORRECTION

    #print('Time diff to IS2 (seconds)', time_diff_to_IS2)
    return time_diff_to_IS2



def main(FILE_PATH, FLIGHT_PATH, HEM_STR, MEASUREMENT_HOUR, MEASUREMENT_MIN, IS2_CROSSOVER_HOUR, \
	IS2_CROSSOVER_MIN, WIND_SPEED_PLANE, WIND_ANGLE_PLANE, PLANE_ALTITUDE, \
	PLANE_LATITUDE, PLANE_LONGITUDE, PRESSURE_LEVS, OUT_SEQUENCE):
	

	#============== 2. Wind data? ======================	

	# Either calculate from raw Ed data or just grab previously generated files..
	all_slope_xy, all_intercept_xy, all_angle_xy, lat_reanalysis, lon_reanalysis = CalcWindData(FILE_PATH, PRESSURE_LEVS, HEM_STR)

	#all_slope_xy, all_intercept_xy, all_angle_xy = GrabWindData(FILE_PATH)

	#============== 3. Surface drift ======================

	time_diff_to_IS2 = IS2timeDiff(MEASUREMENT_HOUR, MEASUREMENT_MIN, IS2_CROSSOVER_MIN, IS2_CROSSOVER_HOUR)
	
	ice_velocity, ice_dir_deg = CalcSurfaceWindsDrift(HEM_STR, PLANE_ALTITUDE, PLANE_LATITUDE, 
		PLANE_LONGITUDE, WIND_ANGLE_PLANE,WIND_SPEED_PLANE, all_intercept_xy, 
		all_angle_xy, all_slope_xy, lon_reanalysis, lat_reanalysis, PRESSURE_LEVS, 
		time_diff_to_IS2)

	distance_traveled = abs(time_diff_to_IS2)*ice_velocity

	
	print('Correction distance (km):', distance_traveled)
	print('Correction bearing (deg N):', ice_dir_deg)

	if OUT_SEQUENCE:
		#============== 4. Waypoint corrections ======================

		original_lons, original_lats, corrected_lats, corrected_lons = CalcSequenceDriftCorrection(distance_traveled, ice_dir_deg, FLIGHT_PATH)

		print('Correction distance (km):', distance_traveled)
		print('Correction bearing (deg N):', ice_dir_deg)
		#============== 5. Plotting ======================
		
		# Drop the last row for plotting!
		#PlotMapCheck(original_lons[0:-1], original_lats[0:-1], corrected_lons[0:-1], corrected_lats[0:-1], savefigure=0)


if __name__ == "__main__":

	MEASUREMENT_HOUR = 13 # Hour (Z) ;valid time for correction (time the winds are measured)
	MEASUREMENT_MIN = 0 # Minutes
	IS2_CROSSOVER_HOUR = 10 # Hour (Z) = Time at center waypoint on IS-2 track
	IS2_CROSSOVER_MIN = 0 #Minutes

	WIND_SPEED_PLANE = 10 # Knots (converted to m/s in code below)
	WIND_ANGLE_PLANE = 270 # Degrees (0 -> 360)
	PLANE_ALTITUDE = 0.0 # Feet, use to calc plane pressure below
	PLANE_LATITUDE = 82.813 # Degrees (-90 -> +90)
	PLANE_LONGITUDE = -119.840 # Degrees
	FILE_PATH='./data/'
	FLIGHT_PATH='./20190419/' # If we want to use directories to store data..Kinda pointless, maybe drop
	PRESSURE_LEVS = np.array([925,950,975,1000]) # pressure levels used to produce reanalysis wind scalings
	#PRESSURE_LEVS = np.array([850, 875, 900, 925,950,975,1000]) # pressure levels used to produce reanalysis wind scalings
	HEM_STR='ARCTIC_SPRING' #'ARCTIC_SPRING', 'ARCTIC_SUMMER', 'ANTARCTIC'
	OUT_SEQUENCE=True # Output a sequence file


	main(FILE_PATH,FLIGHT_PATH, HEM_STR, MEASUREMENT_HOUR, MEASUREMENT_MIN, IS2_CROSSOVER_HOUR, 
		IS2_CROSSOVER_MIN, WIND_SPEED_PLANE, WIND_ANGLE_PLANE, PLANE_ALTITUDE,
		PLANE_LATITUDE, PLANE_LONGITUDE, PRESSURE_LEVS, OUT_SEQUENCE)








