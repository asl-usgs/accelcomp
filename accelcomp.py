#!/usr/bin/env python

####################################
#Code: accelcomp.py
#Adam Ringler
#
#This code compares broadband data to accelerometers by deconvolving all of the data to displacement
####################################



###################################
#Here are the different methods
#getorientation()
#getdip()
#rotatehorizontal()
#choptocommon()
#getlatlon()
#getstalist()
#readcmt()
#getdata()
#getPAZ2()
#getcolor()
#writestats()
###################################

import sys
import os
import glob
import numpy
import matplotlib
import math
import warnings
from matplotlib.pyplot import (figure,axes,plot,xlabel,ylabel,title,subplot,legend,savefig,show,xlim)
from obspy import read, Stream
from obspy.core import UTCDateTime
from obspy.xseed import Parser
from time import gmtime, strftime
from obspy.core.util.geodetics import gps2DistAzimuth
from obspy.signal.cross_correlation import xcorr
from obspy.taup.taup import getTravelTimes

#Here are the various parameters a user might change
#We might want to put these in a config file to avoid them sitting in the code

debug=False
datalessloc = '/APPS/metadata/SEED/'
#Here is the data location use True for xs0 otherwise use false
dataloc = False
userminfre = .05
usermaxfre = .25
lents = 2000
#Use half the value you think you want e.g. 2 gives you a total of 4 poles
filtercornerpoles = 4
#Here is the number of seconds before the P-wave arrival
bfarrival = 120
#Here is the number of seconds after the S-wave arrival
afarrival = 600

manstalist=False
stations=['IU SAML']

def getorientation(net,sta,loc,chan,evetime,xseedval):
#A function to get the orientation of a station at a specific time
#We use the net, sta, loc, and chan to parse the dataless
#we use the eventime time to get the correct metadata

	for cursta in xseedval.stations:
#As we scan through blockettes we need to find blockettes 50 and 52
		for blkt in cursta:
			if blkt.id == 50:
#Pull the station info for blockette 50
				stacall = blkt.station_call_letters.strip()

			if stacall == sta:
				if blkt.id == 52 and blkt.location_identifier == loc and blkt.channel_identifier == chan:
#Okay we are in blockette 52 and we have the correct location and channel
					if type(blkt.end_date) is str:
						curdoy = strftime("%j",gmtime())
						curyear = strftime("%Y",gmtime())
						curtime = UTCDateTime(curyear + "-" + curdoy + "T00:00:00.0") 
						if blkt.start_date <= evetime:
							azimuth = blkt.azimuth
					elif blkt.start_date <= evetime and blkt.end_date >= evetime:
						azimuth = blkt.azimuth
	return azimuth

def getdip(net,sta,loc,chan,evetime,xseedval):
#A function to get the dip of a station at a specific time
#We use net, sta, loc, and chan to isolate the correct blockette
#We use the event time to get the correct epoch
	for cursta in xseedval.stations:
#As we scan through blockettes we need to find blockettes 50 and 52
		for blkt in cursta:
			if blkt.id == 50:
#Pull the station info for blockette 50
				stacall = blkt.station_call_letters.strip()
			if stacall == sta:
				if blkt.id == 52 and blkt.location_identifier == loc and blkt.channel_identifier == chan:
					if type(blkt.end_date) is str:
						curdoy = strftime("%j",gmtime())
						curyear = strftime("%Y",gmtime())
						curtime = UTCDateTime(curyear + "-" + curdoy + "T00:00:00.0") 
						if blkt.start_date <= evetime:
							dip = blkt.dip
					elif blkt.start_date <= evetime and blkt.end_date >= evetime:
						dip = blkt.dip
	return dip

def rotatehorizontal(stream, angle1,angle2):
#This function rotates the stream using the two angles
#This is something currently being debugged with Tyler
	debugrothor = False
	if stream[0].stats.channel in set(['LHE','LHN','LNE','LNN']):
		stream.sort(['channel'],reverse=False)
#We need to check if this is necessary or not and when we should rotate or not
#		stream[1].data = - stream[1].data

#Function to rotate the data by a given angle
        theta_r1 = math.radians(angle1)
	theta_r2 = math.radians(angle2 - 90) 
	if debugrothor:
		print 'Rotating ' + stream[0].stats.channel + ' by: ' + str(angle1) 
		print 'Rotating ' + stream[1].stats.channel + ' by: ' + str(angle2)

# create new trace objects with same info as previous
        rotatedN = stream[0].copy()
        rotatedE = stream[1].copy()

# assign rotated data
        rotatedN.data = stream[0].data*math.cos(theta_r1) - stream[1].data*math.sin(theta_r1)
        rotatedE.data = stream[1].data*math.cos(theta_r2) + stream[0].data*math.sin(theta_r2)
	if debugrothor:
		print 'Rotation matrix: ' + str(math.cos(theta_r1)) + ' -' + str(math.sin(theta_r1))
		print 'Rotation matrix: ' + str(math.cos(theta_r2)) + ' ' + str(math.sin(theta_r2))

#Lets rename the channels to N,E instead of 1 and 2
	if rotatedN.stats.location == '20':
		rotatedN.stats.channel='LNN'
		rotatedE.stats.channel='LNE'
	else:
		rotatedN.stats.channel='LHN'
		rotatedE.stats.channel='LHE'

# return new streams object with rotated traces
        streamsR = Stream(traces = [rotatedN, rotatedE])
        return streamsR

def choptocommon(stream):
#A function to chop the data to a common time window
	debugchoptocommon = False
	stimes = []
	etimes = []

#Lets get the start and end time for each trace in the stream
	for trace in stream:
		stimes.append(trace.stats.starttime)
		etimes.append(trace.stats.endtime)
	newstime = stimes[0]
	newetime = etimes[0]

#Lets now find the latest start time
	for curstime in stimes:
		if debug:
			print(curstime)
		if curstime >= newstime:
			newstime = curstime

#Lets find the earliest end time	
	for curetime in etimes:
		if debug:		
			print(curetime)
		if curetime <= newetime:
			newetime = curetime

	if debugchoptocommon:
		print(newstime)
		print(newetime)
		print(stream)

#Now we trim each trace by our latest start time and earliest end time
	for trace in stream:	
		trace.trim(starttime=newstime,endtime=newetime)
	if debug:
		print(stream)
	return stream

def getlatlon(sta,etime,xseedval):
#A function to get the lat and lon of a station at a given time
#This function uses lat and lon for blockette 50

	for cursta in xseedval.stations:
#As we scan through blockettes we need to find blockettes 50
		for blkt in cursta:
			if blkt.id == 50:
#Pull the station info for blockette 50
				stacall = blkt.station_call_letters.strip()

#Lets check if we have the correct station
				if stacall == sta:
					lat = blkt.latitude
					lon = blkt.longitude	

#Now lets check if the epoch is the correct one for the time we have given
					if type(blkt.end_effective_date) is str:
						curdoy = strftime("%j",gmtime())
						curyear = strftime("%Y",gmtime())
						curtime = UTCDateTime(curyear + "-" + curdoy + "T00:00:00.0") 
						if blkt.start_effective_date <= etime:
							lat = blkt.latitude
							lon = blkt.longitude
					elif blkt.start_effective_date <= etime and blkt.end_effective_date >= etime:
						lat = blkt.latitude
						lon = blkt.longitude	
	return lat,lon

def getstalist(sp,etime,curnet):
#A function to get a station list
	stations = []
	for cursta in sp.stations:
#As we scan through blockettes we need to find blockettes 50 
		for blkt in cursta:
			if blkt.id == 50:
#Pull the station info for blockette 50
				stacall = blkt.station_call_letters.strip()
				if debug:
					print "Here is a station in the dataless" + stacall
				if type(blkt.end_effective_date) is str:
					curdoy = strftime("%j",gmtime())
					curyear = strftime("%Y",gmtime())
					curtime = UTCDateTime(curyear + "-" + curdoy + "T00:00:00.0") 
					
					if blkt.start_effective_date <= etime:
						stations.append(curnet + ' ' + blkt.station_call_letters.strip())
				elif blkt.start_effective_date <= etime and blkt.end_effective_date >= etime:
					stations.append(curnet + ' ' + \
					blkt.station_call_letters.strip())	
	return stations

def readcmt(cmt):
#This function reads the cmt and gets the various important event information from it
	debugreadcmt = False
#Now we can continue like there is no difference between Princeton and our Synthetics
#Lets get the event time from the cmt
	cmtline1 = ' '.join(cmt[0].split())

#Now we want the lat, lon, time shift, half-duration, and depth
	cmtlat = cmt[4].replace('latitude:','').strip()
	cmtlon = cmt[5].replace('longitude:','').strip()
	tshift = float(cmt[2].replace('time shift:','').strip())
	hdur = float(cmt[3].replace('half duration:','').strip())
	dep = float(cmt[6].replace('depth:','').strip())

#Here are some debug statements to make sure we are parsing correctly
	if debugreadcmt:
		print cmtline1
	cmtline1 = cmtline1.split()
	if debugreadcmt:
		print cmtline1[1] + ' ' + cmtline1[2] + ' ' + cmtline1[3] + ' ' + cmtline1[4] + ' ' + cmtline1[5] + ' ' + cmtline1[6]
	eventtime = UTCDateTime(int(cmtline1[1]),int(cmtline1[2]),int(cmtline1[3]),int(cmtline1[4]),int(cmtline1[5]),float(cmtline1[6]))
	if debugreadcmt:
		print 'Year:' + str(eventtime.year)
		print 'Day:' + str(eventtime.julday)
		print 'Hour:' + str(eventtime.hour)
		print 'Minute:' + str(eventtime.minute)
	return cmtlat, cmtlon, eventtime, tshift,hdur,dep

def getdata(net,sta,eventtime,lents,dataloc):
#This function goes to one of the archives and gets the data
	debuggetdata = False

#Here are some hard coded values probably only necessary when messing with this function
	preeventday = eventtime - 24*60*60
	posteventday = eventtime + 24*60*60
	prepostwin= 3000
#If II get off of /tr1 else get the data from /xs0 or /xs1
	if net == 'II':
		dataprefix = '/tr1/telemetry_days/'
	else:
#Not using /tr1 data so lets get it from the archive
		if net in set(['IW','NE','US']):	
			dataprefix = 'xs1'	
		else:
			dataprefix = 'xs0'
		dataprefix = '/' + dataprefix + '/seed/'
	if debug:
		print 'Here is the dataprefix:' + dataprefix
	datastime = eventtime-prepostwin
	dataetime = eventtime+lents+prepostwin
	frstring = dataprefix + net + '_' + sta + '/'

#This data pull is kind of a mess do we want to change this approach
#Maybe we need a function to get the data from /xs0 or /tr1 without this mess
#Here we pull the event data, post event data, and the pre event data for the LH
	st = read( frstring + str(eventtime.year) + \
	'/' + str(eventtime.year) + '_' + str(eventtime.julday).zfill(3) + '*/*LH*.seed', \
	starttime=datastime,endtime=dataetime)
	
	st += read( frstring + str(posteventday.year) + \
	'/' + str(posteventday.year) + '_' + str(posteventday.julday).zfill(3) + '*/*LH*.seed', \
	starttime=datastime,endtime=dataetime)
	
	st += read(frstring + str(preeventday.year) + \
	'/' + str(preeventday.year) + '_' + str(preeventday.julday).zfill(3) + '*/*LH*.seed', \
	starttime=datastime,endtime=dataetime)

#Here we pull the event data, post event data, and the pre event data for the LN
	st += read(frstring + str(eventtime.year) + \
	'/' + str(eventtime.year) + '_' + str(eventtime.julday).zfill(3) + '*/*LN*.seed', \
	starttime=datastime,endtime=dataetime)

	st += read(frstring + str(posteventday.year) + \
	'/' + str(posteventday.year) + '_' + str(posteventday.julday).zfill(3) + '*/*LN*.seed', \
	starttime=datastime,endtime=dataetime)

	st += read(frstring + str(preeventday.year) + \
	'/' + str(preeventday.year) + '_' + str(preeventday.julday).zfill(3) + '*/*LN*.seed', \
	starttime=datastime,endtime=dataetime)

	st.merge(fill_value='latest')
	if debuggetdata:
		print 'We have data'
	return st

def getPAZ2(sp,net,sta,loc,chan,eventtime):
#This function gets the instrument response for the given net, sta, loc, chann at a fixed event time

	debuggetPAZ2 = False
        data = {}
	station_flag = False
	channel_flag = False
	for statemp in sp.stations:
		for blockette in statemp:
#Try to find blockette 50
			if blockette.id == 50:
				station_flag = False
				if net == blockette.network_code and sta == blockette.station_call_letters:
					station_flag = True
					if debuggetPAZ2:
						print 'We found the station blockettes'
			elif blockette.id == 52 and station_flag:
				channel_flag = False
#Okay we are in the right station, location, and chann
				if blockette.location_identifier == loc and blockette.channel_identifier == chan:
					if debuggetPAZ2:
						print 'We are in the location and channel blockette'
						print 'End date: ' + str(blockette.end_date)
						print 'Start date: ' + str(blockette.start_date)
					if type(blockette.end_date) is str:
						curdoy = strftime("%j",gmtime())
						curyear = strftime("%Y",gmtime())
						curtime = UTCDateTime(curyear + "-" + curdoy + "T00:00:00.0") 
						if blockette.start_date <= eventtime:
							channel_flag = True
							if debuggetPAZ2:
								print 'We found the channel blockette'
					elif blockette.start_date <= eventtime and blockette.end_date >= eventtime:
						channel_flag = True
						if debuggetPAZ2:
							print 'We found the channel blockette'
#Now we are in blockette 58 so we can get the various response parameters
			elif blockette.id == 58 and channel_flag and station_flag:
				if blockette.stage_sequence_number == 0:
					data['sensitivity'] = blockette.sensitivity_gain
				elif blockette.stage_sequence_number == 1:
					data['seismometer_gain'] = blockette.sensitivity_gain
				elif blockette.stage_sequence_number == 2:
					data['digitizer_gain'] = blockette.sensitivity_gain
			elif blockette.id == 53 and channel_flag and station_flag:
				data['gain'] = blockette.A0_normalization_factor
				data['poles'] = []
				if not blockette.transfer_function_types == 'A':
					msg = 'Only supporting Laplace transform response ' + \
					'type. Skipping other response information.'
					warnings.warn(msg, UserWarning)
					continue
				for i in range(blockette.number_of_complex_poles):
					p = complex(blockette.real_pole[i], blockette.imaginary_pole[i])
					data['poles'].append(p)
				data['zeros'] = []
				for i in range(blockette.number_of_complex_zeros):
					if debuggetPAZ2:
						print 'Here are the number of zeros: ' + str(blockette.number_of_complex_zeros)
					if blockette.number_of_complex_zeros > 1:
						z = complex(blockette.real_zero[i], blockette.imaginary_zero[i])
#For accels you can get some funny non-array type results so we have this piece for real poles
					else:
						z = complex(blockette.real_zero, blockette.imaginary_zero)
					data['zeros'].append(z)
        return data

def getcolor(chan,loc):
#This function sets the color of the trace by channel and location

#If we have an accelerometer we use black
	if chan in set(['LXN','LXE','LXZ']):
		color = 'k'

#Primary sensors will be green
	elif (loc == '00' or loc ==''):
		color = 'g'

#Secondary sensors will be red
	elif loc == '10':
		color = 'r'

#Third sensors this could be an odd broadband will be cyan
	elif loc == '20':
		color = 'c'
	else:

#Everything else is black
		color = 'b'
	return color



def writestats(statfile,streamin,comp):
#This function does the final stat computations for the accelerometer plots just produced
#This was taken from the synthetics code so the accel plays the role of the synthetic
	debugwstats = False
	try:
		syncomp = "LN" + comp	
		datacomp = "LH" + comp
		
		if debugwstats:
			print(streamin)
			print 'Here is the comp:' + syncomp
		syn = streamin.select(channel = syncomp)
		if debugwstats:
			print(syn)
		for tr in streamin.select(channel = datacomp):	
			if debugwstats:
				print 'Here is the trace value:' + str(numpy.sum(tr.data*syn[0].data))
				print 'Here is the accel value:' + str(numpy.sum(numpy.square(syn[0].data)))
#Here we compute a residual scale factor
			resi = "{0:.2f}".format(numpy.sum(tr.data*syn[0].data)/numpy.sum(numpy.square(syn[0].data)))
			if debugwstats:
				print 'Here is the resi:' + str(resi)
#Lets compute a cross-correlation and lag
			lag, corr = xcorr(tr,syn[0],50)
			corr = "{0:.2f}".format(corr)
			if debugwstats:
				print 'Here is the corr:' + str(corr)
				print 'Here are the results:' + tr.stats.network + "," + tr.stats.station
				print ' Here are more:' + "," + tr.stats.location + "," + tr.stats.channel + "," +  str(resi)
				print 'And more:' + "," + str(lag) + "," + str(corr) + "\n"

#Now we want to write to a file
			statfile.write(tr.stats.network + "," + tr.stats.station)
			statfile.write("," + tr.stats.location + "," + tr.stats.channel + "," +  str(resi))
			statfile.write("," + str(lag) + "," + str(corr) + "\n")
	
	except:	
		if debug:
			print 'No residual for' + cursta + ' ' + 'LH' + comp	
	return



#Start of the main part of the program
if not len(sys.argv) == 4:
	print "Usage: CMT ResultsName Network"
	exit(0)
 
cmtfile = sys.argv[1]

#Read in the CMT solution 
if debug:
	print "We are using local synthetics"
if not os.path.isfile(cmtfile):
	print "No CMT found"
	exit(0)
cmt = tuple(open(cmtfile))
cmtlat, cmtlon, eventtime, tshift, hdur, dep = readcmt(cmt)

#Lets make a local results directory
resultdir = sys.argv[2]
if resultdir[-1] == '/':
	resultdir = resultdir[:-1]
if not os.path.exists(os.getcwd() + '/' + resultdir):
	os.mkdir(os.getcwd() + '/' + resultdir)

#Lets get the current network
curnet = sys.argv[3]

#Lets open the results file to write
statfile = open(os.getcwd() + '/' + resultdir + '/Results' + curnet + '.csv' ,'w')
statfile.write('net,sta,loc,chan,scalefac,lag,corr\n')




#Lets read in the dataless
try:
	sp = Parser(datalessloc + curnet + ".dataless")
except:
	print "Can not read the dataless."
	exit(0)

#If we arent doing manual station lists we need to get one for the network
if not manstalist:
	stations = getstalist(sp,eventtime,curnet)

if debug:
	print "Here are the stations we found"	
	for sta in stations:
		print "Here is a station:" + sta


#Lets start by using a station list and then move to a different approach
for sta in stations:
	cursta = sta.strip()
	if debug:
		print 'Current station:' + cursta
#Time to split the cursta into its network and current station 
	cursta = sta.split()
	net = cursta[0]
	cursta = cursta[1]

#Now we get the data for the event
	try:
		st = getdata(net,cursta,eventtime,lents,dataloc)
	except:
		print('No data for ' + net + ' ' + cursta)
		continue
		
#Lets go through each trace in the stream and deconvolve and filter
	for trace in st:
#Here we get the response and remove it
		
		paz=getPAZ2(sp,net,cursta,trace.stats.location,trace.stats.channel,eventtime)
		try:

			trace.taper(max_percentage=0.05, type='cosine')
#If we have an accelerometer we want an extra zero to go to displacement
			if trace.stats.channel in ('LNZ','LN1','LN2','LNE','LNN'):
				paz['zeros'].append(0+0j)
#Here we remove the response
			trace.simulate(paz_remove=paz)
#Here we filter, integrate, taper, trim, detrend, and filter
			trace.filter("bandpass",freqmin = userminfre,freqmax= usermaxfre, corners=filtercornerpoles)
			trace.integrate()
			trace.taper(max_percentage=0.05, type='cosine')
			trace.trim(starttime=eventtime + tshift/2,endtime=(eventtime+lents + tshift/2))
			trace.detrend()
			trace.filter("bandpass",freqmin = userminfre,freqmax= usermaxfre, corners=filtercornerpoles)
		except:
			print('Can not find the response')
			st.remove(trace)

#Lets check for reverse polarity and fix 
	finalstream=Stream()
	for trace in st.select(component="Z"):
		dipval = getdip(net,cursta,trace.stats.location,trace.stats.channel,eventtime,sp)
		if debug:
			print 'Here is the dip value:' + str(dipval)
		if dipval == 90.0:
			trace.data = -trace.data
		finalstream += trace
		st.remove(trace)

#Now we rotate the horizontals to N/S and E/W
#Should the rotation be put into a function to remove it from the loop?
	locations=[]
	for trace in st:
		locations.append(trace.stats.location)
	locations=set(locations)

	for curloc in locations:
		curlochorizontal = st.select(location=curloc)
		curlochorizontal.sort(['channel'])
		if debug:
			print "Here are the number of traces:" + str(len(curlochorizontal)) + " which should be 2"
			print(curlochorizontal)
		azi1=getorientation(net,cursta,curloc,curlochorizontal[0].stats.channel,eventtime,sp)
		azi2=getorientation(net,cursta,curloc,curlochorizontal[1].stats.channel,eventtime,sp)
		if debug:
			print "Here is the azimuth for " + net + " " + cursta + " " + curloc + " " + curlochorizontal[0].stats.channel + str(azi1)
			print "Here is the azimuth for " + net + " " + cursta + " " + curloc + " " + curlochorizontal[1].stats.channel + str(azi2)
		curlochorizontal = choptocommon(curlochorizontal)
		try:
			finalstream += rotatehorizontal(curlochorizontal,azi1,azi2)	
		except:
			print 'Can not rotate using azi1:' + str(azi1) + ' and azi2:' + str(azi2)
			
	if debug:
		print(finalstream)



#We now need to plot everything and save it
#Lets plot the verticals first
	vertcomps = finalstream.select(component="Z")
	vertcomps = choptocommon(vertcomps)
#We want to get the distance of the event and of the station
#We also want the back-azimuth
	lat,lon = getlatlon(cursta,eventtime,sp)
	dist= gps2DistAzimuth(float(cmtlat),float(cmtlon),lat,lon)
	bazi ="{0:.1f}".format(dist[2])
	dist ="{0:.1f}".format( 0.0089932 * dist[0] / 1000)
	if debug:
		print 'Here is the distance:' + str(dist)
		print 'Here is the depth:' + str(dep)

#Here is the travel time so we can do the final trim
#Should this be in a function to avoid it being in the main loop?
	tt = getTravelTimes(delta=float(dist), depth = dep,model='ak135') 
	firstarrival = tt[0]['time']
	for ttphase in tt:
		phasename = ttphase['phase_name']
		phasename = phasename[:1]
		if phasename == 'S':
			secondarrival = ttphase['time']
			break
#Here we do the trim from the phases	
	for trace in vertcomps:
		newstime = trace.stats.starttime + firstarrival - bfarrival
		newetime = trace.stats.starttime + secondarrival + afarrival
		trace.trim(starttime=newstime,endtime=newetime)


	if debug:
		print 'Here is the first arrival time: ' + str(firstarrival)
		print 'Here is the second arrival time: ' + str(secondarrival)
		print 'Here are the chopped components'
		print(vertcomps)
		
#Here is the mess of plotting info  
#Set the time series
	tz=numpy.arange(0,vertcomps[0].stats.npts / vertcomps[0].stats.sampling_rate, vertcomps[0].stats.delta)

#Get a legend and plot the vertical
	synplot = figure(1)

#Here we setup subplot 1 and do a title
	subplot(311)
	titlelegend = vertcomps[0].stats.network + ' ' + vertcomps[0].stats.station + ' '

#Here is the start time of the plot
	stime = str(vertcomps[0].stats.starttime.year) + ' ' + str(vertcomps[0].stats.starttime.julday) + ' ' + \
	str(vertcomps[0].stats.starttime.hour) + ':' + str(vertcomps[0].stats.starttime.minute) + \
	':' + str("{0:.2f}".format(vertcomps[0].stats.starttime.second))

	titlelegend = titlelegend + stime + ' ' 
	
	if debug:
		print "Latitude:" + str(lat)
		print "Longitude:" + str(lon)	
		print "CMT Latitude:" + str(cmtlat)
		print "CMT Longitude:" + str(cmtlon)
	
#Here we add the distance and the back-azimuth to the legend	
	titlelegend = titlelegend + 'Dist:' + str(dist) 
	titlelegend = titlelegend + ' BAzi:' + str(bazi) 

#Here we add the frequencies
	minper = "{0:.0f}".format(1/usermaxfre)
	maxper = "{0:.0f}".format(1/userminfre)
	titlelegend = titlelegend + ' ' + str(minper) + '-' + str(maxper) + ' s per.'
	title(titlelegend,fontsize=12)
	vertcomps.sort(['location'])
	for comps in vertcomps.select(component="Z"):
		curcolor = getcolor(comps.stats.channel,comps.stats.location)
		plot(tz,(comps.data*(10**3)), curcolor, label=comps.stats.location + ' ' + comps.stats.channel)
	legend(prop={'size':6})
	xlim(0,len(tz))
	finalstream = choptocommon(finalstream)
	finalstream.sort(['location','channel'])
	if debug:
		print "Here is the final stream:"
		print(finalstream)

#We now plot the N/S component of the data
	subplot(312)
	tne=numpy.arange(0,finalstream[0].stats.npts / finalstream[0].stats.sampling_rate, finalstream[0].stats.delta)
	for comps in finalstream.select(component="N"):
		curcolor = getcolor(comps.stats.channel,comps.stats.location)
		plot(tne,(comps.data*(10**3)),curcolor, label=comps.stats.location + ' ' + comps.stats.channel)
	legend(prop={'size':6})
	ylabel('Displacement (mm)')	
	xlim(0,len(tne))
	
#Now we plot the E/W compoent of the data
	subplot(313)
	for comps in finalstream.select(component="E"):
		curcolor = getcolor(comps.stats.channel,comps.stats.location)
		plot(tne,(comps.data*(10**3)),curcolor, label=comps.stats.location + ' ' + comps.stats.channel)
	legend(prop={'size':6})
	xlabel('Time (s)')
	xlim(0,len(tne))

#Finally we need to save the figure
	savefig(os.getcwd() + '/' + resultdir + '/' + cursta + \
	str(vertcomps[0].stats.starttime.year) + str(vertcomps[0].stats.starttime.julday) + \
	str(vertcomps[0].stats.starttime.hour) + str(vertcomps[0].stats.starttime.minute) + '.jpg', format = 'jpeg', dpi=400)

#Lets clear the plot so we have no residual
	synplot.clear()

#Time to write some info into the statfile
#Write the network and the station
	writestats(statfile,vertcomps,'Z')
	writestats(statfile,finalstream,'N')
	writestats(statfile,finalstream,'E')
	
	
#Lets get an RMS from the synthetic and the data
	


statfile.close()
