#!/usr/bin/env python

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import os, sys
import numpy as np
sys.path.extend(map(os.path.abspath, ['/nfs/users/nfs_s/sh16/scripts/modules/']))
import shapefile
from matplotlib.collections import LineCollection
from matplotlib import cm
from optparse import OptionParser, OptionGroup
from math import floor, ceil

##########################
# Error message function #
##########################

def DoError(errorstring):
	print "\nError:", errorstring
	print "\nFor help use -h or --help\n"
	sys.exit()


def revcomp(sequence):
	rev=sequence[::-1]
	revcomp=''
	d={'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'g':'c', 'c':'g', "n":"n", "N":"N"}
	for i in rev:
		if d.has_key(i):
			revcomp=revcomp+d[i]
		else:
			revcomp=revcomp+i
	
	return revcomp


##############################
# Get command line arguments #
##############################

def get_user_options():
	usage = "usage: %prog [options] <list of bam files>"
	version="%prog 1.0. Written by Simon Harris, Wellcome Trust Sanger Institute, 2011"
	parser = OptionParser(usage=usage, version=version)
	
	#do not allow arguments to be interspersed. e.g. -a -b arg1 agr2. MUST be -a arg1 -b arg2.
	parser.disable_interspersed_args()
	
	parser.add_option("-o", "--output", action="store", dest="output", help="Output file prefix", default="")
	parser.add_option("-t", "--type", action="store", dest="filetype", help="Output file type. Choose from png or pdf [Default=%default]", type="choice", choices=["png", "pdf", "screen"], default="pdf")
	parser.add_option("-c", "--continents", action="store_true", dest="continents", help="Colour continents", default=False)
	parser.add_option("-a", "--area", action="store", dest="area", help="Area of world map to show", type="choice", choices=["world", "europe", "africa", "australia", "UK", "USA", "SAmerica", "user"], default="world")
	parser.add_option("-u", "--userarea", action="store", dest="user", help="User-defined area of world map to show. Only valid with user choice of area (-a) option. Must be a comma separated list of lower left latitude, upper right latitude, lower left longitude, upper right longitude", default="")
	parser.add_option("-p", "--projection", action="store", dest="projection", help="Map projection to use.", type="choice", choices=['cyl', 'merc', 'mill', 'gall', 'tmerc'], default="merc")
	parser.add_option("-f", "--fill", action="store", dest="background", help="Background fill to use", type="choice", choices=['white', 'SR', 'etopo', 'bluemarble'], default="white")
	parser.add_option("-C", "--countries", action="store_true", dest="countries", help="Draw countries", default=False)
	parser.add_option("-S", "--states", action="store_true", dest="states", help="Draw states (USA)", default=False)
	parser.add_option("-R", "--rivers", action="store_true", dest="rivers", help="Draw rivers", default=False)
	parser.add_option("-O", "--coastlines", action="store_false", dest="coastlines", help="Do not draw coastlines", default=True)
	parser.add_option("-M", "--meridians", action="store_true", dest="meridians", help="Draw meridians", default=False)
	parser.add_option("-P", "--parallels", action="store_true", dest="parallels", help="Draw parallels", default=False)
	
	
	return parser.parse_args()


################################
# Check command line arguments #
################################

def check_input_validity(options, args):

	if options.output=='':
		options.output="test"
		
		
		
		
		
########
# Main #
########


if __name__ == "__main__":



	(options, args)=get_user_options()
	check_input_validity(options, args)		
		
			
	
	fig = plt.figure(figsize=(11.7,8.3))
	
	
	#Custom adjust of the subplots
	plt.subplots_adjust(left=0.05,right=0.95,top=0.90,bottom=0.05,wspace=0.15,hspace=0.05)
	ax = plt.subplot(111)
	
	#if options.continents or options.area=="world":
	if options.area=="world":
		x1 = -180
		x2 = 180
		y1 = -80
		y2 = 80
		resolution='c'
		#m = Basemap(width=12000000,height=9000000, resolution='c',projection=options.projection, llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2)
		#m = Basemap(resolution='c',projection='tmerc', llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2, lat_0=50, lon_0=0)
		#m = Basemap(projection='robin',lon_0=0,resolution='c')
		#m.fillcontinents()
		
		#m.warpimage(image="HYP_50M_SR_W/HYP_50M_SR_W.tif", scale=0.25)
		#m.drawparallels(np.arange(-90.,120.,30.))
		#m.drawmeridians(np.arange(0.,360.,60.))
	elif options.area=="europe":
		if options.projection=="tmerc":
			x1 = -12
			x2 = 50
			y1 = 32
			y2 = 64
		else:
			x1 = -20
			x2 = 40
			y1 = 32
			y2 = 64
		resolution='i'
		#m = Basemap(resolution='i',projection=options.projection, llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2,lat_ts=(x2-x1)/2)
		#m.shadedrelief()
	elif options.area=="africa":
		if options.projection=="tmerc":
			x1 = -30
			x2 = 65
			y1 = -32
			y2 = 34
		else:
			x1 = -20
			x2 = 54
			y1 = -38
			y2 = 40
		resolution='i'
		#m = Basemap(resolution='i',projection=options.projection, llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2,lat_ts=(x2-x1)/2)
		
	elif options.area=="australia":
		if options.projection=="tmerc":
			x1 = 106
			x2 = 152
			y1 = -45
			y2 = -8
		else:
			x1 = 109
			x2 = 157
			y1 = -45
			y2 = -8
		resolution='i'
		#m = Basemap(resolution='i',projection=options.projection, llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2,lat_ts=(x2-x1)/2)
	elif options.area=="USA":
		if options.projection=="tmerc":
			print "tmerc is broken for USA"
			sys.exit()
			x1 = -100
			x2 = -50
			y1 = 20
			y2 = 72
		else:
			x1 = -170
			x2 = -50
			y1 = 20
			y2 = 72
		resolution='i'
	
	elif options.area=="UK":
		if options.projection=="tmerc":
			x1 = -10.5
			x2 = 3.5
			y1 = 49.5
			y2 = 59.5
			resolution='i'
		else:
			x1 = -11
			x2 = 3.5
			y1 = 49.5
			y2 = 59.5
			resolution='i'
	
	elif options.area=="SAmerica":
		if options.projection=="tmerc":
			x1 = -104
			x2 = -30
			y1 = -52
			y2 = 20
		else:
			x1 = -90
			x2 = -30
			y1 = -60
			y2 = 20
		resolution='i'	
	
	elif options.area=="user":
		try:
			latlongs=map(float,options.user.split(","))
		except StandardError:
			DoError("Illegal user defined map area "+options.user+". Must be a comma-separated list of four numbers")
		x1 = latlongs[2]
		x2 = latlongs[3]
		y1 = latlongs[0]
		y2 = latlongs[1]
		resolution='c'
		if x2-x1>180:
			resolution='c'
		else:
			resolution='i'
	
	if options.projection in ["merc", "mill", "cyl", "gall"]:
		m=Basemap(resolution=resolution, projection=options.projection,llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2,lat_ts=y2-((y2-y1)/2))
	elif options.projection=="tmerc":
		if options.area=="world":
			print "Cannot show whole world with tmerc projection. Use merc instead."
			sys.exit()
		
		
		m=Basemap(resolution=resolution, projection=options.projection,llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2,lat_0=y2-((y2-y1)/2), lon_0=x2-((x2-x1)/2))
		#print y1,y2,x1,x2,y2-((y2-y1)/2), x2-((x2-x1)/2)
#		print m(x1, y1)
#		print m(x1, y2)
#		print m(x2,y1)
#		print m(x2, y2)
		
		
	m.drawmapboundary(fill_color='white')
	if options.background=="SR":
		m.warpimage(image="/lustre/scratch108/bacteria/sh16/maptest/HYP_50M_SR_W/HYP_50M_SR_W.tif", scale=0.25)
	elif options.background=="etopo":
		m.etopo()
	elif options.background=="bluemarble":
		m.bluemarble()
	
	if options.coastlines:
		m.drawcoastlines(linewidth=0.5)
	
	if options.parallels:
		yb1=ceil(y1)
		yb2=floor(y2)
		if yb1%2!=0:
			yb1+=1
		diff=ceil((yb2-yb1)/10)
		
		if diff%2!=0:
			diff+=1
		#print y1, y2, yb1, yb2, diff, np.arange(yb1,yb2+diff*5,diff), np.arange(yb1-(diff*5),yb2+(diff*5),diff)
		m.drawparallels(np.arange(yb1-(diff*5),yb2+(diff*5),diff), labels=[1,1,0,0], color='black', dashes=[1,1], labelstyle='+/-', linewidth=0.2)
		
	if options.meridians:
		xb1=ceil(x1)
		xb2=floor(x2)
		if xb1%2!=0:
			xb1+=1
		diff=ceil((xb2-xb1)/10)
		if diff%2!=0:
			diff+=1
		#print x1, x2, xb1, xb2, diff, np.arange(xb1,xb2,diff), np.arange(xb1-(diff*5),xb2+(diff*5),diff)
		m.drawmeridians(np.arange(xb1-(diff*5),xb2+(diff*5),diff), labels=[0,0,1,1], color='black', dashes=[1,1], labelstyle='+/-', linewidth=0.2)
	
	if options.countries:
		m.drawcountries()
	if options.states:
		m.drawstates()
	if options.rivers:
		m.drawrivers()
	
	
	
#	if options.filetype=="screen":
#		plt.show()	
#	else:
#		plt.savefig(options.output+"."+options.filetype,dpi=300)
#	
#	#print dir(sf)
#	#print dir(sf.record)
#	sys.exit()
	
	
		#m.shadedrelief()
	#m = Basemap(resolution='c',projection='merc', llcrnrlat=y1,urcrnrlat=y2,llcrnrlon=x1,urcrnrlon=x2)
	#m = Basemap(resolution='i',projection='merc')
	#m.drawcountries(linewidth=0.5)
	
	#m.drawparallels(np.arange(y1,y2,2.),labels=[1,0,0,0],color='black',dashes=[1,0],labelstyle='+/-',linewidth=0.2) 
	# draw parallels
	#m.drawmeridians(np.arange(x1,x2,2.),labels=[0,0,0,1],color='black',dashes=[1,0],labelstyle='+/-',linewidth=0.2) 
	# draw meridians
	#m.shadedrelief()
	
	#m.warpimage(image="HYP_50M_SR_W/HYP_50M_SR_W.tif")
	#geography_sf = shapefile.Reader("50m_physical/50m_geography_regions_polys")
	
	#shapes=geography_sf.shapes()
	#print dir(geography_sf)
	#geography_sr=geography_sf.shapeRecords()
	
	#for test in geography_sr:
		#print dir(test)
	#	if test.record[4].lower()=="europe":
	#		print test.record
		#if test.record[1]=="continent" and test.record[2].lower()=="europe":
	#	if test.record[4].lower() in ["europe", "north america"] and test.record[0]<2:
	#		parts= test.shape.parts
	#		points=test.shape.points
	
	
	#		print test.shape.parts
			
	#		partpoints=[]
	#		for x, part in enumerate(parts):
	#			if x<len(parts)-1:
	#				partpoints.append(points[part:parts[x+1]])
	#			else:
	#				partpoints.append(points[part:])
	#		if test.record[4].lower()=="europe":
	#			col=(1,0,0)
	#		else:
	#			col=(0,0,1)
	#		for partpoint in partpoints:
	
	#			shpsegs = []
	#			lats=[]
	#			lons=[]
				#print points
	#			for point in partpoint:
	#				lats.append(point[0])
	#				lons.append(point[1])
	#			x,y = m(lats, lons)
	#			shpsegs.append(zip(x,y))
				#lonpt, latpt = m(xpt,ypt,inverse=True)#m.plot(xpt,ypt,'r*', ms=2)
	#			lines = LineCollection(shpsegs,antialiaseds=(1,))
	#			lines.set_facecolors(col)
	#			lines.set_edgecolors('k')
	#			lines.set_linewidth(0.3)
	#			ax.add_collection(lines)
	
	#plt.show()
	#sys.exit()
	
	
	sf = shapefile.Reader("/lustre/scratch108/bacteria/sh16/maptest/50m_cultural/ne_50m_admin_0_countries")
	shapes=sf.shapes()
	fields=sf.fields
	#print fields
	
	sr=sf.shapeRecords()
	#print dir(sf)
	#print dir(sf.record)
	#print sf.fields
	#sys.exit()
	locations={}
	for x, test in enumerate(sr):
		locations[test.record[16]]=x
		#print test.record
		#print dir(test.shape)
		#sys.exit()
	
	
	
	
	
	
	if options.continents:
		mylocs=locations.keys()
		loc2continent={'Canada': "north america", 'East Timor': "africa", 'Sao Tome and Principe': 201, 'Ashmore and Cartier Is.': 12, 'Turkmenistan': "asia", 'Saint Helena': 192, 'Vatican': "europe", 'Lithuania': "europe", 'N. Korea': "asia", 'Cambodia': "asia", 'Ethiopia': "africa", 'Aruba': 0, 'Swaziland': "africa", 'Grenada': "north america", 'Argentina': "south america", 'Bolivia': "south america", 'Cameroon': "africa", 'Burkina Faso': "africa", 'Aland': 5, 'Bahrain': "asia", 'Saudi Arabia': "asia", 'Fr. Polynesia': 180, 'Japan': "asia", 'Cape Verde': 49, 'W. Sahara': "africa", 'Slovenia': "europe", 'St. Barthelemy': 27, 'Kuwait': "asia", 'Jordan': "asia", 'Taiwan': "asia", 'Dominica': 59, 'Liberia': "africa", 'Congo (Kinshasa)': "africa", 'Marshall Is.': 141, 'Jamaica': "north america", 'Solomon Is.': 193, 'Oman': "asia", 'Tanzania': "africa", 'Mauritania': "africa", 'Greenland': "north america", 'Gabon': "africa", 'Niue': 161, 'Monaco': "europe", 'Wallis and Futuna': 236, 'New Zealand': "australia", 'Yemen': "asia", 'Jersey': "europe", 'Bahamas': "north america", 'Albania': "europe", 'West Bank': "asia", 'Macau': "asia", 'Norfolk Island': 158, 'United Arab Emirates': "asia", 'Guam': 89, 'Uruguay': "south america", 'India': "asia", 'Azerbaijan': "asia", 'St. Vin. and Gren.': "north america", 'Lesotho': "africa", 'Congo (Brazzaville)': "africa", 'Kenya': "africa", 'Tajikistan': "asia", 'Turkey': "europe", 'Afghanistan': "asia", 'Bangladesh': "asia", 'Indian Ocean Ter.': "asia", 'Lebanon': "africa", 'Saint Lucia': "north america", 'Br. Indian Ocean Ter.': 101, 'San Marino': "europe", 'Kyrgyzstan': "asia", 'Mongolia': "asia", 'France': "europe", 'S. Sudan': "africa", 'Bermuda': 30, 'Namibia': "africa", 'Somalia': "africa", 'Peru': "south america", 'Laos': "asia", 'Nauru': "australia", 'Seychelles': 208, 'Norway': "europe", 'Malawi': "africa", 'Benin': "africa", 'Cuba': "north america", 'Turks and Caicos Is.': "north america", 'Montenegro': "europe", 'Cayman Is.': 53, 'Togo': "africa", 'China': "asia", 'Heard I. and McDonald Is.': 92, 'Armenia': "asia", 'Falkland Is.': "south america", 'Ukraine': "europe", 'Ghana': "africa", 'Tonga': "australia", 'Indonesia': "asia", 'Libya': "africa", 'Dominican Rep.': "north america", 'Finland': "europe", 'Mauritius': 152, 'Eq. Guinea': "africa", 'Liechtenstein': "europe", 'Belarus': "europe", 'St. Martin': 134, 'Mali': "africa", 'Poland': "europe", 'Russia': "europe", 'Bulgaria': "europe", 'United States': "north america", 'Romania': "europe", 'Czech Rep.': "europe", 'St. Kitts and Nevis': 118, 'Portugal': "europe", 'Trinidad and Tobago': "north america", 'British Virgin Is.': "north america", 'Cyprus': "europe", 'Sweden': "europe", 'Qatar': "asia", 'Malaysia': "asia", 'Austria': "europe", 'Vietnam': "asia", 'Mozambique': "africa", 'Uganda': "africa", 'Hungary':"europe" , 'Niger': "africa", 'Brazil': "south america", 'Netherlands': "europe", 'Guinea': "africa", 'Panama': "north america", 'Costa Rica': "north america", 'Luxembourg': "europe", 'American Samoa': 10, 'Andorra': 6, 'Chad': "africa", 'Ivory Coast': "africa", 'Pakistan': "asia", 'Palau': "asia", 'Nigeria': "africa", 'Somaliland': "africa", 'Ecuador': "south america", 'U.S. Virgin Is.': "north america", 'Brunei': "asia", 'Australia': "australia", 'Iran': "asia", 'Algeria': "africa", 'El Salvador': "north america", 'Fr. S. and Antarctic Lands': 13, 'Guatemala': "north america", 'Chile': "south america", 'Puerto Rico': "north america", 'Belgium': "europe", 'Kiribati': 117, 'Haiti': "north america", 'Belize': "north america", 'Hong Kong': "asia", 'Sierra Leone': "africa", 'Georgia': "asia", 'Gambia': "africa", 'Philippines': "asia", 'Guinea Bissau': "africa", 'S. Geo. and S. Sandw. Is.': 191, 'Moldova': "europe", 'Morocco': "africa", 'Croatia': "europe", 'Malta': "europe", 'Guernsey': "europe", 'Thailand': "asia", 'Switzerland': "europe", 'Angola': "africa", 'Isle of Man': "europe", 'Myanmar': "asia", 'Estonia': "europe", 'Kosovo': "europe", 'Antigua and Barb.': "north america", 'South Africa': "africa", 'Faroe Is.': "europe", 'Uzbekistan': "asia", 'Tunisia': "africa", 'Djibouti': "africa", 'Rwanda': "africa", 'Spain': "europe", 'Colombia': "south america", 'Burundi': "africa", 'Slovakia': "europe", 'Pitcairn Is.': "europe", 'Fiji': "australia", 'Barbados': "north america", 'Madagascar': "africa", 'Italy': "europe", 'Bhutan': "asia", 'Sudan': "africa", 'Nepal': "asia", 'Micronesia': "aisa", 'Maldives': "asia", 'S. Korea': "asia", 'Suriname': "south america", 'Gaza': "asia", 'Anguilla': "north america", 'Venezuela': "south america", 'Israel': "asia", 'St. Pierre and Miquelon': "north america", 'Central African Rep.': "africa", 'Iceland': "europe", 'Zambia': "africa", 'Samoa': "australia", 'Senegal': "africa", 'Papua New Guinea': "australia", 'Zimbabwe': "africa", 'Germany': "europe", 'Vanuatu': "australia", 'Denmark': "europe", 'Kazakhstan': "asia", 'N. Mariana Is.': "north america", 'Eritrea': "africa", 'Ireland': "europe", 'N. Cyprus': "europe", 'Iraq': "asia", 'Montserrat': "north america", 'New Caledonia': "australia", 'Macedonia': "europe", 'Paraguay': "south america", 'Latvia': "europe", 'Guyana': "south america", 'Syria': "asia", 'Sint Maarten': "europe", 'Cook Is.': "australia", 'Honduras': "north america", 'Bosnia and Herz.': "europe", 'Mexico': "north america", 'Egypt': "africa", 'Nicaragua': "north america", 'Singapore': "acia", 'Serbia': "europe", 'Botswana': "africa", 'United Kingdom': "europe", 'Antarctica': "antarctica", 'Greece': "europe", 'Sri Lanka': "asia", 'Cura\xe7ao': "south america", 'Comoros': "africa", '						  ': "?"}
	
	
		regioncolours={"europe": (1,0,0), "north america": (0,1,1), "asia": (0,0,1), "south america": (0,1,0), "antartica": (0,0,0), "africa": (1,0,1), "australia": (1,1,0)}
	
		mycontinents=["africa", "asia", "australia", "europe", "north america", "south america"]
		regioncolours={}
		if len(mycontinents)==2:
			regioncolours[mycontinents[0]]=(1, 0, 0)
			regioncoloursmycontinents[[1]]=(0, 0, 1)
		else:
			for x, name in enumerate(mycontinents):	
							
				#proportion=(float(x)/(len(colourslist)-1))*1275
				proportion=(float(x)/(len(mycontinents)-1))*1175
						
				red=510-proportion
				blue=proportion-510
				green=proportion
				if red<-510:
					reddiff=510+red
					red=reddiff*-1
				elif red<0:
					red=0
				elif red>255:
					red=255
				if blue<0:
					blue=0
				elif blue>255:
					blue=255
			#			if green>255 and green<765:
			#				green=255
				if green>155 and green<765:
					green=155
				elif green>=765:
					greendiff=765-green
					green=155+greendiff
					if green<0:
						green=0
				regioncolours[name]=(float(red)/255, float(green)/255, float(blue)/255)
		
	else:
		#mylocs=["Argentina", "Australia", "Canada", "China", "Denmark", "Finland", "France", "Hong Kong", "Italy", "Japan", "Kenya", "Netherlands", "Poland", "Russia", "Senegal", "Sweden", "Taiwan", "United Kingdom"]
		#mylocs=["France", "United Kingdom", "Hungary", "Belgium", "Sweden",  "Ireland", "Italy", "Greece", "Spain", "Portugal", "Germany", "Ukraine", "Netherlands", "Poland", "Latvia", "Estonia", "USA", "Canada", "Brazil"]
		mylocs=["Argentina", "Canada", "Denmark", "France", "Gambia", "Italy", "Netherlands", "Russia", "Saudi Arabia", "South Africa", "Sweden", "Taiwan", "Tanzania", "United Kingdom", "United States"]
		regioncolours={}
		if len(mylocs)==2:
			regioncolours[mylocs[0]]=(1, 0, 0)
			regioncolours[mylocs[1]]=(0, 0, 1)
		else:
			for x, name in enumerate(mylocs):	
							
				#proportion=(float(x)/(len(colourslist)-1))*1275
				proportion=(float(x)/(len(mylocs)-1))*1175
						
				red=510-proportion
				blue=proportion-510
				green=proportion
				if red<-510:
					reddiff=510+red
					red=reddiff*-1
				elif red<0:
					red=0
				elif red>255:
					red=255
				if blue<0:
					blue=0
				elif blue>255:
					blue=255
			#			if green>255 and green<765:
			#				green=255
				if green>155 and green<765:
					green=155
				elif green>=765:
					greendiff=765-green
					green=155+greendiff
					if green<0:
						green=0
				regioncolours[name]=(float(red)/255, float(green)/255, float(blue)/255)
	
	
	for myloc in mylocs:
		if myloc in locations:
			test=sr[locations[myloc]]
			points= test.shape.points
			parts= test.shape.parts
			partpoints=[]
	
			for x, part in enumerate(parts):
				if x<len(parts)-1:
					partpoints.append(points[part:parts[x+1]])
				else:
					partpoints.append(points[part:])
			if options.continents and loc2continent[myloc] in regioncolours:
				col=regioncolours[loc2continent[myloc]]
			elif myloc in regioncolours:
				col=regioncolours[myloc]
			else:
				#continue
				col=(0,0,0)
				#col=np.random.rand(1)
			for partpoint in partpoints:
				shpsegs = []
				lats=[]
				lons=[]
				for point in partpoint:
					lats.append(point[0])
					lons.append(point[1])
				x,y = m(lats, lons)
				shpsegs.append(zip(x,y))
				#lonpt, latpt = m(xpt,ypt,inverse=True)
				#m.plot(xpt,ypt,'r*', ms=2)
				lines = LineCollection(shpsegs,antialiaseds=(1,))
				lines.set_facecolors(col)
				lines.set_edgecolors('k')
				lines.set_linewidth(0.3)
				ax.add_collection(lines)
		else:
			print "Could not find", myloc, "in shapefiles"
	
	
	if options.filetype=="screen":
		plt.show()	
	else:
		plt.savefig(options.output+"."+options.filetype,dpi=300)
	
	#print dir(sf)
	#print dir(sf.record)
	sys.exit()
	for x, shape in enumerate(shapes):
		#print dir(shape)
		#print shape.shapeType
		points= shape.points
		for point in points:
			#print point
			xpt,ypt = m(point[0], point[1])
			lonpt, latpt = m(xpt,ypt,inverse=True)
			m.plot(xpt,ypt,'r*', ms=16)
		plt.show()
					  
	
		
		sys.exit()
	
	print shapes
	
	
	plt.savefig('tutorial07.png',dpi=300)
	plt.show()
	
	
	
	
	
	
	
	sys.exit()
	
	fig=plt.figure(figsize=(15,15))
	
	# setup Lambert Conformal basemap.
	m = Basemap(width=12000000,height=9000000,projection='tmerc', resolution='c',llcrnrlon=-30,llcrnrlat=-30,urcrnrlon=65,urcrnrlat=25, lon_0=10, lat_0=10)
	# draw coastlines.
	#m.drawcoastlines()
	# draw a boundary around the map, fill the background.
	# this background will end up being the ocean color, since
	# the continents will be drawn on top.
	#m.drawmapboundary(fill_color='lightblue') 
	# fill continents, set lake color same as ocean color. 
	#m.fillcontinents(color='lightgreen',lake_color='lightblue')
	m.drawcountries(color='gray')
	
	lats= [-33.227792, -32.29684, -29.8167, -29.52984, -28.45411, -26.66386, -26.270759, -15.786111, 6.137778, 6.8166667, 13.211728, 13.245833, 13.31041, 13.47, 13.512668, 13.9, 14.216667, 14.421311, 14.75, 14.758889]
	lons = [21.856859, 26.419389, 30.6167, 31.20209, 26.796785, 25.283759, 28.112268, 35.005833, 1.2125, 3.2, -16.195976, 7.273889, -14.22346, -16.696389, 2.112516, 6.016667, 1.45, 6.044096, -17.333333, 5.774722]
	name = ["Western Cape, South Africa", "Eastern Cape, South Africa", "Mpumalanga, South Africa", "KwaZulu-Natal, South Africa", "Free State, South Africa", "North West, South Africa", "Gauteng, South Africa", "Southern, Malawi", "Lome, Togo", "South West, Nigeria", "Western Division, The Gambia", "Maradi, Niger", "Upper Region Division, The Gambia", "Western Division, The Gambia", "Niamey, Niger", "Dosso, Niger", "Tillabery, Niger", "Tahoua, Niger", "Dakar, Senegal", "Tahoua, Niger"]
	
	x, y = m(lons, lats)
	m.plot(x, y, '.')
	
	m.shadedrelief()
	plt.savefig(option.output+"."+options.filetype)
	sys.exit()
	# meridians on bottom and left
	#parallels = np.arange(0.,81,10.)
	# labels = [left,right,top,bottom]
	#m.drawparallels(parallels,labels=[False,True,True,False])
	#meridians = np.arange(10.,351.,20.)
	#m.drawmeridians(meridians,labels=[True,False,False,True])
	# plot blue dot on Blantyre
	lon, lat = 35.005833, -15.786111 # Location of Blantyre
	# convert to map projection coords. 
	# Note that lon,lat can be scalars, lists or numpy arrays.
	xpt,ypt = m(lon,lat) 
	# convert back to lat/lon
	lonpt, latpt = m(xpt,ypt,inverse=True)
	m.plot(xpt,ypt,'r*', ms=16)  # plot a blue dot there
	# put some text next to the dot, offset a little bit
	# (the offset is in map projection coordinates)
	plt.text(xpt+100000,ypt+100000,'MLW, Blantyre, Malawi', size=14)
	
	lon, lat = 28.047305, -26.204103 # Location of Gauteng
	# convert to map projection coords. 
	# Note that lon,lat can be scalars, lists or numpy arrays.
	xpt,ypt = m(lon,lat) 
	# convert back to lat/lon
	lonpt, latpt = m(xpt,ypt,inverse=True)
	m.plot(xpt,ypt,'r*', ms=16)  # plot a blue dot there
	# put some text next to the dot, offset a little bit
	# (the offset is in map projection coordinates)
	plt.text(xpt+100000,ypt+100000,'NICD, Johannesberg, SA', size=14)
	
	lon, lat = 2.112516, 13.512668 # Location of Niamey
	# convert to map projection coords. 
	# Note that lon,lat can be scalars, lists or numpy arrays.
	xpt,ypt = m(lon,lat) 
	# convert back to lat/lon
	lonpt, latpt = m(xpt,ypt,inverse=True)
	m.plot(xpt,ypt,'r*', ms=16)  # plot a blue dot there
	# put some text next to the dot, offset a little bit
	# (the offset is in map projection coordinates)
	plt.text(xpt+100000,ypt+100000,'CERMES, Niamey, Niger', size=14)
	
	lon, lat = -16.5775, 13.453056 # Location of Banjul
	# convert to map projection coords. 
	# Note that lon,lat can be scalars, lists or numpy arrays.
	xpt,ypt = m(lon,lat) 
	# convert back to lat/lon
	lonpt, latpt = m(xpt,ypt,inverse=True)
	m.plot(xpt,ypt,'r*', ms=16)  # plot a blue dot there
	# put some text next to the dot, offset a little bit
	# (the offset is in map projection coordinates)
	plt.text(xpt+100000,ypt+200000,'MRC, Banjul, The Gambia', size=14)
	
	m.shadedrelief()
	#plt.show()
	plt.savefig('testmap.png')
