import os,shutil,sys,math

class Dye(object):
        # Represents a single dye orientation
        def __init__(self, xin, yin,zin):
            self.x = float(xin)
            self.y = float(yin)
            self.z = float(zin)

        def __repr__(self): #return a point as a string
            return str(self.x)+" "+str(self.y)+" "+str(self.z)


class Sim(object):
	# Represents an entire sampling trajectory
	def __init__(self, coorfile,struct):#Reads data and prints coordinates appended to pdb
		frames = []
                self.ene = []
		fnum = 0
		vect = open(coorfile,'r')
                dyenum=coorfile.split('.')[1]
                if dyenum == 'egfp':
                    dyenum = 999
                if dyenum == 'coords':
                    dyenum = 152
                #ener = open(coorfile.split('.')[0]+'.ener.txt')
                for l in range(14):#discard the header
                        vect.readline()
                        #ener.readline()
                splitv = vect.readline().split()
                #splite = ener.readline().split()
                xsum = 0
                ysum = 0
                zsum = 0
                while len(splitv) != 0:#read all lines
                        xcoor = float(splitv[1])
                        xsum = xsum + xcoor
                        ycoor = float(splitv[2])
                        ysum = ysum + ycoor
                        zcoor = float(splitv[3])
                        zsum = zsum + zcoor
                        if xcoor != 0:
                                struct.write("ATOM  15633  CES CES  9")
                                struct.write("{0:d}     {1:7.3f} {2:7.3f}\
 {3:7.3f}".format(int(dyenum),xcoor,ycoor,zcoor))#mod for gfp
                                struct.write("  1.00  1.70      2G47 \n")
				frame=Dye(xcoor,ycoor,zcoor)
				frames.append(frame)
                                #self.ene.append(splite[1])
                                fnum = fnum + 1
                        splitv = vect.readline().split()
                        #splite = ener.readline().split()
		self.allFrames=frames
                self.ave=Dye(xsum/fnum,ysum/fnum,zsum/fnum)
                self.fnum=fnum

def Dist(coor1,coor2):
	d1 = coor1.x-coor2.x
	d2 = coor1.y-coor2.y
	d3 = coor1.z-coor2.z
	radicand = d1*d1 + d2*d2 + d3*d3
        return math.sqrt(radicand)

def runCharmm(filename, subs, runfile,depend=0):
    template = open(filename,'r').read()
    run = open(runfile,"r").read()
    currtemplate = template.format( **subs )
    subdict = dict(subs)
    fname = getName(filename)+str(subdict['name'])
    currrun = run.format( filename = fname)
    output = open(fname+'.inp','w')
    outname = fname+".run"
    outrun = open(outname,"w")
    output.write(currtemplate)
    outrun.write(currrun)
    output.close()
    outrun.close()
    if depend == 0:
        return runCommand("qsub "+ outname )
    else:
        return runCommand("qsub -W depend=afterok:"+ str(depend) + " " +\
                            outname)

def getName(filename): #parse name.* format
    return filename.split('.')[0]


def runCommand(command):
    child = os.popen(command)
    data = child.read()
    err = child.close()
    if err:
        raise RuntimeError, '%s failed w/ exit code %d' % (command,err)
    return data


##############################################################################
if __name__ == "__main__":
	infile = open(sys.argv[1],'r')
	filenames = infile.readline().split()
	numfiles = len(filenames)
        pdb = sys.argv[2]
	outname = pdb+'_point.pdb'
	#Help
	if (len(sys.argv) == 1):
		print "Usage: CorrelFile PDB_file ; File format:CorrelFile repeat.... \
                       Distance pairs on following lines; Remember to remove last two lines from PDB"

	#shutil.copy(pdb+'_min.pdb', outname)
	#struct = open(outname,'a')
        struct = open(outname,'w')
	darray = [0]*numfiles
	count = 0
        goodframe = [0]*numfiles
        print "Distance, Average Coordinates "
	while count < numfiles:
		filename = filenames[count]
		darray[count] = Sim(filename,struct)
                print filename + '  ' + str(darray[count].ave) 
                #search for structure that best fits the average
                frames = darray[count].allFrames
                #energy = darray[count].ene
                frameD = 1000000000.0 #Terrible
                for i in range(len(frames)):
                    #print Dist(frames[i],darray[count].ave)
                    if Dist(frames[i],darray[count].ave) < frameD:
                        #if float(energy[i]) < float(frameE):
                            goodframe[count] = i+1 
                            #offset by one for frame number in dcd
                            #frameE = energy[i]
                            frameD = Dist(frames[i],darray[count].ave)
                print "Best frame: " + str(goodframe[count])
		count = count + 1
	struct.close()
        
        #Run charmm script to pick out the appropriate dye coordinates
        runfile = "cluster.run"
        subs = {'frm1':goodframe[0],\
                'name':pdb}
        runCharmm("prep_linker3.inp",subs,runfile)


	#Calculate distances
	indist = infile.readline().split()
	while  len(indist) !=0:#read all lines
		d1 = int(indist[0])-1
		d2 = int(indist[1])-1
		outfile = "Distances over frames \n"
		totd = 0
		frmd = 0
		frmm = 0
                cnt = 0
		data = []
                minv = 10000
                maxv = 0 #junk initialize
		#histogram = [0]*500 #bin 1-500 A in 1 A segments
		name = "dye{0}_vs_dye{1}".format(d1,d2)
		outf = open(name,'w')
		while frmd < darray[d1].fnum:
			while frmm < darray[d2].fnum:
				mdist = Dist(darray[d1].allFrames[frmd],darray[d2].allFrames[frmm]) 
				totd = totd + mdist
				data.append(mdist)
				#histogram[binn] = histogram[binn] + 1
                                if float(mdist)< minv:
                                    minv=float(mdist)
                                if float(mdist)>maxv:
                                    maxv=float(mdist)
				frmm = frmm + 1
				cnt = cnt + 1
			frmd = frmd + 1
			frmm = 0
                data = sorted(data) #put in order
                length = len(data)
                mean = sum(data)/float(length)
                dev=0
                for val in data:
                    dev = dev + (val-mean)**2
                stdev = (dev/length)**0.5
                #Determine bin size via Freedman and Diaconis
                iqr = data[3*length/4]-data[length/4]
                width = 2*iqr*length**(-1/3.0)
                rangev = maxv-minv
                nbins = int(math.ceil(rangev/width))
                #hist=[] #if needed for further analysis
                i=0
                temp=0
                midval=width/2.0
                outstring=''
                for n in range(nbins):
                    maxvalue = minv + width*n
                    currcnt=0
                    datapt = data[i]
                    while datapt < maxvalue:
                        currcnt=currcnt+1
                        i=i+1
                        datapt = data[i]
                    #hist.append(currcnt)
                    temp=temp+currcnt/float(length)
                    outstring = outstring + str(maxvalue-midval) + '  ' +\
str((currcnt/float(length*width))) + '\n'
                outf.write(outstring)
                print "count: " + str(length)
                print "mean: " + str(mean)
                print "stdev: " + str(stdev)
                print temp

		#print d1 
		#print d2
		#print totd/cnt 
		indist = infile.readline().split()
