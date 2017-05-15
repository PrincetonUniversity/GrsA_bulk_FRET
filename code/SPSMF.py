import os,shutil,sys
import re
import math
import time
import glob,string
import numpy as np
from scipy.stats import norm


##############################################################################	
	
class SimFrame(object):
	# Represents a single molecular conformation (frame) from the simulation
	def __init__(self, simID, frameID, coordSet, energy):
		self.simID = simID
		self.frameID = frameID
		self.energy = energy

		dx = coordSet[3] - coordSet[0]
		dy = coordSet[4] - coordSet[1]
		dz = coordSet[5] - coordSet[2]
		radicand = dx*dx + dy*dy + dz*dz
		self.d = math.sqrt(radicand)



class Sim(object):
	# Represents the collection of frames from all simulations	
	def __init__(self, coordFilesList, eFilesList):
		allSimFrames = []
		for i in range(len(coordFilesList)):
			coordList = convert(coordFilesList[i])
			eList = convert(eFilesList[i])
			for j in range(len(coordList)):
				# Some frames are missing, so skip over them
				if coordList[j] == []:
					continue
				sID = i 
				fID = j
				coords = coordList[j]
				e = eList[j]				
				newSimFrame = SimFrame(sID, fID, coords, e)
				allSimFrames.append(newSimFrame)
		self.allFrames = allSimFrames
	

##############################################################################	

class MovFrame(object):
	# Represents a frame from the actual trajectory
	def __init__(self, frameID, dt, d, uncertainty):
		self.frameID = frameID
		self.dt = dt
		self.d = d
		self.uncertainty = uncertainty
                self.num = 0
	def getGoodSimFrame(self, allSimFrames):
		# First find all simFrames within the acceptable range of d
		hi = self.d + self.uncertainty
		lo = self.d - self.uncertainty
		goodSimFrames = []
		for simFrame in allSimFrames:
			if simFrame.d >= lo and simFrame.d <= hi:
				goodSimFrames.append(simFrame)

		if len(goodSimFrames) == 0:
			return None
		
                self.num = len(goodSimFrames)
		# Find and return the frame with the lowest total energy
		lowestE = 1
		currentBest = None
		for item in goodSimFrames:
			if item.energy < lowestE:
				lowestE = item.energy
				currentBest = item
		
		return currentBest
			


class Mov(object):
	# Represents all frames in the trajectory movie
	def __init__(self, trajDataFilename):
		trajData = convert(trajDataFilename)
		allMovFrames = []
		for i in range(len(trajData)):
			step = trajData[i]
			dt = step[0]
			d = step[1]
			uncertainty = step[2]
			newMovFrame = MovFrame(i+1, dt, d, uncertainty)
			allMovFrames.append(newMovFrame)
		self.allFrames = allMovFrames

##############################################################################
# Write anything to a .txt file for easy future access
def write(txtFilename, List):
	txtFile = open(txtFilename, "w")
	txtFile.write(str(List))
	txtFile.close()

# convert .txt file back into its original data structure
def convert(filename):
	string = open(filename).read()
	item = eval(string)
	return item

##############################################################################
# Process the file containing the trajectory data
def processXC(xcFilename):
	fin = open(xcFilename)
	dataList = []
	for line in fin:
		a = line.split()
		float_a = [float(P) for P in a]
		dt = float_a[1] - float_a[0]
		dist = float_a[2] * 51		# Angstroms conversion factor
		#uncertainty = rel_uncertainty * dist
                uncertainty = float_a[3] * 51
		data = [dt, dist, uncertainty]
		dataList.append(data)
	return dataList

def runCommand(command):
    child = os.popen(command)
    data = child.read()
    err = child.close()
    if err:
        raise RuntimeError, '%s failed w/ exit code %d' % (command,err)
    return data

def getName(filename): #parse name.* format
    return filename.split('.')[0]

def getNamePath(filename): #parse the file name fron a ../name.* format
    split = filename.split('.')
    return split[2].split('/')[1]

def getCrystal(fname): #function that parses the crystal dimension from *-solvated.pdb
    header = open(fname,'r').readline()
    split = header.split(':')
    return split[1]

def runJob(filename, subs, runfile,depend=0):
    template = open(filename,'r').read()
    run = open(runfile,"r").read()
    currtemplate = template.format( **subs )
    subdict = dict(subs)
    fname = getNamePath(filename)+str(subdict['name'])
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
        return runCommand("qsub -W depend=afterok:"+ str(depend) + " " + outname)

def convertPDB(coordinates,pdb,output):
    #Takes the coordinates from an xyz file and puts them into a pdb file
    outfile = open(output,"w")
    coord = open(coordinates)
    #2 header lines in xyz
    coord.readline()
    coord.readline()
    pdbf = open(pdb)
    #Ignore remark lines in pdb
    pdbl = pdbf.readline()
    pdbls = pdbl.split()
    while pdbls[0] != 'ATOM':
        outfile.write(pdbl)
        pdbl = pdbf.readline()
        pdbls = pdbl.split()
    while pdbls[0] == 'ATOM':
        coordl = coord.readline().split()
        outfile.write(pdbl[0:26]+'      ')
        outfile.write("{0:6.3f}  {1:6.3f}  {2:6.3f}".format(float(coordl[1]),float(coordl[2]),float(coordl[3])))
        outfile.write("  {0} {1:5.2f}      {2}\n".format(pdbls[8],float(pdbls[9]),pdbls[10]))
        pdbl = pdbf.readline()
        pdbls= pdbl.split()
    outfile.close()

def getEner(log):
    #Extract energies with getprop
    ecount=0
    files=[]
    for i in log:
        eset=[]
        os.system("getprop ENER ../"+i)
        efile = open('ENER.DYNA.DAT','r')
        data = efile.readline().split()
        if ecount ==0:#The first line isn't a repeat
            eset.append(float(data[1]))
        data = efile.readline().split()
        while len(data) != 0:#read in the rest of the data
            eset.append(float(data[1]))
            data = efile.readline().split()
        name = "Efile"+str(ecount)+".txt"
        write(name,eset)
        files.append(name)
        ecount=ecount+1
    return files
 

def getCoor(num):
    #Parse correl files
    coorfile = []
    for i in range(num):
        dye1 = open('dye1_'+str(i)+'.txt','r')
        dye2 = open('dye2_'+str(i)+'.txt','r')
        for l in range(14):#discard the header
            dye1.readline()
            dye2.readline()
        coorset=[]
        coor1 = dye1.readline().split()
        coor2 = dye2.readline().split()
        while len(coor1) != 0:#read dye coords
            coorset.append([float(coor1[1]),float(coor1[2]),float(coor1[3]),float(coor2[1]),float(coor2[2]),float(coor2[3])])
            coor1 = dye1.readline().split()
            coor2 = dye2.readline().split()
        if len(coor2) !=0:
            print "Error: dye coordinate files are not the same length"
        name = "Coorfile"+str(i)+".txt"
        write(name,coorset)
        coorfile.append(name)
    return coorfile

def makeColvar(frmname, colvarf, colvarc=0, force=0):
    template = open(colvarf,'r').read()
    currcolv = template.format(colvar=colvarc, name=frmname, force=force)
    fname = getNamePath(colvarf)+'-'+frmname+'.inp'
    output = open(fname,'w')
    output.write(currcolv)


def detSolvent(Coorfiles,Efiles):
	simulation = Sim(Coorfiles, Efiles)
	simulationFrames = simulation.allFrames

	trajectory = Mov("traj.txt")
	trajFrames = trajectory.allFrames
        
	logFile = open("simulation.txt", "w")
	for sFrame in simulationFrames:
		logFile.write(str(sFrame.simID)+" "+str(sFrame.frameID)+" "+str(sFrame.d)+"\n")
       	logFile.close()

	#get info on trajectory
	tmax=0
	tmin=1000
	count=0
	usum=0
	txtFile = open("frameIDs.txt", "w")	
	for mFrame in trajFrames:
		if mFrame.d > tmax:
        		tmax=mFrame.d
                	tmaxFrame = mFrame
        	if mFrame.d < tmin:
                	tmin=mFrame.d
                	tminFrame = mFrame
        	count=count+1
        	usum=usum+mFrame.uncertainty
		good = mFrame.getGoodSimFrame(simulationFrames)
		if good == None:
			txtFile.write(str(mFrame.frameID) + " 0 0 0 " +\
               		str(mFrame.d)+ " " +  str(mFrame.uncertainty) +  "\n")
		else:
			txtFile.write(str(mFrame.frameID)+" "+str(good.simID)+\
                	" " +str(good.frameID) + " " + str(mFrame.num) + " " +\
                	str(mFrame.d)+" " + str(mFrame.uncertainty) + "\n")
	txtFile.close()

	#determine the explicit water simulations to run
        aveunc = usum/count
        step = aveunc/2.0
        start = tmin + step
        end = tmax - step
        span = end - start
        numpts = math.ceil(span/step)

        loc=start
        sampleList = []
        count = 0
        while count < numpts:
        	sampleList.append([0,loc,aveunc])
        	loc = loc + step
        	count=count+1
        write("sampling.txt", sampleList)

        toSimulate = Mov("sampling.txt")
        Samples = toSimulate.allFrames
	summFile = open("simstorun.txt", "w")
	simList = []
        for mFrame in Samples:
        	good = mFrame.getGoodSimFrame(simulationFrames)
        	if good == None:
                	print "ERROR!!!-can't find an appropriate structure in simulation"
            	else:
			written = str(good.simID)+" "+str(good.frameID) +"\n"
			simList.append(good)
                	summFile.write(written) 
	return simList    

def Prep(inputfiles,runfiles):
    #Set up directory structure and moves files up if needed
    if os.path.exists('simulationfiles/') == False:
            os.mkdir('simulationfiles/')
            shutil.copy('inputs/my.stdin','simulationfiles')
            shutil.copy('inputs/top_dye.inp','simulationfiles')
            shutil.copy('inputs/par_dye.inp','simulationfiles')
            #shutil.copy('inputs/water.crd','simulationfiles')
            #shutil.copy('inputs/addions.str','simulationfiles')
    for name in inputfiles:
        if os.path.exists(name) == False:
            shutil.copy('inputs/'+name,'.')
    for name in runfiles:
        if os.path.exists(name) == False:
            shutil.copy('runfiles/'+name,'.')
    os.chdir('simulationfiles')

def prepBins(log,traj,psf,pdb):
        #Get ready for simulations with bins
        shutil.copy('../'+psf,'.')
        shutil.copy('../'+pdb,'.')

    	#Extract energies with getprop
    	Efiles = getEner(log)
	num = len(Efiles)

	if os.path.isfile('correl0.log') == False:
	    #Run correl to get distances
            print "Running Correl"
            ccount=0
            for i in dcd:
                shutil.copy('../'+i,'.')
                subs = {'name':ccount,'dcd':i}
                runJob("../inputs/correl.inp",subs,"../local.run")
                ccount = ccount+1
	    if ccount != num:
                print "Error: number of coordinate and energy files do not agree"
            time.sleep(30)#should be quick

	#Extract coordinate files
	Coorfiles = getCoor(num)
	
	# Parse experimental Trajectory
	write("traj.txt",processXC('../'+traj))
	
	# Determine simulations to run
        return detSolvent(Coorfiles,Efiles)

def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[result.size/2:]

def writePlot(time,array,filename):
    outstr = ''
    for i in range(len(array)):
        outstr = outstr + str(time[i]) + "   " + str(array[i]) + "\n"
    outfi = open(filename,"w")
    #print filename
    outfi.write(outstr)
    outfi.close()

def writeGnuplot(fits,name):
    outstr = ''
    allTheLetters = string.ascii_lowercase
    for i in range(len(fits)): #write fits for gnuplot
        #outstr = outstr + allTheLetters[i] + "(x) =   " + str(fits[i][0])\
        #    + "* exp ( -" + str(fits[i][1]) + " * x ) \n"
        outstr = outstr + allTheLetters[i] + "(x) =   x**"+ str(fits[i][0])\
            + " * " + str(np.exp(fits[i][1])) + "  ) \n"
    outfi = open(name,"w")
    outfi.write(outstr)
    outfi.close()

def getPDB(fname,frnum,simdcd):
	#pull a pdb file out of a dcd simulation trajectory	
	convertCommand = "catdcd -o %s.pdb -otype pdb \
                                 -s feck.psf -stype psf \
                                 -first %d -last %d %s > outprodcd.log" \
                                 %(fname,frnum,frnum,simdcd)
        os.system(convertCommand)

#write histogram from correl distance file
def writeHist(dist):
        minv = 10000
        maxv = 0 #junk initialize
        #read from correl distance file
        infile=open('../'+str(d[0])+'_'+str(d[1])+'.txt','r')
        for l in range(14):#discard the header
            infile.readline()
        inline=infile.readline().split()
        data=[]
        while len(inline) != 0:
            inval = float(inline[1])#import from file
            data.append(inval)#import from file
            if float(inline[1])< minv:
                minv=float(inline[1])
            if float(inline[1])>maxv:
                maxv=float(inline[1])
            inline=infile.readline().split()
        data = sorted(data) #put in order
        length = len(data)
        #Determine bin size via Freedman and Diaconis
        iqr = data[3*length/4]-data[length/4]
        width = 2*iqr*length**(-1/3.0)
        rangev = maxv-minv
        nbins = int(math.ceil(rangev/width))
        #hist=[] #if needed for further analysis
        i=0
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
            outstring = outstring + str(maxvalue-midval) + '  ' +\
str(currcnt) + '\n'
        outfile=open('../'+str(d[0])+'_'+str(d[1])+'.hist','w')
        outfile.write(outstring)

##############################################################################
if __name__ == "__main__":
    # Bad hard coded variables
    #Should fix so program will check whether step1 ran and auto assign for
    #most cases
    #dcd = ["feck.dcd", "feck2.dcd"]
    #log = ["high_frict.log", "low_frict.log"]
    #traj = "apo_ak/ak.apo.43.m1.1323.xc"
    #psf = "feck.psf"
    #pdb = "feck_min.pdb"

    #Help
    if (len(sys.argv) == 1):
        print "Usage: SPSMF.py solvate protein_name replicas \
#Solvate a number of structures \"protein_name.replica.crd\""
        print "SPSMF.py neutralize protein_name replicas \
#Neutralize a number of simulations (step 2)"
        print "SPSMF.py md protein_name replicas \
#Run NAMD dynamics for a number of solvated replicas (step 3)"
        print "SPSMF.py restart protein_name replicas #Restart NAMD (step 4)"
        print "SPSMF.py combine_colvars protein_namd replicas \
#Combine colvar traj files in /colvars"
        print "SPSMF.py initial_colvar protein_name replicas \
#Print out the initial value of colvars for a set of structure"
        print "SPSMF.py setup protein_name (pdb file name without extension)\
offset (optional-first resi in pdb file)\
#Sets up a pdb file for charmm simulation, generates protein_name_min.pdb"
        print "SPSMF.py setup_dye protein_name (pdb file name without\
extension) offset (optional-first resi in pdb file)\
#Sets up a pdb file for charmm simulation with dyes, generates\
protein_name_min.pdb"
        print "SPSMF.py sgld protein_name offset #Runs an initial SGLD\
simulation (step 2)"
        print "SPSMF.py restart_sgld protein_name offset \
#Restarts a SGLD Simulation(step 3)"
        print "SPSMF.py dist protein_name offset \
#Gets dye-dye distance"
        print "SPSMF.py sampling protein_name offset (required) dyes \
(list of all dyes)a#Run sampling for dyes"
        print "SPSMF.py bins solvate #Solvate all frames of trajectory"
        print "SPSMF.py bins md #Run md simulations (step 2)"
        print "SPSMF.py bins restart #Restart md simulation (step 3 if\
desired)"
        print "SPSMF.py bins plot #Plot distnces (step 4 or 3)"
        print "SPSMF.py bins restraint #Run SGLD simulations with restraints h"
        print "SPSMF.py bins dist #gets distnce measurements"
        print "SPSMF.py convert xyz pdb output #Puts xyz coordinates into a \
pdb template-Files must be located in execuation directory"
        exit(1)


    if (sys.argv[1] == 'convert'):
        if (len(sys.argv) != 5):
            print "You need to provide the xyz, pdb template, and output file name"
            exit(1)        
        xyz = sys.argv[2]
        pdb = sys.argv[3]
        out = sys.argv[4]
        convertPDB(xyz,pdb,out)

    elif (sys.argv[1] == 'setup'):
        Prep(["prepprot.inp"],["local.run"])
        if (sys.argv[2] == ''):
            print "The pdb file name (without extension) is required after setup"
            exit(1)
        else:
            if (sys.argv[3] == ''):#get offset
                offset = 0
            else:
                offset = sys.argv[3]
            fname = sys.argv[2]
            runfile = "../local.run"
            subs = {'name':fname}
            runJob("../prepprot.inp",subs,runfile)
            print "Submitted"

    elif (sys.argv[1] == 'setup_dye'):
        Prep(["step1.inp"],["cluster.run"])
        if (sys.argv[2] == ''):
            print "The pdb file name (without extension) is required \
                    after setup_dye"
            exit(1)
        else:
            if (len(sys.argv)== 3):
                offset = 0
            else:
                offset = sys.argv[3]
            fname = sys.argv[2]
            shutil.copy('../'+fname+'.pdb','.')
            runfile = "../cluster.run"
            subs = {'name':fname,'offset':offset}
            runJob("../step1.inp",subs,runfile)

    elif (sys.argv[1] == 'sampling'):
        Prep(["sampling.inp"],["cluster.run"])
        if (len(sys.argv) < 5):
            print "The pdb file name (without extension), offest \
                    ,and at least one dye  number is required"
            exit(1)
        else:
            offset = sys.argv[3]
            count = 4#get dye listing
            dyes = []
            while count < len(sys.argv):
                dye = sys.argv[count]
                dyes.append(dye)
                count = count + 1
            segi = sys.argv[2]
            prot = sys.argv[2]
            runfile = "../cluster.run"
            for f in range(len(dyes)):
                fname = prot+str(dyes[f])
                subs =\
{'name':fname,'offset':offset,'dye':dyes[f],'prot':prot,'segi':segi}
                runJob("../sampling.inp",subs,runfile)
                print "Sampling for dye " + str(dyes[f]) + " submitted!"

    elif (sys.argv[1] == 'sgld'):
        Prep(["start.inp"],["cluster.run"])
        if (sys.argv[2] == ''):
            print "The pdb file name (without extension) is required after sgld"
            exit(1)
        else:
            fname = sys.argv[2]
            if (os.path.isfile('../'+fname+'_min.pdb') == False):
                print "The sgld pdb file must be set up first"
                exit(1)
            else:
                if (sys.argv[3] == ''):#get offset
                    offset = 0
                else:
                    offset = sys.argv[3]
                runfile = "../cluster.run"
                subs = {'name':fname}
                output = runJob("../start.inp",subs,runfile)
                depend = getName(output)
                print "Submitted job for " + fname

    elif (sys.argv[1] == 'restart_sgld'):
        Prep(["restart3.inp"],["cluster.run"])
        #need to automate numbering
        if (sys.argv[3] == ''):#get offset
            offset = 0
        else:
            offset = sys.argv[3]
        fname=sys.argv[2]
        runfile = "../cluster.run"
        subs = {'name':fname}
        runJob("../restart3.inp",subs,runfile)
        print "Submitted job for " + fname

    elif (sys.argv[1] == 'dist'):
        Prep(["dist.inp"],["local.run"])
        fname=sys.argv[2]
        runfile = "../local.run"
        subs = {'name':fname}
        #runJob("../dist.inp",subs,runfile)
        #print "Submitted job for " + fname
        #time.sleep(15)
        dists = [[152,492],[152,508]] #could be generalized
        for d in dists:
            writeHist(d)

    elif ( sys.argv[1] == 'solvate'):
        Prep()
        name = sys.argv[2]
        num = int(sys.argv[3])
        i = 1
        if (sys.argv[4] == ''):#get offset
            offset = 0
        else:
            offset = sys.argv[4]
        while i < num+1:
            #Run the charmm calculation
            fname = name + '.' + str(i)
            runfile = "../cluster.run"
            subs = {'name':fname,'psf':name}
            output = runJob("../solvate.inp",subs,runfile)
            depend = getName(output)
            output = runJob("../neutralize.inp",subs,runfile,depend)
            depend = getName(output)
            output = runJob("../final-minimization.inp",subs,runfile,depend)
            print "Submitted" + str(i)
            i = i + 1

    elif ( sys.argv[1] == 'neutralize'):
        Prep()
        name = sys.argv[2]
        num = int(sys.argv[3])
        if (sys.argv[4] == ''):#get offset
            offset = 0
        else:
            offset = sys.argv[4]
        i = 1
        while i < num+1:
            fname = name + '.' + str(i)
            if os.path.isfile(fname+'-solvated.pdb') == False:
                print "Solvation must finish before continuing"
                exit(1)
            #crystal = getCrystal(fname+'-solvated.pdb')
            #Run the charmm calculation
            runfile = "../cluster.run"
            subs = {'name':fname,'psf':name}
            output = runJob("../neutralize.inp",subs,runfile)
            depend = getName(output)
            output = runJob("../final-minimization.inp",subs,runfile,depend)
            print "Submitted" + str(i)
            i = i + 1

    elif (sys.argv[1] == 'md'):
        Prep()
        name = sys.argv[2]
        num = int(sys.argv[3])
        if (sys.argv[4] == ''):#get offset
            offset = 0
        else:
            offset = sys.argv[4]
        i = 1
        if os.path.isfile('../colvar_input') == False:
                print "You need to make an input file for colvar with"
                print "center, jump, restraint weight as colvar_input"
                exit(1)
        infile =  open('../colvar_input','r')
        startc = float(infile.readline().split()[1])
        jump = float(infile.readline().split()[1])
        force = float(infile.readline().split()[1])
        #could be improved via checking of column one labels


        currc = startc
        colvarf = "../colvars.inp"
        while i < num+1:
            print currc
            fname = name + '.' + str(i)
            runfile = "../namd.run"
            makeColvar(fname,colvarf,currc,force)
            subs = {'name':fname,'psf':name}
            output = runJob("../namd_heat.inp",subs,runfile)
            depend = getName(output)
            output = runJob("../namd.inp",subs,runfile,depend)
            depend = getName(output)
            output = runJob("../namd_restart.inp",subs,runfile,depend)
            print "Submitted" + str(i)
            currc = currc + jump
            i = i + 1

    elif (sys.argv[1] == 'restart'):
        Prep()
        name = sys.argv[2]
        num = int(sys.argv[3])
        if (sys.argv[4] == ''):#get offset
            offset = 0
        else:
            offset = sys.argv[4]
        i = 1
        while i < num+1:
            fname = name + '.' + str(i)
            runfile = "../namd.run"
            subs = {'name':fname,'psf':name}
            output = runJob("../namd_restart9.inp",subs,runfile)
            print "Submitted" + str(i)
            i = i + 1

    elif (sys.argv[1] == 'pmf'):
        name = sys.argv[2]
        num = int(sys.argv[3])
        #num = 1
        if os.path.isfile('colvar_input') == False:
                print "You need to make an input file for colvar with"
                print "center, jump, restraint weight as colvar_input"
                exit(1)
        infile =  open('colvar_input','r')
        startc = float(infile.readline().split()[1])
        jump = float(infile.readline().split()[1])
        force = float(infile.readline().split()[1])
        #could be improved via checking of column one labels

        os.chdir('colvars')
        i = 0
        forces = []
        restrval = startc

        while i < num:#for every replica
            #prep rep
            forces.append([])
            maxval = 0
            minval = 10000
            pmfsum = 0
            count = 0
            filein = open('{0}.{1}.colvars.traj'.format(name,i+1),"r")
            #read in data
            linein = filein.readline()
            while linein != '':#for every frame
                val = float(linein.split()[1])
                if val > maxval:
                    maxval = val
                if val < minval:
                    minval = val
                mf =  -force*(val-restrval)
                forces[i].append(mf)
                pmfsum = pmfsum + mf
                count = count + 1
                linein = filein.readline()
            restrval = restrval + jump
            print i
            i=i+1

        print "data in"

        #avgforce = []#contains each force in the bin
        meanf = []#contains the mean of each bin
        uncertainty = []#contains the sigma for each bin
        mus = []#contains mu for each bina
        muvals = []
        sigvals = []
        #sigvals.append([]*num)
        i = 0
        binlens =\
np.array([50,500,1000,5000,10000,50000,100000,150000,200000,250000,300000,350000,400000,450000,500000])
        for f in range(len(binlens)):
            histf = []
            avgforce = []#contains each force in the bin
            for i in range(num):#for each replica
                sigvals.append([])
                histf.append([0]*800)
                avgforce.append([])
                begin = 0
                end = binlens[f]
                sumtot = 0
                count = 0
                #numbins = len(forces[i])/timavg
                numbins = int(math.floor(len(forces[i]) / binlens[f]))
                print numbins
                for j in range(numbins): #for each bin calc average force
                #while count < len(forces[i]): #for entire trajectory
                    sumint = 0
                    #print end
                    #count = begin
                    while count < end:
                        sumint = sumint + forces[i][count]
                        count = count + 1
                    binavg = sumint/binlens[f]
                    avgforce[i].append(binavg)
                    sumtot = sumtot + binavg
                    point = int(math.floor(binavg * 10.0) + 400)
                    histf[i][point] = histf[i][point] + 1
                    end = end + binlens[f]
                    #begin = begin + 1
                xvals = np.arange(-40,40,0.1)
                writePlot(xvals,histf[i],"Hist_Rep"+str(i))
                meanf.append(sumtot/numbins)
                mu,sigma = norm.fit(avgforce[i])
                #muvals.append(mu)
                sigvals[i].append(sigma)
                uncertainty.append(sigma/(numbins**0.5))
        fits = []
        extunc = [] #extrapolated uncertainty in force
        for i in range(num):
            outstr = ''
            fitx = np.array([])
            fity = np.array([])
            for f in range(len(binlens)):
                outstr = outstr + str(binlens[f]) + "  " + str(sigvals[i][f])\
+" \n"
            outfi = open("convergence"+str(i),"w")
            outfi.write(outstr)
            outfi.close()
            fity =\
np.array([sigvals[i][0],sigvals[i][1],sigvals[i][2],sigvals[i][3],sigvals[i][4],sigvals[i][5],sigvals[i][6],\
sigvals[i][7],sigvals[i][8],sigvals[i][9],sigvals[i][10],sigvals[i][11],sigvals[i][12],sigvals[i][13],sigvals[i][14]])
            #print binlens
            #print fity
            #A = np.vstack([binlens**0.5, np.ones(len(binlens))]).T
            #m,c = np.linalg.lstsq(A, fity)[0]
            popt = np.polyfit(np.log(binlens),np.log(fity),1)
            #popt,pcov =curve_fit(fitFunct,binlens,fity))
            #print popt
            fits.append(popt)
            extunc.append(3500000**popt[0] * np.exp(popt[1]))
            #!!!!Hard coded 7 ns extrapolation
        restrval = startc
        print "Replica PMF"
        print str(restrval) + "       0"
        #now integrate and get profile
        i = 1
        pmf = 0
        errors = 0
        histp = []
        pmfe = []##contains binned pmf values
        while i < num:#for every replica
            pmfe.append([])
            histp.append([0]*800)
            restrval = restrval + jump
            pmfv = 0.5 * (meanf[i-1] + meanf[i]) * jump
            pmf = pmf + pmfv
            properror = 0.5 * (extunc[i-1]**2 + extunc[i]**2)**0.5 * jump
            #print extunc[i-1], extunc[i]
            #print str(properror) + "delta error"
            errors = (errors**2 + properror**2)**0.5
            print str(restrval) + '    ' +  str(pmf) + '    ' + str(errors)
            #writePlot(xvals,histp[i-2],"HistPmf_Rep"+str(i-2))
            i = i + 1
        #writePlot(muvals,sigvals,"FinalHistFitting")

    elif (sys.argv[1] == 'combine_trajectory'):
        name = sys.argv[2]
        num = int(sys.argv[3])
        i = 1
        outstr = ''
        count = 0
        while i < num+1:
            infi =\
open('simulationfiles/{0}.{1}-namd-2.colvars.traj'.format(name,i),'r')
            line = infi.readline()
            while line != '':
                split = line.split()
                if split[0] == '#': #trash comments
                    line = infi.readline()
                else:
                    count = count + 1
                    if count%1000 == 0:
                        outstr = outstr + str(count) + "   " + str(split[1]) \
+ "\n"
                    line = infi.readline()
            print i
            i = i + 1
        outfi = open('summary_trajectory',"w")
        outfi.write(outstr)
        outfi.close()

    elif (sys.argv[1] == 'combine_colvars'):
        name = sys.argv[2]
        num = int(sys.argv[3])
        os.chdir('colvars')
        i = 1
        while i < num+1:
            count = 0
            outstr = ''#find files that match the colvar format
            files = glob.glob('{0}.{1}-namd*.colvars.traj'.format(name,i))
            files = sorted(files,key=lambda\
x:int(x.split('.')[1].split('-')[0]))
#because I couldn't use the same punctuation...
            for j in range(len(files)):
                infi = open(files[j],"r")
                line = infi.readline() #trash header
                #line = infi.readline() #trash first line 
                #(repeated from previous colvar)
                line = infi.readline()
                while line != '':
                    split = line.split()
                    if split[0] == '#': #trash comments
                        line = infi.readline()
                    else:
                        count = count +1
                        outstr = outstr + line
                        line = infi.readline()
            outfi = open('{0}.{1}.colvars.traj'.format(name,i),"w")
            outfi.write(outstr)
            outfi.close()
            print i, count
            i = i+1

    elif (sys.argv[1] == 'initial_colvar'):
        Prep()
        name = sys.argv[2]
        num = int(sys.argv[3])
        i = 1
        colvarf = "../colvars-ini.inp"
        while i < num+1:
            fname = name + '.' + str(i)
            runfile = "../namd.run"
            makeColvar(fname,colvarf)
            subs = {'name':fname,'psf':name}
            output = runJob("../namd_ini.inp",subs,runfile)
            print "Submitted" + str(i)
            i = i + 1

    elif (sys.argv[1] == 'get_initial'):
        name = sys.argv[2]
        num = int(sys.argv[3])
        i = 1
        outstr1 = ''
        outstr2 = ''
        points = []
        while i < num+1:
            infi =\
open('simulationfiles/{0}.{1}-namd-ini.colvars.traj'.format(name,i),"r")
            infi.readline()#kill comment
            indat = infi.readline().split()
            outstr1 = outstr1 + str(i) + "   " + indat[1] + "\n"
            i = i + 1
        outfi = open("colvars_summ","w")
        outfi.write(outstr1)
        outfi.close()


    elif (sys.argv[1] == 'bins'):
        Prep()
	# Run simulations in bins
        #Parse input files and figure out what simulations we're running
        toRun = prepBins(log,traj,psf,pdb)
        #Get results from prep
        toSimulate = Mov("sampling.txt")
        uncertainty = toSimulate.allFrames[0].uncertainty #average uncertainty

        maxnum = 1 # control variables
        start = 0
        end = len(toRun)
        if end > maxnum:
            end = start+maxnum
        sim = start
        while sim < end:
            #Get info on the simulation to run
            simnum=toRun[sim].simID
            frnum=toRun[sim].frameID
            simdcd=dcd[simnum]
            name = getName(simdcd)
            fname = name+'_'+str(frnum)
            distance=toRun[sim].d
            #Check for duplicates with diff distance
            temp=fname
            duplicates=0
            for i in range(sim):
                if simnum==toRun[i].simID and frnum==toRun[i].frameID:
                    duplicates=duplicates+1
                    fname = temp + '_' + str(duplicates)
	    #run actual simulations, depending on command line args 
	    if ( sys.argv[2] == 'solvate'):
            	#Pull the starting structure out of the dcd file
            	getPDB(fname,frnum,simdcd)
            	#Run the charmm calculations
                runfile = "../local.run"
            	subs = {'name':fname}
            	runJob("../solvate.inp",subs,runfile)
        
	    if ( sys.argv[2] == 'md'):
            	#Solvation must be finished before next step
            	if os.path.isfile(fname+'-solvated.pdb') == False:
                	print "Solvation must finish before continuing"
                	exit(1)
            	runfile = "../cluster.run"
            	crystal = getCrystal(fname+'-solvated.pdb')
            	subs ={'name':fname,'crystal':crystal}
            	output = runJob("../neutralize.inp",subs,runfile)
            	depend = getName(output)
            	subs = {'name':fname,'distance':distance,'resval':'0.5','crystal':crystal}
            	output = runJob("../final-minimization.inp",subs,runfile,depend)
            	depend = getName(output)
            	output = runJob("../md-heat.inp",subs,runfile,depend)
            	depend = getName(output)
            	output = runJob("../md-dyna.inp",subs,runfile,depend)
            	depend = getName(output)
            	output = runJob("../md-dyna2.inp",subs,runfile,depend)
                subs = {'name':fname,'crystal':crystal,'dcd':fname+'-dyna-2.dcd'}
            	depend = getName(output)
            	runJob("../get_dist.inp",subs,runfile,depend)

            if ( sys.argv[2] == 'restart'):
                runfile = "../cluster.run"
                crystal = getCrystal(fname+'-solvated.pdb')
                subs = {'name':fname,'distance':distance,'crystal':crystal}
                runJob("../md-dyna3.inp",subs,runfile)

	    if ( sys.argv[2] == 'restraint'):
		getPDB(fname,frnum,simdcd)
		runfile = "../cluster.run"
		subs = {'name':fname,'distance':distance,'resval':'0.5'}
		runJob("../sgld_restraint.inp",subs,runfile)
            
            if (sys.argv[2] == 'bin_dist'):
                runfile = "../local.run"
                subs = {'name':fname,'psf':psf,'dcd':fname+'.dcd'}
                runJob("../get_dist.inp",subs,runfile)
    
            if (sys.argv[2] == 'plot'):
                dname = 'dist_'+fname+'.txt'
                subs={'file':dname,'output':'dist_'+fname+'.pdf','stdev':uncertainty/2.0,'mean':distance}
                writeHist()#!!!Needs to be changed to new implementation
            else:
                print "Huh?"
            sim= sim+1
            print "Simulation " + fname + " submitted!"
    else:
        print "I don't understand what you want me to do"
