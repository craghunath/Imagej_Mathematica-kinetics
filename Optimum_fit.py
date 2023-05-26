from dataclasses import dataclass, field
from typing import List #
# import matplotlib as mpl
# mpl.use("pgf")
import matplotlib.pyplot as plt
# plt.style.use('ggplot')
# plt.rcParams.update({'axes.facecolor':'azure'})
plt.rcParams.update({ 'pgf.texsystem':"pdflatex"})
from numpy import array, mean, std, arange

# from re import sub, split, findall
# from os import listdir
# import matplotlib.pyplot as plt
# from numpy import arange, array, diff, sqrt, append, zeros,  mean as nmean  #, abs, where ,convolve, ones, split as nsplit,
# from scipy.integrate import odeint
# from scipy.optimize import curve_fit 

@dataclass
class Tbd:
	"""
	These are essentials required for computing Two boat dynamics
		ppcm : pixcels per centimeter, so the entire experimental data of distance in cm and time in s
		fpath : file or folder path
		frames : number of frames per second at which the dynamics video is recorded
	"""
	cordNvel: tuple = field(init=False)
	ppcm: float = 10
	frames: float = 10
	boatMass: float = 5.6*10**-3
	fpath: str = '000'
	# folderpath: str = '000'
	def __post_init__(self):

		self.cordNvel = self.dataexp(self.fpath)
		# self.optm = ()

	def dataexp(self, filepath:str):
		from re import sub, split
		from numpy import array, diff, sqrt, append, zeros  #,  mean as nmean
		# filepath = self.fpath
		info = open(filepath).read()
		datStr = list(filter(None, split ("\n|,",sub(".+Y\n|.+,\d{2,3},","",info))))
		end = len(datStr)
		
		btxy = array( list(map(float, datStr)) )
		
		btxyDiff = btxy[int(end/2):] - btxy[0:int(end/2)]
		
		btx = btxyDiff[0:end-1:2]
		bty = btxyDiff[1:end:2]
		btr = sqrt(btx**2 + bty**2)
		ppcm = self.ppcm
		cords = btr/ppcm #__in cm #btr*0.01/ppcm __in m   #diff(btr)[0:rend-2:2]*0.01/47
		
		tSlice = 1/self.frames
		vel = diff(cords)/tSlice
		vels = append(zeros(1),vel)
		return cords,vels


	
	def modelFun(self, t, y, mbot:float, km2:float, n1:float, n2:float):
		"""
		'odeint' function for numerical solution needs the above kind of structure of the function to process
		t : independent cordinate array for the solution
		y : tuple or list of all the intial conditions, here r0, and rdot0
		"""
		r, z = y
		zdot = [z, ( km2/(r**4)-( n1*z + n2*z**2) )/mbot]
		return zdot

	def modelFun3by2(self, t, y, mbot:float, km2:float, n2:float, n3:float):
		"""
		'odeint' function for numerical solution needs the above kind of structure of the function to process
		t : independent cordinate array for the solution
		y : tuple or list of all the intial conditions, here r0, and rdot0
		"""
		r, z = y
		zdot = [z, ( km2/(r**4)-(  n2*z**2 + n3*z**1.5) )/mbot]
		return zdot
	
	def modelSolver(self, tim, km2:float, n1:float, n2:float,  filepath:str=fpath, flag:str='fit'):
		"""
		'flag' has two options default 'fit' or 'normal', either to output only velocities or n x 2 numpy array of [cordinates, velocities]
			-- output of only velocities is the necesscessity for the 'curve_fit' function
		cm is distance unit, but the initial r0 for the solver need to be fed in 'm'
		OUTPUT :: also in cm, therefore 100 multiplied in num_soln
		'xcords' in the function is dependent on the outside evaluation of experimental data 
		"""
		from scipy.integrate import odeint
		flp:str = filepath
		if filepath == '000':
			xcords, vels = self.cordNvel
		else:
			xcords, vels = self.dataexp(filepath=flp)
		

		r0, z0 = xcords[0]*10**-2, 0.0
		y0 = [r0, z0]
		mbot = self.boatMass   #7.16*10**-3 #5.6  7.16
		num_soln = 100*odeint(self.modelFun, t=tim, y0=y0,  args=(mbot, km2*10**-9, n1, n2), tfirst=True) # , tfirst=True
		rtn = {'normal': num_soln, 'fit' : num_soln[:,1]}
		return rtn[flag]

	def dataFit(self, filepath:str=fpath, fit:str='actual', duration:float=0):
		"""
		'fit' has two options default 'actual' or 'movavg', either to fit actual or moving averaged data
		use 'duration'= 5 s for fitting about 15 cm of the path covered (avg ~3 cm/s is the spped of the swimmer)
		"""		
		from scipy.optimize import curve_fit
		from numpy import arange, convolve, ones

		flp:str = filepath
		if filepath == '000':
			_, vels = self.cordNvel
		else:		
			_, vels = self.dataexp(filepath=flp)
		
		tInc = 1/self.frames
		if duration != 0:
			tStop: float = duration
		else:
			tStop: float = len(vels)*tInc
		"""
		-- check the validity of all point solver with tStop
		-- it is better to evaluate fit for first 15 seconds
		"""
		pts = int(tStop*self.frames) + 1
		tdata = arange(start= 0.0, stop=tStop, step=tInc)

		def mav():
			vavg = convolve(vels, ones(10)/10, mode='same')
			return vavg[:pts]

		vdat = {'actual':vels[:pts], 'movavg':mav()}

		popt, _ = curve_fit(self.modelSolver, xdata=tdata, ydata=vdat[fit], p0=(0.2,0.0001,0.0001), bounds=((0,0,0), (2,0.1,0.1)))

		return popt

	def movAvg(self, filepath:str=fpath):
		from numpy import arange, convolve, ones

		flp:str = filepath
		if filepath == '000':
			_, v = self.cordNvel
		else:
			_, v = self.dataexp(filepath=flp)

		vavg = convolve(v, ones(10)/10, mode='same')
		return vavg
	
	def plotEach(self, savefolder:str, optink, filepath:str=fpath):

		plt.rc('figure', max_open_warning = 0)
		plt.style.use('ggplot')
		plt.rcParams.update({'axes.facecolor':'oldlace'})
		from numpy import arange,  where , array
		from re import sub, findall

		flp:str = filepath
		if filepath == '000':
			cords, vels = self.cordNvel
		else:
			cords, vels = self.dataexp(filepath=flp)

		bm = round(10**3 * self.boatMass, 3)

		match = findall(pattern=r"/\d+_", string=self.fpath)
		fno = sub(pattern=r"/",repl="",string=match[0], count=1)


		tInc = 1/self.frames
		tStop: float = len(vels)*tInc

		"""
		-- check the validity of all point solver with tStop
		-- it is better to evaluate fit for first 5 seconds
		"""
		pts = int(tStop*self.frames) +1
		tdata = arange(start= 0.0, stop=tStop, step=tInc)
		soln = self.modelSolver(tim=tdata, km2=optink[0], n1=optink[1], n2=optink[2],  flag='normal')
		# fratio = list(map(lambda r,z :   (optink[2]*10**-9 )/ (r**4)/ (optink[0]*z + optink[1]*z**2 )  ,soln[:,0]*10**-2, soln[:,1]*10**-2 ))
		
		#  (optink[2]*10**-9/(r**4) )
		"""Individual fit block"""
		plt.figure()
		plt.plot(tdata, vels, "+", color='indigo', markersize=0.8) #arange(start= 0.0, stop=tStop, step=tInc)
		plt.plot(tdata,soln[:,1], "b-", linewidth=1.1) #tdata
		# plt.plot(tdata[:240],fratio[:240]) #(tdata[25:49],fratio[24:48])
		plt.grid(color='midnightblue', linestyle='-', linewidth=1, alpha=0.85)
		fln = "".join((savefolder,'n12k1.4chk_',fno,"ppcm",str(self.ppcm),"_mbot",str(bm),".jpg"))
		plt.xlabel('time (s)')
		plt.ylabel('Speed (cm/s)')
		plt.savefig(fln, dpi=300)


		return fno, tdata[where(max(vels)==array(vels) )[0][0]] 




def fileList(path:str)-> List[str]:

	from os import listdir
	from re import  findall

	flnames = list(filter(lambda strn: findall('(?=.*)(?=.*\.csv)', strn), listdir(path)));
	return list(map(lambda fl: ''.join([path, fl]), flnames))

def chisq(instance):

	from numpy import arange, sum, where
	from re import findall, sub
	_, vels = instance.cordNvel
	optm  = instance.dataFit()
	fps = instance.frames
	tstop:float = len(vels)/fps
	tdat = arange(start=0.0,stop=tstop,step=(1/fps))
	fitvels = instance.modelSolver(tim=tdat, n2=optm[0], n3=optm[1], km2=optm[2])
	dfr = fitvels[1:len(vels)] - vels[1:]
	chis = dfr*dfr/fitvels[1:len(vels)]
		
	match = findall(pattern=r"/\d+_", string=instance.fpath)
	fno = sub(pattern=r"/",repl="",string=match[0], count=1)

	return fno,round(sum(chis),3)  

def timcodrvel(instance, filepath:str):
	from numpy import arange
	from re import findall, sub

	"""
	Wirting time, cordinates, speeds to a csv file
	"""

	match = findall(pattern=r"/\d+_", string=filepath)
	fno = sub(pattern=r"/",repl="",string=match[0], count=1)

	cords, vels = instance.cordNvel
	fps = instance.frames
	tstop:float = len(vels)/fps
	tdat = arange(start=0.0,stop=tstop,step=(1/fps))
	# list(zip(tdat, cords, vels))
	stringFormat = "{},{},{}\n"
	strdat = ''.join( list(map(lambda a,b,c: stringFormat.format(a, b, c) , tdat, cords, vels))  )
	path = "/dwitiya/magintpy/cordNvels/rri_" + fno + ".csv"
	flWrite = open(path,"w")
	flWrite.write(strdat)
	flWrite.close()
	return None

def aggplts(instance, ax):

	cords, vels = instance.cordNvel
	ax.plot(cords, vels, '+', markersize=0.6, alpha=0.5)

def aggpltsT(instance, ax):
	_, vels = instance.cordNvel
	tInc = 1/instance.frames
	tStop: float = len(vels)*tInc
	tdata = arange(start= 0.0, stop=tStop, step=tInc)
	
	ax.plot(tdata, vels, '+', markersize=0.6, alpha=0.5)

def modelFun123(t, y, mbot:float, km2:float, n1:float, n2:float, n3:float):
	"""
	'odeint' function for numerical solution needs the above kind of structure of the function to process
	t : independent cordinate array for the solution
	y : tuple or list of all the intial conditions, here r0, and rdot0
	"""
	r, z = y
	zdot = [z, ( km2/(r**4)-( n1*z + n2*z**2 + n3*z**1.5) )/mbot]
	return zdot

def modelFun(t, y, mbot:float, km2:float, n1:float, n2:float):
	"""
	'odeint' function for numerical solution needs the above kind of structure of the function to process
	t : independent cordinate array for the solution
	y : tuple or list of all the intial conditions, here r0, and rdot0
	"""
	r, z = y
	zdot = [z, ( km2/(r**4)-( n1*z + n2*z**2 ) )/mbot]
	return zdot

def modelSolver(tim, km2:float, n1:float, n2:float, flag:str='fit'): #  km2:float, n1:float, n2:float, flag:str='fit'
	"""
	'flag' has two options default 'fit' or 'normal', either to output only velocities or n x 2 numpy array of [cordinates, velocities]
	-- output of only velocities is the necesscessity for the 'curve_fit' function
	cm is distance unit, but the initial r0 for the solver need to be fed in 'm'
	OUTPUT :: also in cm, therefore 100 multiplied in num_soln
	'xcords' in the function is dependent on the outside evaluation of experimental data 
	"""
	from scipy.integrate import odeint
	r0, z0 = mcord[0]*10**-2, mvel[0]*10**-2

	y0 = [r0, z0]
	mbot = 5.6*10**-3   #7.16*10**-3 #5.6  7.16
	num_soln = 100*odeint(modelFun, t=tim, y0=y0,  args=(mbot, km2*10**-9, n1, n2), tfirst=True) # , tfirst=True
	rtn = {'normal': num_soln, 'fit' : num_soln[:,1]}
	return rtn[flag]


def dataFit(frames:int, Vels, fit:str='actual'):
	"""
	'fit' has two options default 'actual' or 'movavg', either to fit actual or moving averaged data
	use 'duration'= 5 s for fitting about 15 cm of the path covered (avg ~3 cm/s is the spped of the swimmer)
	"""	
	from scipy.optimize import curve_fit
	from numpy import arange, convolve, ones

	tInc: float = 1/frames
	tStop: float = len(Vels)*tInc
	"""
	-- check the validity of all point solver with tStop
	-- it is better to evaluate fit for first 15 seconds
	"""
	# pts = int(tStop*self.frames) + 1
	tdata = arange(start= 0.0, stop=tStop, step=tInc)

	def mav():
		vavg = convolve(Vels, ones(10)/10, mode='same')
		return vavg

	vdat = {'actual':Vels, 'movavg':mav()}


	popt, _ = curve_fit(modelSolver, xdata=tdata, ydata=vdat[fit], p0=(0.1,0.0001,0.0001), bounds=((0,0,0), (2,0.1,0.1)))


	return popt



if __name__ == "__main__":

	# listf = fileList(path="/dwitiya/bigtank/prjct_finalise/data_PPE_Repulsion/")
	# tbd_instances = [Tbd(ppcm=42,frames=24,boatMass=5.6*10**-3,fpath=fls) for fls in listf]
	# alfits = [inst.dataFit(fit='actual') for inst in tbd_instances]
	# for inst in tbd_instances:
	# 	optimals = inst.dataFit()
	# 	flno, tpeak = inst.plotEach(savefolder="/dwitiya/magintpy/aggplots/nonMatch/",optink=optimals)
	# 	#print(flno,"  n1=",optimals[0],"  n2=",optimals[1],"  km2=",optimals[2])
	# 	print(flno,"  tpeak=",round(tpeak,3),"  n1=",round(optimals[1],5),"  n2=",round(optimals[2],5),"  km2=",round(optimals[0],3))
	# # print(listf) # round(optimals[2],3)


	# listf = fileList(path="/dwitiya/bigtank/prjct_finalise/data_PPE_Repulsion/")
	# tbd_instances = [Tbd(ppcm=42,frames=24,boatMass=5.6*10**-3,fpath=fls) for fls in listf]
	# alchisqs = [chisq(inst) for inst in tbd_instances]
	# #print([abs(a-b)/a for a, b in alchisqs]) 
	# for ii in alchisqs: print(ii)


	listf = fileList(path="/folder_path/")
	tbd_instances = [Tbd(ppcm=53,frames=60,boatMass=7.16*10**-3,fpath=fls) for fls in listf]
	alfits = [timcodrvel(tbd_instances[ll],listf[ll]) for ll in range(len(listf))]

	listf = fileList(path="/folder_path/")
	tbd_instances = [Tbd(ppcm=42,frames=24,boatMass=5.6*10**-3,fpath=fls) for fls in listf]
	crdNvl = [inst.cordNvel for inst in tbd_instances]
	cvl = list(sum(crdNvl,()))
	lf = len(cvl)
	cord, vel = cvl[:lf-1:2], cvl[1:lf:2]
			# print(len(listf), len(cord))
			# for co in range(0,len(listf)): print(len(cord[co]), listf[co], co )
	minlen = min((len(co) for co in cord))

	acord = array([crd[:minlen] for crd in cord])
	avel = array([vl[:minlen] for vl in vel])
	mcord = array([mean(acord[:,sl]) for sl in range(0,minlen)])

	mvel = array([mean(avel[:,sl]) for sl in range(0,minlen)])
	errvel = array([std(avel[:,sl]) for sl in range(0,minlen)])

	frames = 24
	tInc = 1/frames
	tStop: float = len(mvel)*tInc
	tdata = arange(start= 0.0, stop=tStop, step=tInc)


			# stringFormat = "{},{},{}\n"
			# strdat = ''.join( list(map(lambda a,b,c: stringFormat.format(a, b, c) , tdata, mcord, mvel))  )
			# path = "/dwitiya/magintpy/cordNvels/avg_" + ".csv"
			# flWrite = open(path,"w")
			# flWrite.write(strdat)
			# flWrite.close()


	optimals = dataFit(frames=frames, Vels=mvel)
	mdlsoln = modelSolver(tim=tdata,km2=optimals[0], n1=optimals[1],n2=optimals[2], flag='normal' )
	# km2=optimals[0], n1=optimals[1],n2=optimals[2], flag='normal'

	fig, (dis, tim) = plt.subplots(nrows=2, ncols=1, sharey=True)

	# plt.errorbar(x=mcord,y=mvel, yerr=errvel, elinewidth=1.2, capsize=1.5, ecolor='maroon', color='lime', linewidth=1,label="Error:1 stdev")
	
	alfits = [aggplts(instance=inst,  ax=dis) for inst in tbd_instances]
	alfits = [aggpltsT(instance=inst, ax=tim) for inst in tbd_instances]

	dis.fill_between(mcord, mvel - errvel, mvel + errvel, color = 'maroon',label="Error:1 stdev")
	dis.plot(mcord, mvel, '-', color='lime', linewidth=0.7, label="Average Speed")

	tim.fill_between(tdata, mvel - errvel, mvel + errvel, color = 'maroon',label="Error:1 stdev")
	tim.plot(tdata, mvel, '-', color='lime', linewidth=0.7, label="Average Speed")
	lbl = ''.join([' n1=',str(round(optimals[1],5)),' n2=',str(round(optimals[2],4)),' km2=',str(round(optimals[0],3)),r'$\times 10^{-9}$'])
	# r'  km2 =1.4 $\times 10^{-9}$', ' n3 =',str(round(optimals[2],4))
	tim.plot(tdata,mdlsoln[:,1], 'b--',  label='Fit' )
	dis.plot(mdlsoln[:,0],mdlsoln[:,1], 'b--',  label='Fit' )
	dis.set_title(lbl)

	dis.grid(True)
	tim.grid(True)
	dis.set_xlabel('Distance (cm)')
	dis.set_ylabel('Speed (cm s$^{-1}$)')
	tim.set_xlabel('Time (s)')
	tim.set_ylabel('Speed (cm s$^{-1}$)')
	dis.legend(loc='upper right')
	tim.legend(loc='upper right')
	fig.tight_layout()	
	# plt.rcParams['pgf.texsystem']
	plt.savefig("/path/nam.pdf", dpi=400)






