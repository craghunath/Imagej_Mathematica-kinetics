from dataclasses import dataclass, field
from typing import List


# from re import sub, split, findall
# from os import listdir
# import matplotlib.pyplot as plt
# from numpy import arange, array, diff, sqrt, append, zeros,  mean as nmean  #, abs, where ,convolve, ones, split as nsplit,
# from scipy.integrate import odeint
# from scipy.optimize import curve_fit 

@dataclass
class UTbd:
	"""
	These are essentials required for computing Two boat dynamics
		ppcm : pixcels per centimeter, so the entire experimental data of distance in cm and time in s
		fpath : Strictly folder path. used for averaging the all the data
		frames : number of frames per second at which the dynamics video is recorded
	"""
	cordNvel: tuple = field(init=False)
	ppcm: float = 10
	frames: float = 10
	boatMass: float = 5.6*10**-3
	fpath: str = '000'
	# folderpath: str = '000'
	def __post_init__(self):

		self.decel = int(3*self.frames)
		self.endt = int(15*self.frames)
		self.cordNvel = self.dataexpAvg(self.fpath)

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
		cords = btr/ppcm 
		
		#tSlice = 1/self.frames
		vel =  diff(cords)*self.frames   #diff(cords)/tSlice
		vels = append(zeros(1),vel)
		"""initial velocity is non zero check this """
		return [cords[self.decel:self.endt],vels[self.decel:self.endt]]



	def dataexpAvg(self,folderpath):
		#from os.path import isdir, isfile
		from os import listdir
		from re import  findall
		from numpy import mean, array
		"""			
		OUTPUT :: If fpath is a folder/directory, output is list of tuples with ndarray pair. 
			:: If fpath is a file, output is a tuple of ndarray pair.
		"""
	
		foldr = folderpath
		#(?=.*15\d{4}_)(?=.*23Aug)ll)/" #replace the folder path
		flnames = list(filter(lambda strn: findall('(?=.*)(?=.*\.csv)', strn), listdir(foldr)));
		filepaths = list(map(lambda fl: ''.join([foldr, fl]), flnames));
		#print(filepaths)(?=.*PPE_Repul)

		cordsNvels =  array(list(map(lambda flps: self.dataexp(filepath=flps), filepaths)))

		return mean(cordsNvels[:,0], axis=0), mean(cordsNvels[:,1], axis=0)

	def u1(self,t,c:float):
		from numpy import exp
		u0 = self.cordNvel[1][0]
		return u0 * exp(-c*t) 	

	def u2(self,t,c:float):
		u0 = self.cordNvel[1][0]
		return 1/(1/u0 + c*t) 

	def u3by2(self,t,c:float):
		from numpy import sqrt
		u0 = self.cordNvel[1][0]
		return 1/(1/sqrt(u0) + c*t)**2

	def u12(self,t,c1:float, c2:float):
		from numpy import exp
		u0 = self.cordNvel[1][0]
		return u0 * exp(-c1*t) + 1/(1/u0 + c2*t)

	def modelFun(self, t, y, mbot:float, km2:float, n:float ):
		"""
		'odeint' function for numerical solution needs the above kind of structure of the function to process
		t : independent cordinate array for the solution
		y : tuple or list of all the intial conditions, here r0, and rdot0
		"""
		r, z = y
		zdot = [z, ( km2/(r**4)-( n*z**(2) ) )/mbot]
		return zdot

	def modelFunn3(self, t, y, mbot:float, km2:float, n:float ):
		"""
		'odeint' function for numerical solution needs the above kind of structure of the function to process
		t : independent cordinate array for the solution
		y : tuple or list of all the intial conditions, here r0, and rdot0
		"""
		r, z = y
		zdot = [z, ( km2/(r**4)-( n*z**(1.5) ) )/mbot]
		return zdot

	def modelFunn1(self, t, y, mbot:float, km2:float, n:float ):
		"""
		'odeint' function for numerical solution needs the above kind of structure of the function to process
		t : independent cordinate array for the solution
		y : tuple or list of all the intial conditions, here r0, and rdot0
		"""
		r, z = y
		zdot = [z, ( km2/(r**4)-( n*z ) )/mbot]
		return zdot

	def modelFunn12(self, t, y, mbot:float, km2:float, n1:float, n2:float ):
		"""
		'odeint' function for numerical solution needs the above kind of structure of the function to process
		t : independent cordinate array for the solution
		y : tuple or list of all the intial conditions, here r0, and rdot0
		"""
		r, z = y
		zdot = [z, ( km2/(r**4)-( n1*z + n2*z**2 ) )/mbot]
		return zdot

	def modelSolver(self, tim, n:float, flag:str='fit'):
		from scipy.integrate import odeint
		"""
		'flag' has two options default 'fit' or 'normal', either to output [cordinates, velocities] or only velocities
			-- output of only velocities is the necesscessity for the 'curve_fit' function
		cm is distance unit, but the initial r0 for the solver need to be fed in 'm'
		OUTPUT :: also in cm, therefore 100 multiplied in num_soln
		'xcords' in the function is dependent on the outside evaluation of experimental data 
		"""

		xcords, vels = self.cordNvel

		r0, z0 = xcords[0]*10**-2, vels[0]*10**-2
		y0 = [r0, z0]
		mbot = self.boatMass   #7.16*10**-3 #5.6  7.16
		num_soln = 100*odeint(self.modelFun, t=tim, y0=y0,  args=(mbot, 0, n), tfirst=True) # , tfirst=True
		rtn = {'normal': num_soln, 'fit' : num_soln[:,1]}
		return rtn[flag]

	def modelSolvern3(self, tim, n:float, flag:str='fit'):
		from scipy.integrate import odeint
		"""
		'flag' has two options default 'fit' or 'normal', either to output [cordinates, velocities] or only velocities
			-- output of only velocities is the necesscessity for the 'curve_fit' function
		cm is distance unit, but the initial r0 for the solver need to be fed in 'm'
		OUTPUT :: also in cm, therefore 100 multiplied in num_soln
		'xcords' in the function is dependent on the outside evaluation of experimental data 
		"""

		xcords, vels = self.cordNvel
		

		r0, z0 = xcords[0]*10**-2, vels[0]*10**-2
		y0 = [r0, z0]
		mbot = self.boatMass   #7.16*10**-3 #5.6  7.16
		num_soln = 100*odeint(self.modelFunn3, t=tim, y0=y0,  args=(mbot, 0, n), tfirst=True) # , tfirst=True
		rtn = {'normal': num_soln, 'fit' : num_soln[:,1]}
		return rtn[flag]

	def modelSolvern1(self, tim, n:float, flag:str='fit'):
		from scipy.integrate import odeint
		"""
		'flag' has two options default 'fit' or 'normal', either to output [cordinates, velocities] or only velocities
			-- output of only velocities is the necesscessity for the 'curve_fit' function
		cm is distance unit, but the initial r0 for the solver need to be fed in 'm'
		OUTPUT :: also in cm, therefore 100 multiplied in num_soln
		'xcords' in the function is dependent on the outside evaluation of experimental data 
		"""

		xcords, vels = self.cordNvel
		

		r0, z0 = xcords[0]*10**-2, vels[0]*10**-2
		y0 = [r0, z0]
		mbot = self.boatMass   #7.16*10**-3 #5.6  7.16
		num_soln = 100*odeint(self.modelFunn1, t=tim, y0=y0,  args=(mbot, 0, n), tfirst=True) # , tfirst=True
		rtn = {'normal': num_soln, 'fit' : num_soln[:,1]}
		return rtn[flag]
	
	def modelSolvern12(self, tim, n1:float, n2:float, flag:str='fit'):
		from scipy.integrate import odeint
		"""
		'flag' has two options default 'fit' or 'normal', either to output [cordinates, velocities] or only velocities
			-- output of only velocities is the necesscessity for the 'curve_fit' function
		cm is distance unit, but the initial r0 for the solver need to be fed in 'm'
		OUTPUT :: also in cm, therefore 100 multiplied in num_soln
		'xcords' in the function is dependent on the outside evaluation of experimental data 
		"""

		xcords, vels = self.cordNvel
		

		r0, z0 = xcords[0]*10**-2, vels[0]*10**-2
		y0 = [r0, z0]
		mbot = self.boatMass   #7.16*10**-3 #5.6  7.16
		num_soln = 100*odeint(self.modelFunn12, t=tim, y0=y0,  args=(mbot, 0, n1, n2), tfirst=True) # , tfirst=True
		rtn = {'normal': num_soln, 'fit' : num_soln[:,1]}
		return rtn[flag]
	


	def dataFit(self, fit:str='actual'):
		"""
		'fit' has two options default 'actual' or 'movavg', either to fit actual or moving averaged data
		use 'duration'= 5 s for fitting about 15 cm of the path covered (avg ~3 cm/s is the speed of the swimmer)
		"""		
		from scipy.optimize import curve_fit
		from numpy import arange, convolve, ones

		_, vels = self.cordNvel
		
		tInc = 1/self.frames
		tStop = float((self.endt - self.decel)*tInc)

		"""
		-- check the validity of all point solver with tStop
		"""
		#pts = int(tStop*self.frames)
		tdata = arange(start=0.0, stop=tStop, step=tInc)

		def mav():
			vavg = convolve(vels, ones(10)/10, mode='same')
			return vavg

		vdat = {'actual':vels, 'movavg':mav()}

		# optRes = {
		# 'n2':curve_fit(self.modelSolver, xdata=tdata, ydata=vdat[fit], p0=(0.0001), bounds=((0), (1)), method='trf'),
		# 'n3':curve_fit(self.modelSolvern3, xdata=tdata, ydata=vdat[fit], p0=(0.0001), bounds=((0), (1)), method='trf'),
		# 'u2':curve_fit(self.u2,xdata=tdata, ydata=vdat[fit], p0=(0.1), bounds=((0), (5)), method='trf'),
		# 'u3by2':curve_fit(self.u3by2,xdata=tdata, ydata=vdat[fit], p0=(0.1), bounds=((0), (5)), method='trf')
		# }
		optRes = [
		curve_fit(self.modelSolver, xdata=tdata, ydata=vdat[fit], p0=(0.0001), bounds=((0), (1)), method='trf')[0][0],
		curve_fit(self.modelSolvern3, xdata=tdata, ydata=vdat[fit], p0=(0.0001), bounds=((0), (1)), method='trf')[0][0],
		curve_fit(self.u2,xdata=tdata, ydata=vdat[fit], p0=(0.1), bounds=((0), (5)), method='trf')[0][0],
		curve_fit(self.u3by2,xdata=tdata, ydata=vdat[fit], p0=(0.1), bounds=((0), (10)), method='trf')[0][0],
		curve_fit(self.u1,xdata=tdata, ydata=vdat[fit], p0=(0.1), bounds=((0), (10)), method='trf')[0][0],
		curve_fit(self.modelSolvern1,xdata=tdata, ydata=vdat[fit], p0=(0.0001), bounds=((0), (1)), method='trf')[0][0],
		curve_fit(self.modelSolvern12,xdata=tdata, ydata=vdat[fit], p0=(0.0001,0.0001), bounds=((0,0), (1,1)), method='trf')[0],
		# curve_fit(self.u12,xdata=tdata, ydata=vdat[fit], p0=(0.1,0.1), bounds=((0,0), (2,2)), method='trf')[0]
		]
		
		# popt, pcov
		return optRes

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
		import matplotlib.pyplot as plt
		# plt.rcParams.update({'figure.max_open_warning': 0}) # this is another method for max warning silence
		plt.rc('figure', max_open_warning = 0)
		plt.rcParams.update({'axes.facecolor':'whitesmoke'})

		from numpy import arange
		from re import sub, findall

		flp:str = filepath
		if filepath == '000':
			cords, vels = self.cordNvel
		else:
			cords, vels = self.dataexp(filepath=flp)

		bm = round(10**3 * self.boatMass, 3)

		# match = findall(pattern=r"/\d+_", string=self.fpath)
		# fno = sub(pattern=r"/",repl="",string=match[0], count=1)
		
		# "_5sor15cm",


		tInc = 1/self.frames
		tStop = float((self.endt - self.decel)*tInc)
		"""
		-- check the validity of all point solver with tStop
		-- it is better to evaluate fit for first 5 seconds
		"""
		#pts = int(tStop*self.frames)
		tdata = arange(start= 0.0, stop=tStop, step=tInc)
		soln2 = self.modelSolver(tim=tdata, n=optink[0], flag='normal')
		soln3 = self.modelSolvern3(tim=tdata, n=optink[1], flag='normal')
		solnu2 = self.u2(t=tdata, c=optink[2])
		solnu3by2 = self.u3by2(t=tdata, c=optink[3])
		solnu1 = self.u1(t=tdata, c=optink[4])
		soln1 = self.modelSolvern1(tim=tdata, n=optink[5], flag='normal')
		soln12 = self.modelSolvern12(tim=tdata, n1=optink[6][0], n2=optink[6][1], flag='normal')
		# solnu12 = self.u12(t=tdata, c1=optink[7][0], c2=optink[7][1])
		plt.figure()
		
		plt.plot(tdata,soln2[:,1], "b--", linewidth=2, alpha=0.6, label='Fit: n2 = %.2e'%optink[0])
		# plt.plot(tdata,soln3[:,1], "-.",color='navy', linewidth=0.8,  label='Fit: n3 = %.2e'%optink[1])
		plt.plot(tdata,solnu2, "k:", linewidth=1, label=r'Fit:$1/(1/U_{0}+C t)$, C = %.2e'%optink[2])
		plt.plot(tdata,solnu3by2, "-",color='lime', linewidth=2.5, label=r'Fit: $1/(1/\sqrt{U_{0}}+C t)^{2}$, C = %.2e'%optink[3])
		plt.plot(tdata,soln3[:,1], "-.",color='navy', linewidth=0.8,  label='Fit: n3 = %.2e'%optink[1])
		plt.plot(tdata,soln1[:,1], "-",color='cyan', linewidth=2, label='Fit: n1 = %.2e'%optink[5])
		plt.plot(tdata,solnu1, color='darkgreen',linestyle=(0, (3, 1, 1, 1)), linewidth=1, label=r'Fit: $U_{0} e^{-C t}$, C = %.2e'%optink[4])
		plt.plot(tdata, vels, "r-", linewidth=2, alpha=0.75, label='data')
		plt.plot(tdata,soln12[:,1], "-",color='black', linewidth=1, label='Fit: n1 = %.2e, n2 = %.2e'%tuple(optink[6]))

		# plt.plot(tdata,solnu12, color='indigo',linestyle=(0, (5, 1)), linewidth=0.7, label=r'Fit:C1 = %.2e, C2 = %.2e'%tuple(optink[7]))
		# $U_{0} e^{-C_{1} t} + 1/(1/U_{0}+C_{2} t)$
		# plt.grid(color='k', linestyle='-', linewidth=1, alpha=0.4)
		plt.title( 'Deceleration Fits')
		# fln = "".join((savefolder,"Uallk0_",fno,str(self.ppcm),"_ppe.jpg"))
		fln = "".join((savefolder,"Uall3sk0_",str(self.ppcm),"_ppe.jpg"))
		plt.xlabel('time (s)')
		plt.ylabel('Speed (cm/s)')
		plt.legend()
		# plt.ticklabel_format(axis='y', style='scientific',scilimits=(0,0))
		plt.savefig(fln, dpi=400)

		return None #fno


def fileList(path:str)-> List[str]:
	#from os.path import isdir, isfile
	from os import listdir
	from re import  findall

	#(?=.*15\d{4}_)(?=.*23Aug)ll)/" #replace the folder path
	flnames = list(filter(lambda strn: findall('(?=.*)(?=.*\.csv)', strn), listdir(path)));
	return list(map(lambda fl: ''.join([path, fl]), flnames))


from numpy import shape
if __name__ == "__main__":

	fls = "/folder_path/"
	tbdavg = UTbd(ppcm=42,frames=24,boatMass=5.6*10**-3,fpath=fls)
	tbdavg.plotEach(savefolder="/to_folder/",optink=tbdavg.dataFit())

