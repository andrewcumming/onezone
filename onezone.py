# onezone.py
# One zone model for Type I X-ray bursts

from scipy.integrate import ode
from scipy import optimize
import numpy
import matplotlib.pyplot as plt
import sys
from math import exp

class OneZone:
	def __init__(self):
		# parameters
		self.arad=7.5657e-15
		self.clight=3e10
		self.gravity = 1.9e14
		self.mdot_Edd = 8.8e4
		self.mproton = 1.67e-24
		self.kB = 1.38e-16
		self.Q3a = 7.275
		# adjust fluxes into per gram units
		Qconv = 1.6e-6/self.mproton   # 1 MeV per nucleon in erg/g
		self.Q3a*=Qconv
		self.Qbase = 0.15*Qconv   # heating from below in MeV/nucleon
		# He mass fraction
		self.XHe=0.3
		self.Ye=1.0-self.XHe/2.0
		self.Yi=1.0-3.0*self.XHe/4.0

	def pressure(self,rho,T,Ye,Yi):
		Ped=9.91e12*(Ye*rho)**(5.0/3.0)   # non-relativistic degenerate electron pressure
		Pend=Ye*rho*1.38e-16*T/1.67e-24   # non-degenerate electron pressure 
		Pe=(Ped**2+Pend**2)**0.5   # fitting formula for electron pressure (Pacynski 1983)
		Pion=Yi*rho*1.38e-16*T/1.67e-24	# ion pressure
		Prad=self.arad*T**4/3.0	# radiation pressure
		return Pion+Pe+Prad

	def find_rho_eqn(self,rho,P,Ye,Yi,T):
		return self.pressure(rho,T,Ye,Yi)-P
		
	def find_rho(self,P,T,Ye,Yi):
		# inverts the equation of state to find the density
		rho = optimize.brentq(self.find_rho_eqn,1.0,1e8,xtol=1e-6,args=(P,Ye,Yi,T))
		return rho
	
	def kappa(self,rho5,T8,Ye):
		# radiative opacity
		kff = 0.4*0.753*Ye*rho5/T8**3.5    # free-free Gaunt factor set to 0.4
		kes = 0.4*Ye/((1+2.7*rho5/T8**2)*(1.0+(T8/4.5)**0.86))
		# Thompson scattering with corrections from Paczynski 1983
		return kes + kff   # note that strictly Rosseland mean opacities don't add (<30% error)
	
	def derivs(self,t,z):
		T = z[0]
		T8=1e-8*T
		ycolumn = z[1]
		rho5 = 1e-5 * self.find_rho(self.gravity*ycolumn,T,self.Ye,self.Yi)

		# nuclear burning rate for triple alpha
		eps3a=5.3e21*rho5**2*self.XHe**3*exp(-44.0/T8)/T8**3
		# boost the triple alpha energy release to burn to iron group
		eps3a*=(1.6+(1.0-self.XHe)*4.9)/0.606
		# and include hot CNO burning and the base flux
		eps = eps3a + 5.8e15*0.01 + self.Qbase*self.mdot/ycolumn

		# assume ideal gas for heat capacity of ions and electrons
		cp = 2.5*(self.Yi+self.Ye)*self.kB/self.mproton
		cp += 4.0*self.arad*1e19*T8**3/rho5;
    
		# cooling rate
		eps_cool = self.clight * self.arad * T**4 / (3.0 * self.kappa(rho5,T8,self.Ye) * ycolumn**2)

		# calculate derivatives
		dTdt = (eps - eps_cool) / cp
		dydt = self.mdot - 12.0 * eps3a * ycolumn/self.Q3a
	
		return [dTdt,dydt]
	
	def evolve(self,mdot=1.0,tend=1e4,rtol=1e-3,atol=1e-3):
		self.mdot=mdot*self.mdot_Edd	
		self.tend=tend
		# Integrate
		solver=ode(self.derivs)
		t=[]
		result=[]
		inic = [2e8,2e8]  # initial temperature and column depth
		solver.set_initial_value(inic,0.0)
		solver.set_integrator('vode',atol=atol,rtol=rtol)
		while solver.successful() and solver.t<tend:
			solver.integrate(self.tend,step=True)
			t.append(solver.t)
			result.append(solver.y)
		t=numpy.array(t)
		result=numpy.array(result)
		
		# Calculate the lightcurve
		lum = numpy.array([])
		llast = -10.0
		tlast=0.0
		goingup=True
		deltat=0.0
		for counter, T in enumerate(result[:,0]):
			ycolumn = result[counter,1]

			rho5 = 1e-5 * self.find_rho(self.gravity*ycolumn,T,self.Ye,self.Yi)
			T8=1e-8*T

			eps_cool = self.clight * self.arad * T**4 / (3.0 * self.kappa(rho5,T8,self.Ye) * ycolumn**2)
			luminosity = 4.0*3.1415*1e12 * ycolumn * eps_cool
			if (luminosity > llast):
				goingup=True
			else:
				if (goingup==True):
					deltat = t[counter]-tlast
					tlast=t[counter]
				goingup=False
			llast=luminosity
			lum=numpy.append(lum, luminosity)
		
		return t, result[:,1], result[:,0], lum, deltat
		
		
		
if __name__ == '__main__':

	# accretion rate and time to run from the command line
	if (len(sys.argv)<3):
		print 'Syntax: onezone.py <mdot> <time>'
		exit()
	mdot=float(sys.argv[1])
	tend=float(sys.argv[2])

	onezone = OneZone()
	time, y, T, L, deltat = onezone.evolve(mdot,tend)

	print 'Time between flashes=',deltat
	n=len(L)
	print 'Luminosity range = ',max(L[0.5*n:])-min(L[0.5*n:])

	# PLOT THE RESULTS

	# Plot the temperature against time
	#plt.plot(t,result[:,0])
	#plt.show()

	# Plot the lightcurve
	plt.plot(time,L)
	plt.xlabel('Time (s)')
	plt.ylabel('Luminosity (erg/s)')
	plt.show()
	#plt.savefig('lc.pdf')

	# Plot the temperature against column
	plt.plot(y,T)
	plt.xlabel('Column depth (g/cm^2)')
	plt.ylabel('Temperature (K)')
	plt.show()
	#plt.savefig('yT.pdf')
