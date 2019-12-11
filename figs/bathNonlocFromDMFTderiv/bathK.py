import numpy as np

t = 1.0
U = 2.0
mu = 0*U/2
nw = 200
nk = 100
wmin = -10.0
wmax = +10.0
delta = 0.1

def eps(k):
	return -2*t*np.cos(k)

def sigma(k):
	return mu  - 0.0*np.cos(k)  - 0.0j * w**2 * np.exp(-0.1*w**2)  * (1.0 - 0.0*np.cos(2*k) )

dw = (wmax-wmin)/nw
outf = open('bathhybK.dat','w')
for n in range(nw):
	w = wmin + n*dw + 1.0j*delta

	G00 = 0.0j
	Ga0 = 0.0j
	G0a = 0.0j
	Gb0 = 0.0j
	G0b = 0.0j
	Gab = 0.0j
	Gba = 0.0j
	sloc = 0.0j
	for ik in range(nk):
		k = ik*2*np.pi/nk

		sloc += sigma(k)/nk

		G = 1.0/( w+mu - eps(k) - sigma(k) )

		G00 += G / nk 
		Ga0 += G * np.exp(+1.0j*k) / nk 
		G0a += G * np.exp(-1.0j*k) / nk 
		Gab += G * np.exp(+2.0j*k) / nk 
		Gba += G * np.exp(-2.0j*k) / nk 

	Gb0 = G0a
	G0b = Ga0

	hyb = 0.0j
	# i = a j = a
	hyb += t**2 * ( G00 - Ga0*G0a / G00  )

	# i = a j = b
	hyb += t**2 * ( Gab - Ga0*G0b / G00  )

	# i = b j = a
	hyb += t**2 * ( Gba - Gb0*G0a / G00  )

	# i = b j = b
	hyb += t**2 * ( G00 - Gb0*G0b / G00  )

	bath = 1.0/( w+mu*0 - hyb )

	################################

	bath_dyson = 1.0/( 1/G00 + sloc )
	hyb_dyson = w+mu - 1.0/bath_dyson

	#################################
	
	Gsnonloc = 0.0j
	for ik in range(nk):
		k = ik*2*np.pi/nk

		Gsnonloc += 1.0/( w+mu - eps(k) - sigma(k) + sloc ) / nk

	#################################

	outf.write(str(w.real) + '\t')	

#	outf.write(str(-G00.real) + '\t' )	
	outf.write(str(-Gab.real) + '\t' )	
	outf.write(str(-Gba.real) + '\t' )	
#	outf.write(str(-Gsnonloc.imag) + '\t' )	

	outf.write(str( (1.0/bath).real) + '\t')	
	outf.write(str(hyb.real) + '\t')	

	outf.write(str(-bath.imag) + '\t')	
	outf.write(str(-hyb.imag) + '\t')	
	outf.write(str(-bath_dyson.imag) + '\t')	
	outf.write(str(-hyb_dyson.imag) + '\t')	

	outf.write('\n')	

outf.close()
