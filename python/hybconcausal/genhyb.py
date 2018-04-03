import numpy as np

readsigma = 0
nk = 200
nw = 200
emin = -10.0
emax = 10.0
de = (emax-emin)/nw
ntau = 500
beta = 40.0

# k runs from 0 to 1

def eps(k):
	return -np.cos(k*np.pi*2)

def gauss(x,mu,s):
	return np.exp(-0.5*(x-mu)**2/s**2)/np.sqrt(s**2 * 2*np.pi)

def ImSig(k,w):
	fac = 25
	disp = 0.45
	return -gauss(w,-5 -disp*np.cos(k*np.pi*2),1)*fac - gauss(w,+5 -disp*np.cos(k*np.pi*2),1)*fac
	#return -w**2 * np.exp(-0.04*( w + disp*np.cos(k*np.pi*2) )**2)

Sigma = np.zeros((nw,nk),dtype=complex)
Glat  = np.zeros((nw,nk),dtype=complex)
Gloc  = np.zeros((nw),dtype=complex)

# generate complete Selfenergy from Kramers Kronig
print 'generate sigma...'
if (readsigma==1):
	inf = open('sigma.dat','r')
	for iw in range(nw):
		data = [float(x) for x in inf.readline().split()]
		for ik in range(nk):
			Sigma[iw,ik] = data[1+2*ik+0] + data[1+2*ik+1]*1.0j
	inf.close()
outf = open('sigma.dat','w')
outfg0 = open('gf0.dat','w')
outfg = open('gf.dat','w')
outfg0l = open('g0loc.dat','w')
outfgl = open('gloc.dat','w')
outfgb = open('gbath.dat','w')
outfsl = open('sigma_loc.dat','w')
outfhyb = open('hybrid.dat','w')
hybw = []
hybpos = 0.0
hybneg = 0.0
for iw in range(nw):
	print iw*100.0/nw,'% done'
	w = emin + iw*de
	outf.write( str(w) + '\t' )
	outfg.write( str(w) + '\t' )
	outfg0.write( str(w) + '\t' )
	outfgl.write( str(w) + '\t' )
	outfg0l.write( str(w) + '\t' )
	outfsl.write( str(w) + '\t' )
	outfgb.write( str(w) + '\t' )
	outfhyb.write( str(w) + '\t' )
	gloc = 0.0j
	gloc0 = 0.0j
	sloc = 0.0j
	for ik in range(nk):
		k = ik*1.0/nk
#######################################
		if (readsigma!=1):
			Sigma[iw,ik] = 1.0j*ImSig(k,w)
		
			for iw2 in range(nw):
				if (iw != iw2 ):
					w2 = emin + iw2*de
					Sigma[iw,ik] += ImSig(k,w2)*de/(np.pi*( w2-w ))
#######################################
		Sigma[iw,ik] += -np.cos(k*2*np.pi)*0.2
		outf.write( str(Sigma[iw,ik].real) + '\t' +  str(Sigma[iw,ik].imag) + '\t' )

		# green's function
		gk0 = 1.0/( w + 0.01j - eps(k) )
		gk = 1.0/( w + 0.01j - eps(k) - Sigma[iw,ik] )
		outfg.write(  str(gk.real) + '\t' +  str(gk.imag) + '\t' )
		outfg0.write( str(gk0.real) + '\t' +  str(gk0.imag) + '\t' )

		gloc += gk/nk
		gloc0 += gk0/nk
		sloc += Sigma[iw,ik]/nk

	outf.write('\n')
	outfg.write('\n')
	outfg0.write('\n')

	outfgl.write( str(gloc.real) + '\t' +  str(gloc.imag) + '\n')
	outfg0l.write(str(gloc0.real) + '\t' +  str(gloc0.imag) +'\n')
	outfsl.write(str(sloc.real) + '\t' +  str(sloc.imag) +'\n')

	gbath = 1.0/( 1.0/gloc + sloc )
	outfgb.write(str(gbath.real) + '\t' +  str(gbath.imag) +'\n')

	hyb = w - 1.0/gbath
	hybw.append(hyb)
	outfhyb.write(str(hyb.real) + '\t' +  str(hyb.imag) +'\n')

	if (hyb.imag>0):
		hybpos += hyb.imag*de
	else:
		hybneg += hyb.imag*de
outf.close()
outfg0.close()
outfg.close()
outfg0l.close()
outfgl.close()
outfgb.close()
outfsl.close()
outfhyb.close()

print 'Hybridization positiv:',hybpos
print 'Hybridization negativ:',hybneg
print np.round(-hybpos*100.0/hybneg,1),',% positive!'

outf = open('hyb_tau.dat','w')
for it in range(ntau+1):
	tau = it*beta/ntau
	hybt = 0.0
	hybtc = 0.0
	outf.write( str(tau) + '\t' )
	for iw in range(nw):
		w = emin + iw*de
		hybt += hybw[iw].imag*np.exp(-tau*w)/( np.exp(-beta*w)+1 ) * de/beta

		if (hybw[iw].imag<0):
			hybtc += hybw[iw].imag*np.exp(-tau*w)/( np.exp(-beta*w)+1 ) * de/beta
	outf.write( str(hybt) + '\t' + str(hybtc) + '\n' )
outf.close()

outf = open('hyb_mats.dat','w')
for n in range(1000):
	iwn = (2*n+1)*np.pi/beta
	hybiw = 0.0
	hybiwc = 0.0
	outf.write( str(iwn) + '\t' )
	for iw in range(nw):
		w = emin + iw*de
		hybiw += hybw[iw].imag/( w - 1.0j*iwn ) * de

		if (hybw[iw].imag<0):
			hybiwc += hybw[iw].imag/( w - 1.0j*iwn ) * de

	outf.write( str(hybiw.real) + '\t' + str(hybiw.imag) + '\t' + str(hybiwc.real) + '\t' + str(hybiwc.imag) +  '\n' )
outf.close()
