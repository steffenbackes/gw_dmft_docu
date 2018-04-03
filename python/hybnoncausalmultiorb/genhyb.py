import numpy as np

readsigma = 0
nk = 50
nw = 200
emin = -10.0
emax = 10.0
de = (emax-emin)/nw
ntau = 500
beta = 40.0
norb = 2

# k runs from 0 to 1

def eps(k):
	return -np.cos(k*np.pi*2)

def gauss(x,mu,s):
	return np.exp(-0.5*(x-mu)**2/s**2)/np.sqrt(s**2 * 2*np.pi)

def ImSig(k,w):
	fac = 25
	disp = 0.0
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
hybpos = np.zeros((norb))
hybneg = np.zeros((norb))
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
	gloc = np.zeros((norb,norb),dtype=complex)
	gloc0 = np.zeros((norb,norb),dtype=complex)
	sloc = np.zeros((norb,norb),dtype=complex)
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
#		Sigma[iw,ik] += -np.cos(k*2*np.pi)*0.2
		outf.write( str(Sigma[iw,ik].real) + '\t' +  str(Sigma[iw,ik].imag) + '\t' )

		# green's function
		#gk0 = 1.0/( w + 0.01j - eps(k) )
		#gk = 1.0/( w + 0.01j - eps(k) - Sigma[iw,ik] )

		epsmat = np.array([ [eps(k), 0], [0, -0.8*eps(k) ] ])
		hybmat = np.array([ [0, 0.1*np.cos(k*2*np.pi) ],[0.1*np.cos(k*2*np.pi),0] ] )
		sigmamat = np.array([ [Sigma[iw,ik], 0], [0, 0.4*Sigma[iw,ik]] ] )

		gk0 = np.linalg.inv( (w + 0.01j)*np.identity(norb) -epsmat - hybmat  )
		gk = np.linalg.inv(  (w + 0.01j)*np.identity(norb) -epsmat - hybmat - sigmamat  )

		for m in range(norb):
			outfg.write(  str(gk[m,m].real) + '\t' +  str(gk[m,m].imag) + '\t' )
			outfg0.write( str(gk0[m,m].real) + '\t' +  str(gk0[m,m].imag) + '\t' )

		gloc += gk/nk
		gloc0 += gk0/nk
		sloc += sigmamat/nk

	outf.write('\n')
	outfg.write('\n')
	outfg0.write('\n')

	gbath = np.linalg.inv( np.linalg.inv(gloc) + sloc )
	hyb = w*np.identity(norb) - np.linalg.inv(gbath)
	hybw.append(hyb)

	for m in range(norb):
		outfgl.write( str(gloc[m,m].real) + '\t' +  str(gloc[m,m].imag) + '\t')
		outfg0l.write(str(gloc0[m,m].real) + '\t' +  str(gloc0[m,m].imag) +'\t')
		outfsl.write(str(sloc[m,m].real) + '\t' +  str(sloc[m,m].imag) +'\t')

		outfgb.write(str(gbath[m,m].real) + '\t' +  str(gbath[m,m].imag) +'\t')

		outfhyb.write(str(hyb[m,m].real) + '\t' +  str(hyb[m,m].imag) +'\t')

		if (hyb[m,m].imag>0):
			hybpos[m] += hyb[m,m].imag*de
		else:
			hybneg[m] += hyb[m,m].imag*de

	outfgl.write('\n' )
	outfg0l.write('\n' )
	outfsl.write('\n' )
	outfgb.write('\n' )
	outfhyb.write( '\n' )


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
	outf.write( str(tau) + '\t' )
	for m in range(norb):
		hybt = 0.0
		hybtc = 0.0
		for iw in range(nw):
			w = emin + iw*de
			hybt += hybw[iw][m,m].imag*np.exp(-tau*w)/( np.exp(-beta*w)+1 ) * de/beta

			if (hybw[iw][m,m].imag<0):
				hybtc += hybw[iw][m,m].imag*np.exp(-tau*w)/( np.exp(-beta*w)+1 ) * de/beta
		outf.write( str(hybt) + '\t' + str(hybtc) + '\t' )
	outf.write( '\n' )
outf.close()

outf = open('hyb_mats.dat','w')
for n in range(1000):
	iwn = (2*n+1)*np.pi/beta
	outf.write( str(iwn) + '\t' )
	for m in range(norb):
		hybiw = 0.0
		hybiwc = 0.0
		for iw in range(nw):
			w = emin + iw*de
			hybiw += hybw[iw][m,m].imag/( w - 1.0j*iwn ) * de

			if (hybw[iw][m,m].imag<0):
				hybiwc += hybw[iw][m,m].imag/( w - 1.0j*iwn ) * de

		outf.write( str(hybiw.real) + '\t' + str(hybiw.imag) + '\t' + str(hybiwc.real) + '\t' + str(hybiwc.imag) +  '\t' )
	outf.write( '\n' )
outf.close()

