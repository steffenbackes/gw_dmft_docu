import numpy as np
np.set_printoptions(precision=4, suppress=True)

t = 1.0
tp = 0.0
U = 2.0
corrOrb = [0]
uncorrOrb = [1]
nk = 100
mu = 0*U/2
nw = 200
wmin = -10.0
wmax = +10.0
delta = 0.1
dw = (wmax-wmin)/nw

Rmax = (nk-1)/2
orb = corrOrb + uncorrOrb
ncorr = len(corrOrb)
nuncorr = len(uncorrOrb)
norb = len(orb)

def H(k):
	return np.array([ [-2*t*np.cos(k), tp], [tp, -2*t*np.cos(k) - 3.0] ])

def sigma(k,w):
	#return np.zeros((norb,norb))
	return np.array( [ [mu  - 0.5*np.cos(k)  - 0.0j * w**2 * np.exp(-0.1*w**2)  * (1.0 - 0.0*np.cos(2*k) ),  0.0],[0.0, 0.0] ] )

def getSubMat(mat,i,j):
	sub = np.zeros(( len(i),len(j) ),dtype=complex)
	for m1 in range(len(i)):
		for m2 in range(len(j)):
			sub[m1,m2] = mat[i[m1],j[m2]]
	return sub

thop = np.zeros((2*Rmax+1, norb,norb),dtype=complex)
for i in range(2*Rmax+1):
	R = -Rmax+i
	for ik in range(nk):
		k = ik*2*np.pi/nk

		thop[i,:,:] += H(k) * np.exp(+1.0j*k*R) / nk 
	

	if ( np.linalg.norm(thop[i,:,:]) > 1.0E-5  ):
		print 'R=',R,', thop:\n',thop[i,:,:]


Gk = np.zeros((norb,norb,nk,nw),dtype=complex)
G00w = np.zeros((norb,norb,nw),dtype=complex)
sloc = np.zeros((norb,norb,nw),dtype=complex)
for ik in range(nk):
	k = ik*2*np.pi/nk
	for n in range(nw):
		w = wmin + n*dw + 1.0j*delta
		Gk[:,:,ik,n] = np.linalg.inv( (w+mu)*np.identity(norb) - H(k) - sigma(k,w)  )

		G00w[:,:,n]  += Gk[:,:,ik,n] / nk
		sloc[:,:,n] += sigma(k,w) / nk

Hybrid = np.zeros((ncorr,ncorr,nw),dtype=complex)
for i in range(2*Rmax+1):
	for j in range(2*Rmax+1):
		Ri = -Rmax+i
		Rj = -Rmax+j

		t0i = getSubMat(thop[i,:,:], corrOrb, orb )
		tj0 = getSubMat(thop[j,:,:], orb, corrOrb )

		if ( np.linalg.norm(t0i) < 1.0E-5 or np.linalg.norm(tj0) < 1.0E-5  ):
			continue

		Gijw = np.zeros((norb,norb,nw),dtype=complex)
		Gi0w = np.zeros((norb,norb,nw),dtype=complex)
		G0jw = np.zeros((norb,norb,nw),dtype=complex)

		for ik in range(nk):
			k = ik*2*np.pi/nk

			Gijw += Gk[:,:,ik,:] * np.exp(+1.0j*k*(Ri-Rj) ) / nk
			Gi0w += Gk[:,:,ik,:] * np.exp(+1.0j*k*(Ri   ) ) / nk
			G0jw += Gk[:,:,ik,:] * np.exp(+1.0j*k*(  -Rj) ) / nk


		if Ri==0:
			t0i = getSubMat(thop[i,:,:], corrOrb, uncorrOrb )
		if Rj==0:
			tj0 = getSubMat(thop[j,:,:], uncorrOrb, corrOrb )

		# Now we have all Green's fucntions
		for n in range(nw):
			w = wmin + n*dw + 1.0j*delta

			Gi0 = getSubMat(Gi0w[:,:,n], orb, corrOrb)
			G0j = getSubMat(G0jw[:,:,n], corrOrb, orb)
			Gij = getSubMat(Gijw[:,:,n], orb, orb)

			G00inv = np.linalg.inv( getSubMat( G00w[:,:,n], corrOrb,corrOrb ) )

			if Ri==0:
				Gi0 = getSubMat(Gi0w[:,:,n], uncorrOrb, corrOrb)
				Gij = getSubMat(Gijw[:,:,n], uncorrOrb, orb)
			if Rj==0:
				G0j = getSubMat(G0jw[:,:,n], corrOrb, uncorrOrb)
				Gij = getSubMat(Gijw[:,:,n], orb, uncorrOrb)
			if Ri==0 and Rj==0:
				Gij = getSubMat(Gijw[:,:,n], uncorrOrb, uncorrOrb)
	
			G0ij = Gij - np.dot( np.dot( Gi0 , G00inv ) , G0j )

			Hybrid[:,:,n] +=  np.dot( np.dot( t0i , G0ij ),  tj0 ) 


outf = open('bathhybKmultiorb.dat','w')
Hybrid_dyson = np.zeros((ncorr,ncorr,nw),dtype=complex)
for n in range(nw):
	w = wmin + n*dw + 1.0j*delta

	hyb = Hybrid[:,:,n]
	bath = np.linalg.inv( (w+mu)*np.identity(ncorr) - hyb )

	################################

	bath_dyson = np.linalg.inv( np.linalg.inv( getSubMat( G00w[:,:,n], corrOrb,corrOrb ) ) +  getSubMat( sloc[:,:,n], corrOrb,corrOrb ) )
	hyb_dyson = (w+mu)*np.identity(ncorr) - np.linalg.inv(bath_dyson)

	Hybrid_dyson[:,:,n] = hyb_dyson

	#################################
	
	Gsnonloc = 0.0j
	for ik in range(nk):
		k = ik*2*np.pi/nk

		Gsnonloc += getSubMat( np.linalg.inv( (w+mu)*np.identity(norb) - H(k) - sigma(k,w) + sloc[:,:,n]  ) , corrOrb,corrOrb ) / nk

	#################################

	outf.write(str(w.real) + '\t')	
	m = 0

#	outf.write(str(-G00.real) + '\t' )	
#	outf.write(str(-Gab.real) + '\t' )	
#	outf.write(str(-Gba.real) + '\t' )	
#	outf.write(str(-Gsnonloc.imag) + '\t' )	
#	outf.write(str( (1.0/bath).real) + '\t')	
#	outf.write(str(hyb.real) + '\t')	

	outf.write(str(-bath[m,m].imag) + '\t')	
	outf.write(str(-hyb[m,m].imag) + '\t')	
	outf.write(str(-bath_dyson[m,m].imag) + '\t')	
	outf.write(str(-hyb_dyson[m,m].imag) + '\t')	
	outf.write(str(-Gsnonloc[m,m].imag) + '\t' )	
	outf.write('\n')	
outf.close()

outf = open('delta0.dat','w')
beta = 40.0
ntau = 1000
for it in range(ntau+1):
	tau = it*beta/ntau
	
	hybt = 0.0
	for n in range(nw):
		w = wmin + n*dw 
	
		hybt += ( Hybrid[0,0,n].imag/( np.exp(w*(tau-beta)) + np.exp(w*tau) ) )*dw/np.pi

	outf.write( str(tau) + '\t' + str(hybt) + '\t' + str(hybt) + '\n' )
outf.close()
