import numpy as np

print 'Read DFT Greens function...'
print 'Only read a single orbital...'

norb = 3
nspin = 2
Nw = 3000
beta = 40.0

############################################################################################

############################################################################################

Gf = np.zeros((Nw),dtype=complex)
P = np.zeros((norb,norb,Nw))

inf = open('gloc.dat','r')
for n in range(Nw):
	data = [ float(x) for x in inf.readline().split() ]
	Gf[n] = data[1] + data[2]*1.0j
inf.close()

def getgf(n):
	if (n<-Nw):
		return +1.0j/( (2*n+1)*np.pi/beta )
	if (n<0):
		return np.conj(Gf[-n-1])
	if (n<Nw):
		return Gf[n]
	else:
		return -1.0j/( (2*n+1)*np.pi/beta )

print 'Calculate polarization...'

outf = open('P.dat','w')
Nwmax = 500
for n in range(Nwmax):
#	print n*100.0/Nwmax,'% done...'
	Ptmp = 0.0
	for m in range(-Nw+1,Nw):
		Ptmp += ( getgf(m)*getgf(m-n) ).real
	Ptmp = Ptmp*nspin/beta

	outf.write( str(2*n*np.pi/beta) + '\t'  )
	for m1 in range(norb):
		for m2 in range(norb):
			P[m1,m2,n] = Ptmp - 0.0017
			
			outf.write( str(P[m1,m2,n]) + '\t' )

	outf.write( '\n' )
outf.close()

#################

print 'Read U and J matrices...'

Uscr = np.zeros((norb,norb,Nw))
Jscr = np.zeros((norb,norb,Nw))
Ubare = np.zeros((norb,norb,Nw))
Jbare = np.zeros((norb,norb,Nw))

UscrIn = open('Uscr.dat','r')
JscrIn = open('Jscr.dat','r')
UbareIn = open('Ubare.dat','r')
JbareIn = open('Jbare.dat','r')

for n in range(Nwmax):
	tmp1 = [float(x) for x in UscrIn.readline().split()[1:] ]
	tmp2 = [float(x) for x in JscrIn.readline().split()[1:] ]
	tmp3 = [float(x) for x in UbareIn.readline().split()[1:] ]
	tmp4 = [float(x) for x in JbareIn.readline().split()[1:] ]
	for m1 in range(norb):
		for m2 in range(norb):
			Uscr[m1,m2,n] = tmp1[m1*norb+m2]
			Jscr[m1,m2,n] = tmp2[m1*norb+m2]
			Ubare[m1,m2,n] = tmp3[m1*norb+m2]
			Jbare[m1,m2,n] = tmp4[m1*norb+m2]

UscrIn.close()
JscrIn.close()
UbareIn.close()
JbareIn.close()

############################################################
print 'Calculate screened parameters from bare ones...'

outf1 = open('Uscr_loc.dat','w')
outf2 = open('Jscr_loc.dat','w')
for n in range(Nwmax):

	Pdiag = np.zeros((norb,norb))
	for m in range(norb):
		Pdiag[m,m] = P[m,m,n]

	Uscrloc = np.dot( Ubare[:,:,n] , np.linalg.inv( np.identity(norb) - np.dot(Pdiag[:,:], Ubare[:,:,n])  )  )

	outf1.write( str(2*n*np.pi/beta) + '\t'  )
	outf2.write( str(2*n*np.pi/beta) + '\t'  )
	for m1 in range(norb):
		for m2 in range(norb):
	
			Jscrloc = Jbare[m1,m2,n] / ( 1 - P[m1,m2,n]*Jbare[m1,m2,n] )

			outf1.write( str(Uscrloc[m1,m2]) + '\t' )
			outf2.write( str(Jscrloc) + '\t' )

	outf1.write( '\n' )
	outf2.write( '\n' )
outf1.close()
outf2.close()

###########################################################

print 'Read DMFT bare Umatrix and F0(ivn) and calculate the effective screened interaction...'

F0 = np.zeros((Nw))
inf = open('F0_F2F4_pade.dat','r')
for n in range(Nw):
	F0[n] = float( inf.readline().split()[1] )
inf.close()

Ubare = np.zeros((norb,norb,Nw))
Jbare = np.zeros((norb,norb,Nw))
inf = open('umatrix.dat','r')
for m1 in range(norb):	
	m1up = [float(x) for x in inf.readline().split() ]
	m1dn = [float(x) for x in inf.readline().split() ]

	for n in range(Nw):
		#Ubare[m1,m1,n] = m1up[m1*2+1]
		#Jbare[m1,m1,n] = m1up[m1*2+1]

		for m2 in range(norb):	
			Ubare[m1,m2,n] = m1dn[m2*2]                  #orb1dn <-> orb2up
			Jbare[m1,m2,n] = m1up[m2*2+1] - m1up[m2*2+0] 

			Ubare[m1,m2,n] += -15.4536 + F0[n]
			if (m1==m2):
				Jbare[m1,m2,n] += -15.4536 + F0[n]
inf.close()

print 'Ubare[inf] obtained from umatrix.dat: \n',Ubare[:,:,-1]
print 'Jbare[inf] obtained from umatrix.dat: \n',Jbare[:,:,-1]
print 'Ubare[0] obtained from umatrix.dat: \n',Ubare[:,:,0]
print 'Jbare[0] obtained from umatrix.dat: \n',Jbare[:,:,0]

outf1 = open('Uscr_loc_F0scr.dat','w')
outf2 = open('Jscr_loc_F0scr.dat','w')
for n in range(Nwmax):

	Pdiag = np.zeros((norb,norb))
	for m in range(norb):
		Pdiag[m,m] = P[m,m,n]

	Uscrloc = np.dot( Ubare[:,:,n] , np.linalg.inv( np.identity(norb) - np.dot(Pdiag[:,:], Ubare[:,:,n])  )  )

	outf1.write( str(2*n*np.pi/beta) + '\t'  )
	outf2.write( str(2*n*np.pi/beta) + '\t'  )
	for m1 in range(norb):
		for m2 in range(norb):
	
			Jscrloc = Jbare[m1,m2,n] / ( 1 - P[m1,m2,n]*Jbare[m1,m2,n] )

			outf1.write( str(Uscrloc[m1,m2]) + '\t' )
			outf2.write( str(Jscrloc) + '\t' )

	outf1.write( '\n' )
	outf2.write( '\n' )
outf1.close()
outf2.close()

