import numpy as np

mu = 0.0
nw = 100
nk = 25
idelta = 0.05j
wf = []

def eps(kx,ky):
	return -2*( np.cos(kx*2*np.pi/nk) + np.cos(ky*2*np.pi/nk) )

Sigma = np.zeros((nk,nk,nw),dtype=complex)
SigmaHF = np.zeros((nk,nk),dtype=complex)
SigmaW0 = np.zeros((nk,nk),dtype=complex)
Sloc = np.zeros((nw),dtype=complex)
inf = open('sigma_k_n0.5U4.dat','r')
for n in range(nw):
	for kx in range(nk):
		for ky in range(nk):
			data = [float(x) for x in inf.readline().split()]

			if (kx==0 and ky==0):
				wf.append( data[2] )

			Sigma[kx,ky,n] = data[19] + 1.0j*data[20]
			Sloc[n] += Sigma[kx,ky,n]/nk**2

			if (n==0):
				for i in range(6): # Sum all freq-indep. parts
					SigmaHF[kx,ky] += data[3 + 2*i]

		inf.readline()
inf.close()

iw0 = int(-wf[0]/(wf[1]-wf[0])+1)
print wf[iw0]
slocre = 0.0
for kx in range(nk):
	for ky in range(nk):
		slocre += ( Sigma[kx,ky,iw0] ).real/nk**2

		SigmaW0[kx,ky] = (Sigma[kx,ky,iw0]).real

mu = mu + slocre + 0.11

outf = open('gloc.dat','w')
Gloc = np.zeros((nw),dtype=complex)
Gloc2 = np.zeros((nw),dtype=complex)
Gloc3 = np.zeros((nw),dtype=complex)
Gloc4 = np.zeros((nw),dtype=complex)
for n in range(nw):
	for kx in range(nk):
		for ky in range(nk):
			Gloc[n]  += 1.0/( wf[n] + idelta + mu - eps(kx,ky) - Sigma[kx,ky,n] ) / nk**2  # full Sigma
			Gloc2[n] += 1.0/( wf[n] + idelta + mu - eps(kx,ky) - Sloc[n]        ) / nk**2  # only local part
			Gloc3[n] += 1.0/( wf[n] + idelta + mu - eps(kx,ky) - Sloc[n] - SigmaW0[kx,ky] + 2.47 ) / nk**2 # local + kdepSigma at Ef
			Gloc4[n] += 1.0/( wf[n] + idelta + mu - eps(kx,ky) - Sloc[n] - SigmaHF[kx,ky] + 2.7 ) / nk**2 # local + k-dep hartree

	outf.write( str(wf[n]) + '\t' + str(Gloc[n].real) + '\t' + str(Gloc[n].imag) + '\t' + str(Gloc2[n].real) + '\t' + str(Gloc2[n].imag) + '\t' + str(Gloc3[n].real) + '\t' + str(Gloc3[n].imag) + '\t' + str(Gloc4[n].real) + '\t' + str(Gloc4[n].imag) + '\n' )
outf.close()

outf = open('bathhyb.dat','w')
for n in range(nw):
	gbath = 1.0/( 1.0/Gloc[n] + Sloc[n] )
	gbath2 = 1.0/( 1.0/Gloc2[n] + Sloc[n] )
	gbath3 = 1.0/( 1.0/Gloc3[n] + Sloc[n] )

	hyb = wf[n] + idelta - 1.0/gbath
	hyb2 = wf[n] + idelta - 1.0/gbath2
	hyb3 = wf[n] + idelta - 1.0/gbath3

	#outf.write(str(wf[n]) + '\t' + str(gbath.real) + '\t' + str(gbath.imag) + '\t' + str(gbath2.real) + '\t' + str(gbath2.imag) + '\t' + str(hyb.real) + '\t' + str(hyb.imag) + '\t' + str(hyb2.real) + '\t' + str(hyb2.imag) + '\n')
	outf.write(str(wf[n]) + '\t' + str(hyb.real) + '\t' + str(hyb.imag) + '\t' + str(hyb2.real) + '\t' + str(hyb2.imag) + '\t' + str(hyb3.real) + '\t' + str(hyb3.imag) + '\n')
outf.close()

outf = open('bathhyb_general.dat','w')
outf2 = open('bathhyb_df.dat','w')
for n in range(nw):

	epsG = 0.0j
	epsGeps = 0.0j
	for kx in range(nk):
		for ky in range(nk):
			epsG +=  eps(kx,ky)     /( wf[n] + idelta + mu - eps(kx,ky) - Sigma[kx,ky,n] ) / nk**2
			epsGeps +=  eps(kx,ky)**2 /( wf[n] + idelta + mu - eps(kx,ky) - Sigma[kx,ky,n] ) / nk**2

	hyb = epsGeps - epsG**2/Gloc[n]
	gbath = 1.0/( wf[n] + idelta + mu - hyb )

	outf.write(str(wf[n]) + '\t' + str(gbath.real) + '\t' + str(gbath.imag) + '\t' + str(hyb.real) + '\t' + str(hyb.imag) +  '\n')

	p = -( 1 + 2*epsG )/Gloc[n]
	q = epsGeps/Gloc[n]

	hyb1 = 0.5*( -p + np.sqrt( p**2 - 4*q ) )
	hyb2 = 0.5*( -p - np.sqrt( p**2 - 4*q ) )

	if (hyb1.imag>0):
		tmp = hyb1
		hyb1 = hyb2
		hyb2 = tmp
	
	outf2.write(str(wf[n]) + '\t' + str(hyb1.real) + '\t' + str(hyb1.imag) + '\t' + str(hyb2.real) + '\t' + str(hyb2.imag) +  '\n')
outf.close()

