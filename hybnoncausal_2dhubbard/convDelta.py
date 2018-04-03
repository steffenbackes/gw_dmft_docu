import numpy as np

norb = 1
ntau = 1000
beta = 40.0

##############################################

inf = open('Hyb_fullk.dat','r')
erg = []
delta = []

#tmp = inf.readline()
tmp = inf.readline()
while ( tmp!=''):
	data = [float(x) for x in inf.readline().split()]

	erg.append( data[0] )
	
	delta.append( [0.0 for m in range(norb)] )
	for m in range(norb):
		delta[-1][m] = data[ 2*m +1 ] + data[2*m+1 +1]*1.0j

	#	if delta[-1][m].imag > 0.0:
	#		delta[-1][m] = delta[-1][m].real

	tmp = inf.readline()
inf.close()

outfi = open('delta_tau.dat','w')
nw = len(delta)
for it in range(ntau+1):
	#if ( np.mod(int(it*10.0/ntau),2)==0  ):
#		print it*100.0/ntau,'% done'
	
	tau = it*beta/ntau
	
	outfi.write( str(tau) + '\t' )
	for m in range(norb):
		val = 0.0
		for n in range(1,nw):
			dw = erg[n]-erg[n-1]

			val += ( delta[n][m].imag * np.exp(-tau*erg[n])/(np.exp(-beta*erg[n])+1) ) * dw/np.pi
		outfi.write( str(val) + '\t' + str(val) + '\t')
	outfi.write( '\n' )

outfi.close()



