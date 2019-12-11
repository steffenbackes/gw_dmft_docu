import numpy as np

t = 1
U = 2.0
nw = 1000
wmin = -6.0
wmax = +6.0
delta = 0.03

dw = (wmax-wmin)/nw
outf = open('bathhyb.dat','w')
for n in range(nw):
	w = wmin + n*dw + 1.0j*delta

	Ham = np.array([[0,-t],[-t,0]])

	sig00 = (U**2/4) * w/(w**2-9*t**2)
	sig01 = (U**2/4) * 3*t/(w**2-9*t**2)
	sig = np.array([[sig00,sig01],[sig01,sig00]])

	G = np.linalg.inv( w*np.identity(2) - Ham - sig )
	Gdmft = np.linalg.inv( w*np.identity(2) - Ham - np.diag(np.diag(sig)) )

	hyb = t**2 * ( G[1,1] - G[1,0]*G[0,1]/G[0,0] )
	hyb_dmft = t**2 * ( Gdmft[1,1] - Gdmft[1,0]*Gdmft[0,1]/Gdmft[0,0] )

	bath = 1.0/( w - hyb )
	bath_dmft = 1.0/( w - hyb_dmft )

	###########
	bath_dyson = 1.0/( 1.0/G[0,0] + sig00 )
	hyb_dyson = w - 1.0/bath_dyson

	outf.write(str(w.real) + '\t')	

	outf.write(str(bath.real) + '\t' + str(bath.imag) + '\t')	
	outf.write(str(hyb.real) + '\t' + str(hyb.imag) + '\t')	
	outf.write(str(bath_dyson.real) + '\t' + str(bath_dyson.imag) + '\t')	
	outf.write(str(hyb_dyson.real) + '\t' + str(hyb_dyson.imag) + '\t')	
	outf.write(str(bath_dmft.real) + '\t' + str(bath_dmft.imag) + '\t')	
	outf.write(str(hyb_dmft.real) + '\t' + str(hyb_dmft.imag) + '\t')	

	outf.write('\n')	

outf.close()
