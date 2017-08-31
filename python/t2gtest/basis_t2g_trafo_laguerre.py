import numpy as np
import scipy.special as spsp

#################################
# prefactor for generalized Laguerre polynomial, same for all orbs
# L^alpha_n = 1 for 3d orbs
LagPref = np.sqrt((2.0/3)**3 /( 6*120.0 ))*10

###################################
# normalized one-particle basis
###################################
def psi1(r,t,p):
	return np.sin(2*p)*np.sin(t)**2 / ( 0.25*np.pi*np.sqrt(6.0) )      * LagPref*np.exp(-r/3.0)*(r*2/3.0)**2 # 3dxy
def psi2(r,t,p):
	return np.cos(p)*np.sin(t)*np.cos(t) / ( 0.25*np.pi*np.sqrt(2.0) ) * LagPref*np.exp(-r/3.0)*(r*2/3.0)**2 #3dyz
def psi3(r,t,p):
	return np.sin(p)*np.sin(t)*np.cos(t) / ( 0.25*np.pi*np.sqrt(2.0) ) * LagPref*np.exp(-r/3.0)*(r*2/3.0)**2 #3dxz

##########################################
# Product Basis stuff
##########################################
def Bf1(r,t,p):
	return psi1(r,t,p)*psi1(r,t,p)
def Bf2(r,t,p):
	return psi2(r,t,p)*psi2(r,t,p)
def Bf3(r,t,p):
	return psi3(r,t,p)*psi3(r,t,p)
##########################
def Bf4(r,t,p):
	return psi1(r,t,p)*psi2(r,t,p)
def Bf5(r,t,p):
	return psi1(r,t,p)*psi3(r,t,p)
def Bf6(r,t,p):
	return psi2(r,t,p)*psi3(r,t,p)
############################
def Bf7(r,t,p):
	return psi2(r,t,p)*psi1(r,t,p)
def Bf8(r,t,p):
	return psi3(r,t,p)*psi1(r,t,p)
def Bf9(r,t,p):
	return psi3(r,t,p)*psi2(r,t,p)


def rotbasis(r,t,p, i,basis,rot):
	ret = 0.0
	for j in range(len(basis)):
		ret += rot[j,i] * basis[j](r,t,p)
	return ret

def getOmatrix(basis,nb):
	omat = np.zeros((nb,nb))
	for i in range(nb):
		print np.round(i*100.0/nb,1),'% done'
		for j in range(nb):
			for it in range(Nx):
				for ip in range(Nx):
					for ir in range(Nr):
						t = theta_start +it*dtheta
						p = phi_start +ip*dphi
						r = ir*dr

						omat[i,j] += basis[i](r,t,p)*basis[j](r,t,p)
			print omat[i,j]
	
	return omat*dphi*dtheta*dr

def getOmatrixRot(basis,rot,nb):
	omat = np.zeros((nb,nb))
	for i in range(nb):
		print np.round(i*100.0/nb,1),'% done'
		for j in range(nb):
			for it in range(Nx):
				for ip in range(Nx):
					for ir in range(Nr):
						t = theta_start +it*dtheta
						p = phi_start +ip*dphi
						r = ir*dr

						omat[i,j] += rotbasis(r,t,p,i,basis,rot) * rotbasis(r,t,p,j,basis,rot)
	return omat*dphi*dtheta*dr


def getTrafoElement(i,j,a, basisProd, basis2Part, rot):
	elem = 0.0
	for it in range(Nx):
		for ip in range(Nx):
			for ir in range(Nr):
				t = theta_start +it*dtheta
				p = phi_start +ip*dphi
				r = ir*dr

				elem += basis2Part[i](r,t,p) * rotbasis(r,t,p,a,basisProd,rot) * basis2Part[j](r,t,p)

	return elem*dphi*dtheta*dr


########################################################
Nr = 30
Rend = 25.0
dr = Rend/Nr

Nx = 30
phi_start = 0.0
phi_end = 2*np.pi
dphi = (phi_end-phi_start)/(Nx-1)

theta_start = 0.0
theta_end = np.pi
dtheta = (theta_end-theta_start)/(Nx-1)

#basis1 = np.array([ Bf1,Bf2,Bf3,Bf4,Bf5,Bf6,Bf7,Bf8,Bf9 ])
basis1 = np.array([ Bf1,Bf2,Bf3,Bf4,Bf5,Bf6 ])
basis2 = np.array([ psi1,psi2,psi3 ])
##################################################################
# now product basis

Omatrix = getOmatrix(basis1, len(basis1))
print 'Omatrix: \n',Omatrix.round(4)

w,v = np.linalg.eigh(Omatrix)
print 'Eigenvalues:\n',w.round(4)
print 'Eigenvectors:\n',v.round(4)

delitems = []
for i in range(len(w)):
	if abs(w[i]) < 0.0001:
		delitems.append(i)
		print 'Remove Eigenvector ',i,' with Eigenvalue ',w[i]
w = np.delete(w, delitems)
v = np.delete(v, delitems, 1)

nobasis = len(w)
D=np.diag((1.0/np.sqrt(w)))
print 'We have reduced the basis to ',nobasis,' elements'
print 'Eigenvalues 1/sqrt:\n',D.round(5)
print 'Remaining Eigenvectors:\n',v.round(5)

rot = np.dot(v,D)
print 'Rotation matrix: \n', rot.round(5)

Omatrix = getOmatrixRot(basis1,rot, nobasis)
print 'Omatrix2: \n',Omatrix.round(5)

Pmat = np.zeros(( len(basis2)**2, len(basis1)  ))
for i in range(len(basis2)):
	for j in range(len(basis2)):
		print np.round((i*len(basis2)+j)*100.0/len(basis2)**2,1),'% done'
		for a in range(len(basis1)):
			Pmat[i*len(basis2)+j,a] = getTrafoElement(i,j,a, basis1, basis2, rot)
print 'Final Pmat: \n', Pmat.round(5)

##############################################

print 'sum_a P[ii,a]*P[jj,a]:'
for i in range(len(basis2)):
	for j in range(len(basis2)):
		tmp = 0.0
		for a in range(len(basis1)):
			tmp += Pmat[i*len(basis2)+j,a]*Pmat[i*len(basis2)+j,a]

	print 'i:',i,' j:',j,'  val = ',tmp


################################################

print 'Read DFT Greens function and calculate polarization...'
print 'Only read a single orbital...'
Nw = 3000
beta = 40.0
norb = len(basis2)

Gf = np.zeros((Nw),dtype=complex)
P = np.zeros((Nw))

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

outf = open('P.dat','w')
Nwmax = 500
for n in range(Nwmax):
#	print n*100.0/Nwmax,'% done...'
	for m in range(-Nw+1,Nw):
		P[n] += ( getgf(m)*getgf(m-n) ).real
	P[n] /= beta
	outf.write( str(2*n*np.pi/beta) + '\t' + str(P[n]) + '\n' )
outf.close()

#################
print 'Read F0 and F2+F4/14...'

F = np.zeros(( Nwmax,2  ))
inf = open('F0_F2F4_pade.dat','r')
for n in range(Nwmax):
	data = [float(x) for x in inf.readline().split()]
	F[n,0] = data[1]
	F[n,1] = data[2]
inf.close()

Wmat = np.zeros(( Nwmax, norb,norb,norb,norb  ))
Polmat = np.zeros(( Nwmax, norb,norb,norb,norb  ))
print 'Construct Wmatrix and Polarization matrix...'

for n in range(Nwmax):
	for i in range(norb):
		for j in range(norb):
			for k in range(norb):
				for l in range(norb):
			
					# Wmatrix
					if (i==j and j==k and k==l):
						Wmat[n,i,j,k,l] = F[n,0] + 0.3           # U0
					elif (i==k and j==l ):
						Wmat[n,i,j,k,l] = F[n,0]+0.3 - 2*F[n,1]  # Up
					elif (i==l and j==k ):
						Wmat[n,i,j,k,l] = F[n,1]                 # J

					# Pmatrix
					if (i==k and j==l ):
						Polmat[n,i,j,k,l] = P[n]

#######################################################

print 'Now convert Wmat and Polmat into product basis'
WmatProd = np.zeros(( Nwmax, nobasis,nobasis  ))
PolmatProd = np.zeros(( Nwmax, nobasis,nobasis  ))

for a in range(nobasis):
	for b in range(nobasis):
		for n in range(Nwmax):
			for i in range(norb):
				for j in range(norb):
					for k in range(norb):
						for l in range(norb):

							WmatProd[n,a,b] += Pmat[i*norb+k,a] * Wmat[n,i,j,k,l] * Pmat[l*norb+j,b]

							PolmatProd[n,a,b] += Pmat[i*norb+k,a] * Polmat[n,i,j,k,l] * Pmat[l*norb+j,b]

print 'WmatProd[0]:'
print np.round(WmatProd[0,:,:],4)

print 'PolmatProd[0]:'
print np.round(PolmatProd[0,:,:],4)

#####################################################

print 'Do Dyson to get U(iwn)'

UmatProd = np.zeros(( Nwmax, nobasis,nobasis  ))

for n in range(Nwmax):
	UmatProd[n,:,:] = np.linalg.inv(  np.linalg.inv(WmatProd[n,:,:]) + PolmatProd[n,:,:]   )
print 'UmatProd[0]:'
print np.round(UmatProd[0,:,:],4)

###################################################

print 'Transform WmatProd and UmatProd into ijkl basis...'

print ' W[0]:', np.round( Wmat[0,0,0,0,0], 4 )
print 'Wp[0]:', np.round( Wmat[0,0,1,0,1], 4 )
print 'WJ[0]:', np.round( Wmat[0,0,1,1,0], 4 )
print ''

Umat = np.zeros(( Nwmax, norb,norb,norb,norb  ))
Wmat = np.zeros(( Nwmax, norb,norb,norb,norb  ))
for i in range(norb):
	for j in range(norb):
		for k in range(norb):
			for l in range(norb):
				for a in range(nobasis):
					for b in range(nobasis):
						for n in range(Nwmax):
							Umat[n,i,j,k,l] += Pmat[i*norb+k,a] * UmatProd[n,a,b] * Pmat[l*norb+j,b]
							Wmat[n,i,j,k,l] += Pmat[i*norb+k,a] * WmatProd[n,a,b] * Pmat[l*norb+j,b]

print ' W[0]:', np.round( Wmat[0,0,0,0,0], 4 )
print 'Wp[0]:', np.round( Wmat[0,0,1,0,1], 4 )
print 'WJ[0]:', np.round( Wmat[0,0,1,1,0], 4 )
print ''

print ' U[0]:', np.round( Umat[0,0,0,0,0], 4 )
print 'Up[0]:', np.round( Umat[0,0,1,0,1], 4 )
print ' J[0]:', np.round( Umat[0,0,1,1,0], 4 )



