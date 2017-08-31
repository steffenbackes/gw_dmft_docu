import numpy as np

###################################
# normalized one-particle basis
###################################
def psi1(t,p):
	return np.sin(2*p)*np.sin(t)**2 / ( 0.25*np.pi*np.sqrt(6.0) ) # 3dxy
def psi2(t,p):
	return np.cos(p)*np.sin(t)*np.cos(t) / ( 0.25*np.pi*np.sqrt(2.0) ) #3dyz
def psi3(t,p):
	return np.sin(p)*np.sin(t)*np.cos(t) / ( 0.25*np.pi*np.sqrt(2.0) )#3dxz

##########################################
# Product Basis stuff
##########################################
def Bf1(t,p):
	return psi1(t,p)*psi1(t,p)
def Bf2(t,p):
	return psi2(t,p)*psi2(t,p)
def Bf3(t,p):
	return psi3(t,p)*psi3(t,p)
##########################
def Bf4(t,p):
	return psi1(t,p)*psi2(t,p)
def Bf5(t,p):
	return psi1(t,p)*psi3(t,p)
def Bf6(t,p):
	return psi2(t,p)*psi3(t,p)
############################
def Bf7(t,p):
	return psi2(t,p)*psi1(t,p)
def Bf8(t,p):
	return psi3(t,p)*psi1(t,p)
def Bf9(t,p):
	return psi3(t,p)*psi2(t,p)


def rotbasis(t,p, i,basis,rot):
	ret = 0.0
	for j in range(len(basis)):
		ret += rot[j,i] * basis[j](t,p)
	return ret

def getOmatrix(basis,nb):
	omat = np.zeros((nb,nb))
	for i in range(nb):
		print np.round(i*100.0/nb,1),'% done'
		for j in range(nb):
			for it in range(Nx):
				for ip in range(Nx):
					t = theta_start +it*dtheta
					p = phi_start +ip*dphi

					omat[i,j] += basis[i](t,p)*basis[j](t,p)
	return omat*dphi*dtheta

def getOmatrixRot(basis,rot,nb):
	omat = np.zeros((nb,nb))
	for i in range(nb):
		print np.round(i*100.0/nb,1),'% done'
		for j in range(nb):
			for it in range(Nx):
				for ip in range(Nx):
					t = theta_start +it*dtheta
					p = phi_start +ip*dphi

					omat[i,j] += rotbasis(t,p,i,basis,rot) * rotbasis(t,p,j,basis,rot)
	return omat*dphi*dtheta


def getTrafoElement(i,j,a, basisProd, basis2Part, rot):
	elem = 0.0
	for it in range(Nx):
		for ip in range(Nx):
			t = theta_start +it*dtheta
			p = phi_start +ip*dphi

			elem += basis2Part[i](t,p) * rotbasis(t,p,a,basisProd,rot) * basis2Part[j](t,p)

	return elem*dphi*dtheta


########################################################
Nx = 40
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

################################################

exit()

A = np.zeros((len(basis2),len(basis2),len(basis2),len(basis2)))
for i in range(len(basis2)):
	for j in range(len(basis2)):
		A[i,j,i,j] = 100.0
		A[i,j,j,i] = 100.0
#A = np.ones((len(basis2),len(basis2),len(basis2),len(basis2)))

Aprod = np.zeros(( len(basis1), len(basis1)  ))
for a in range(nobasis):
	for b in range(nobasis):
		for i in range(len(basis2)):
			for j in range(len(basis2)):
				for k in range(len(basis2)):
					for l in range(len(basis2)):
						Aprod[a,b] += Pmat[i*len(basis2)+k,a] * A[i,j,k,l] * Pmat[l*len(basis2)+j,b]
print 'Aprod: \n'
for a in range(nobasis):
	for b in range(nobasis):
		print '{:+08.4f}'.format(Aprod[a,b]),'  ',
	print ''


print 'A back to 2Part:'
A2part = np.zeros(( len(basis2), len(basis2),len(basis2),len(basis2)  ))
for i in range(len(basis2)):
	for j in range(len(basis2)):
		for k in range(len(basis2)):
			for l in range(len(basis2)):
				for a in range(nobasis):
					for b in range(nobasis):
						A2part[i,j,k,l] += Pmat[i*len(basis2)+k,a] * Aprod[a,b] * Pmat[l*len(basis2)+j,b]
				print '{:1d} {:1d} {:1d} {:1d} :'.format(i,j,k,l),
				print '{:+08.4f}'.format(A2part[i,j,k,l]),'  ',
				print ''

#################################
Aprod = np.zeros(( len(basis1), len(basis1)  ))
for a in range(nobasis):
	for b in range(nobasis):
		for i in range(len(basis2)):
			for j in range(len(basis2)):
				for k in range(len(basis2)):
					for l in range(len(basis2)):
						Aprod[a,b] += Pmat[i*len(basis2)+k,a] * A2part[i,j,k,l] * Pmat[l*len(basis2)+j,b]
print 'And back to Aprod: \n'
for a in range(nobasis):
	for b in range(nobasis):
		print '{:+08.4f}'.format(Aprod[a,b]),'  ',
	print ''

