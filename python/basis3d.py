import numpy as np

def A(r1,r2):
#	ret= 0.0
#	if ( abs(r1-r2)<0.0001):
#		ret = 1.0/dx
#	return r1*r2 + r1+r2 - r1**2
#	return np.exp((-r1**2-r2**2)*2)
#	return 1/(0.5+abs(r1-r2))
	return 1.0

def psi1(t,p):
	return np.sin(2*p)*np.sin(t)**2
def psi2(t,p):
	return np.cos(p)*np.sin(t)*np.cos(t)
def psi3(t,p):
	return np.sin(p)*np.sin(t)*np.cos(t)


##########################################
##########################################
# 2-particle basis stuff
##########################################

def getOpElement2pb(i,j,k,l,Op,basis):
	ret = 0.0
	for x in range(Nx):
		r1 = a + x*dx
		for y in range(Nx):
			r2 = a + y*dx
		
			ret +=   basis[i](r1)*basis[j](r2) * Op(r1,r2) * basis[l](r2)*basis[k](r1)
	return ret*dx*dx

def evalOp2pb(r1,r2,Optensor,basis,nb):
	ret = 0.0
	for i in range(nb):
		for j in range(nb):
			for k in range(nb):
				for l in range(nb):
					ret += basis[i](r1)*basis[j](r2) * Optensor[i,j,k,l] * basis[l](r2)*basis[k](r1)
	return ret

##########################################
##########################################
# Product Basis stuff
##########################################
def Bf1(t,p):
	return psi1(t,p)*psi1(t,p)/np.sqrt(2.024)
def Bf2(t,p):
	return psi1(t,p)*psi2(t,p)/np.sqrt(0.1928)
def Bf3(t,p):
	return psi2(t,p)*psi1(t,p)/np.sqrt(0.1928)
def Bf4(t,p):
	return psi2(t,p)*psi2(t,p)/np.sqrt(0.175)
##########################
def Bf5(t,p):
	return psi3(t,p)*psi1(t,p)/np.sqrt(0.1928)
def Bf6(t,p):
	return psi3(t,p)*psi2(t,p)/np.sqrt(0.0578)
def Bf7(t,p):
	return psi3(t,p)*psi3(t,p)/np.sqrt(0.1735)
def Bf8(t,p):
	return psi2(t,p)*psi3(t,p)/np.sqrt(0.0578)
def Bf9(t,p):
	return psi1(t,p)*psi3(t,p)/np.sqrt(0.1928)

def getOpElement(i,j,Op,basis,rot,nb):
	ret = 0.0
	for x in range(Nx):
		r1 = a + x*dx
		for y in range(Nx):
			r2 = a + y*dx
		
			ret +=  rotbasis(r1,i,basis,rot) * Op(r1,r2) * rotbasis(r2,j,basis,rot)
	return ret*dx*dx


def evalOp(r1,r2,Opmat,basis,rot,nb):
	ret = 0.0
	for i in range(nb):
		for j in range(nb):
			ret += rotbasis(r1,i,basis,rot) * Opmat[i,j] * rotbasis(r2,j,basis,rot)
	return ret

def rotbasis(x,i,basis,rot):
	ret = 0.0
	for j in range(len(basis)):
		ret += rot[j,i] * basis[j](x)
	return ret

def getOmatrix(basis,nb):
	omat = np.zeros((nb,nb))
	for i in range(nb):
		print np.round(i*100.0/nb,1),'% done'
		for j in range(nb):
		#if 1==1:
		#	j = i
			for it in range(Nx):
				for ip in range(Nx):
					t = theta_start +it*dtheta
					p = phi_start +ip*dphi

					omat[i,j] += basis[i](t,p)*basis[j](t,p)
	return omat*dphi*dtheta

def getOmatrixRot(basis,rot,nb):
	omat = np.zeros((nb,nb))
	for i in range(nb):
		for j in range(nb):

			for k in range(Nx):
				x = a +k*dx
				omat[i,j] += rotbasis(x,i,basis,rot)*rotbasis(x,j,basis,rot)
	return omat*dx

def getIdentityrrp(r1,r2,basis,rot,nb):
	ret = 0.0
	for i in range(nb):
		ret += rotbasis(r1,i,basis,rot)*rotbasis(r2,i,basis,rot)
	return ret
##########################################
##########################################

########################################################
Nx = 100
phi_start = 0.0
phi_end = 2*np.pi
dphi = (phi_end-phi_start)/(Nx-1)

theta_start = 0.0
theta_end = np.pi
dtheta = (theta_end-theta_start)/(Nx-1)

basis1 = np.array([ Bf1,Bf2,Bf3,Bf4,Bf5,Bf6,Bf7,Bf8,Bf9 ])
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
exit()
Omatrix = getOmatrixRot(basis1,rot, nobasis)

print 'Omatrix2: \n',Omatrix.round(5)

Opmatrix = np.zeros((nobasis,nobasis))
print 'Calculate Operator in product basis representation...'
for i in range(nobasis):
	for j in range(nobasis):
		Opmatrix[i,j] = getOpElement(i,j,A,basis1,rot,nobasis)
	print np.round((i+1)*100.0/nobasis,1),'% done'

print 'Opmatrix: \n',Opmatrix.round(5)
#print 'inverse Opmatrix: \n', np.linalg.inv(Opmatrix)

exit()

##################################################################
########################################################

print 'Calculate 2-particle tensor...'
nb2p = len(basis2)
Opmatrix2p = np.zeros((nb2p,nb2p,nb2p,nb2p))
for i in range(nb2p):
	print np.round(i*100.0/nb2p,1),'% done'
	for j in range(nb2p):
		for k in range(nb2p):
			for l in range(nb2p):
				Opmatrix2p[i,j,k,l] = getOpElement2pb(i,j,k,l,A,basis2)
np.set_printoptions(precision=3,suppress=True)
print '2-part Tensor: \n', Opmatrix2p

##################################################################
##################################################################

outf = open('Ar.dat','w')
fac = 3
sqrtdiff = 0.0
maxdiff = 0.0
reldiff = 0.0
print 'Writing A(r,r\') on a ',Nx/fac,'x',Nx/fac,' Grid...'
for i in range(Nx/fac):
	#print np.round(i*100.0*3/Nx,1),'% done writing to file'
	for j in range(Nx/fac):
		x = a+i*dx*fac
		y = a+j*dx*fac

		opvalue = evalOp( x,y, Opmatrix,basis1,rot,nobasis)
		opvalue2pb = evalOp2pb(x,y,Opmatrix2p,basis2,nb2p)
		sqrtdiff += (opvalue-A(x,y))**2
		if maxdiff<abs(opvalue-A(x,y)):
			maxdiff = abs(opvalue-A(x,y))
			reldiff = abs(opvalue-A(x,y))/abs(opvalue)

		outf.write( str(x) + '\t' + str(y) + '\t' + str(A(x,y)) + '\t' \
                + str(opvalue) + '\t'                                \
                + str(opvalue2pb) + '\t'                             \
                + str(getIdentityrrp(x,y,basis1,rot,nobasis)) + '\n' )
	outf.write('\n')
outf.close()

print 'Avg.sq.dev: ',sqrtdiff/(Nx/fac)**2
print 'Max.dev: ',maxdiff,' (', np.round(reldiff*100.0,3) ,'% rel.dev.)'

