import numpy as np

def A(r1,r2):
#	ret= 0.0
#	if ( abs(r1-r2)<0.0001):
#		ret = 2.0/dx
#	return ret
#	return 2.0
#	return 2.0
	return r1**2 + r2**2+1
#	return np.exp((-r1**2-r2**2)*2)
#	return r1+2*r2

def psi1(x):
	return np.sqrt(0.5)
def psi2(x):
	return x / np.sqrt(2.0/3)


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
def Bf1(x):
	return psi1(x)*psi1(x)
def Bf2(x):
	return psi1(x)*psi2(x)
def Bf3(x):
	return psi2(x)*psi1(x)
def Bf4(x):
	return psi2(x)*psi2(x)

def getOpElement(i,j,Op,basis,rot,nb):
	ret = 0.0
	for x in range(Nx):
		r1 = a + x*dx
		for y in range(Nx):
			r2 = a + y*dx
		
			ret +=  rotbasis(r1,i,basis,rot) * Op(r1,r2) * rotbasis(r2,j,basis,rot)
	return ret*dx*dx

def conv2partToProd(aa,bb,Optensor,prodbasis,rot,p2basis,nb):
	ret = 0.0
	for x in range(Nx):
		r1 = a + x*dx
		for y in range(Nx):
			r2 = a + y*dx
				
			for i in range(len(p2basis)):
				for j in range(len(p2basis)):
					for k in range(len(p2basis)):
						for l in range(len(p2basis)):
							ret +=  rotbasis(r1,aa,prodbasis,rot) \
                      * p2basis[i](r1)*p2basis[j](r2) * Optensor[i,j,k,l] * p2basis[l](r2)*p2basis[k](r1) \
                           * rotbasis(r2,bb,prodbasis,rot)

	return ret*dx*dx

def convProdToTensor(i,j,k,l,Opmatrix,prodbasis,rot,p2basis,nb):
	ret = 0.0
	for x in range(Nx):
		r1 = a + x*dx
		for y in range(Nx):
			r2 = a + y*dx
				
			for aa in range(nb):
				for bb in range(nb):
					ret +=  rotbasis(r1,aa,prodbasis,rot) \
                * p2basis[i](r1)*p2basis[j](r2) * Opmatrix[aa,bb] * p2basis[l](r2)*p2basis[k](r1) \
                           * rotbasis(r2,bb,prodbasis,rot)

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

def getOmatrix(basis,rot,nb):
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
Nx = 50
a = -1.0
b = +1.0 #2*np.pi
dx = (b-a)/(Nx-1)
#basis1 = np.array([ Bf1,Bf2,Bf3,Bf4 ])
basis1 = np.array([ Bf1,Bf2,Bf3,Bf4 ])
basis2 = np.array([ psi1,psi2 ])
########################################################

print 'Calculate 2-particle tensor...'
nb2p = len(basis2)
Optensor = np.zeros((nb2p,nb2p,nb2p,nb2p))
for i in range(nb2p):
	print np.round(i*100.0/nb2p,1),'% done'
	for j in range(nb2p):
		for k in range(nb2p):
			for l in range(nb2p):
				Optensor[i,j,k,l] = getOpElement2pb(i,j,k,l,A,basis2)
		
#				Optensor[i,j,k,l] = 0.0
#				if (i==k and j==l):
#					Optensor[i,j,k,l] = 1.0

np.set_printoptions(precision=3,suppress=True)
print 'This is the 2-part Tensor: \n', Optensor
print 'Now combine left two indices and right two indices and set up a matrix like the DMFT community does...'
print 'I.e. A_ijkl =^= A_(ij)(kl)'
Optensor2matrix = np.zeros((nb2p*nb2p,nb2p*nb2p))
for i in range(nb2p):
	for j in range(nb2p):
		for k in range(nb2p):
			for l in range(nb2p):
				x = i*nb2p + j
				y = k*nb2p + l
				Optensor2matrix[x,y] = Optensor[i,j,k,l]
print 'This is the 2-part Tensor in matrix form: \n', Optensor2matrix
print '\n \n'
##################################################################
##################################################################
# now product basis
print 'Now go to product basis:'
print 'Calculate the Overlap matrix of the product basis...'

Omatrix = getOmatrix(basis1,np.eye(len(basis1)), len(basis1))
print 'Overlap matrix: \n',Omatrix.round(4)

w,v = np.linalg.eigh(Omatrix)
print 'Eigenvalues:\n',w.round(4)
print 'Eigenvectors:\n',v.round(4)

delitems = []
for i in range(len(w)):
	if abs(w[i]) < 0.001:
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
Omatrix = getOmatrix(basis1,rot, nobasis)

print 'Overlap matrix of the new basis: \n',Omatrix.round(5)
print 'Now we have a good product basis!'
##############################################################
##############################################################

Opmatrix = np.zeros((nobasis,nobasis))
print 'Calculate Operator in product basis representation in the standard way from A(r1,r2)...'
for i in range(nobasis):
	for j in range(nobasis):
		Opmatrix[i,j] = getOpElement(i,j,A,basis1,rot,nobasis)
	print np.round((i+1)*100.0/nobasis,1),'% done'

print 'Opmatrix in product basis: \n',Opmatrix.round(5)
#print 'inverse Opmatrix: \n', np.linalg.inv(Opmatrix)

##################################################################
##################################################################

print 'Compare to <=> Convert the Operator tensor in the 2particle basis to product basis...'
OpmatrixConv = np.zeros((nobasis,nobasis))
for i in range(nobasis):
	for j in range(nobasis):
		print np.round((i*nobasis+j)*100.0/(nobasis**2),1),'% done'
		OpmatrixConv[i,j] = conv2partToProd(i,j,Optensor,basis1,rot,basis2,nobasis)

print 'This is the Opmatrix in the product basis generated from A_ijkl \n',OpmatrixConv
##############################################################
##############################################################

print 'Now we have three representations: '
print '1) The operator in the combined index notation: \n', Optensor2matrix
print '2) The operator in the product basis: \n', Opmatrix
print '3) The operator in the product basis converted from 1): \n', OpmatrixConv
print '\n'

print 'Now invert them:'
try:
	OpInv2Part = np.linalg.inv( Optensor2matrix )
	print '1) Inverse in the combined index notation: \n', OpInv2Part
	print 'Convert inverse 1) into tensor and then to product basis...'
	InvTensor = np.zeros((nb2p,nb2p,nb2p,nb2p))
	for i in range(nb2p):
		for j in range(nb2p):
			for k in range(nb2p):
				for l in range(nb2p):
					x = i*nb2p + j
					y = k*nb2p + l
					InvTensor[i,j,k,l] = OpInv2Part[x,y]
	InvOpmatrixConv = np.zeros((nobasis,nobasis))
	for i in range(nobasis):
		for j in range(nobasis):
			print np.round((i*nobasis+j)*100.0/(nobasis**2),1),'% done'
			InvOpmatrixConv[i,j] = conv2partToProd(i,j,InvTensor,basis1,rot,basis2,nobasis)
	print '1) Inverse Tensor in product basis from 1): \n', InvOpmatrixConv
except:
	print '1) Inverse in the combined index notation is singular!'

print '\n'
try:
	OpInvProd = np.linalg.inv( Opmatrix )
	print '2) Inverse in the product basis: \n', OpInvProd
	print 'Convert inverse 2) from product basis to 2particle basis...'
	Inv2Pmatrix = np.zeros((nb2p*nb2p,nb2p*nb2p))
	for i in range(nb2p):
		for j in range(nb2p):
			print np.round( (i*nb2p+j)*100.0/nb2p**2, 1),'% done'
			for k in range(nb2p):
				for l in range(nb2p):
					x = i*nb2p + j
					y = k*nb2p + l
					Inv2Pmatrix[x,y]     = convProdToTensor(i,j,k,l,OpInvProd,basis1,rot,basis2,nobasis)
	print '2) Inverse in 2particle matrix: \n',Inv2Pmatrix
except:
	print '2) Inverse in the product basis is singular!'

print '\n'
try:
	OpInvProdConv = np.linalg.inv( OpmatrixConv )
	print '3) Inverse in the product basis converted from 1): \n', OpInvProdConv
	print 'Convert inverse 3) from product basis to 2particle basis...'
	Inv2PmatrixConv = np.zeros((nb2p*nb2p,nb2p*nb2p))
	for i in range(nb2p):
		for j in range(nb2p):
			print np.round( (i*nb2p+j)*100.0/nb2p**2, 1),'% done'
			for k in range(nb2p):
				for l in range(nb2p):
					x = i*nb2p + j
					y = k*nb2p + l
					Inv2PmatrixConv[x,y] = convProdToTensor(i,j,k,l,OpInvProdConv,basis1,rot,basis2,nobasis)
	print '3) Inverse in 2particle matrix: \n',Inv2PmatrixConv
except:
	print '3) Inverse in the product basis converted from 1) is singular!'

print '\n'

##############################################################
##############################################################


##############################################################
##############################################################
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
		invopvalue = evalOp( x,y, OpInvProd,basis1,rot,nobasis)
		opvalue2pb = evalOp2pb(x,y,Optensor,basis2,nb2p)
		sqrtdiff += (opvalue-A(x,y))**2
		if maxdiff<abs(opvalue-A(x,y)):
			maxdiff = abs(opvalue-A(x,y))
			reldiff = abs(opvalue-A(x,y))/abs(opvalue)

		outf.write( str(x) + '\t' + str(y) + '\t' + str(A(x,y)) + '\t' \
                + str(opvalue) + '\t'                                \
                + str(opvalue2pb) + '\t'                             \
                + str(invopvalue) + '\n' )
#                + str(getIdentityrrp(x,y,basis1,rot,nobasis)) + '\n' )
	outf.write('\n')
outf.close()

print 'Avg.sq.dev: ',sqrtdiff/(Nx/fac)**2
print 'Max.dev: ',maxdiff,' (', np.round(reldiff*100.0,3) ,'% rel.dev.)'

