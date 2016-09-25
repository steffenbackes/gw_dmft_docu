import numpy as np

# use 2-dim one-particle basis with np.sin and np.cos
Nx = 100
a=0
b=2*np.pi
dx = (b-a)/Nx

def skalarprod(f,g):
	res = 0.0
	for i in range(Nx):
		x = a+dx*i
		#res += dx*np.conjugate( f(x) )*g(x)
		res += dx*f(x)*g(x)
	return res

def psi1(x):
	return np.sin(x)
def psi2(x):
	return np.cos(x)

def Bf1(x):
	return psi1(x)*psi1(x)
def Bf2(x):
	return psi1(x)*psi2(x)
def Bf3(x):
	return psi2(x)*psi1(x)
def Bf4(x):
	return psi2(x)*psi2(x)

def Bff1(x):
	return psi1(x)*psi1(x)*(-0.70711) + 0.70711*psi2(x)*psi2(x)
def Bff2(x):
	return psi1(x)*psi2(x)*(0.70711) + 0.70711*psi2(x)*psi1(x)
def Bff3(x):
	return psi1(x)*psi2(x)*(-0.70711) + 0.70711*psi2(x)*psi1(x)
def Bff4(x):
	return psi1(x)*psi1(x)*(0.70711) + 0.70711*psi2(x)*psi2(x)

def fillOmatrix(Omatrix):
	norms = np.zeros((4))
	i=0
	for b in [Bf1,Bf2,Bf3,Bf4]:
		norms[i] = skalarprod(b,b)
		i+=1

	i=0
	for b1 in [Bf1,Bf2,Bf3,Bf4]:
		j = 0
		for b2 in [Bf1,Bf2,Bf3,Bf4]:
			Omatrix[i,j] = skalarprod(b1,b2)/np.sqrt(norms[i]*norms[j])
			j+=1
		i+=1

def fillOmatrix2(Omatrix):
	i=0
	for b1 in [Bff1,Bff2,Bff3,Bff4]:
		j = 0
		for b2 in [Bff1,Bff2,Bff3,Bff4]:
			Omatrix[i,j] = skalarprod(b1,b2)
			j+=1
		i+=1

Omatrix = np.zeros((4,4))
fillOmatrix(Omatrix)

print 'Omatrix: \n',Omatrix.round(5)

w,v = np.linalg.eig(Omatrix)
print 'Eigenvalues:\n',w.round(5)
print 'Eigenvectors:\n',v.round(5)

Omatrix2 = np.zeros((4,4))
fillOmatrix2(Omatrix2)

print 'Omatrix2: \n',Omatrix2.round(5)

#print np.sqrt( np.diag( np.append( w[:-1], [0.0] ) ) )
