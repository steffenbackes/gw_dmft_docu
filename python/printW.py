
def W(a,b,c,d):
	if ( (a==c and b==d) or (a==d and b==c) or (a==b and c==d) ):
		return 1
	else:
		return 0

orbs = ['a','b','c']
for i in orbs:
	for j in orbs:
		for k in orbs:
			for l in orbs:
				m1 = 'a'
				m2 = 'a'
				m3 = 'b'
				m4 = 'b'
				if ( W(m1,i,m3,j)==1 and W(k,m2,l,m4)==1 ):
					print 'W'+str(m1)+str(i)+str(m3)+str(j)+' P'+str(i)+str(j)+str(k)+str(l)+' v'+str(k)+str(m2)+str(l)+str(m4)
