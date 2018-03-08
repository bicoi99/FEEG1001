import pylab as p
import numpy as np
from aclabtools import choose
        
def bezier(a):
    n=np.shape(a)[0]-1
    b=np.zeros([101,2])
    terms=np.zeros([n+1,2])
    t=np.linspace(0,1,101)
    for i in range(0,101):
        for j in range(0,n+1):
            terms[j,:]=a[j,:]*choose(n,j)*t[i]**j*(1-t[i])**(n-j)
        b[i,:]=sum(terms,0)
    # plot cubic Bezier
    p.plot(b[:,0],b[:,1])
    # plot control points
    p.plot(a[:,0],a[:,1],'ko')
    # plot control polygon
    p.plot(a[:,0],a[:,1],'k')
    return b[:,0:2]
    
def rational_bezier(a,z):
    n=np.shape(a)[0]-1
    b=np.zeros([101,2])
    terms=np.zeros([n+1,2])
    rat_terms=np.zeros([n+1,2])
    t=np.linspace(0,1,101)
    for i in range(0,101):
        for j in range(0,n+1):
            terms[j,:]=a[j,:]*z[j]*choose(n,j)*t[i]**j*(1-t[i])**(n-j)
            rat_terms[j,:]=z[j]*choose(n,j)*t[i]**j*(1-t[i])**(n-j)
        b[i,:]=sum(terms,0)/sum(rat_terms,0)
        #b[i,:]=b[i,:]/b[i,2]
    # plot cubic Bezier
    p.plot(b[:,0],b[:,1])
    # plot control points
    p.plot(a[:,0],a[:,1],'ko')
    # plot control polygon
    p.plot(a[:,0],a[:,1],'k')
    return b