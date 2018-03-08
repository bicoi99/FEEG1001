import matplotlib.pyplot as pl
import numpy as np
from scipy.optimize import fsolve

def choose(n, k):
    if 0 <= k <= n:
        p = 1
        for t in range(0,min(k, n - k),1):
            p = (p * (n - t)) // (t + 1)
        return p
    else:
        return 0

def mls(x, X, y, sigma): #1D quadratic moving least squares
    N = max(np.shape(y))
    weights = np.zeros(N)
    A = np.zeros([N, 3])
    A[:, 0] = X**2
    A[:, 1] = X
    A[:, 2] = np.ones([1, N])
    for i in range(0, N):
        weights[i] = np.exp(-np.sum((x - X[i])**2) / (2 * sigma))
    W = np.diag(weights)
    a = np.linalg.lstsq(np.dot(np.dot(A.conj().T, W), A), np.dot(np.dot(A.conj().T, W), y) )
    f = a[0][0] * x**2 + a[0][1] * x + a[0][2]
    return f

def mls_error(sigma, X, y): #1D quadratic moving least squares cross-validation
    y_test = np.zeros(len(y))
    error = np.zeros(len(y))
    for i in range(0,len(y)):
        y_test[i] = mls(X[i], np.append(X[0: i], X[i+1:-1]), np.append(y[0:i], y[i+1:-1]), sigma)
        error[i] = (y[i]-y_test[i])**2
    sum_error = sum(error)
    return sum_error

def mls_curve_fit(w_array,cd,x):
    sigma_best = fsolve(mls_error,0.5,args=(w_array,cd)) #fit moving least squares
    w_fine = np.linspace(np.min(w_array) ,np.max(w_array) ,101)
    y_pred= np.zeros(101)
    for i in range(0,101):
        y_pred[i] = mls(w_fine[i], w_array, cd, sigma_best)   
    pl.plot(w_array, cd, 'o', label='data points')
    pl.plot(w_fine, y_pred, label='MLS fit')
    pl.legend(bbox_to_anchor=(1.05, 1), loc=2)
    f = mls(x, w_array, cd, sigma_best)
    return f

def vortex_panel(pointsDef,alpha_deg,plot):
    xb=pointsDef[:,0]
    yb=pointsDef[:,1]
    npanel=np.shape(xb)[0]-1
    xc=np.zeros(npanel)
    yc=np.zeros(npanel)
    ds=np.zeros(npanel)
    dx=np.zeros(npanel)
    dy=np.zeros(npanel)
    theta=np.zeros(npanel)
    InfluenceMat=np.zeros((npanel+1,npanel+1))
    TangentialMat=np.zeros((npanel,npanel+1))
    R=np.zeros((npanel+1))
    psi=np.zeros((npanel+1))
    U=1;
    alpha=alpha_deg*np.pi/180.0;
    for i in range(0,npanel):
        xc[i]=(xb[i]+xb[i+1])/2.0
        yc[i]=(yb[i]+yb[i+1])/2.0
        ds[i]=((xb[i+1]-xb[i])**2+(yb[i+1]-yb[i])**2)**0.5
        dx[i]=xb[i+1]-xb[i]
        dy[i]=yb[i+1]-yb[i]
        if xc[i]<xb[i]:
            theta[i]=-np.arcsin((yb[i+1]-yb[i])/ds[i])
        else:
            theta[i]=np.arcsin((yb[i+1]-yb[i])/ds[i])-np.pi
    for i in range(0,npanel):
        for j in range(0,npanel):
            S=ds[j];
            sinij=np.sin(theta[i]-theta[j])
            cosij=np.cos(theta[i]-theta[j])
            sinji=np.sin(theta[j]-theta[i])
            cosji=np.cos(theta[j]-theta[i])
            # work in panel frame of ref
            xt=xc[i]-xb[j]
            yt=yc[i]-yb[j]
            # rotate
            xpc=xt*np.cos(theta[j])+yt*np.sin(theta[j])
            ypc=-xt*np.sin(theta[j])+yt*np.cos(theta[j])
        
            xt=xb[j+1]-xb[j]
            yt=yb[j+1]-yb[j]
            # rotate
            xpc2=xt*np.cos(theta[j])+yt*np.sin(theta[j])
            
            R1=(xpc**2+ypc**2)**0.5
            R2=((xpc-xpc2)**2+ypc**2)**0.5
            B1=np.arctan(ypc/xpc)
            B2=np.arctan(ypc/(xpc-xpc2))
            
            if ypc<0 and xpc<0:
                B1=B1-np.pi
            elif ypc>0 and xpc<0:
                B1=B1+np.pi
            if ypc<0 and (xpc-xpc2)<0:
                B2=B2-np.pi
            elif ypc>0 and (xpc-xpc2)<0:
                B2=B2+np.pi           
            B=B2-B1
            Ustar=-np.log(R2/R1)*(1/(2.0*np.pi))
            Vstar=B*(1/(2.0*np.pi))
            Ustar_v=Vstar;
            Vstar_v=-Ustar;
            if i==j:
                Ustar=0.0
                Vstar=0.5
                Ustar_v=0.5;
                Vstar_v=0;
            InfluenceMat[i,j]=-sinij*Ustar+cosij*Vstar
            InfluenceMat[i,npanel]=InfluenceMat[i,npanel]-sinij*Ustar_v+cosij*Vstar_v
            if i==npanel-1 or i==0:
                InfluenceMat[npanel,j]=InfluenceMat[npanel,j]+cosji*Ustar-sinji*Vstar
                InfluenceMat[npanel,npanel]=InfluenceMat[npanel,npanel]+cosji*Ustar_v-sinji*Vstar_v             
           
            TangentialMat[i,j]=cosji*Ustar-sinji*Vstar
            TangentialMat[i,npanel]=TangentialMat[i,npanel]+cosji*Ustar_v-sinji*Vstar_v
            
    R[0:-1]=U*np.sin(theta-alpha)
    R[-1]=-U*np.cos(theta[0]-alpha) - U*np.cos(theta[-1]-alpha)
    q, r=np.linalg.qr(InfluenceMat)
    p=np.dot(q.T,R)
    psi=np.dot(np.linalg.inv(r), p)
    vt=(U*np.cos(theta-alpha))+np.dot(TangentialMat,psi)
    cp=1.0-(vt/U)**2
    cx=sum(cp*dy)
    cy=sum(cp*dx)
    cl=-(cy*np.cos(alpha)-cx*np.sin(alpha))
    if plot==1:
        pl.plot(xc,-cp,'k')
        pl.xlabel('x/c')
        pl.ylabel(r'$-C_P$')
        pl.xlim(-0.1, 1.1)
    return cl,cp,xc