from pylab import *
from matplotlib import patches
def crlbUsblAlone(x,y,xb,yb,dtUsbl,varr,vara,varv,vel,varhdg):
    T = dtUsbl
    return crlbUsblOdo(T,x,y,xb,yb,dtUsbl,varr,vara,varv,vel,varhdg)
def crlbUsblOdoConverge(x,y,xb,yb,dtUsbl,varr,vara,varv,vel,varhdg,
                        tol=0.01):
    T = dtUsbl
    crlb,cep,srss0 = crlbUsblOdo(T,x,y,xb,yb,
                                 dtUsbl,varr,vara,varv,vel,varhdg)
    dSrss = 1
    N = 1
    while (dSrss > tol) and (N < 100):
        T += dtUsbl
        crlb,cep,srss = crlbUsblOdo(T,x,y,xb,yb,
                                    dtUsbl,varr,vara,varv,vel,varhdg)
        dSrss = abs( (srss0-srss)/srss0 )
        srss0 = srss
        N+=1
    return (crlb,cep,srss,N)
        
def crlbUsblOdo(T,x,y,xb,yb,dtUsbl,varr,vara,varv,vel,varhdg):
    '''
    INPUTS
    T - Duration of observation window. Number of obs. determined by dtUsbl
    x,y - location of observations 
    xb,yb - location of USBL origin
    dtUsbl - sampling time of USBL
    varr - USBL range variance
    vara - USBL angle variance
    varv - velocity variance
    vel - mobile velocity
    varhdg - heading velocity
    
    OUTPUTS:
    CRLBxy - 2x2 array, cov. of x and y
    CEP - circular error probable
    SRSS - sqrt sum of sqares for diag of CRLBxy
    '''
    # Number of USBL measurement
    NUsbl = int(floor(T/dtUsbl))
    Ns = 2*NUsbl      # number of states
    # USBL Range measurements
    Cr = zeros((NUsbl,Ns))
    r = sqrt((x-xb)**2+(y-yb)**2)
    dy = y-yb
    dx = x-xb
    for ii in range(NUsbl):
        Cr[ii,ii*2] = dx/r
        Cr[ii,ii*2+1] = dy/r
    # USBL Angle
    Ca = zeros((NUsbl,Ns))
    for ii in range(NUsbl):
        Ca[ii,ii*2]=-dy/(r**2)
        Ca[ii,ii*2+1]=dx/(r**2)
    # DVL and Hdg
    No = 2*(NUsbl-1)  # number of odo. observations
    #Co = zeros((NUsbl-1,2*NAsbl))
    C1 = hstack((zeros((No,2)),eye(No)))
    C2 = hstack((-1*eye(No),zeros((No,2))))
    Co = C1+C2
    ddist = dtUsbl*vel
    varoy = varv*dtUsbl
    varox = varv*dtUsbl+varhdg*ddist
    ro = zeros(No)
    for ii in range(NUsbl-1):
        ro[ii*2]=varox
        ro[ii*2+1]=varoy

    # Combine
    C = vstack((Cr,Ca,Co))
    R = diag(hstack((ones(NUsbl)*varr,
                     ones(NUsbl)*vara,
                     ro)))

    crlb = inv(dot(dot(C.T,inv(R)),C))

    # Extract middle of the cov. matrix
    NN = crlb.shape[0]/2
    if NUsbl%2:  # odd
        ll = NUsbl-1
    else:
        ll = NUsbl-2
    
    # metrics
    crlbxy = crlb[ll:ll+2,ll:ll+2]
    cep = sum(sqrt(eigvals(crlbxy)))*0.59
    srss = sqrt(sum(diag(crlbxy)))

    return (crlbxy,cep,srss)

def plotErrEllipse(ax,covm,ecenter=(0,0),chi=4.61):
    '''
    Chi-square distribution for 2 dof
    chi = 4.61 90%  Default
    chi = 5.99 95% 
    chi = 6.17 two sigma
    '''

    d,V = eig(inv(covm))
    majMag = sqrt(abs(chi/d[1]))
    minMag = sqrt(abs(chi/d[0]))
    th = math.atan2(V[1,1],V[0,1])
    e1 = patches.Ellipse(ecenter,width=majMag,height=minMag,angle=th*180.0/pi,linewidth=2,fill=False)
    ax.add_patch(e1)
    #ax.set_aspect('equal')
    ax.autoscale()


