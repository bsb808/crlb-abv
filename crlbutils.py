from pylab import *
from matplotlib import patches


def crlbLbl(Xb,X,sigR):
    # Cmatrix
    Cr = zeros((len(Xb),2))
    for xb,ii in zip(Xb,range(len(Xb))):
        r = sqrt((X[0]-xb[0])**2+(X[1]-xb[1])**2)
        r = max((r,0.1))
        for jj in range(2):
            d = X[jj]-xb[jj]
            Cr[ii,jj] = d/r

    # Rmatrix
    R = eye(Cr.shape[0])*(sigR**2)
    lam = dot(dot(Cr.T,inv(R)),Cr)
    try:
        crlb = inv(lam)
    except LinAlgError:
        crlb = diag([100,100])*sigR**2
    cep = sum(sqrt(eigvals(crlb)))*0.59
    srss = sqrt(sum(diag(crlb)))
    return crlb, cep, srss

def crlbUsblAlone(x,y,xb,yb,dtUsbl,varr,vara,varv,vel,varhdg):
    T = dtUsbl
    return crlbUsblOdo(T,x,y,xb,yb,dtUsbl,varr,vara,varv,vel,varhdg)
def crlbUsblOdoConverge(x,y,xb,yb,dtUsbl,varr,vara,varv,vel,varhdg,
                        tol=0.01):
    T = dtUsbl
    crlb,cep,srss0 = crlbUsblOdo(T,x,y,xb,yb,
                                 dtUsbl,varr,vara,varv,vel,varhdg)
    dSrss0 = 1
    dSrss1 = 1
    srss1 = 1e6
    N = 1
    while (max((dSrss0,dSrss1)) > tol) and (N < 100):
        T += dtUsbl
        crlb,cep,srss = crlbUsblOdo(T,x,y,xb,yb,
                                    dtUsbl,varr,vara,varv,vel,varhdg)
        dSrss0 = abs( (srss1-srss)/srss1 )
        dSrss1 = abs( (srss0-srss1)/srss0)

        srss0 = srss1
        srss1 = srss
        N+=1
    if N > 99:
        print("WARNING: crlbutils - covergence req'd %d interations"%N)
    return (crlb,cep,srss,N)

def crlbLblOdoConverge(X,Xb,dtLbl,sigR,sigV,vel,sigHdg,tol=0.01):
    T = dtLbl
    crlb,cep,srss0 = crlbLblOdo(T,X,Xb,dtLbl,sigR,sigV,vel,sigHdg)
    dSrss0 = 1
    dSrss1 = 1
    srss1 = 1e6
    N = 1
    while (max((dSrss0,dSrss1)) > tol) and (N < 100):
        T += dtLbl
        crlb,cep,srss = crlbLblOdo(T,X,Xb,dtLbl,sigR,sigV,vel,sigHdg)

        dSrss0 = abs( (srss1-srss)/srss1 )
        dSrss1 = abs( (srss0-srss1)/srss0)

        srss0 = srss1
        srss1 = srss
        N+=1
    if N > 99:
        print("WARNING: crlbutils - covergence req'd %d interations"%N)
    return (crlb,cep,srss,N)

def crlbLblOdo(T,X,Xb,dtLbl,sigR,sigV,vel,sigHdg):

    '''
    LBL and Odometry
    '''
    Nlbl = int(floor(T/dtLbl))
    Ns = 2*Nlbl      # number of states
    # Lbl Range measurements
    Cr = zeros((len(Xb)*Nlbl,Ns))
    for kk in range(Nlbl):
        for xb,ii in zip(Xb,range(len(Xb))):
            r = sqrt((X[0]-xb[0])**2+(X[1]-xb[1])**2)
            r = max((r,0.1))
            for jj in range(2):
                d = X[jj]-xb[jj]
                Cr[ii+(kk*len(Xb)),jj+(kk*2)] = d/r

    # DVL and Hdg
    No = 2*(Nlbl-1)  # number of odo. observations
    #Co = zeros((NUsbl-1,2*NAsbl))
    C1 = hstack((zeros((No,2)),eye(No)))
    C2 = hstack((-1*eye(No),zeros((No,2))))
    Co = C1+C2
    ddist = dtLbl*vel
    varoy = (sigV**2)*dtLbl
    varox = (sigV**2)*dtLbl+(sigHdg**2)*ddist
    ro = zeros(No)
    for ii in range(Nlbl-1):
        ro[ii*2]=varox
        ro[ii*2+1]=varoy

    # Combine
    C = vstack((Cr,Co))
    R = diag(hstack((ones(Nlbl*len(Xb))*(sigR**2),
                     ro)))

    lam = dot(dot(C.T,inv(R)),C)
    try:
        crlb = inv(lam)
    except LinAlgError:
        crlb = eye(R.shape[0])*100*sigR**2

    # Extract middle of the cov. matrix
    NN = crlb.shape[0]/2
    if Nlbl%2:  # odd
        ll = Nlbl-1
    else:
        ll = Nlbl-2

    # metrics
    crlbxy = crlb[ll:ll+2,ll:ll+2]
    cep = sum(sqrt(eigvals(crlbxy)))*0.59
    srss = sqrt(sum(diag(crlbxy)))

    return (crlbxy,cep,srss)

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
    majMag = 2*sqrt(abs(chi/d[1]))
    minMag = 2*sqrt(abs(chi/d[0]))
    th = math.atan2(V[1,1],V[0,1])
    e1 = patches.Ellipse(ecenter,width=majMag,height=minMag,angle=th*180.0/pi,linewidth=2,fill=False)
    ax.add_patch(e1)
    #ax.set_aspect('equal')
    ax.autoscale()



def annoteParams(ax,x,y,dtUsbl,varr,vara,varv,vel,varhdg,loc=(0.1,0.5),usbl=True,odo=True):
    if usbl and odo:
        textstr = ("x=%.1f m, y=%.1f m \n"
                   "USBL: $\delta t= %.1f s, \, \sigma_r= %.2f m, \, "
                   "\sigma_{\\theta}=%.2f deg$ \n"
                   "DVL: $ \sigma_v=%.1f mm/s, \, \sigma_{hdg}=%.2f deg $\n"
                   "vel=%.1f m/s"%
                   (x,y,dtUsbl,sqrt(varr),sqrt(vara)*180/pi,
                    sqrt(varv)*1000,sqrt(varhdg)*180/pi,vel))
    elif usbl and (not odo):
        textstr = ("x=%.1f m, y=%.1f m \n"
                   "USBL: $\delta t= %.1f s, \, \sigma_r= %.2f m, \, "
                   "\sigma_{\\theta}=%.2f deg$"%
                   (x,y,dtUsbl,sqrt(varr),sqrt(vara)*180/pi))
    props = dict(boxstyle='round', facecolor='wheat', alpha=1)
    ax.text(loc[0],loc[1],textstr,transform=ax.transAxes,verticalalignment='top',bbox=props)


def annoteLblParams(ax,Xb,dtLbl,sigR,sigV,vel,sigHdg,
                    loc=(0.1,0.5)):
    textstr = ("LBL: $\delta t= %.1f s, \, \sigma_r= %.2f m$ \n"
               "DVL: $ \sigma_v=%.1f mm/s, \, \sigma_{hdg}=%.2f deg $\n"
               "vel=%.1f m/s"%
               (dtLbl,sigR,
                sigV*1000,sigHdg*180/pi,vel))
    props = dict(boxstyle='round', facecolor='wheat', alpha=1)
    ax.text(loc[0],loc[1],textstr,transform=ax.transAxes,verticalalignment='top',bbox=props)
