

# SBL + DVL/Model + Heading

D = 1100.0;  # Track length [m]
vel = 0.2;  # vehicle speed in m/s
# Assume 10 second update rate - that means every 2 meters
Updaterate = 10 
Nstations = round(D/vel/Updaterate);
#Nstations = 10
Nusbl = Nstations-1;  # Number of USBL fixes


# SBL Uncertaint
sig_r = 0.1;  # 1 sigma in range [m]
BL = 6.0 # baseline [m]


# Interstitial
dX = D/Nusbl;
dT = dX/vel;

# DVL Uncertainty
sig_v = 0.003;  # velocity [m/s]
sig_x = sqrt((sig_v**2)*dT)
# model uncertainty as percent of distance traveled
percentE = 0.1
sig_x = dX*percentE
sig_y = sig_x

# Heading 
sig_h = 0.1*pi/180;  # heading [rad]


# Positions along X
XX = linspace(dX,D,Nusbl);

def CmatrixSbl(DistList):
    '''
    Build a C matrix for standalone SBL
    DistList is the list of distances along the line
    Assumes motion along the X axis
    '''
    Nobs = len(DistList)
    return eye(2*Nobs)

def RmatrixSbl(DistList,sig_r,BL):
    '''
    Build a C matrix for standalone SBL
    DistList is the list of distances along the line
    Assumes motion along the X axis
    '''
    Nobs = len(DistList)
    R = zeros((2*Nobs,2*Nobs))
    for d,ii in zip(DistList,range(len(DistList))):
        R[ii*2,ii*2]=(sig_r**2)/2.0  # unc. in X direction
        r = sqrt(d**2+(BL/2.0)**2)
        R[ii*2+1,ii*2+1]=(sig_r**2)/(2.0*(BL/(2.0*r))**2)  # unc. in Y direction
    return R
        

def CmatrixDvlCompass(DistList):
    Nobs = len(DistList)
    C = eye(2*Nobs)
    for ii in range(2,2*Nobs):
        C[ii,ii-2]=-1
    return C
    
def RmatrixDvlCompass(DistList,sigx,sigy):
    Nobs = len(DistList)
    R = sigx*eye(2*Nobs)
    for ii in range(1,Nobs,2):
        R[ii,ii]=sigy
    return R


Csbl = CmatrixSbl(XX)
Rsbl = RmatrixSbl(XX,sig_r,BL)

def CRLB(C,R):
    F = dot(dot(C.transpose(),(inv(R))),C)
    lam = inv(F)
    CC = sqrt(diag(lam))
    alongT = CC[::2];
    acrossT = CC[1::2];
    return (alongT,acrossT)

alongTsbl,acrossTsbl = CRLB(Csbl,Rsbl)

# CEP
def CEP(sigx,sigy):
    sigl = max((sigx,sigy));
    sigs = min((sigx,sigy));
    w = sigs/sigl;
    if w < 0.5:
        return sigl*(0.67+0.8*w**2);
    else:
        return 0.59*sigl*(1+w);

def xy2cep(alongT,acrossT):
    cep = zeros(len(alongT));
    for ii in range(len(alongT)):
        sigx = alongT[ii]
        sigy = acrossT[ii]
        cep[ii]=CEP(sigx,sigy)
    return cep

cepsbl = xy2cep(alongTsbl,acrossTsbl)
    
    
Cdvl = CmatrixDvlCompass(XX)
C = vstack((Csbl,Cdvl))
Rdvl = RmatrixDvlCompass(XX,sig_x,sig_y)
R = hstack((diag(Rsbl),diag(Rdvl)))
R = diag(R)

alongT,acrossT = CRLB(C,R)
cep = xy2cep(alongT,acrossT)


# Plot
figure(4)
XXX = hstack((0,XX))
clf()
plot(XXX,hstack((0,alongTsbl)),'b--',linewidth=2,label='Along: SBL');
plot(XXX,hstack((0,alongT)),'b-',linewidth=2,label='Along: SBL+Model');
hold(True)
plot(XXX,hstack((0,acrossTsbl)),'r--',linewidth=2,label='Across: SBL');
plot(XXX,hstack((0,acrossT)),'r-',linewidth=2,label='Across: SBL+Model');
plot(XXX,hstack((0,cepsbl)),'k--',linewidth=2,label='CEP: SBL');
plot(XXX,hstack((0,cep)),'k-',linewidth=2,label='CEP: SBL+Model');
grid('on')
ylabel('Position Uncertainty [m]')
xlabel('Distance Traveled [m]')
xlim([0,1000])
legend(loc=2)
title('SBL BL=%2.0f m, Unc. %2.0f cm  Update %3.1f s; Model Unc. %2.0f%%;'%(BL,sig_r*100, D/vel/Nstations,percentE*100))
xlim([0,1000])
show()
show()

figure(5)
clf()
plot(XXX,hstack((0,cepsbl)),'k--',linewidth=2,label='SBL');
plot(XXX,hstack((0,cep)),'k-',linewidth=2,label='SBL+Model');
legend(loc=2)
grid('on')
ylabel('Position Uncertainty (CEP) [m]')
xlabel('Distance Traveled [m]')
title('SBL BL=%2.0f m, Unc. %2.0f cm  Update %3.1f s; Model Unc. %2.0f%%;'%(BL,sig_r*100, D/vel/Nstations,percentE*100))
xlim([0,1000])
show()

def printall(fname):
    figure(4)
    savefig(fname+'.png',dpi=300)
    figure(5)
    savefig(fname+'2.png',dpi=300)
    
