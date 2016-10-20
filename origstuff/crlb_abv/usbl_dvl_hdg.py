# USBL + DVL + Heading

# Assume path is straight line along the X-axis

# Setup
D = 1000;  # Track length [m]
# Assume 10 second update rate - that means every 2 meters
Nstations = 500;
#Nstations = 3
Nusbl = Nstations-1;  # Number of USBL fixes
vel = 0.2;  # vehicle speed in m/s
# USBL Uncertaint
sig_r = 0.2;  # 1 sigma in range [m]
sig_b = 0.25*pi/180;  # 1 sigma in bearing [rad]
# DVL Uncertainty
sig_v = 0.003;  # velocity [m/s]
# Heading 
sig_h = 0.1*pi/180;  # heading [rad]

# Interstitial
dX = D/Nusbl;
dT = dX/vel;
# X positions

# Build C matrix one row at a time
XX = linspace(dX,D,Nusbl);
for ii in range(Nusbl):
    # Range obs.
    cr = zeros(2*Nusbl)
    cr[ii*2]=1
    if ii == 0:
        C = cr
    else:
        C = vstack((C,cr))
    # Bearing obs.
    br = zeros(2*Nusbl)
    br[2*ii-1]=1/XX[ii];
    C = vstack((C,br))


# Now DVL and Hdg
cx = zeros(2*Nusbl)
cx[0] = 1;
C = vstack((C,cx))
cx = zeros(2*Nusbl)
cx[1] = 1;
C = vstack((C,cx))

for ii in range(1,Nusbl):
    cx = zeros(2*Nusbl)
    cx[(ii-1)*2]=-1;
    cx[(ii-1)*2+2]=1
    C = vstack((C,cx));
    cy = zeros(2*Nusbl)
    cy[(ii-1)*2+1]=-1;
    cy[(ii-1)*2+3]=1;
    C = vstack((C,cy));


# Build R matrix 

Nobs = shape(C)[0]
R = eye(Nobs)
for ii in range(Nusbl):
    iii = (ii)*2;
    R[iii,iii] = sig_r**2;
    jjj = 2*ii+1;
    R[jjj,jjj] = sig_b**2;

for ii in range(1,Nusbl):
    iii = Nusbl*2+(ii-1)*2
    R[iii,iii] = dT*(sig_v**2);
    jjj =  Nusbl*2+(ii-1)*2+1
    R[jjj,jjj] = dT*(sig_v**2)+(dX*sig_h)**2;


F = dot(dot(C.transpose(),(inv(R))),C)
lam = inv(F)
CC = sqrt(diag(lam));
alongT = CC[::2];
acrossT = CC[1::2];
# CEP
cep = zeros(len(alongT));
for ii in range(len(alongT)):
    sigx = alongT[ii]
    sigy = acrossT[ii]
    sigl = max((sigx,sigy));
    sigs = min((sigx,sigy));
    w = sigs/sigl;
    if w < 0.5:
        cep[ii]= sigl*(0.67+0.8*w**2);
    else:
        cep[ii]=0.59*sigl*(1+w);

    
        

# Plot
figure(3)
XXX = hstack((0,XX))
clf()
plot(XXX,hstack((0,alongT)),'b-',linewidth=2,label='Along-Track');
hold(True)
plot(XXX,hstack((0,acrossT)),'r-',linewidth=2,label='Across-Track');
plot(XXX,hstack((0,cep)),'k-',linewidth=2,label='CEP');
legend()
grid('on')
ylabel('Position Uncertainty [m]')
xlabel('Distance Traveled [m]')
show()
