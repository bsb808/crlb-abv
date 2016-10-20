    

# Build C matrix for multiple measurements
x = 10
y = 0
xb = 0
yb = 0

TT = linspace(1,50,50)
CRx = []
CRy = []
tsig = []
Cep = []
#TT = [2.0]
for II in range(len(TT)):
# How many USBL measurmeents
    T = TT[II] #1.0  #s
    dtUsbl = 1.0
    NUsbl = int(floor(T/dtUsbl))

    # How many DVL measurements
    dtDvl = 0.1
    NDvl = int(floor(T/dtDvl))

    # Range measurements
    Cr = zeros((NUsbl,2*NUsbl))
    r = sqrt((x-xb)**2+(y-yb)**2)
    dy = y-yb
    dx = x-xb
    for ii in range(NUsbl):
        Cr[ii,ii*2] = dx/r
        Cr[ii,ii*2+1] = dy/r
    varr = (1.0)**2

    # Angle
    Ca = zeros((NUsbl,2*NUsbl))
    for ii in range(NUsbl):
        Ca[ii,ii*2]=-dy/(r**2)
        Ca[ii,ii*2+1]=dx/(r**2)
    vara = (1*pi/180.0)**2

    # DVL and Hdg
    No = 2*(NUsbl-1)
    Ns = 2*NUsbl
    #Co = zeros((NUsbl-1,2*NAsbl))
    C1 = hstack((zeros((No,2)),eye(No)))
    C2 = hstack((-1*eye(No),zeros((No,2))))
    Co = C1+C2
    varv = (0.01)**2   # velocity variance
    vel = 1  # velocity
    ddist = dtUsbl*vel
    varox = varv*dtUsbl
    varhdg = (10*pi/180.0)**2  # heading variance
    varoy = varv*dtUsbl+varhdg*ddist
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
    NN = crlb.shape[0]/2
    CRx.append(sqrt(crlb[0,0]))
    CRy.append(sqrt(crlb[1,1]))
    # Extract middle of the cov. matrix
    if NUsbl%2:  # odd
        ll = NUsbl-1
    else:
        ll = NUsbl-2
    crlbxy = crlb[ll:ll+2,ll:ll+2]
    Cep.append(sum(sqrt(eigvals(crlbxy)))*0.59)
    tsig.append(sqrt(sum(diag(crlbxy))))
    #print("NUsbl %d, ll %d"%(NUsbl,ll))
    #print crlbxy
    #for jj in range(crlb.shape[0]/2):
    #    crlbxy = crlb[jj*2:jj*2+2,jj*2:jj*2+2]
        #print crlbxy

figure(1)
clf()
plot(TT,CRx,label='CRx')
plot(TT,CRy,label='CRy')
plot(TT,tsig,label='sqrtss')
plot(TT,Cep,label='CEP')
legend()
grid(True)
show()
