


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


# Transponder locations
# Equilateral triangle
xb1 = (10,10)
xb2 = (190,10)
bl = xb2[0]-xb1[0]
xb3 = (xb1[0]+bl/2.0,xb1[1]+bl/2.0*tan(60*pi/180))
Xb = [xb1,xb2,xb3]
# oblique
xb1a = (100, xb1[1]+0.5*bl/2.0*tan(60*pi/180))
#Xb =[xb1a,xb2,xb3]

# 2 line
#Xb = [xb1,(190,190)]

sigR = 0.23  #0.1  # range measurment stdev


# Mesh
X,Y = (200.0,200.0)  # extend of area
Nx, Ny = (100,100)
xx = linspace(0,X,Nx)
yy = linspace(0,Y,Ny)

xv,yv = meshgrid(xx,yy)

# Iterate through mesh
zCep = zeros(xv.shape)
zSrss = zeros(xv.shape)

for x,ii in zip(xx,range(len(xx))):
    for y,jj in zip(yy,range(len(yy))):
        crlb,cep,srss = crlbLbl(Xb,(x,y),sigR)
        zCep[jj,ii] = cep
        zSrss[jj,ii] = srss

for z,ii,tit in zip([zCep, zSrss],range(2),
                    ['Circular Error Probable (CEP)',
                     'Sqrt Sum of Squares']):
    figure(ii)
    clf()
    #CP = contour(xv,yv,z,levels=array([1., 2., 5., 10. , 30.])*sigR)
    CP = contour(xv,yv,z)
    clabel(CP,inline=1,fontsize=10)
    for xb in Xb:
        plot(xb[0],xb[1],'ro',markersize=12)
    title("%s: $\sigma_r = %.3f \, m$"%(tit,sigR))
    xlabel('X [m]')
    ylabel('Y [m]')

show()





 
