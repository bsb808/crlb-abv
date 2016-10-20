# Meant to be run thorugh ipython --pylab


def Cusbl(x,y,xb,yb):
    ''' Evaluate C matrix for position x and beacon at xb,yb '''
    C = zeros((2,2))
    # Range
    r = sqrt((x-xb)**2+(y-yb)**2)
    dy = y-yb
    dx = x-xb
    C[0,0] = dx/r
    C[0,1] = dy/r
    # Angle
    C[1,0] = -dy/(r**2)
    C[1,1] = dx/(r**2)
    return C

X,Y = (100.0,100.0)
Nx, Ny = (100,100)
xx = linspace(0,X,Nx)
yy = linspace(0,Y,Ny)

xv,yv = meshgrid(xx,yy)

zv = zeros(xv.shape)
zvx = zeros(xv.shape)
zvy = zeros(xv.shape)

R = zeros((2,2))
R[0,0] = 0.0001  # range variance
R[1,1] = (pi/180.0)**2  # angle variance

xb = X/2.0
yb = -1.0

xb = -1
yb = -1


for x,ii in zip(xx,range(len(xx))):
    for y,jj in zip(yy,range(len(yy))):
        C = Cusbl(x,y,xb,yb)
        crlb = inv(dot(dot(C.T,inv(R)),C))
        zv[jj,ii] = sqrt(crlb[0,0]+crlb[1,1])
        zvx[jj,ii] = sqrt(crlb[0,0])
        zvy[jj,ii] = sqrt(crlb[1,1])
        #zv[jj,ii] = Cusbl(x,y,xb,yb)[1,0]

for z,ii in zip([zv, zvx, zvy],range(3)):
    figure(ii)
    clf()
    CP = contour(xv,yv,z)
    clabel(CP,inline=1,fontsize=10)
    
show()


