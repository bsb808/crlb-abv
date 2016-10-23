
import crlbutils as cu
reload(cu)


# Transponder locations
# Equilateral triangle
xb1 = (10,10)
xb2 = (190,10)
bl = xb2[0]-xb1[0]
xb3 = (xb1[0]+bl/2.0,xb1[1]+bl/2.0*tan(60*pi/180))
Xb = [xb1,xb2,xb3]
constLabel='equil'
# oblique
xb1a = (100, xb1[1]+0.5*bl/2.0*tan(60*pi/180))
Xb =[xb1a,xb2,xb3]
constLabel='oblique'

# 2 line
Xb = [xb1,(190,190)]
constLabel='line'
sigR = 0.23  #0.1  # range measurment stdev
dtLbl = 1.0
sigV = 0.003
vel = 1
sigHdg = 2.0*pi/180.0

# Mesh
X,Y = (200.0,200.0)  # extend of area
Nx, Ny = (100,100)
xx = linspace(0,X,Nx)
yy = linspace(0,Y,Ny)

xv,yv = meshgrid(xx,yy)

# Iterate through mesh
zCep = zeros(xv.shape)
zSrss = zeros(xv.shape)
zN = zeros(xv.shape)


for x,ii in zip(xx,range(len(xx))):
    for y,jj in zip(yy,range(len(yy))):
        #crlb,cep,srss = cu.crlbLbl(Xb,(x,y),sigR)
        #N = 1
        #crlb,cep,srss = cu.crlbLblOdo(10,(x,y),Xb,1,sigR,1,1,1)
        crlb,cep,srss,N = cu.crlbLblOdoConverge((x,y),Xb,dtLbl,sigR,sigV,vel,sigHdg)
        zCep[jj,ii] = cep
        zSrss[jj,ii] = srss
        zN[jj,ii]=N

for z,ii,tit in zip([zCep, zSrss, zN],range(3),
                    ['Circular Error Probable (CEP)',
                     'Sqrt Sum of Squares',
                     'N interations']):
    figure(ii)
    clf()
    #CP = contour(xv,yv,z,levels=array([1., 2., 5., 10. , 30.])*sigR)
    CP = contour(xv,yv,z)
    clabel(CP,inline=1,fontsize=10)
    for xb in Xb:
        plot(xb[0],xb[1],'ro',markersize=12)
    title("%s: $\sigma_r = %.3f \, m$"%(tit,sigR))
    ax = gca()
    cu.annoteLblParams(ax,Xb,dtLbl,sigR,sigV,vel,sigHdg,loc=(0.01,0.95))
    xlabel('X [m]')
    ylabel('Y [m]')

    savefig('lbl-odo-%s-%d.png'%(constLabel,ii),dpi=300)
    #savefig('lbl-standalone-%s-%d.png'%(constLabel,ii),dpi=300)

show()





 
