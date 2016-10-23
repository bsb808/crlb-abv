import crlbutils as cu
reload(cu)

dtLbl = 1
T = 2*dtLbl
xb1 = (10,10)
xb2 = (190,10)
bl = xb2[0]-xb1[0]
xb3 = (xb1[0]+bl/2.0,xb1[1]+bl/2.0*tan(60*pi/180))
Xb = [xb1,xb2,xb3]
X = (50,50)
sigR = 1
sigV = 0.003
sigHdg = 1

dtLbl = 1
NN = linspace(1,100,100)
Cep = []
Srss = []
for N in NN:
    crlb,cep,srss = cu.crlbLblOdo(N*dtLbl,X,Xb,dtLbl,sigR,sigV,vel,sigHdg)
    Cep.append(cep)
    Srss.append(srss)

crlb,cep,srss,N = cu.crlbLblOdoConverge(X,Xb,dtLbl,sigR,sigV,vel,sigHdg,tol=0.01)

figure(2)
clf()
plot(NN,Cep,'o',label='cep')
plot(NN,Srss,'o',label='srss')
legend()
show()
