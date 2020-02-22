
import crlbutils as cu
reload(cu)

#close('all')

figure(1)
clf()

VV = [5, 50, 500] #, 200, 800]
Labels = ['Slow','Med','Fast']
for V,L in zip(VV,Labels):
    varv = (V*0.003)**2   # velocity variance

    # Parameters
    # GPS
    dtGps = 5*60.0
    varxy = (1852.0)**2 #0.3)**2 #(0.14)**2  # USBL range variance

    vel = 0.5  # vehicle velocity [m/s]
    varhdg = (20.0*pi/180.0)**2  # heading variance

    Cepo = []
    Srsso = []
    Sigxo = []
    Sigyo = []
    CIo = []
    Crlbo = []
    No = []
    
    NN = arange(2,100,1)
    dtMins = arange(5,50,5)
    dtMins = linspace(.01,1000,50)
    dtMins = logspace(0,4,50)
    for dtMin in dtMins:
        N = 50
        # With odo
        dtGps = dtMin*60.0
        T = dtGps*N
        crlb, cep, srss = cu.crlbGpsOdo(T,dtGps,varxy,varv,vel,varhdg)
        ###Yy        crlb, cep, srss, N = cu.crlbGpsOdoConverge(x,y,xb,yb,dtUsbl,varr,vara,varv,vel,varhdg)
        Cepo.append(cep)
        Srsso.append(srss)
        Sigxo.append(sqrt(crlb[0,0]))
        Sigyo.append(sqrt(crlb[1,1]))

        CIo.append(sqrt(5.99*max(crlb[0,0],crlb[1,1])))
        #CIo.append(sqrt(5.99)*(sqrt(Sigxo[-1]*Sigyo[-1])))
        Crlbo.append(crlb)
        No.append(N)


    #semilogx(dtMins,Cepo,label=L)
    semilogx(dtMins/100,Cepo/max(Cepo),label=L)
    grid(True)

ylabel(r'$\overline{\mathrm{CEP}}$ [dimensionless]')
xlabel('Normalized Update Rate  [dimensionless]')
legend(loc=2)
show()
