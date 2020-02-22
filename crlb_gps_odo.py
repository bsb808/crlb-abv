
import crlbutils as cu
reload(cu)

#close('all')

figure(1)
clf()

VV = [50, 200, 800]
Labels = ['High','Med','Low']
for V,L in zip(VV,Labels):
    varv = (V*0.003)**2   # velocity variance

    # Parameters
    # GPS
    dtGps = 5*60.0
    varxy = (1852.0)**2 #0.3)**2 #(0.14)**2  # USBL range variance

    vel = 5.0  # vehicle velocity [m/s]
    varhdg = (20.0*pi/180.0)**2  # heading variance

    Cepo = []
    Srsso = []
    Sigxo = []
    Sigyo = []
    CIo = []
    Crlbo = []
    No = []

    NN = arange(2,100,1)
    for N in NN:

        # With odo
        T = dtGps*N
        crlb, cep, srss = cu.crlbGpsOdo(T,dtGps,varxy,varv,vel,varhdg)
        N = 1
        ###Yy        crlb, cep, srss, N = cu.crlbGpsOdoConverge(x,y,xb,yb,dtUsbl,varr,vara,varv,vel,varhdg)
        Cepo.append(cep)
        Srsso.append(srss)
        Sigxo.append(sqrt(crlb[0,0]))
        Sigyo.append(sqrt(crlb[1,1]))

        CIo.append(sqrt(5.99*max(crlb[0,0],crlb[1,1])))
        #CIo.append(sqrt(5.99)*(sqrt(Sigxo[-1]*Sigyo[-1])))
        Crlbo.append(crlb)
        No.append(N)


    plot(NN,Cepo/max(Cepo),label=L)
    grid(True)

ylabel(r'$\overline{\mathrm{CEP}}$ [dimenstionless]')
xlabel('Position Updates [N]')
legend()
show()
