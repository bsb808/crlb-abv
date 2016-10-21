
import crlbutils as cu
reload(cu)
import copy
    
# Position
x = 200
y = 0
xb = 0
yb = 0

# Parameters
# USBL
sens = {}
sens['dtUsbl'] = 1.0
sens['varr'] = (0.1)**2  # USBL range variance
sens['vara'] = (2.24*pi/180.0)**2  # USBL angle variance
sens['varv'] = (0.003)**2   # velocity variance
sens['vel'] = 1  # vehicle velocity
sens['varhdg'] = (2.0*pi/180.0)**2  # heading variance

osens = copy.deepcopy(sens)
close('all')
for kk in sens.keys():
    sens = copy.deepcopy(osens)
    Npts = 20
    if 'var' in kk:
        pvec = logspace(-4,4,Npts)*osens[kk]
    elif 'dtUsbl' in kk:
        pvec = logspace(-1,10,Npts)*osens[kk]
    else:
        pvec = logspace(-2,2,Npts)*osens[kk]
    Cep = []
    Srss = []
    Sigx = []
    Sigy = []
    Crlb = []
    NN = []
    for p in pvec:
        sens[kk]=p
        crlb, cep, srss, N = cu.crlbUsblOdoConverge(x,y,xb,yb,
                                                    sens['dtUsbl'],
                                                    sens['varr'],
                                                    sens['vara'],
                                                    sens['varv'],
                                                    sens['vel'],
                                                    sens['varhdg'])
        Cep.append(cep)
        Srss.append(srss)
        Sigx.append(sqrt(crlb[0,0]))
        Sigy.append(sqrt(crlb[1,1]))
        Crlb.append(crlb)
        NN.append(N)
    figure()
    clf()
    subplot(111)
    if kk=='dtUsbl':
        xlabel('USBL: $\delta t$ [s]')
    elif kk=='varr':
        xlabel('USBL: $\sigma_r$ [m]')
        pvec = sqrt(pvec)
    elif kk=='vara':
        xlabel('USBL: $\sigma_{\\theta}$ [deg]')
        pvec = sqrt(pvec)*180/pi
    elif kk=='varv':
        xlabel('DVL: $\sigma_v$ [mm/s]')
        pvec = sqrt(pvec)*1000
    elif kk=='varhdg':
        xlabel('DVL: $\sigma_{hdg}$ [deg]')
        pvec = sqrt(pvec)*180/pi
    elif kk=='vel':
        xlabel('Mobile Velocity [m/s]')
        pvec = pvec
    else:
        xlabel(kk)
        pvec=pvec
    pfun = semilogx  # semilog or just plot?
    pfun(pvec,Srss,label='SRSS')
    pfun(pvec,Sigx,label='SigX')
    pfun(pvec,Sigy,label='SigY')
    title('Sensitivity: %s'%kk)
    legend(loc='center right')
    ax = gca()
    cu.annoteParams(ax,x,y,
                 osens['dtUsbl'],
                 osens['varr'],
                 osens['vara'],
                 osens['varv'],
                 osens['vel'],
                 osens['varhdg'])

show()

