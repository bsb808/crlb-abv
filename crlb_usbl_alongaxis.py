
import crlbutils as cu
reload(cu)


# Position
XX = linspace(1,10,10)
y = 0
xb = 0
yb = 0

# Parameters
# USBL
dtUsbl = 5.0
varr = (0.3)**2 #(0.14)**2  # USBL range variance
vara = (2.0*pi/180.0)**2  # USBL angle variance
varv = (0.003)**2   # velocity variance
vel = 1  # vehicle velocity
varhdg = (2.0*pi/180.0)**2  # heading variance
# 
Cep = []
Srss = []
Sigx = []
Sigy = []
Crlb = []

Cepo = []
Srsso = []
Sigxo = []
Sigyo = []
Crlbo = []
No = []
for x in XX:
    # Stand Alone
    crlb, cep, srss = cu.crlbUsblAlone(x,y,xb,yb,dtUsbl,varr,vara,varv,vel,varhdg)
    Cep.append(cep)
    Srss.append(srss)
    Sigx.append(sqrt(crlb[0,0]))
    Sigy.append(sqrt(crlb[1,1]))
    Crlb.append(crlb)

    # With odo
    crlb, cep, srss, N = cu.crlbUsblOdoConverge(x,y,xb,yb,dtUsbl,varr,vara,varv,vel,varhdg)
    Cepo.append(cep)
    Srsso.append(srss)
    Sigxo.append(sqrt(crlb[0,0]))
    Sigyo.append(sqrt(crlb[1,1]))
    Crlbo.append(crlb)
    No.append(N)

figure(1)
clf()
ax=subplot(211)
plot(XX,Srss,label='$\sqrt{\sigma_x^2+\sigma_y^2}$')
plot(XX,Sigx,label='$\sigma_x$')
plot(XX,Sigy,label='$\sigma_y$')
#title('USBL Standalone: $\sigma_{range}$=%.2f m, $\sigma_{angle}$=%.2f deg'%(sqrt(varr),sqrt(vara)*180.0/pi))
title('USBL Standalone')
legend(loc='lower right')
grid(True)
ylabel('Uncertainty [m]')
cu.annoteParams(ax,x,y,
                 dtUsbl,
                 varr,
                 vara,
                 varv,
                 vel,
                 varhdg,
                loc=(0.01,0.95),odo=False)
ax = subplot(212,sharex=ax)
for crlb,x in zip(Crlb,XX):
    cu.plotErrEllipse(ax,crlb,ecenter=(x,y))
xlabel('Distance [m]')
ylabel('Uncertainty, 90% [m]')
grid(True)
show()

figure(2)
clf()
ax=subplot(211)
plot(XX,Srsso,label='$\sqrt{\sigma_x^2+\sigma_y^2}$')
plot(XX,Sigxo,label='$\sigma_x$')
plot(XX,Sigyo,label='$\sigma_y$')
#plot(XX,No,label='N')
#title('USBL with DVL: sigR=%.2f m, sigT=%.2f deg'%(sqrt(varr),sqrt(vara)*180.0/pi))
title('USBL with DVL')
legend(loc='upper right')
grid(True)
ylabel('Uncertainty [m]')
cu.annoteParams(ax,x,y,
                 dtUsbl,
                 varr,
                 vara,
                 varv,
                 vel,
                 varhdg,
                loc=(0.01,0.95))
ax = subplot(212,sharex=ax)
for crlb,x in zip(Crlbo,XX):
    cu.plotErrEllipse(ax,crlb,ecenter=(x,y))
xlabel('Distance [m]')
ylabel('Uncertainty, 90% [m]')
grid(True)
show()


