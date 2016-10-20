
import crlbutils as cu
reload(cu)


# Position
XX = linspace(1,100,10)
y = 0
xb = 0
yb = 0

# Parameters
# USBL
dtUsbl = 1.0
varr = (1.0)**2  # USBL range variance
vara = (1*pi/180.0)**2  # USBL angle variance
varv = 1000**2  
varv = (0.01)**2   # velocity variance
vel = 1  # vehicle velocity
varhdg = 10**2  
varhdg = (10*pi/180.0)**2  # heading variance
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
subplot(211)
plot(XX,Cep,label='CEP')
plot(XX,Sigx,label='SigX')
plot(XX,Sigy,label='SigY')
title('USBL Standalone: sigR=%.2f m, sigT=%.2f deg'%(sqrt(varr),sqrt(vara)*180.0/pi))
legend(loc='lower right')
grid(True)
ylabel('Uncertainty [m]')
ax = subplot(212)
for crlb,x in zip(Crlb,XX):
    cu.plotErrEllipse(ax,crlb,ecenter=(x,y))
xlabel('Distance [m]')
ylabel('Uncertainty, 90% [m]')
grid(True)
show()

figure(2)
clf()
subplot(211)
plot(XX,Cepo,label='CEP')
plot(XX,Sigxo,label='SigX')
plot(XX,Sigyo,label='SigY')
#plot(XX,No,label='N')
title('USBL with Odo.: sigR=%.2f m, sigT=%.2f deg'%(sqrt(varr),sqrt(vara)*180.0/pi))
legend(loc='lower right')
grid(True)
ylabel('Uncertainty [m]')
ax = subplot(212)
for crlb,x in zip(Crlbo,XX):
    cu.plotErrEllipse(ax,crlb,ecenter=(x,y))
xlabel('Distance [m]')
ylabel('Uncertainty, 90% [m]')
grid(True)
show()


