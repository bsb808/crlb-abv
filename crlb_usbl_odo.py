
import crlbutils as cu
reload(cu)


# Position
x = 10
y = 10
xb = 0
yb = 0

# Parameters
# USBL
dtUsbl = 1.0
varr = (1.0)**2  # USBL range variance
vara = (1*pi/180.0)**2  # USBL angle variance
varv = 1000**2  
#varv = (0.01)**2   # velocity variance
vel = 1  # vehicle velocity
varhdg = 10**2  
#varhdg = (10*pi/180.0)**2  # heading variance
# 

CRx = []
CRy = []

Cep = []

figure(3)
clf()
figure(4)
clf()
figure(5)
clf()

vals = linspace(1,100,10)
vals = [1.0]
SRSS = []
NN = []
for val in vals:
    dtUsbl = val
    Srss = []
    TT = []
    TT.append(dtUsbl)
    crlb,cep,srss = cu.crlbUsblOdo(TT[-1],x,y,xb,yb,
                                dtUsbl,varr,vara,varv,vel,varhdg)
    Srss.append(srss)
    dSrss = 1
    # determine convergence within some tolerance
    while dSrss > 0.01:
        TT.append(TT[-1]+dtUsbl)
        crlb,cep,srss = cu.crlbUsblOdo(TT[-1],x,y,xb,yb,
                                    dtUsbl,varr,vara,varv,vel,varhdg)
        Srss.append(srss)
        dSrss = abs( (Srss[-2]-Srss[-1])/Srss[-2] )
        N = int(TT[-1]/dtUsbl)
    
    SRSS.append(srss)
    NN.append(N)

    figure(3)
    plot(TT,Srss,label='SRSS %f'%val)
    figure(4)
    plot(Srss,label='SRSS %f'%val)
    figure(5)
    plot(diff(Srss)/Srss[:-1],label='SRSS %f'%val)
legend()
grid(True)

figure(6)
clf()
subplot(211)
plot(vals,SRSS)
ylabel('SRSS')
subplot(212)
plot(vals,NN)
ylabel('N')
xlabel('vals')



show()
