
import crlbutils as cu
reload(cu)

close('all')

# Position
L=300
L=10.0
LL = [300.0,10.0]
for L,ll in zip(LL,range(len(LL))):

    XX = linspace(1,L,10)
    y = 0
    xb = 0
    yb = 0

    # Parameters
    # USBL
    dtUsbl = 5.0
    varr = (0.3)**2 #(0.14)**2  # USBL range variance
    vara = (0.3*pi/180.0)**2  # USBL angle variance
    varv = (0.003)**2   # velocity variance
    vel = 0.5  # vehicle velocity [m/s]
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
    ax=subplot(2,1,ll+1)
    plot(XX,Srss,label='$\sqrt{\sigma_x^2+\sigma_y^2}$')
    plot(XX,Sigx,label='$\sigma_x$')
    plot(XX,Sigy,label='$\sigma_y$')
    if ll < 1:
        title('USBL Standalone')
        legend(loc='lower right')

        ylabel('Uncertainty [m]')
        cu.annoteParams(ax,x,y,
                        dtUsbl,
                        varr,
                        vara,
                        varv,
                        vel,
                        varhdg,
                        loc=(0.01,0.95),odo=False)
    else:
        xlabel('Along-Track Distance [m]')
    grid(True)



    for jj in range(2):
        figure()
        clf()
        if jj > 0:
            ax = subplot(111)
        else:
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
        if jj > 0:
            pass
        else:
            ax = subplot(212,sharex=ax)
            for crlb,x in zip(Crlb,XX):
                cu.plotErrEllipse(ax,crlb,ecenter=(x,y),chi=1.0)

                ylabel('Across-Track Distance [m]')
                grid(True)
        xlabel('Along-Track Distance [m]')

        figure()
        clf()
        if jj > 0:
            ax=subplot(111)
        else:
            ax=subplot(211)
        plot(XX,Srsso,label='$\sqrt{\sigma_x^2+\sigma_y^2}$')
        plot(XX,Sigxo,label='$\sigma_x$')
        plot(XX,Sigyo,label='$\sigma_y$')
        #plot(XX,No,label='N')
        #title('USBL with DVL: sigR=%.2f m, sigT=%.2f deg'%(sqrt(varr),sqrt(vara)*180.0/pi))
        title('USBL with DVL')
        if jj > 0:
            #legend(loc='lower right')
            legend(loc='upper right')
            ylim([0,max(Srsso)*1.5])
        else:
            legend(loc='upper right')
            ylim([0,max(Srsso)*2.5])
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
        if jj > 0:
            pass
        else: 
            ax = subplot(212,sharex=ax)
            for crlb,x in zip(Crlbo,XX):
                cu.plotErrEllipse(ax,crlb,ecenter=(x,y),chi=1.0)

            ylabel('Across-Track Distance [m]')
            grid(True)
        xlabel('Along-Track Distance [m]')

    figure(20)
    ax=subplot(2,1,ll+1)
    plot(XX,Srsso,label='$\sqrt{\sigma_x^2+\sigma_y^2}$')
    plot(XX,Sigxo,label='$\sigma_x$')
    plot(XX,Sigyo,label='$\sigma_y$')
    if ll < 1:
        title('USBL with DVL')
        legend(loc='upper right')
        cu.annoteParams(ax,x,y,
                        dtUsbl,
                        varr,
                        vara,
                        varv,
                        vel,
                        varhdg,
                        loc=(0.01,0.95))
    else:
        xlabel('Along-Track Distance [m]')
    ylabel('Uncertainty [m]')
    grid(True)

show()



figure(1)
savefig('crlb_alongaxis_usblonly_twoscales.png',dpi=300)
figure(2)
savefig('crlb_alongaxis_usblonly_subplots_%.0f.png'%LL[0],dpi=300)
figure(3)
savefig('crlb_alongaxis_usbldvl_subplots_%.0f.png'%LL[0],dpi=300)
figure(4)
savefig('crlb_alongaxis_usblonly_1plot_%.0f.png'%LL[0],dpi=300)
figure(5)
savefig('crlb_alongaxis_usbldvl_1plot_%.0f.png'%LL[0],dpi=300)
figure(20)
savefig('crlb_alongaxis_usbldvl_twoscales.png',dpi=300)
figure(21)
savefig('crlb_alongaxis_usblonly_subplots_%.0f.png'%LL[1],dpi=300)
figure(22)
savefig('crlb_alongaxis_usbldvl_subplots_%.0f.png'%LL[1],dpi=300)
figure(23)
savefig('crlb_alongaxis_usblonly_1plot_%.0f.png'%LL[1],dpi=300)
figure(24)
savefig('crlb_alongaxis_usbldvl_1plot_%.0f.png'%LL[1],dpi=300)

