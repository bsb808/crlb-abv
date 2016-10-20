% USBL + DVL + Heading

% Assume path is straight line along the X-axis

% Setup
D = 1000;  % Track length [m]
% Assume 10 second update rate - that means every 2 meters
Nstations = 500;
Nusbl = Nstations-1;  % Number of USBL fixes
vel = 0.2;  % vehicle speed in m/s
% USBL Uncertaint
sig_r = 0.2;  % 1 sigma in range [m]
sig_b = 0.25*pi/180;  % 1 sigma in bearing [rad]
% DVL Uncertainty
sig_v = 0.003;  % velocity [m/s]
% Heading 
sig_h = 1.0*pi/180;  % heading [rad]

% Interstitial
dX = D/Nusbl;
dT = dX/vel;
% X positions

% Build C matrix one row at a time
C = [];
XX = linspace(dX,D,Nusbl);
for ii = 1:Nusbl
    % Range obs.
    cr = zeros(1,2*Nusbl);
    cr((ii-1)*2+1)=1;
    C = [C;cr];
    % Bearing obs.
    br = zeros(1,2*Nusbl);
    br(2*ii)=1/XX(ii);
    C = [C;br];
end

% Now DVL and Hdg
cx = zeros(1,2*Nusbl);
cx(1) = 1;
C = [C;cx];
cx = zeros(1,2*Nusbl);
cx(2) = 1;
C = [C;cx];

for ii = 2:Nusbl
    cx = zeros(1,2*Nusbl);
    cx((ii-2)*2+1)=-1;
    cx((ii-2)*2+3)=1;that 
    C = [C;cx];
    cy = zeros(1,2*Nusbl);
    cy((ii-2)*2+2)=-1;
    cy((ii-2)*2+4)=1;
    C = [C;cy];
end

% Build R matrix 
Nobs = size(C,1);
R = eye(Nobs);
for ii = 1:Nusbl
    iii = (ii-1)*2+1;
    R(iii,iii) = sig_r^2;
    jjj = (2*ii);
    R(jjj,jjj) = sig_b^2;
end
for ii = 2:Nusbl+1
    iii = Nusbl*2+(ii-2)*2+1;
    R(iii,iii) = dT*(sig_v^2);
    jjj =  Nusbl*2+(ii-2)*2+2;
    R(jjj,jjj) = dT*(sig_v^2)+(dX*sig_h)^2;
end

F = C'*(inv(R))*C;
lam = inv(F);
CC = sqrt(diag(lam));
alongT = CC(1:2:end);
acrossT = CC(2:2:end);
% CEP
cep = zeros(length(alongT),1);
for ii=1:length(alongT);
    sigx = alongT(ii);
    sigy = acrossT(ii);
    sigl = max([sigx,sigy]);
    sigs = min([sigx,sigy]);
    w = sigs/sigl;
    if w < 0.5
        cep(ii)= sigl*(0.67+0.8*w^2);
    else
        cep(ii)=0.59*sigl*(1+w);
    end
end
    
        

% Plot
figure(3)
XXX = [0,XX];
clf()
plot(XXX,[0;alongT],'b-','linewidth',2);
hold on
plot(XXX,[0;acrossT],'r-','linewidth',2);
plot(XXX,[0;cep],'k-','linewidth',2);
legend('Along-Track','Across-Track','CEP','location','northwest')
grid('on')
ylabel('Position Uncertainty [m]')
xlabel('Distance Traveled [m]')







