%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Script name:        Magnetics wireless charging
%   Author:             N. van Rooijen
%   Initial date:       24-04-2017
%   Version:            1.0
%   Comments:           script describing the transformer model for the
%                       coil of the wireless charging project.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%% electric constants\
f=50;       %frequency
omega_e=      f*2*pi;%electrical speed
rho=    1.75e-8;%electrical resistivity
Ip=     15;%current through primary coil.
Vp=70;%source voltage
%% geometric constants
%mode=   1;%1=vertical stacking, 0 is horizontal stacking
Np=     10;%number of turns primary coil
Ns=     10;%number of turns secondary coil
n=     90 ;%number of strands
ds=    0.001 ;%diameter of individual strand
D=     sqrt(n)*1.154*ds ;%diameter of individual wire
a= Np/Ns; %turn ratio
rp=    0.5 ;%radius primary coil
rs=     0.4;%radius secondary coil
lp= 2*pi*(rp+(Np/2)*D)   ;%length of each turn on primary coil
ls= 2*pi*(rs+(Ns/2)*D)     ;%length of each turn on secondary coil
h=     0.1 ;%height between prim. and sec. coil.
NumSt=  10;%number of steps in fluxdens function.
%% magnetic constants
Uo=4*pi*10^(-7);%magnetic constant of air
Ur=1.2566290e-6;%magnetic permeability copper
NumEl=10; %number of horizontal filaments in SELF function
Rp_outp=rp+(D*Np); Rp_inp=rp;%radius primary coil
Rp_outs=rs+(D*Ns); Rp_ins=rs;%radius secondary coil
Lp=(Np^2)*SELF(Rp_outp,Rp_inp,NumEl);%total self inductance primary coil; %self inductance coil 1
Ls=(Ns^2)*SELF(Rp_outs,Rp_ins,NumEl);%total self inductance secondary coil; %self inductance coil 2

%% total flux density from primary coil
flx_p=0;
for i=1:Np %horizontal stacking
    rs=rs+(D); %radius secondary coil
    flx_p=flx_p+FLUXDENS(rp,rs,h,Ip,NumSt);
end
flux_linkage=0;
for j=1:Ns
    rsec=rs+j*D;
    set=1;
    for p=1:NumSt
        Dis(p)=(rp/NumSt)*p;
        if rsec > Dis(p) && set==1
            test=flx_p(p);          
            if p==1
                flux_linkage=flux_linkage+(flx_p(p)*pi*(Dis(p))^2);
            else
                flux_linkage=flux_linkage+((flx_p(p)*pi*(Dis(p))^2)-(flx_p(p)*pi*(Dis(p-1))^2));
            end
        elseif rsec <= Dis
            flux_linkage=flux_linkage+flx(p)*pi*rsec^2;
            set=0;
        end
    end
    Ssn=pi*((rs+j*D)^2); %enclosed area of secondary turn.
end
%mutual flux
M=flux_linkage/Ip;
%% calculations
%% circuit parameter model
Llks= Lp-a*M; %leakage inductance primary coil
Llkp= Ls - (1/a)*M; %leakage inductance secondary coil
Lm=a*M; %magnetizing inductance
%% skin effect:
delta=sqrt(rho/(pi*omega_e*Uo*Ur));%skin depth
if(delta>=(ds/2))
    ks=1; %skin effect correction factor
else
    ks=delta/(rs);
end
%% calculation of resistances
Rp_dc= (rho * 4*lp)/(pi*n*((ds)^2)*ks);%copper dc resistance primary coil
xs_sqr_p=9*pi*omega_e*(10^-7)/Rp_dc;
ys_p=(xs_sqr_p^2)/(192+(xs_sqr_p^2));
Rp_ac=Rp_dc*(1+ys_p); %ac resistance primary coil

Rs_dc= (rho * 4*ls)/(pi*n*((ds)^2)*ks);%copper dc resistance secondary coil
xs_sqr_s=(9*pi*omega_e*(10^-7))/Rs_dc;
ys_s=(xs_sqr_s^2)/(192+(xs_sqr_s^2));
Rs_ac=Rs_dc*(1+ys_s); %ac resistance secondary coil
%total resistance of both primary and secondary coils.
Rp=Rp_dc+Rp_ac;
Rs=Rs_dc+Rs_ac;
%% equivalent resistances Z1, Z2 and Z3:
Z1=(Rp+1i*omega_e*Llkp);
Z2=(1i*omega_e*Lm);
Z3=Rs+1i*omega_e*Llks;
%% power flow within circuit
Vpc= Vp-(Ip*Z1); %voltage through primary coil.
Vsc=Vpc/a; % voltage in secondary coil
Im=Vpc./Z2; %magnetizing current
Ipc=Ip-Im; %current primary coil, after magnetizing losses
Isec=Ipc*a; %current in secondary coil
Vs=Vsc-Isec.*Z3; %output voltage.
Sout=Vs.*Isec;
Pin=Vp.*Ip;
pf=cos(angle(Vs)-angle(Isec));
Pout=Sout*pf;
eff=abs(Pout)/Pin;
% %% linked flux
% flux=zeros(NumSt);
% flux_avg=zeros(NumSt);
% DN=D*
% dM=D*Np/M;
% for i=1:Nf
% for j=1:Mf
% for k=1:Pf
% for l=1:Qf
%     Rp=Rp+(j*(D*Np))/Mf;
%     Rs=Rs+(l*(D*Np))/Qf;
%     Rp=Rp+(j*(D*Np))/Mf;
%     Rp=Rp+(j*(D*Np))/Mf;
%
%     flux= Flux+(FLUXDENS(Rp,Rs,h,I,NumSt)*Ssn);
% end
% end
% end
% end
% flux_avg=flux/N*M*P*Q;
%

%Mutual inductance

%% self inductance
    function[Lp]=SELF(Rp_out,Rp_in,NumEL)
        R=1e-3;%assumewireradiusequalto1 mm;
        mu0= 4*pi*1.0e-7;
        da=linspace(Rp_in,Rp_out,NumEL);
        db=linspace(Rp_in,Rp_out,NumEL);
        lp=0;
        for aa=1:length(da)
            for bb= 1:length(db)
                a=da(aa);
                b=db(bb);
                if a~=b%mutualinductancebetweentwoturns
                    k=sqrt((4*a*b)/((a+b)^2));
                    [EK EE]=ellipke(k^2);
                    mut= mu0*sqrt(a*b)*(((2/k)-k)*EK-(2/k)*EE);
                    lp=lp+mut;
                else%self-inductanceofa singleturn
                    lp=lp+mu0*a*(log((8*a)/R)-2);
                end
            end
        end
        Lp= lp/(length(da)*length(db));
    end
    
%% calculation of (mutual)flux:
%function fluxdensity,
function[Bz,rs]=FLUXDENS(Rp,Rs,h,I,NumSt)
    mu0= 4*pi*1.0e-7;
    rs=linspace(0,Rs,NumSt);%Accuracycanbeimprovedbyincreasingthenumberofelements
    for nn=1:length(rs)
        alpha=rs(nn)/Rp;
        beta=h/Rp;
        Q = ((1+alpha)^2+beta^2);
        k =sqrt((4*alpha)/Q);
        [KK EK]=ellipke(k^2);
        B0=(I*mu0)/(2*Rp);
        Bz(nn)= B0*(1/(pi*sqrt(Q)))*(EK*((1-alpha^2-beta^2)/(Q-4*alpha))+KK);
    end
end