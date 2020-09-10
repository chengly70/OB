function [fr_MC,fr_PGC,v_mcAvg,v_pgcAvg]=getTrans_PGC(tEnd,Ibg,wghts)
% [fr_MC,fr_PGC,v_mcAvg,v_pgcAvg]=getTrans_PGC(tEnd,Ibg,wghts);
% only saving counts, no GC cells; getting MC Transfer function coupled to PGC
% tEnd=time in ms; Ibg=[ PGC ; MC ];
% Ibg/ev in muAmps/cm^2 (smaller), dc current and random AMPA bg for all cells
% wghts=[wM2P;wM2Pn;wP2M]
% wghts=[100;100;100];
% SUBnet of only PGC/MC (Li & Cleland); Ibg:(PGC;MC)
% !!! t is in ms  !!! return spike times (spTimes) 
% meant for SHORT times to show spikes AND voltage (not 200,000ms or crazyness)
% using Crank-Nicholson/ImplicTrapezoidalRule and RK4 time-stepping

rng('shuffle'); %random seed

Nrlz=10; %number of realizations

dt=0.05; %assuming equally spaced
Lt=round(tEnd/dt)+1;

F=9.64853329e4; %Coulombs per Mole
T=273.15+35; %Kelvin
RTovF=8314.472*T/F; %R*1000 so in mV
thick=1; %thickeness of membrane shell, 1micron for Mitral Cell
thick_g=0.2; %thickeness of membrane shell, 0.2micron for GC/PGC

% PARAMS 
C=1.2; %micro-F/cm^2 for MC & PGC
Rm=30; %kilo-Ohm*cm^2 for MC & GC
Rm_p=20; %kilo-Ohm*cm^2 for PGC
gL=1/Rm; %leak conduct (mS/cm^2)
gL_p=1/Rm_p; %leak conduct (mS/cm^2)
tau_m=C*Rm; %membrane time constant (ms)
diam=20; %microns
diam_p=8; %microns
Ra=70; %Ohm*cm
Ra_p=80; %Ohm*cm
taup_m=C*Rm_p; %membrane time constant (ms)

%PGC: Spine(1), Dend/Soma(2)
%MC: Tuft(3=2+1), PDend(4=2+2), Soma(5=2+3), Dend(6=2+4)
Len_p=[20;2]; %lengths of MC pieces (microns)
Len_m=[4;74;5;100]; 
dx = .25; %1 / (#pts per micron)
Spcv=(0:dx:sum(Len_m)-dx)'; %true length (microns)

Ntt=length(Spcv); %total # compartments in MC
Nttp=length( (0:dx:sum(Len_p)-dx)' );
Ncm=Len_m./dx; %# compartments IF dx divides all entries of Len
Ncmp=Len_p./dx;
Ncmsum=cumsum(Ncm);
Ncmsump=cumsum(Ncmp);
if(Ncmsum(end)~=Ntt || Ncmsump(end)~=Nttp)
    disp('index mismatch!');
    return
end

%record all midpoints
idRec_p=round(Ncmp./2); %mid-point
idRec_p(2)=cumsum(Ncmp(1))+idRec_p(2); %true index
idRec=round(Ncm./2); %mid-point
idRec(2:end)=cumsum(Ncm(1:end-1))+idRec(2:end); %true index
%___ ONLY record at soma for each cell ___
idRec_p=idRec_p(2);
idRec=idRec(3);
lenRec=length(idRec);
lenRec_p=length(idRec_p);
% specify where put gating vars
idGv=(1:Ntt)';
nmGv=length(idGv); %# gating variables
idSegTyp=zeros(nmGv,1); %identity of segment type: Tuft(1), PDend(2), Soma(3), Dend(4)
idSegTyp(1:Ncmsum(1))=1; %putting them everywhere
for j=2:4
    idSegTyp(Ncmsum(j-1)+1:Ncmsum(j))=j;
end
idGv_p=(1:Nttp)';
nmGv_p=length(idGv_p); %# gating variables
idSegTyp_p=zeros(nmGv_p,1); %identity of segment type: Spine(1), Dend(2)/Soma
idSegTyp_p(1:Ncmsump(1))=1; %putting them everywhere
j=2;
    idSegTyp_p(Ncmsump(j-1)+1:Ncmsump(j))=j; %only 2 compartms

%...specify injected current here..
IdcP=sparse(zeros(nmGv_p,1)); %index corresponds to idGv_p(j) value
IdcP(1:round(Ncmp(1)/4))=Ibg(1); %originally Ibg in micro-Amps/cm^2
IdcM=sparse(zeros(nmGv,1));
IdcM(1:Ncm(1))=Ibg(2);

%maximal conduct
% -- for PGC --
gNa_3=[20;50]./gL_p; %(mS/cm^2) dimLess, Spine, Dend, Soma
gA_3=[30;10]./gL_p;
gM_3=[0;1]./gL_p; %only in soma
gDR_3=[5;20]./gL_p; 
gH_3=[0.2;0]./gL_p;  %only in PGC dendrite
gCaPN_p_3=[1;0]./gL_p; %only in dendrite
gCaT_p_3=[3;0]./gL_p;
gCAN_3=[1.5;0]./gL_p; %only in dendrite
gKCa_p_3=[2;0]./gL_p; %only in soma
% -- for MC --
gNa_1=[20;20;50;25]./gL; %(mS/cm^2) dimLess, Tuft, PDend, Soma, Dend
gNaP_1=[0.1;0.1;0.2;0.15]./gL; 
gDR_1=[10;10;10;20]./gL; 
gA_1=[0;0;1;0]./gL; %only in soma
gKS_1=[18;18;20;20]./gL;
gCaL_1=[0.2;0.2;0.4;.05]./gL;
gKCa_1=[0;0;1;0]./gL; %only in soma
%-- augment to lineup with idGv & segment identity (idSegTyp) ---
gNa_p=zeros(nmGv_p,1);
gDR_p=zeros(nmGv_p,1);
gA_p=zeros(nmGv_p,1);
gM_p=zeros(nmGv_p,1);
gH=zeros(nmGv_p,1);
gCaPN_p=zeros(nmGv_p,1);
gCaT_p=zeros(nmGv_p,1);
gCAN_p=zeros(nmGv_p,1);
gKCa_p=zeros(nmGv_p,1);
for j=1:2
    gNa_p(idSegTyp_p==j)=gNa_3(j);
    gDR_p(idSegTyp_p==j)=gDR_3(j);
    gA_p(idSegTyp_p==j)=gA_3(j);
    gM_p(idSegTyp_p==j)=gM_3(j);
    gH(idSegTyp_p==j)=gH_3(j);
    gCaPN_p(idSegTyp_p==j)=gCaPN_p_3(j);
    gCaT_p(idSegTyp_p==j)=gCaT_p_3(j);
    gCAN_p(idSegTyp_p==j)=gCAN_3(j);
    gKCa_p(idSegTyp_p==j)=gKCa_p_3(j);
end
% -- repeat for MC --
gNa=zeros(nmGv,1);
gNaP=zeros(nmGv,1);
gDR=zeros(nmGv,1);
gA=zeros(nmGv,1);
gKS=zeros(nmGv,1);
gCaL=zeros(nmGv,1);
gKCa=zeros(nmGv,1);
for j=1:4
    gNa(idSegTyp==j)=gNa_1(j);
    gNaP(idSegTyp==j)=gNaP_1(j);
    gDR(idSegTyp==j)=gDR_1(j);
    gA(idSegTyp==j)=gA_1(j);
    gKS(idSegTyp==j)=gKS_1(j);
    gCaL(idSegTyp==j)=gCaL_1(j);
    gKCa(idSegTyp==j)=gKCa_1(j);
end
%reversal Poten
El=-60; %mV
Ek=-80;
Ena=45;
El_p=-65; %mV
Eh=0; 
Ecan=0;
Excs=0; %excit synapses (AMPA & NMDA same)
Einh=-80; %inhib synapses (GABAa)

lambd_p=500*.5*sqrt(diam_p*Rm_p*1000/Ra_p)/100; %space-length in microns
lambd=500*.5*sqrt(diam*Rm*1000/Ra)/100; %space-length in microns
%matrices for diffusion part
MD=-2*diag(ones(Nttp,1))+diag(ones(Nttp-1,1),1)+diag(ones(Nttp-1,1),-1);
MD(1,1)=-1; MD(Nttp,Nttp)=-1; %no-flux BC
Al_p = sparse( eye(Nttp) - .5*dt/taup_m*(lambd_p/dx)^2*MD );
Ar_p = sparse( eye(Nttp) + .5*dt/taup_m*(lambd_p/dx)^2*MD );
bv_p = zeros(Nttp,1);
MD=-2*diag(ones(Ntt,1))+diag(ones(Ntt-1,1),1)+diag(ones(Ntt-1,1),-1);
MD(1,1)=-1; MD(Ntt,Ntt)=-1; %no-flux BC
Al = sparse( eye(Ntt) - .5*dt/tau_m*(lambd/dx)^2*MD );
Ar = sparse( eye(Ntt) + .5*dt/tau_m*(lambd/dx)^2*MD );
bv = zeros(Ntt,1);
% -- output -- 
if(lenRec_p==1)
    fr_PGC=0;
else
    fr_PGC=zeros(lenRec_p,1); 
end
if(lenRec==1)
    fr_MC=0;
else
    fr_MC=zeros(lenRec,1);
end
v_mcAvg=0; %average MC voltage
v_pgcAvg=0; %average PGC voltage

vAll_pgc=-75.8*ones(Nttp,1); %all voltages
vAll_mc=-70*ones(Ntt,1); %all voltages

minTspk=3; %min time between spikes (ms)
indLst_p=zeros(lenRec_p,1); %index of time of last spike
indLst=zeros(lenRec,1); %index of time of last spike
vltThres=-5; %(mV) voltage threshold for spike

%#rows=#places with gating variables
caCon=0.05*ones(nmGv,1); %only in soma, micrMol/l?
Eca=zeros(nmGv,1); %dyn variable
xVr=zeros(nmGv,11); %gating variables 
xiR1=zeros(nmGv,11); %aux for RK4
xiR2=zeros(nmGv,11); %aux for RK4
xiR3=zeros(nmGv,11); %aux for RK4
xiR4=zeros(nmGv,11); %aux for RK4
caCr1=zeros(nmGv,1);
caCr2=zeros(nmGv,1);
caCr3=zeros(nmGv,1);
caCr4=zeros(nmGv,1);
rhsV1=zeros(nmGv,1);
rhsV2=zeros(nmGv,1);
rhsV3=zeros(nmGv,1);
rhsV4=zeros(nmGv,1);

caCon_p=0.05*ones(nmGv_p,1); %not in soma, micrMol/l?
Eca_p=zeros(nmGv_p,1); %dyn variable
pVr=zeros(nmGv_p,13); %gating variables 
piR1=zeros(nmGv_p,13); %aux for RK4
piR2=zeros(nmGv_p,13); %aux for RK4
piR3=zeros(nmGv_p,13); %aux for RK4
piR4=zeros(nmGv_p,13); %aux for RK4
caCr1p=zeros(nmGv_p,1);
caCr2p=zeros(nmGv_p,1);
caCr3p=zeros(nmGv_p,1);
caCr4p=zeros(nmGv_p,1);
rhsV1p=zeros(nmGv_p,1);
rhsV2p=zeros(nmGv_p,1);
rhsV3p=zeros(nmGv_p,1);
rhsV4p=zeros(nmGv_p,1);

%synapses
synP2M=zeros(Ncm(1),1);  %vect of #inputs to MC from PGC
    sidP2M=( idGv <= Ncm(1) ); %start of MC Compart (Tuft)  
    numP2M=length(synP2M); %#inputs to MC from PGC
synM2P=zeros(round(Ncmp(1)/4),1); %vect of #inputs to PGC (AMPA)
synM2Pn=zeros(round(Ncmp(1)/4),1); %vect of #inputs to PGC (NMDA)
    sidM2P=( idGv_p <= round(Ncmp(1)/4) ); %start of PGC (Spine)
    numM2P=length(synM2P); %#inputs to PGC from MC
%weights for synapses
wM2P=wghts(1)/length(synM2P);
wM2Pn=wghts(2)/length(synM2P);
wP2M=wghts(3)/length(synP2M);
%background ampa synapses
bgAmp_p=zeros(round(Ncmp(1)/4),1); %background AMPA inputs modeled by low-pass filter
bgAmp_m=zeros(Ncm(1),1); 
wP_bg=2;  %weights for background synapses
wM_bg=2;

% temp vars for presyn volt
tmpPreVp=zeros(numP2M,1);
tmpPreVmp=zeros(numM2P,1);
    
%only used for for Presv Interp 
    idPre_M2P=(1:numM2P)';%presyVol P->M
    idPre_P2M=(1:numP2M)';
%make interp biophys, stayin bounds of presyn
    idP_P2M=(1 : (numM2P-1)/(numP2M-1) : numM2P)';
    idP_M2P=(1 : (numP2M-1)/(numM2P-1) : numP2M)';
%params for background noisy inputs
lam=1/25; %kHz
bgJmp=1; %size of Background Jump

%run initial loop to get rid of transients
ppIn=(rand(round(Ncmp(1)/4),round(1000/dt)) < dt*lam);
ppInM=(rand(Ncm(1),round(1000/dt)) < dt*lam);
for j=1:round(1000/dt)
    bgAmp_p = bgAmp_p + dt/5*(-bgAmp_p) + bgJmp*ppIn(:,j);
    bgAmp_m=bgAmp_m + dt/5*(-bgAmp_m) + bgJmp*ppInM(:,j);
    %get presynaptic voltages
    tmpPreVp=interp1(idPre_M2P,vAll_pgc(sidM2P),idP_P2M,'nearest');
    tmpPreVmp=interp1(idPre_P2M,vAll_mc(sidP2M),idP_M2P,'nearest');
    %synapses
    synP2M=synP2M+dt*(1./(1.25*(1+exp(-((tmpPreVp+40)/2)))).*(1-synP2M)-1/18*synP2M);
    synM2P=synM2P+dt*(1./(1+exp(-(tmpPreVmp)/0.2)).*(1-synM2P)-1/5.5*synM2P);
    synM2Pn=synM2Pn+dt*(1/52./(1+exp(-(tmpPreVmp)/0.2)).*(1-synM2Pn)-1/343*synM2Pn);
    % --- step 1 ---
    rhsV1p=(vAll_pgc-El_p)+gNa_p.*pVr(:,1).^3.*pVr(:,2).*(vAll_pgc-Ena)+gA_p.*pVr(:,3).*pVr(:,4).*(vAll_pgc-Ek)+...
        gM_p.*pVr(:,5).*(vAll_pgc-Ek)+gH.*pVr(:,6).*(vAll_pgc-Eh)+...
        gCaPN_p.*pVr(:,7).^2.*pVr(:,8).*(vAll_pgc-Eca_p)+gKCa_p.*pVr(:,9).*(vAll_pgc-Ek)+...
        gDR_p.*pVr(:,10).*(vAll_pgc-Ek)+gCaT_p.*(pVr(:,11)).^2.*pVr(:,12).*(vAll_pgc-Eca_p)+...
        gCAN_p.*caCon_p./(200+caCon_p).*pVr(:,13).*(vAll_pgc-Ecan);
    rhsV1p=dt/taup_m*(IdcP-rhsV1p);
    piR1(:,1)=dt*2.1*(ana_g(vAll_pgc)./(ana_g(vAll_pgc)+bna_g(vAll_pgc))-pVr(:,1))./tauM(vAll_pgc);
    piR1(:,2)=dt*2.1*(hInf(vAll_pgc)-pVr(:,2))./tauH(vAll_pgc);
    piR1(:,3)=dt*3.3./atau(vAll_pgc).*(aminf_g(vAll_pgc)-pVr(:,3)); %I_A
    piR1(:,4)=dt*3.3./hatau_g(vAll_pgc).*(ahinf_g(vAll_pgc)-pVr(:,4)); %I_A
    piR1(:,5)=dt*(minfMusc(vAll_pgc)-pVr(:,5))./tauMusc(vAll_pgc); %I_M
    piR1(:,6)=dt*2.1./hcurmtau(vAll_pgc).*(hcurminf(vAll_pgc)-pVr(:,6));
    piR1(:,7)=dt*(capminf(vAll_pgc)-pVr(:,7))./capmtau(vAll_pgc); %I_CaPN
    piR1(:,8)=dt*(caphinf(vAll_pgc)-pVr(:,8))./caphtau(vAll_pgc);
    piR1(:,9)=dt*(kcaa(vAll_pgc,caCon_p).*(1-pVr(:,9)) - 0.05*pVr(:,9));
    piR1(:,10)=dt*3.3*(mss_dr(vAll_pgc)-pVr(:,10))./taum_dr(vAll_pgc);    %I_DR
    piR1(:,11)=dt*(catminfP(vAll_pgc)-pVr(:,11))./cattaumP(vAll_pgc); %I_CaT
    piR1(:,12)=dt*(cathinf(vAll_pgc)-pVr(:,12))./cattauh(vAll_pgc); 
    piR1(:,13)=dt*(canminf(vAll_pgc)-pVr(:,13))./cantau(vAll_pgc); %I_CAN
    Ica=gCaPN_p.*pVr(:,7).^2.*pVr(:,8).*(vAll_pgc-Eca_p)+gCaT_p.*(pVr(:,11)).^2.*pVr(:,12).*(vAll_pgc-Eca_p); 
    caCr1p=dt*(-Ica./(2*F*thick_g) + (.05-caCon_p)/10 );
    rhsV1=(vAll_mc-El)+gNa.*xVr(:,1).^3.*xVr(:,2).*(vAll_mc-Ena)+gNaP.*minf(vAll_mc).*(vAll_mc-Ena)+...
        gA.*xVr(:,3).*xVr(:,4).*(vAll_mc-Ek)+gKS.*xVr(:,5).*xVr(:,6).*(vAll_mc-Ek)+...
        gCaL.*xVr(:,7).*xVr(:,8).*(vAll_mc-Eca)+gKCa.*xVr(:,9).*(vAll_mc-Ek)+gDR.*xVr(:,10).^2.*xVr(:,11).*(vAll_mc-Ek);
    rhsV1=dt/tau_m*(IdcM-rhsV1);
    xiR1(:,1)=dt*(ana(vAll_mc).*(1-xVr(:,1)) - bna(vAll_mc).*xVr(:,1));
    xiR1(:,2)=dt*(ahna(vAll_mc).*(1-xVr(:,2)) - bhna(vAll_mc).*xVr(:,2));
    xiR1(:,3)=dt*3.3./atau(vAll_mc).*(aminf(vAll_mc)-xVr(:,3)); %I_A
    xiR1(:,4)=dt*3.3./hatau(vAll_mc).*(ahinf(vAll_mc)-xVr(:,4)); %I_A
    xiR1(:,5)=dt/10*(ksminf(vAll_mc)-xVr(:,5)); %I_KS
    xiR1(:,6)=dt./kshtau(vAll_mc).*(kshinf(vAll_mc)-xVr(:,6));
    xiR1(:,7)=dt*(cala(vAll_mc).*(1-xVr(:,7)) - calb(vAll_mc).*xVr(:,7)); %I_CaL
    xiR1(:,8)=dt*(calha(vAll_mc).*(1-xVr(:,8)) - calhb(vAll_mc).*xVr(:,8));
    xiR1(:,9)=dt*(kcaa(vAll_mc,caCon).*(1-xVr(:,9)) - 0.05*xVr(:,9));
    xiR1(:,10)=dt*(nss_dr(vAll_mc)-xVr(:,10))./taun_dr(vAll_mc);    %I_DR
    xiR1(:,11)=dt*(kss_dr(vAll_mc)-xVr(:,11))./50;
    Ica=gCaL.*xVr(:,7).*xVr(:,8).*(vAll_mc-Eca);
    caCr1=dt*(-Ica./(2*F*thick) + (.05-caCon)/10 );
    % --- step 2 ---
    rhsV2p=(vAll_pgc+.5*rhsV1p-El_p)+gNa_p.*(pVr(:,1)+.5*piR1(:,1)).^3.*(pVr(:,2)+.5*piR1(:,2)).*(vAll_pgc+.5*rhsV1p-Ena)+...
        gA_p.*(pVr(:,3)+.5*piR1(:,3)).*(pVr(:,4)+.5*piR1(:,4)).*(vAll_pgc+.5*rhsV1p-Ek)+...
        gM_p.*(pVr(:,5)+.5*piR1(:,5)).*(vAll_pgc+.5*rhsV1p-Ek)+gH.*(pVr(:,6)+.5*piR1(:,6)).*(vAll_pgc+.5*rhsV1p-Eh)+...
        gCaPN_p.*(pVr(:,7)+.5*piR1(:,7)).^2.*(pVr(:,8)+.5*piR1(:,8)).*(vAll_pgc+.5*rhsV1p-Eca_p)+gKCa_p.*(pVr(:,9)+...
        .5*piR1(:,9)).*(vAll_pgc+.5*rhsV1p-Ek)+gDR_p.*(pVr(:,10)+.5*piR1(:,10)).*(vAll_pgc+.5*rhsV1p-Ek)+...
        gCaT_p.*(pVr(:,11)+.5*piR1(:,11)).^2.*(pVr(:,12)+.5*piR1(:,12)).*(vAll_pgc+.5*rhsV1p-Eca_p)+...
        gCAN_p.*(caCon_p+.5*caCr1p)./(200+(caCon_p+.5*caCr1p)).*(pVr(:,13)+.5*piR1(:,13)).*(vAll_pgc+.5*rhsV1p-Ecan);
    rhsV2p=dt/taup_m*(IdcP-rhsV2p);
    piR2(:,1)=dt*2.1*(ana_g(vAll_pgc+.5*rhsV1p)./(ana_g(vAll_pgc+.5*rhsV1p)+bna_g(vAll_pgc+.5*rhsV1p)) -(pVr(:,1)+.5*piR1(:,1)) )./tauM(vAll_pgc+.5*rhsV1p);
    piR2(:,2)=dt*2.1*(hInf(vAll_pgc+.5*rhsV1p)-(pVr(:,2)+.5*piR1(:,2)))./tauH(vAll_pgc+.5*rhsV1p);
    piR2(:,3)=dt*3.3./atau(vAll_pgc+.5*rhsV1p).*(aminf_g(vAll_pgc+.5*rhsV1p)-(pVr(:,3)+.5*piR1(:,3))); %I_A
    piR2(:,4)=dt*3.3./hatau_g(vAll_pgc+.5*rhsV1p).*(ahinf_g(vAll_pgc+.5*rhsV1p)-(pVr(:,4)+.5*piR1(:,4))); %I_A
    piR2(:,5)=dt*(minfMusc(vAll_pgc+.5*rhsV1p)-(pVr(:,5)+.5*piR1(:,5)))./tauMusc(vAll_pgc+.5*rhsV1p); %I_M
    piR2(:,6)=dt*2.1./hcurmtau(vAll_pgc+.5*rhsV1p).*(hcurminf(vAll_pgc+.5*rhsV1p)-(pVr(:,6)+.5*piR1(:,6)));
    piR2(:,7)=dt*(capminf(vAll_pgc+.5*rhsV1p)-(pVr(:,7)+.5*piR1(:,7)))./capmtau(vAll_pgc+.5*rhsV1p); %I_CaPN
    piR2(:,8)=dt*(caphinf(vAll_pgc+.5*rhsV1p)-(pVr(:,8)+.5*piR1(:,8)))./caphtau(vAll_pgc+.5*rhsV1p);
    piR2(:,9)=dt*(kcaa(vAll_pgc+.5*rhsV1p,caCon_p+.5*caCr1p).*(1-(pVr(:,9)+.5*piR1(:,9))) - 0.05*(pVr(:,9)+.5*piR1(:,9)));
    piR2(:,10)=dt*3.3*(mss_dr(vAll_pgc+.5*rhsV1p)-(pVr(:,10)+.5*piR1(:,10)))./taum_dr(vAll_pgc+.5*rhsV1p);    %I_DR
    piR2(:,11)=dt*(catminfP(vAll_pgc+.5*rhsV1p)-(pVr(:,11)+.5*piR1(:,11)))./cattaumP(vAll_pgc+.5*rhsV1p); %I_CaT
    piR2(:,12)=dt*(cathinf(vAll_pgc+.5*rhsV1p)-(pVr(:,12)+.5*piR1(:,12)))./cattauh(vAll_pgc+.5*rhsV1p); 
    piR2(:,13)=dt*(canminf(vAll_pgc+.5*rhsV1p)-(pVr(:,13)+.5*piR1(:,13)))./cantau(vAll_pgc+.5*rhsV1p); %I_CAN
    Ica=gCaPN_p.*(pVr(:,7)+.5*piR1(:,7)).^2.*(pVr(:,8)+.5*piR1(:,8)).*(vAll_pgc+.5*rhsV1p-Eca_p)+...
        gCaT_p.*(pVr(:,11)+.5*piR1(:,11)).^2.*(pVr(:,12)+.5*piR1(:,12)).*(vAll_pgc+.5*rhsV1p-Eca_p);
    caCr2p=dt*(-Ica./(2*F*thick_g) + (.05-(caCon_p+.5*caCr1p))/10 );
    rhsV2=(vAll_mc+.5*rhsV1-El)+gNa.*(xVr(:,1)+.5*xiR1(:,1)).^3.*(xVr(:,2)+.5*xiR1(:,2)).*(vAll_mc+.5*rhsV1-Ena)+...
        gNaP.*minf(vAll_mc+.5*rhsV1).*(vAll_mc+.5*rhsV1-Ena)+gA.*(xVr(:,3)+.5*xiR1(:,3)).*(xVr(:,4)+...
        .5*xiR1(:,4)).*(vAll_mc+.5*rhsV1-Ek)+gKS.*(xVr(:,5)+.5*xiR1(:,5)).*(xVr(:,6)+.5*xiR1(:,6)).*(vAll_mc+.5*rhsV1-Ek)+...
        gCaL.*(xVr(:,7)+.5*xiR1(:,7)).*(xVr(:,8)+.5*xiR1(:,8)).*(vAll_mc+.5*rhsV1-Eca)+gKCa.*(xVr(:,9)+...
        .5*xiR1(:,9)).*(vAll_mc+.5*rhsV1-Ek)+gDR.*(xVr(:,10)+.5*xiR1(:,10)).^2.*(xVr(:,11)+.5*xiR1(:,11)).*(vAll_mc+.5*rhsV1-Ek);
    rhsV2=dt/tau_m*(IdcM-rhsV2);
    xiR2(:,1)=dt*(ana(vAll_mc+.5*rhsV1).*(1-(xVr(:,1)+.5*xiR1(:,1))) - bna(vAll_mc+.5*rhsV1).*(xVr(:,1)+.5*xiR1(:,1)));
    xiR2(:,2)=dt*(ahna(vAll_mc+.5*rhsV1).*(1-(xVr(:,2)+.5*xiR1(:,2))) - bhna(vAll_mc+.5*rhsV1).*(xVr(:,2)+.5*xiR1(:,2)));
    xiR2(:,3)=dt*3.3./atau(vAll_mc+.5*rhsV1).*(aminf(vAll_mc+.5*rhsV1)-(xVr(:,3)+.5*xiR1(:,3))); %I_A
    xiR2(:,4)=dt*3.3./hatau(vAll_mc+.5*rhsV1).*(ahinf(vAll_mc+.5*rhsV1)-(xVr(:,4)+.5*xiR1(:,4))); %I_A
    xiR2(:,5)=dt/10*(ksminf(vAll_mc+.5*rhsV1)-(xVr(:,5)+.5*xiR1(:,5))); %I_KS
    xiR2(:,6)=dt./kshtau(vAll_mc+.5*rhsV1).*(kshinf(vAll_mc+.5*rhsV1)-(xVr(:,6)+.5*xiR1(:,6)));
    xiR2(:,7)=dt*(cala(vAll_mc+.5*rhsV1).*(1-(xVr(:,7)+.5*xiR1(:,7))) - calb(vAll_mc+.5*rhsV1).*(xVr(:,7)+.5*xiR1(:,7))); %I_CaL
    xiR2(:,8)=dt*(calha(vAll_mc+.5*rhsV1).*(1-(xVr(:,8)+.5*xiR1(:,8))) - calhb(vAll_mc+.5*rhsV1).*(xVr(:,8)+.5*xiR1(:,8)));
    xiR2(:,9)=dt*(kcaa(vAll_mc+.5*rhsV1,caCon+.5*caCr1).*(1-(xVr(:,9)+.5*xiR1(:,9))) - 0.05*(xVr(:,9)+.5*xiR1(:,9)));
    xiR2(:,10)=dt*(nss_dr(vAll_mc+.5*rhsV1)-(xVr(:,10)+.5*xiR1(:,10)))./taun_dr(vAll_mc+.5*rhsV1);    %I_DR
    xiR2(:,11)=dt*(kss_dr(vAll_mc+.5*rhsV1)-(xVr(:,11)+.5*xiR1(:,11)))./50;
    Ica=gCaL.*(xVr(:,7)+.5*xiR1(:,7)).*(xVr(:,8)+.5*xiR1(:,8)).*(vAll_mc+.5*rhsV1-Eca);
    caCr2=dt*(-Ica./(2*F*thick) + (.05-(caCon+.5*caCr1))/10 );
    % --- step 3 ---
    rhsV3p=(vAll_pgc+.5*rhsV2p-El_p)+gNa_p.*(pVr(:,1)+.5*piR2(:,1)).^3.*(pVr(:,2)+.5*piR2(:,2)).*(vAll_pgc+.5*rhsV2p-Ena)+...
        gA_p.*(pVr(:,3)+.5*piR2(:,3)).*(pVr(:,4)+.5*piR2(:,4)).*(vAll_pgc+.5*rhsV2p-Ek)+...
        gM_p.*(pVr(:,5)+.5*piR2(:,5)).*(vAll_pgc+.5*rhsV2p-Ek)+gH.*(pVr(:,6)+.5*piR2(:,6)).*(vAll_pgc+.5*rhsV2p-Eh)+...
        gCaPN_p.*(pVr(:,7)+.5*piR2(:,7)).^2.*(pVr(:,8)+.5*piR2(:,8)).*(vAll_pgc+.5*rhsV2p-Eca_p)+gKCa_p.*(pVr(:,9)+...
        .5*piR2(:,9)).*(vAll_pgc+.5*rhsV2p-Ek)+gDR_p.*(pVr(:,10)+.5*piR2(:,10)).*(vAll_pgc+.5*rhsV2p-Ek)+...
        gCaT_p.*(pVr(:,11)+.5*piR2(:,11)).^2.*(pVr(:,12)+.5*piR2(:,12)).*(vAll_pgc+.5*rhsV2p-Eca_p)+...
        gCAN_p.*(caCon_p+.5*caCr2p)./(200+(caCon_p+.5*caCr2p)).*(pVr(:,13)+.5*piR2(:,13)).*(vAll_pgc+.5*rhsV2p-Ecan);
    rhsV3p=dt/taup_m*(IdcP-rhsV3p);
    piR3(:,1)=dt*2.1*(ana_g(vAll_pgc+.5*rhsV2p)./(ana_g(vAll_pgc+.5*rhsV2p)+bna_g(vAll_pgc+.5*rhsV2p)) -(pVr(:,1)+.5*piR2(:,1)) )./tauM(vAll_pgc+.5*rhsV2p);
    piR3(:,2)=dt*2.1*(hInf(vAll_pgc+.5*rhsV2p)-(pVr(:,2)+.5*piR2(:,2)))./tauH(vAll_pgc+.5*rhsV2p);
    piR3(:,3)=dt*3.3./atau(vAll_pgc+.5*rhsV2p).*(aminf_g(vAll_pgc+.5*rhsV2p)-(pVr(:,3)+.5*piR2(:,3))); %I_A
    piR3(:,4)=dt*3.3./hatau_g(vAll_pgc+.5*rhsV2p).*(ahinf_g(vAll_pgc+.5*rhsV2p)-(pVr(:,4)+.5*piR2(:,4))); %I_A
    piR3(:,5)=dt*(minfMusc(vAll_pgc+.5*rhsV2p)-(pVr(:,5)+.5*piR2(:,5)))./tauMusc(vAll_pgc+.5*rhsV2p); %I_M
    piR3(:,6)=dt*2.1./hcurmtau(vAll_pgc+.5*rhsV2p).*(hcurminf(vAll_pgc+.5*rhsV2p)-(pVr(:,6)+.5*piR2(:,6)));
    piR3(:,7)=dt*(capminf(vAll_pgc+.5*rhsV2p)-(pVr(:,7)+.5*piR2(:,7)))./capmtau(vAll_pgc+.5*rhsV2p); %I_CaPN
    piR3(:,8)=dt*(caphinf(vAll_pgc+.5*rhsV2p)-(pVr(:,8)+.5*piR2(:,8)))./caphtau(vAll_pgc+.5*rhsV2p);
    piR3(:,9)=dt*(kcaa(vAll_pgc+.5*rhsV2p,caCon_p+.5*caCr2p).*(1-(pVr(:,9)+.5*piR2(:,9))) - 0.05*(pVr(:,9)+.5*piR2(:,9)));
    piR3(:,10)=dt*3.3*(mss_dr(vAll_pgc+.5*rhsV2p)-(pVr(:,10)+.5*piR2(:,10)))./taum_dr(vAll_pgc+.5*rhsV2p);    %I_DR
    piR3(:,11)=dt*(catminfP(vAll_pgc+.5*rhsV2p)-(pVr(:,11)+.5*piR2(:,11)))./cattaumP(vAll_pgc+.5*rhsV2p); %I_CaT
    piR3(:,12)=dt*(cathinf(vAll_pgc+.5*rhsV2p)-(pVr(:,12)+.5*piR2(:,12)))./cattauh(vAll_pgc+.5*rhsV2p); 
    piR3(:,13)=dt*(canminf(vAll_pgc+.5*rhsV2p)-(pVr(:,13)+.5*piR2(:,13)))./cantau(vAll_pgc+.5*rhsV2p); %I_CAN
    Ica=gCaPN_p.*(pVr(:,7)+.5*piR2(:,7)).^2.*(pVr(:,8)+.5*piR2(:,8)).*(vAll_pgc+.5*rhsV2p-Eca_p)+...
        gCaT_p.*(pVr(:,11)+.5*piR2(:,11)).^2.*(pVr(:,12)+.5*piR2(:,12)).*(vAll_pgc+.5*rhsV2p-Eca_p);
    caCr3p=dt*(-Ica./(2*F*thick_g) + (.05-(caCon_p+.5*caCr2p))/10 );
    rhsV3=(vAll_mc+.5*rhsV2-El)+gNa.*(xVr(:,1)+.5*xiR2(:,1)).^3.*(xVr(:,2)+.5*xiR2(:,2)).*(vAll_mc+.5*rhsV2-Ena)+...
        gNaP.*minf(vAll_mc+.5*rhsV2).*(vAll_mc+.5*rhsV2-Ena)+gA.*(xVr(:,3)+.5*xiR2(:,3)).*(xVr(:,4)+...
        .5*xiR2(:,4)).*(vAll_mc+.5*rhsV2-Ek)+gKS.*(xVr(:,5)+.5*xiR2(:,5)).*(xVr(:,6)+.5*xiR2(:,6)).*(vAll_mc+.5*rhsV2-Ek)+...
        gCaL.*(xVr(:,7)+.5*xiR2(:,7)).*(xVr(:,8)+.5*xiR2(:,8)).*(vAll_mc+.5*rhsV2-Eca)+gKCa.*(xVr(:,9)+...
        .5*xiR2(:,9)).*(vAll_mc+.5*rhsV2-Ek)+gDR.*(xVr(:,10)+.5*xiR2(:,10)).^2.*(xVr(:,11)+.5*xiR2(:,11)).*(vAll_mc+.5*rhsV2-Ek);
    rhsV3=dt/tau_m*(IdcM-rhsV3);
    xiR3(:,1)=dt*(ana(vAll_mc+.5*rhsV2).*(1-(xVr(:,1)+.5*xiR2(:,1))) - bna(vAll_mc+.5*rhsV2).*(xVr(:,1)+.5*xiR2(:,1)));
    xiR3(:,2)=dt*(ahna(vAll_mc+.5*rhsV2).*(1-(xVr(:,2)+.5*xiR2(:,2))) - bhna(vAll_mc+.5*rhsV2).*(xVr(:,2)+.5*xiR2(:,2)));
    xiR3(:,3)=dt*3.3./atau(vAll_mc+.5*rhsV2).*(aminf(vAll_mc+.5*rhsV2)-(xVr(:,3)+.5*xiR2(:,3))); %I_A
    xiR3(:,4)=dt*3.3./hatau(vAll_mc+.5*rhsV2).*(ahinf(vAll_mc+.5*rhsV2)-(xVr(:,4)+.5*xiR2(:,4))); %I_A
    xiR3(:,5)=dt/10*(ksminf(vAll_mc+.5*rhsV2)-(xVr(:,5)+.5*xiR2(:,5))); %I_KS
    xiR3(:,6)=dt./kshtau(vAll_mc+.5*rhsV2).*(kshinf(vAll_mc+.5*rhsV2)-(xVr(:,6)+.5*xiR2(:,6)));
    xiR3(:,7)=dt*(cala(vAll_mc+.5*rhsV2).*(1-(xVr(:,7)+.5*xiR2(:,7))) - calb(vAll_mc+.5*rhsV2).*(xVr(:,7)+.5*xiR2(:,7))); %I_CaL
    xiR3(:,8)=dt*(calha(vAll_mc+.5*rhsV2).*(1-(xVr(:,8)+.5*xiR2(:,8))) - calhb(vAll_mc+.5*rhsV2).*(xVr(:,8)+.5*xiR2(:,8)));
    xiR3(:,9)=dt*(kcaa(vAll_mc+.5*rhsV2,caCon+.5*caCr2).*(1-(xVr(:,9)+.5*xiR2(:,9))) - 0.05*(xVr(:,9)+.5*xiR2(:,9)));
    xiR3(:,10)=dt*(nss_dr(vAll_mc+.5*rhsV2)-(xVr(:,10)+.5*xiR2(:,10)))./taun_dr(vAll_mc+.5*rhsV2);    %I_DR
    xiR3(:,11)=dt*(kss_dr(vAll_mc+.5*rhsV2)-(xVr(:,11)+.5*xiR2(:,11)))./50;
    Ica=gCaL.*(xVr(:,7)+.5*xiR2(:,7)).*(xVr(:,8)+.5*xiR2(:,8)).*(vAll_mc+.5*rhsV2-Eca);
    caCr3=dt*(-Ica./(2*F*thick) + (.05-(caCon+.5*caCr2))/10 );
    % --- step 4 ---
    rhsV4p=(vAll_pgc+rhsV3p-El_p)+gNa_p.*(pVr(:,1)+piR3(:,1)).^3.*(pVr(:,2)+piR3(:,2)).*(vAll_pgc+rhsV3p-Ena)+...
        gA_p.*(pVr(:,3)+piR3(:,3)).*(pVr(:,4)+piR3(:,4)).*(vAll_pgc+rhsV3p-Ek)+...
        gM_p.*(pVr(:,5)+piR3(:,5)).*(vAll_pgc+rhsV3p-Ek)+gH.*(pVr(:,6)+piR3(:,6)).*(vAll_pgc+rhsV3p-Eh)+...
        gCaPN_p.*(pVr(:,7)+piR3(:,7)).^2.*(pVr(:,8)+piR3(:,8)).*(vAll_pgc+rhsV3p-Eca_p)+gKCa_p.*(pVr(:,9)+...
        piR3(:,9)).*(vAll_pgc+rhsV3p-Ek)+gDR_p.*(pVr(:,10)+piR3(:,10)).*(vAll_pgc+rhsV3p-Ek)+...
        gCaT_p.*(pVr(:,11)+piR3(:,11)).^2.*(pVr(:,12)+piR3(:,12)).*(vAll_pgc+rhsV3p-Eca_p)+...
        gCAN_p.*(caCon_p+caCr3p)./(200+(caCon_p+caCr3p)).*(pVr(:,13)+piR3(:,13)).*(vAll_pgc+rhsV3p-Ecan);
    rhsV4p=dt/taup_m*(IdcP-rhsV4p);
    piR4(:,1)=dt*2.1*(ana_g(vAll_pgc+rhsV3p)./(ana_g(vAll_pgc+rhsV3p)+bna_g(vAll_pgc+rhsV3p)) -(pVr(:,1)+piR3(:,1)) )./tauM(vAll_pgc+rhsV3p);
    piR4(:,2)=dt*2.1*(hInf(vAll_pgc+rhsV3p)-(pVr(:,2)+piR3(:,2)))./tauH(vAll_pgc+rhsV3p);
    piR4(:,3)=dt*3.3./atau(vAll_pgc+rhsV3p).*(aminf_g(vAll_pgc+rhsV3p)-(pVr(:,3)+piR3(:,3))); %I_A
    piR4(:,4)=dt*3.3./hatau_g(vAll_pgc+rhsV3p).*(ahinf_g(vAll_pgc+rhsV3p)-(pVr(:,4)+piR3(:,4))); %I_A
    piR4(:,5)=dt*(minfMusc(vAll_pgc+rhsV3p)-(pVr(:,5)+piR3(:,5)))./tauMusc(vAll_pgc+rhsV3p); %I_M
    piR4(:,6)=dt*2.1./hcurmtau(vAll_pgc+rhsV3p).*(hcurminf(vAll_pgc+rhsV3p)-(pVr(:,6)+piR3(:,6)));
    piR4(:,7)=dt*(capminf(vAll_pgc+rhsV3p)-(pVr(:,7)+piR3(:,7)))./capmtau(vAll_pgc+rhsV3p); %I_CaPN
    piR4(:,8)=dt*(caphinf(vAll_pgc+rhsV3p)-(pVr(:,8)+piR3(:,8)))./caphtau(vAll_pgc+rhsV3p);
    piR4(:,9)=dt*(kcaa(vAll_pgc+rhsV3p,caCon_p+caCr3p).*(1-(pVr(:,9)+piR3(:,9))) - 0.05*(pVr(:,9)+piR3(:,9)));
    piR4(:,10)=dt*3.3*(mss_dr(vAll_pgc+rhsV3p)-(pVr(:,10)+piR3(:,10)))./taum_dr(vAll_pgc+rhsV3p);    %I_DR
    piR4(:,11)=dt*(catminfP(vAll_pgc+rhsV3p)-(pVr(:,11)+piR3(:,11)))./cattaumP(vAll_pgc+rhsV3p); %I_CaT
    piR4(:,12)=dt*(cathinf(vAll_pgc+rhsV3p)-(pVr(:,12)+piR3(:,12)))./cattauh(vAll_pgc+rhsV3p); 
    piR4(:,13)=dt*(canminf(vAll_pgc+rhsV3p)-(pVr(:,13)+piR3(:,13)))./cantau(vAll_pgc+rhsV3p); %I_CAN
    Ica=gCaPN_p.*(pVr(:,7)+piR3(:,7)).^2.*(pVr(:,8)+piR3(:,8)).*(vAll_pgc+rhsV3p-Eca_p)+...
        gCaT_p.*(pVr(:,11)+piR3(:,11)).^2.*(pVr(:,12)+piR3(:,12)).*(vAll_pgc+rhsV3p-Eca_p);
    caCr4p=dt*(-Ica./(2*F*thick_g) + (.05-(caCon_p+caCr3p))/10 );
    rhsV4=(vAll_mc+rhsV3-El)+gNa.*(xVr(:,1)+xiR3(:,1)).^3.*(xVr(:,2)+xiR3(:,2)).*(vAll_mc+rhsV3-Ena)+...
        gNaP.*minf(vAll_mc+rhsV3).*(vAll_mc+rhsV3-Ena)+gA.*(xVr(:,3)+xiR3(:,3)).*(xVr(:,4)+...
        xiR3(:,4)).*(vAll_mc+rhsV3-Ek)+gKS.*(xVr(:,5)+xiR3(:,5)).*(xVr(:,6)+xiR3(:,6)).*(vAll_mc+rhsV3-Ek)+...
        gCaL.*(xVr(:,7)+xiR3(:,7)).*(xVr(:,8)+xiR3(:,8)).*(vAll_mc+rhsV3-Eca)+gKCa.*(xVr(:,9)+...
        xiR3(:,9)).*(vAll_mc+rhsV3-Ek)+gDR.*(xVr(:,10)+xiR3(:,10)).^2.*(xVr(:,11)+xiR3(:,11)).*(vAll_mc+rhsV3-Ek);
    rhsV4=dt/tau_m*(IdcM-rhsV4);
    xiR4(:,1)=dt*(ana(vAll_mc+rhsV3).*(1-(xVr(:,1)+xiR3(:,1))) - bna(vAll_mc+rhsV3).*(xVr(:,1)+xiR3(:,1)));
    xiR4(:,2)=dt*(ahna(vAll_mc+rhsV3).*(1-(xVr(:,2)+xiR3(:,2))) - bhna(vAll_mc+rhsV3).*(xVr(:,2)+xiR3(:,2)));
    xiR4(:,3)=dt*3.3./atau(vAll_mc+rhsV3).*(aminf(vAll_mc+rhsV3)-(xVr(:,3)+xiR3(:,3))); %I_A
    xiR4(:,4)=dt*3.3./hatau(vAll_mc+rhsV3).*(ahinf(vAll_mc+rhsV3)-(xVr(:,4)+xiR3(:,4))); %I_A
    xiR4(:,5)=dt/10*(ksminf(vAll_mc+rhsV3)-(xVr(:,5)+xiR3(:,5))); %I_KS
    xiR4(:,6)=dt./kshtau(vAll_mc+rhsV3).*(kshinf(vAll_mc+rhsV3)-(xVr(:,6)+xiR3(:,6)));
    xiR4(:,7)=dt*(cala(vAll_mc+rhsV3).*(1-(xVr(:,7)+xiR3(:,7))) - calb(vAll_mc+rhsV3).*(xVr(:,7)+xiR3(:,7))); %I_CaL
    xiR4(:,8)=dt*(calha(vAll_mc+rhsV3).*(1-(xVr(:,8)+xiR3(:,8))) - calhb(vAll_mc+rhsV3).*(xVr(:,8)+xiR3(:,8)));
    xiR4(:,9)=dt*(kcaa(vAll_mc+rhsV3,caCon+caCr3).*(1-(xVr(:,9)+xiR3(:,9))) - 0.05*(xVr(:,9)+xiR3(:,9)));
    xiR4(:,10)=dt*(nss_dr(vAll_mc+rhsV3)-(xVr(:,10)+xiR3(:,10)))./taun_dr(vAll_mc+rhsV3);    %I_DR
    xiR4(:,11)=dt*(kss_dr(vAll_mc+rhsV3)-(xVr(:,11)+xiR3(:,11)))./50;
    %update calcium concentration
    Ica=gCaL.*(xVr(:,7)+xiR3(:,7)).*(xVr(:,8)+xiR3(:,8)).*(vAll_mc+rhsV3-Eca);
    caCr4=dt*(-Ica./(2*F*thick) + (.05-(caCon+caCr3))/10 );
    % FINAL Step... update voltage; add implicit trap (Crank-nicholson) diffusion
    bv_p(idGv_p)=1/6*(rhsV1p+2*rhsV2p+2*rhsV3p+rhsV4p);
    bv_p(sidM2P)=bv_p(sidM2P)+dt/C*(-(wM2P.*synM2P+wM2Pn.*synM2Pn+wP_bg*bgAmp_p).*(vAll_pgc(sidM2P)-Excs)); %synap coupling
    vAll_pgc = Al_p\( Ar_p*vAll_pgc + bv_p );
    bv(idGv)=1/6*(rhsV1+2*rhsV2+2*rhsV3+rhsV4);
    bv(sidP2M)=bv(sidP2M)+dt/C*(-wP2M.*synP2M.*(vAll_mc(sidP2M)-Einh)); %synap coupling
    bv(sidP2M)=bv(sidP2M)+dt/C*(-wM_bg*bgAmp_m.*(vAll_mc(sidP2M)-Excs));
    vAll_mc = Al\( Ar*vAll_mc + bv );
    %update gating variables with RK4
    pVr=pVr+1/6*(piR1+2*piR2+2*piR3+piR4);
    xVr=xVr+1/6*(xiR1+2*xiR2+2*xiR3+xiR4);
    %update calcium concentration
    caCon_p=caCon_p+1/6*(caCr1p+2*caCr2p+2*caCr3p+caCr4p);
    caCon=caCon+1/6*(caCr1+2*caCr2+2*caCr3+caCr4);
    %update calcium reversal potential, Nernst eqn
    Eca_p=RTovF/2*log(10./caCon_p); %get 70mV when caCon=0.05
    Eca=RTovF/2*log(10./caCon); 
    if(lenRec_p==1)
        if(vAll_pgc(idRec_p)>vltThres && (dt*(j-indLst_p )>minTspk))
            indLst_p=j;
        end
    else
        %check if spike in PGC
        indCurSpk=find((vAll_pgc(idRec_p)>vltThres).*(dt*(j-indLst_p )>minTspk)); 
        % won't evaluate for-loop if no spikes
        for jCurSpk=1:length(indCurSpk)
            indLst_p(indCurSpk(jCurSpk))=j; %update index of last spike (could have done this with a vector)..
        end
    end
    if(lenRec==1)
        if(vAll_mc(idRec)>vltThres && (dt*(j-indLst)>minTspk))
            indLst=j;
        end
    else
        indCurSpk=find((vAll_mc(idRec)>vltThres).*(dt*(j-indLst)>minTspk));  %in MC
        for jCurSpk=1:length(indCurSpk)
            indLst(indCurSpk(jCurSpk))=j; %update index of last spike (could have done this with a vector)..
        end
    end
end

for k=1:Nrlz
% --- reset all before new time loop ---
ppIn=(rand(round(Ncmp(1)/4),Lt) < dt*lam);
ppInM=(rand(Ncm(1),Lt) < dt*lam);
if(k==1)
    indLst_p=indLst_p-round(1000/dt); %index of time of last spike, sub prior trans time
    indLst=indLst-round(1000/dt);     %index of time of last spike
else
    indLst_p=indLst_p-Lt; %index of time of last spike, sub prior trans time
    indLst=indLst-Lt;     %index of time of last spike
end
%--- main time-loop --
for j=1:Lt
    bgAmp_p = bgAmp_p + dt/5*(-bgAmp_p) + bgJmp*ppIn(:,j);
    bgAmp_m=bgAmp_m + dt/5*(-bgAmp_m) + bgJmp*ppInM(:,j);
    %get presynaptic voltages
    tmpPreVp=interp1(idPre_M2P,vAll_pgc(sidM2P),idP_P2M,'nearest');
    tmpPreVmp=interp1(idPre_P2M,vAll_mc(sidP2M),idP_M2P,'nearest');
    %synapses
    synP2M=synP2M+dt*(1./(1.25*(1+exp(-((tmpPreVp+40)/2)))).*(1-synP2M)-1/18*synP2M);
    synM2P=synM2P+dt*(1./(1+exp(-(tmpPreVmp)/0.2)).*(1-synM2P)-1/5.5*synM2P);
    synM2Pn=synM2Pn+dt*(1/52./(1+exp(-(tmpPreVmp)/0.2)).*(1-synM2Pn)-1/343*synM2Pn);
   
    % --- step 1 ---
    rhsV1p=(vAll_pgc-El_p)+gNa_p.*pVr(:,1).^3.*pVr(:,2).*(vAll_pgc-Ena)+gA_p.*pVr(:,3).*pVr(:,4).*(vAll_pgc-Ek)+...
        gM_p.*pVr(:,5).*(vAll_pgc-Ek)+gH.*pVr(:,6).*(vAll_pgc-Eh)+...
        gCaPN_p.*pVr(:,7).^2.*pVr(:,8).*(vAll_pgc-Eca_p)+gKCa_p.*pVr(:,9).*(vAll_pgc-Ek)+...
        gDR_p.*pVr(:,10).*(vAll_pgc-Ek)+gCaT_p.*(pVr(:,11)).^2.*pVr(:,12).*(vAll_pgc-Eca_p)+...
        gCAN_p.*caCon_p./(200+caCon_p).*pVr(:,13).*(vAll_pgc-Ecan);
    rhsV1p=dt/taup_m*(IdcP-rhsV1p);
    piR1(:,1)=dt*2.1*(ana_g(vAll_pgc)./(ana_g(vAll_pgc)+bna_g(vAll_pgc))-pVr(:,1))./tauM(vAll_pgc);
    piR1(:,2)=dt*2.1*(hInf(vAll_pgc)-pVr(:,2))./tauH(vAll_pgc);
    piR1(:,3)=dt*3.3./atau(vAll_pgc).*(aminf_g(vAll_pgc)-pVr(:,3)); %I_A
    piR1(:,4)=dt*3.3./hatau_g(vAll_pgc).*(ahinf_g(vAll_pgc)-pVr(:,4)); %I_A
    piR1(:,5)=dt*(minfMusc(vAll_pgc)-pVr(:,5))./tauMusc(vAll_pgc); %I_M
    piR1(:,6)=dt*2.1./hcurmtau(vAll_pgc).*(hcurminf(vAll_pgc)-pVr(:,6));
    piR1(:,7)=dt*(capminf(vAll_pgc)-pVr(:,7))./capmtau(vAll_pgc); %I_CaPN
    piR1(:,8)=dt*(caphinf(vAll_pgc)-pVr(:,8))./caphtau(vAll_pgc);
    piR1(:,9)=dt*(kcaa(vAll_pgc,caCon_p).*(1-pVr(:,9)) - 0.05*pVr(:,9));
    piR1(:,10)=dt*3.3*(mss_dr(vAll_pgc)-pVr(:,10))./taum_dr(vAll_pgc);    %I_DR
    piR1(:,11)=dt*(catminfP(vAll_pgc)-pVr(:,11))./cattaumP(vAll_pgc); %I_CaT
    piR1(:,12)=dt*(cathinf(vAll_pgc)-pVr(:,12))./cattauh(vAll_pgc); 
    piR1(:,13)=dt*(canminf(vAll_pgc)-pVr(:,13))./cantau(vAll_pgc); %I_CAN
    %update calcium concentration; use both Ca-currents
    Ica=gCaPN_p.*pVr(:,7).^2.*pVr(:,8).*(vAll_pgc-Eca_p)+gCaT_p.*(pVr(:,11)).^2.*pVr(:,12).*(vAll_pgc-Eca_p); 
    caCr1p=dt*(-Ica./(2*F*thick_g) + (.05-caCon_p)/10 );
    
    rhsV1=(vAll_mc-El)+gNa.*xVr(:,1).^3.*xVr(:,2).*(vAll_mc-Ena)+gNaP.*minf(vAll_mc).*(vAll_mc-Ena)+...
        gA.*xVr(:,3).*xVr(:,4).*(vAll_mc-Ek)+gKS.*xVr(:,5).*xVr(:,6).*(vAll_mc-Ek)+...
        gCaL.*xVr(:,7).*xVr(:,8).*(vAll_mc-Eca)+gKCa.*xVr(:,9).*(vAll_mc-Ek)+gDR.*xVr(:,10).^2.*xVr(:,11).*(vAll_mc-Ek);
    rhsV1=dt/tau_m*(IdcM-rhsV1);
    xiR1(:,1)=dt*(ana(vAll_mc).*(1-xVr(:,1)) - bna(vAll_mc).*xVr(:,1));
    xiR1(:,2)=dt*(ahna(vAll_mc).*(1-xVr(:,2)) - bhna(vAll_mc).*xVr(:,2));
    xiR1(:,3)=dt*3.3./atau(vAll_mc).*(aminf(vAll_mc)-xVr(:,3)); %I_A
    xiR1(:,4)=dt*3.3./hatau(vAll_mc).*(ahinf(vAll_mc)-xVr(:,4)); %I_A
    xiR1(:,5)=dt/10*(ksminf(vAll_mc)-xVr(:,5)); %I_KS
    xiR1(:,6)=dt./kshtau(vAll_mc).*(kshinf(vAll_mc)-xVr(:,6));
    xiR1(:,7)=dt*(cala(vAll_mc).*(1-xVr(:,7)) - calb(vAll_mc).*xVr(:,7)); %I_CaL
    xiR1(:,8)=dt*(calha(vAll_mc).*(1-xVr(:,8)) - calhb(vAll_mc).*xVr(:,8));
    xiR1(:,9)=dt*(kcaa(vAll_mc,caCon).*(1-xVr(:,9)) - 0.05*xVr(:,9));
    xiR1(:,10)=dt*(nss_dr(vAll_mc)-xVr(:,10))./taun_dr(vAll_mc);    %I_DR
    xiR1(:,11)=dt*(kss_dr(vAll_mc)-xVr(:,11))./50;
    %update calcium concentration
    Ica=gCaL.*xVr(:,7).*xVr(:,8).*(vAll_mc-Eca);
    caCr1=dt*(-Ica./(2*F*thick) + (.05-caCon)/10 );
    
    % --- step 2 ---
    rhsV2p=(vAll_pgc+.5*rhsV1p-El_p)+gNa_p.*(pVr(:,1)+.5*piR1(:,1)).^3.*(pVr(:,2)+.5*piR1(:,2)).*(vAll_pgc+.5*rhsV1p-Ena)+...
        gA_p.*(pVr(:,3)+.5*piR1(:,3)).*(pVr(:,4)+.5*piR1(:,4)).*(vAll_pgc+.5*rhsV1p-Ek)+...
        gM_p.*(pVr(:,5)+.5*piR1(:,5)).*(vAll_pgc+.5*rhsV1p-Ek)+gH.*(pVr(:,6)+.5*piR1(:,6)).*(vAll_pgc+.5*rhsV1p-Eh)+...
        gCaPN_p.*(pVr(:,7)+.5*piR1(:,7)).^2.*(pVr(:,8)+.5*piR1(:,8)).*(vAll_pgc+.5*rhsV1p-Eca_p)+gKCa_p.*(pVr(:,9)+...
        .5*piR1(:,9)).*(vAll_pgc+.5*rhsV1p-Ek)+gDR_p.*(pVr(:,10)+.5*piR1(:,10)).*(vAll_pgc+.5*rhsV1p-Ek)+...
        gCaT_p.*(pVr(:,11)+.5*piR1(:,11)).^2.*(pVr(:,12)+.5*piR1(:,12)).*(vAll_pgc+.5*rhsV1p-Eca_p)+...
        gCAN_p.*(caCon_p+.5*caCr1p)./(200+(caCon_p+.5*caCr1p)).*(pVr(:,13)+.5*piR1(:,13)).*(vAll_pgc+.5*rhsV1p-Ecan);
    rhsV2p=dt/taup_m*(IdcP-rhsV2p);
    piR2(:,1)=dt*2.1*(ana_g(vAll_pgc+.5*rhsV1p)./(ana_g(vAll_pgc+.5*rhsV1p)+bna_g(vAll_pgc+.5*rhsV1p)) -(pVr(:,1)+.5*piR1(:,1)) )./tauM(vAll_pgc+.5*rhsV1p);
    piR2(:,2)=dt*2.1*(hInf(vAll_pgc+.5*rhsV1p)-(pVr(:,2)+.5*piR1(:,2)))./tauH(vAll_pgc+.5*rhsV1p);
    piR2(:,3)=dt*3.3./atau(vAll_pgc+.5*rhsV1p).*(aminf_g(vAll_pgc+.5*rhsV1p)-(pVr(:,3)+.5*piR1(:,3))); %I_A
    piR2(:,4)=dt*3.3./hatau_g(vAll_pgc+.5*rhsV1p).*(ahinf_g(vAll_pgc+.5*rhsV1p)-(pVr(:,4)+.5*piR1(:,4))); %I_A
    piR2(:,5)=dt*(minfMusc(vAll_pgc+.5*rhsV1p)-(pVr(:,5)+.5*piR1(:,5)))./tauMusc(vAll_pgc+.5*rhsV1p); %I_M
    piR2(:,6)=dt*2.1./hcurmtau(vAll_pgc+.5*rhsV1p).*(hcurminf(vAll_pgc+.5*rhsV1p)-(pVr(:,6)+.5*piR1(:,6)));
    piR2(:,7)=dt*(capminf(vAll_pgc+.5*rhsV1p)-(pVr(:,7)+.5*piR1(:,7)))./capmtau(vAll_pgc+.5*rhsV1p); %I_CaPN
    piR2(:,8)=dt*(caphinf(vAll_pgc+.5*rhsV1p)-(pVr(:,8)+.5*piR1(:,8)))./caphtau(vAll_pgc+.5*rhsV1p);
    piR2(:,9)=dt*(kcaa(vAll_pgc+.5*rhsV1p,caCon_p+.5*caCr1p).*(1-(pVr(:,9)+.5*piR1(:,9))) - 0.05*(pVr(:,9)+.5*piR1(:,9)));
    piR2(:,10)=dt*3.3*(mss_dr(vAll_pgc+.5*rhsV1p)-(pVr(:,10)+.5*piR1(:,10)))./taum_dr(vAll_pgc+.5*rhsV1p);    %I_DR
    piR2(:,11)=dt*(catminfP(vAll_pgc+.5*rhsV1p)-(pVr(:,11)+.5*piR1(:,11)))./cattaumP(vAll_pgc+.5*rhsV1p); %I_CaT
    piR2(:,12)=dt*(cathinf(vAll_pgc+.5*rhsV1p)-(pVr(:,12)+.5*piR1(:,12)))./cattauh(vAll_pgc+.5*rhsV1p); 
    piR2(:,13)=dt*(canminf(vAll_pgc+.5*rhsV1p)-(pVr(:,13)+.5*piR1(:,13)))./cantau(vAll_pgc+.5*rhsV1p); %I_CAN
    %update calcium concentration
    Ica=gCaPN_p.*(pVr(:,7)+.5*piR1(:,7)).^2.*(pVr(:,8)+.5*piR1(:,8)).*(vAll_pgc+.5*rhsV1p-Eca_p)+...
        gCaT_p.*(pVr(:,11)+.5*piR1(:,11)).^2.*(pVr(:,12)+.5*piR1(:,12)).*(vAll_pgc+.5*rhsV1p-Eca_p);
    caCr2p=dt*(-Ica./(2*F*thick_g) + (.05-(caCon_p+.5*caCr1p))/10 );
    
    rhsV2=(vAll_mc+.5*rhsV1-El)+gNa.*(xVr(:,1)+.5*xiR1(:,1)).^3.*(xVr(:,2)+.5*xiR1(:,2)).*(vAll_mc+.5*rhsV1-Ena)+...
        gNaP.*minf(vAll_mc+.5*rhsV1).*(vAll_mc+.5*rhsV1-Ena)+gA.*(xVr(:,3)+.5*xiR1(:,3)).*(xVr(:,4)+...
        .5*xiR1(:,4)).*(vAll_mc+.5*rhsV1-Ek)+gKS.*(xVr(:,5)+.5*xiR1(:,5)).*(xVr(:,6)+.5*xiR1(:,6)).*(vAll_mc+.5*rhsV1-Ek)+...
        gCaL.*(xVr(:,7)+.5*xiR1(:,7)).*(xVr(:,8)+.5*xiR1(:,8)).*(vAll_mc+.5*rhsV1-Eca)+gKCa.*(xVr(:,9)+...
        .5*xiR1(:,9)).*(vAll_mc+.5*rhsV1-Ek)+gDR.*(xVr(:,10)+.5*xiR1(:,10)).^2.*(xVr(:,11)+.5*xiR1(:,11)).*(vAll_mc+.5*rhsV1-Ek);
    rhsV2=dt/tau_m*(IdcM-rhsV2);
    xiR2(:,1)=dt*(ana(vAll_mc+.5*rhsV1).*(1-(xVr(:,1)+.5*xiR1(:,1))) - bna(vAll_mc+.5*rhsV1).*(xVr(:,1)+.5*xiR1(:,1)));
    xiR2(:,2)=dt*(ahna(vAll_mc+.5*rhsV1).*(1-(xVr(:,2)+.5*xiR1(:,2))) - bhna(vAll_mc+.5*rhsV1).*(xVr(:,2)+.5*xiR1(:,2)));
    xiR2(:,3)=dt*3.3./atau(vAll_mc+.5*rhsV1).*(aminf(vAll_mc+.5*rhsV1)-(xVr(:,3)+.5*xiR1(:,3))); %I_A
    xiR2(:,4)=dt*3.3./hatau(vAll_mc+.5*rhsV1).*(ahinf(vAll_mc+.5*rhsV1)-(xVr(:,4)+.5*xiR1(:,4))); %I_A
    xiR2(:,5)=dt/10*(ksminf(vAll_mc+.5*rhsV1)-(xVr(:,5)+.5*xiR1(:,5))); %I_KS
    xiR2(:,6)=dt./kshtau(vAll_mc+.5*rhsV1).*(kshinf(vAll_mc+.5*rhsV1)-(xVr(:,6)+.5*xiR1(:,6)));
    xiR2(:,7)=dt*(cala(vAll_mc+.5*rhsV1).*(1-(xVr(:,7)+.5*xiR1(:,7))) - calb(vAll_mc+.5*rhsV1).*(xVr(:,7)+.5*xiR1(:,7))); %I_CaL
    xiR2(:,8)=dt*(calha(vAll_mc+.5*rhsV1).*(1-(xVr(:,8)+.5*xiR1(:,8))) - calhb(vAll_mc+.5*rhsV1).*(xVr(:,8)+.5*xiR1(:,8)));
    xiR2(:,9)=dt*(kcaa(vAll_mc+.5*rhsV1,caCon+.5*caCr1).*(1-(xVr(:,9)+.5*xiR1(:,9))) - 0.05*(xVr(:,9)+.5*xiR1(:,9)));
    xiR2(:,10)=dt*(nss_dr(vAll_mc+.5*rhsV1)-(xVr(:,10)+.5*xiR1(:,10)))./taun_dr(vAll_mc+.5*rhsV1);    %I_DR
    xiR2(:,11)=dt*(kss_dr(vAll_mc+.5*rhsV1)-(xVr(:,11)+.5*xiR1(:,11)))./50;
    %update calcium concentration
    Ica=gCaL.*(xVr(:,7)+.5*xiR1(:,7)).*(xVr(:,8)+.5*xiR1(:,8)).*(vAll_mc+.5*rhsV1-Eca);
    caCr2=dt*(-Ica./(2*F*thick) + (.05-(caCon+.5*caCr1))/10 );
    
    % --- step 3 ---
    rhsV3p=(vAll_pgc+.5*rhsV2p-El_p)+gNa_p.*(pVr(:,1)+.5*piR2(:,1)).^3.*(pVr(:,2)+.5*piR2(:,2)).*(vAll_pgc+.5*rhsV2p-Ena)+...
        gA_p.*(pVr(:,3)+.5*piR2(:,3)).*(pVr(:,4)+.5*piR2(:,4)).*(vAll_pgc+.5*rhsV2p-Ek)+...
        gM_p.*(pVr(:,5)+.5*piR2(:,5)).*(vAll_pgc+.5*rhsV2p-Ek)+gH.*(pVr(:,6)+.5*piR2(:,6)).*(vAll_pgc+.5*rhsV2p-Eh)+...
        gCaPN_p.*(pVr(:,7)+.5*piR2(:,7)).^2.*(pVr(:,8)+.5*piR2(:,8)).*(vAll_pgc+.5*rhsV2p-Eca_p)+gKCa_p.*(pVr(:,9)+...
        .5*piR2(:,9)).*(vAll_pgc+.5*rhsV2p-Ek)+gDR_p.*(pVr(:,10)+.5*piR2(:,10)).*(vAll_pgc+.5*rhsV2p-Ek)+...
        gCaT_p.*(pVr(:,11)+.5*piR2(:,11)).^2.*(pVr(:,12)+.5*piR2(:,12)).*(vAll_pgc+.5*rhsV2p-Eca_p)+...
        gCAN_p.*(caCon_p+.5*caCr2p)./(200+(caCon_p+.5*caCr2p)).*(pVr(:,13)+.5*piR2(:,13)).*(vAll_pgc+.5*rhsV2p-Ecan);
    rhsV3p=dt/taup_m*(IdcP-rhsV3p);
    piR3(:,1)=dt*2.1*(ana_g(vAll_pgc+.5*rhsV2p)./(ana_g(vAll_pgc+.5*rhsV2p)+bna_g(vAll_pgc+.5*rhsV2p)) -(pVr(:,1)+.5*piR2(:,1)) )./tauM(vAll_pgc+.5*rhsV2p);
    piR3(:,2)=dt*2.1*(hInf(vAll_pgc+.5*rhsV2p)-(pVr(:,2)+.5*piR2(:,2)))./tauH(vAll_pgc+.5*rhsV2p);
    piR3(:,3)=dt*3.3./atau(vAll_pgc+.5*rhsV2p).*(aminf_g(vAll_pgc+.5*rhsV2p)-(pVr(:,3)+.5*piR2(:,3))); %I_A
    piR3(:,4)=dt*3.3./hatau_g(vAll_pgc+.5*rhsV2p).*(ahinf_g(vAll_pgc+.5*rhsV2p)-(pVr(:,4)+.5*piR2(:,4))); %I_A
    piR3(:,5)=dt*(minfMusc(vAll_pgc+.5*rhsV2p)-(pVr(:,5)+.5*piR2(:,5)))./tauMusc(vAll_pgc+.5*rhsV2p); %I_M
    piR3(:,6)=dt*2.1./hcurmtau(vAll_pgc+.5*rhsV2p).*(hcurminf(vAll_pgc+.5*rhsV2p)-(pVr(:,6)+.5*piR2(:,6)));
    piR3(:,7)=dt*(capminf(vAll_pgc+.5*rhsV2p)-(pVr(:,7)+.5*piR2(:,7)))./capmtau(vAll_pgc+.5*rhsV2p); %I_CaPN
    piR3(:,8)=dt*(caphinf(vAll_pgc+.5*rhsV2p)-(pVr(:,8)+.5*piR2(:,8)))./caphtau(vAll_pgc+.5*rhsV2p);
    piR3(:,9)=dt*(kcaa(vAll_pgc+.5*rhsV2p,caCon_p+.5*caCr2p).*(1-(pVr(:,9)+.5*piR2(:,9))) - 0.05*(pVr(:,9)+.5*piR2(:,9)));
    piR3(:,10)=dt*3.3*(mss_dr(vAll_pgc+.5*rhsV2p)-(pVr(:,10)+.5*piR2(:,10)))./taum_dr(vAll_pgc+.5*rhsV2p);    %I_DR
    piR3(:,11)=dt*(catminfP(vAll_pgc+.5*rhsV2p)-(pVr(:,11)+.5*piR2(:,11)))./cattaumP(vAll_pgc+.5*rhsV2p); %I_CaT
    piR3(:,12)=dt*(cathinf(vAll_pgc+.5*rhsV2p)-(pVr(:,12)+.5*piR2(:,12)))./cattauh(vAll_pgc+.5*rhsV2p); 
    piR3(:,13)=dt*(canminf(vAll_pgc+.5*rhsV2p)-(pVr(:,13)+.5*piR2(:,13)))./cantau(vAll_pgc+.5*rhsV2p); %I_CAN
    %update calcium concentration
    Ica=gCaPN_p.*(pVr(:,7)+.5*piR2(:,7)).^2.*(pVr(:,8)+.5*piR2(:,8)).*(vAll_pgc+.5*rhsV2p-Eca_p)+...
        gCaT_p.*(pVr(:,11)+.5*piR2(:,11)).^2.*(pVr(:,12)+.5*piR2(:,12)).*(vAll_pgc+.5*rhsV2p-Eca_p);
    caCr3p=dt*(-Ica./(2*F*thick_g) + (.05-(caCon_p+.5*caCr2p))/10 );
    
    rhsV3=(vAll_mc+.5*rhsV2-El)+gNa.*(xVr(:,1)+.5*xiR2(:,1)).^3.*(xVr(:,2)+.5*xiR2(:,2)).*(vAll_mc+.5*rhsV2-Ena)+...
        gNaP.*minf(vAll_mc+.5*rhsV2).*(vAll_mc+.5*rhsV2-Ena)+gA.*(xVr(:,3)+.5*xiR2(:,3)).*(xVr(:,4)+...
        .5*xiR2(:,4)).*(vAll_mc+.5*rhsV2-Ek)+gKS.*(xVr(:,5)+.5*xiR2(:,5)).*(xVr(:,6)+.5*xiR2(:,6)).*(vAll_mc+.5*rhsV2-Ek)+...
        gCaL.*(xVr(:,7)+.5*xiR2(:,7)).*(xVr(:,8)+.5*xiR2(:,8)).*(vAll_mc+.5*rhsV2-Eca)+gKCa.*(xVr(:,9)+...
        .5*xiR2(:,9)).*(vAll_mc+.5*rhsV2-Ek)+gDR.*(xVr(:,10)+.5*xiR2(:,10)).^2.*(xVr(:,11)+.5*xiR2(:,11)).*(vAll_mc+.5*rhsV2-Ek);
    rhsV3=dt/tau_m*(IdcM-rhsV3);
    xiR3(:,1)=dt*(ana(vAll_mc+.5*rhsV2).*(1-(xVr(:,1)+.5*xiR2(:,1))) - bna(vAll_mc+.5*rhsV2).*(xVr(:,1)+.5*xiR2(:,1)));
    xiR3(:,2)=dt*(ahna(vAll_mc+.5*rhsV2).*(1-(xVr(:,2)+.5*xiR2(:,2))) - bhna(vAll_mc+.5*rhsV2).*(xVr(:,2)+.5*xiR2(:,2)));
    xiR3(:,3)=dt*3.3./atau(vAll_mc+.5*rhsV2).*(aminf(vAll_mc+.5*rhsV2)-(xVr(:,3)+.5*xiR2(:,3))); %I_A
    xiR3(:,4)=dt*3.3./hatau(vAll_mc+.5*rhsV2).*(ahinf(vAll_mc+.5*rhsV2)-(xVr(:,4)+.5*xiR2(:,4))); %I_A
    xiR3(:,5)=dt/10*(ksminf(vAll_mc+.5*rhsV2)-(xVr(:,5)+.5*xiR2(:,5))); %I_KS
    xiR3(:,6)=dt./kshtau(vAll_mc+.5*rhsV2).*(kshinf(vAll_mc+.5*rhsV2)-(xVr(:,6)+.5*xiR2(:,6)));
    xiR3(:,7)=dt*(cala(vAll_mc+.5*rhsV2).*(1-(xVr(:,7)+.5*xiR2(:,7))) - calb(vAll_mc+.5*rhsV2).*(xVr(:,7)+.5*xiR2(:,7))); %I_CaL
    xiR3(:,8)=dt*(calha(vAll_mc+.5*rhsV2).*(1-(xVr(:,8)+.5*xiR2(:,8))) - calhb(vAll_mc+.5*rhsV2).*(xVr(:,8)+.5*xiR2(:,8)));
    xiR3(:,9)=dt*(kcaa(vAll_mc+.5*rhsV2,caCon+.5*caCr2).*(1-(xVr(:,9)+.5*xiR2(:,9))) - 0.05*(xVr(:,9)+.5*xiR2(:,9)));
    xiR3(:,10)=dt*(nss_dr(vAll_mc+.5*rhsV2)-(xVr(:,10)+.5*xiR2(:,10)))./taun_dr(vAll_mc+.5*rhsV2);    %I_DR
    xiR3(:,11)=dt*(kss_dr(vAll_mc+.5*rhsV2)-(xVr(:,11)+.5*xiR2(:,11)))./50;
    %update calcium concentration
    Ica=gCaL.*(xVr(:,7)+.5*xiR2(:,7)).*(xVr(:,8)+.5*xiR2(:,8)).*(vAll_mc+.5*rhsV2-Eca);
    caCr3=dt*(-Ica./(2*F*thick) + (.05-(caCon+.5*caCr2))/10 );
    
    % --- step 4 ---
    rhsV4p=(vAll_pgc+rhsV3p-El_p)+gNa_p.*(pVr(:,1)+piR3(:,1)).^3.*(pVr(:,2)+piR3(:,2)).*(vAll_pgc+rhsV3p-Ena)+...
        gA_p.*(pVr(:,3)+piR3(:,3)).*(pVr(:,4)+piR3(:,4)).*(vAll_pgc+rhsV3p-Ek)+...
        gM_p.*(pVr(:,5)+piR3(:,5)).*(vAll_pgc+rhsV3p-Ek)+gH.*(pVr(:,6)+piR3(:,6)).*(vAll_pgc+rhsV3p-Eh)+...
        gCaPN_p.*(pVr(:,7)+piR3(:,7)).^2.*(pVr(:,8)+piR3(:,8)).*(vAll_pgc+rhsV3p-Eca_p)+gKCa_p.*(pVr(:,9)+...
        piR3(:,9)).*(vAll_pgc+rhsV3p-Ek)+gDR_p.*(pVr(:,10)+piR3(:,10)).*(vAll_pgc+rhsV3p-Ek)+...
        gCaT_p.*(pVr(:,11)+piR3(:,11)).^2.*(pVr(:,12)+piR3(:,12)).*(vAll_pgc+rhsV3p-Eca_p)+...
        gCAN_p.*(caCon_p+caCr3p)./(200+(caCon_p+caCr3p)).*(pVr(:,13)+piR3(:,13)).*(vAll_pgc+rhsV3p-Ecan);
    rhsV4p=dt/taup_m*(IdcP-rhsV4p);
    piR4(:,1)=dt*2.1*(ana_g(vAll_pgc+rhsV3p)./(ana_g(vAll_pgc+rhsV3p)+bna_g(vAll_pgc+rhsV3p)) -(pVr(:,1)+piR3(:,1)) )./tauM(vAll_pgc+rhsV3p);
    piR4(:,2)=dt*2.1*(hInf(vAll_pgc+rhsV3p)-(pVr(:,2)+piR3(:,2)))./tauH(vAll_pgc+rhsV3p);
    piR4(:,3)=dt*3.3./atau(vAll_pgc+rhsV3p).*(aminf_g(vAll_pgc+rhsV3p)-(pVr(:,3)+piR3(:,3))); %I_A
    piR4(:,4)=dt*3.3./hatau_g(vAll_pgc+rhsV3p).*(ahinf_g(vAll_pgc+rhsV3p)-(pVr(:,4)+piR3(:,4))); %I_A
    piR4(:,5)=dt*(minfMusc(vAll_pgc+rhsV3p)-(pVr(:,5)+piR3(:,5)))./tauMusc(vAll_pgc+rhsV3p); %I_M
    piR4(:,6)=dt*2.1./hcurmtau(vAll_pgc+rhsV3p).*(hcurminf(vAll_pgc+rhsV3p)-(pVr(:,6)+piR3(:,6)));
    piR4(:,7)=dt*(capminf(vAll_pgc+rhsV3p)-(pVr(:,7)+piR3(:,7)))./capmtau(vAll_pgc+rhsV3p); %I_CaPN
    piR4(:,8)=dt*(caphinf(vAll_pgc+rhsV3p)-(pVr(:,8)+piR3(:,8)))./caphtau(vAll_pgc+rhsV3p);
    piR4(:,9)=dt*(kcaa(vAll_pgc+rhsV3p,caCon_p+caCr3p).*(1-(pVr(:,9)+piR3(:,9))) - 0.05*(pVr(:,9)+piR3(:,9)));
    piR4(:,10)=dt*3.3*(mss_dr(vAll_pgc+rhsV3p)-(pVr(:,10)+piR3(:,10)))./taum_dr(vAll_pgc+rhsV3p);    %I_DR
    piR4(:,11)=dt*(catminfP(vAll_pgc+rhsV3p)-(pVr(:,11)+piR3(:,11)))./cattaumP(vAll_pgc+rhsV3p); %I_CaT
    piR4(:,12)=dt*(cathinf(vAll_pgc+rhsV3p)-(pVr(:,12)+piR3(:,12)))./cattauh(vAll_pgc+rhsV3p); 
    piR4(:,13)=dt*(canminf(vAll_pgc+rhsV3p)-(pVr(:,13)+piR3(:,13)))./cantau(vAll_pgc+rhsV3p); %I_CAN
    %update calcium concentration
    Ica=gCaPN_p.*(pVr(:,7)+piR3(:,7)).^2.*(pVr(:,8)+piR3(:,8)).*(vAll_pgc+rhsV3p-Eca_p)+...
        gCaT_p.*(pVr(:,11)+piR3(:,11)).^2.*(pVr(:,12)+piR3(:,12)).*(vAll_pgc+rhsV3p-Eca_p);
    caCr4p=dt*(-Ica./(2*F*thick_g) + (.05-(caCon_p+caCr3p))/10 );
    
    rhsV4=(vAll_mc+rhsV3-El)+gNa.*(xVr(:,1)+xiR3(:,1)).^3.*(xVr(:,2)+xiR3(:,2)).*(vAll_mc+rhsV3-Ena)+...
        gNaP.*minf(vAll_mc+rhsV3).*(vAll_mc+rhsV3-Ena)+gA.*(xVr(:,3)+xiR3(:,3)).*(xVr(:,4)+...
        xiR3(:,4)).*(vAll_mc+rhsV3-Ek)+gKS.*(xVr(:,5)+xiR3(:,5)).*(xVr(:,6)+xiR3(:,6)).*(vAll_mc+rhsV3-Ek)+...
        gCaL.*(xVr(:,7)+xiR3(:,7)).*(xVr(:,8)+xiR3(:,8)).*(vAll_mc+rhsV3-Eca)+gKCa.*(xVr(:,9)+...
        xiR3(:,9)).*(vAll_mc+rhsV3-Ek)+gDR.*(xVr(:,10)+xiR3(:,10)).^2.*(xVr(:,11)+xiR3(:,11)).*(vAll_mc+rhsV3-Ek);
    rhsV4=dt/tau_m*(IdcM-rhsV4);
    xiR4(:,1)=dt*(ana(vAll_mc+rhsV3).*(1-(xVr(:,1)+xiR3(:,1))) - bna(vAll_mc+rhsV3).*(xVr(:,1)+xiR3(:,1)));
    xiR4(:,2)=dt*(ahna(vAll_mc+rhsV3).*(1-(xVr(:,2)+xiR3(:,2))) - bhna(vAll_mc+rhsV3).*(xVr(:,2)+xiR3(:,2)));
    xiR4(:,3)=dt*3.3./atau(vAll_mc+rhsV3).*(aminf(vAll_mc+rhsV3)-(xVr(:,3)+xiR3(:,3))); %I_A
    xiR4(:,4)=dt*3.3./hatau(vAll_mc+rhsV3).*(ahinf(vAll_mc+rhsV3)-(xVr(:,4)+xiR3(:,4))); %I_A
    xiR4(:,5)=dt/10*(ksminf(vAll_mc+rhsV3)-(xVr(:,5)+xiR3(:,5))); %I_KS
    xiR4(:,6)=dt./kshtau(vAll_mc+rhsV3).*(kshinf(vAll_mc+rhsV3)-(xVr(:,6)+xiR3(:,6)));
    xiR4(:,7)=dt*(cala(vAll_mc+rhsV3).*(1-(xVr(:,7)+xiR3(:,7))) - calb(vAll_mc+rhsV3).*(xVr(:,7)+xiR3(:,7))); %I_CaL
    xiR4(:,8)=dt*(calha(vAll_mc+rhsV3).*(1-(xVr(:,8)+xiR3(:,8))) - calhb(vAll_mc+rhsV3).*(xVr(:,8)+xiR3(:,8)));
    xiR4(:,9)=dt*(kcaa(vAll_mc+rhsV3,caCon+caCr3).*(1-(xVr(:,9)+xiR3(:,9))) - 0.05*(xVr(:,9)+xiR3(:,9)));
    xiR4(:,10)=dt*(nss_dr(vAll_mc+rhsV3)-(xVr(:,10)+xiR3(:,10)))./taun_dr(vAll_mc+rhsV3);    %I_DR
    xiR4(:,11)=dt*(kss_dr(vAll_mc+rhsV3)-(xVr(:,11)+xiR3(:,11)))./50;
    %update calcium concentration
    Ica=gCaL.*(xVr(:,7)+xiR3(:,7)).*(xVr(:,8)+xiR3(:,8)).*(vAll_mc+rhsV3-Eca);
    caCr4=dt*(-Ica./(2*F*thick) + (.05-(caCon+caCr3))/10 );
    
    % FINAL Step... update voltage; add implicit trap (Crank-nicholson) diffusion
    bv_p(idGv_p)=1/6*(rhsV1p+2*rhsV2p+2*rhsV3p+rhsV4p);
    bv_p(sidM2P)=bv_p(sidM2P)+dt/C*(-(wM2P.*synM2P+wM2Pn.*synM2Pn+wP_bg*bgAmp_p).*(vAll_pgc(sidM2P)-Excs)); %synap coupling
    vAll_pgc = Al_p\( Ar_p*vAll_pgc + bv_p );
    
    bv(idGv)=1/6*(rhsV1+2*rhsV2+2*rhsV3+rhsV4);
    bv(sidP2M)=bv(sidP2M)+dt/C*(-wP2M.*synP2M.*(vAll_mc(sidP2M)-Einh)); %synap coupling
    bv(sidP2M)=bv(sidP2M)+dt/C*(-wM_bg*bgAmp_m.*(vAll_mc(sidP2M)-Excs));   
    vAll_mc = Al\( Ar*vAll_mc + bv );
    
    %update gating variables with RK4
    pVr=pVr+1/6*(piR1+2*piR2+2*piR3+piR4);
    xVr=xVr+1/6*(xiR1+2*xiR2+2*xiR3+xiR4);
    
    %update calcium concentration
    %Ica=gCaL*xVr(:,7).*xVr(:,8).*(V(:,j)-Eca);
    caCon_p=caCon_p+1/6*(caCr1p+2*caCr2p+2*caCr3p+caCr4p);
    caCon=caCon+1/6*(caCr1+2*caCr2+2*caCr3+caCr4);
    %update calcium reversal potential, Nernst eqn
    Eca_p=RTovF/2*log(10./caCon_p); %get 70mV when caCon=0.05
    Eca=RTovF/2*log(10./caCon); 
    
    %diff depending on # recording sites (faster w/ 1?)
    if(lenRec_p==1) %PGC spike?
        if(vAll_pgc(idRec_p)>vltThres && (dt*(j-indLst_p )>minTspk))
            fr_PGC=fr_PGC+1; %counts
            indLst_p=j;
        end
    else
        %check if spike in PGC
        indCurSpk=find((vAll_pgc(idRec_p)>vltThres).*(dt*(j-indLst_p )>minTspk)); 
        % won't evaluate for-loop if no spikes
        for jCurSpk=1:length(indCurSpk)
            fr_PGC(indCurSpk(jCurSpk))=fr_PGC(indCurSpk(jCurSpk))+1; %counts
            indLst_p(indCurSpk(jCurSpk))=j; %update index of last spike (could have done this with a vector)..
        end
    end
    if(lenRec==1) %MC spike?
        if(vAll_mc(idRec)>vltThres && (dt*(j-indLst)>minTspk))
            fr_MC=fr_MC+1; %only cunt spikes
            indLst=j;
        end
    else
        indCurSpk=find((vAll_mc(idRec)>vltThres).*(dt*(j-indLst)>minTspk));  %in MC
        for jCurSpk=1:length(indCurSpk)
            fr_MC(indCurSpk(jCurSpk))=fr_MC(indCurSpk(jCurSpk))+1; %counts
            indLst(indCurSpk(jCurSpk))=j; %update index of last spike (could have done this with a vector)..
        end
    end
   %store average voltage
   v_mcAvg=v_mcAvg+vAll_mc(idRec)/(Lt*Nrlz); %running sum
   v_pgcAvg=v_pgcAvg+vAll_pgc(idRec_p)/(Lt*Nrlz); %running sum
   
end %time-loop j=1:Lt
end %realization loop k=1:Nrlz

%scale firing rates, in Hz 
fr_MC=fr_MC/(tEnd*Nrlz)*1000;
fr_PGC=fr_PGC/(tEnd*Nrlz)*1000;

%--- subfunctions for MC---
    %Ina m alph & bet
    function a_na=ana(v)
        a_na=.32*(v+45)./(1-exp(-(v+45)./4));
    end
    function b_na=bna(v)
        b_na=-.28*(v+18)./(1-exp((v+18)./5));
    end
    %Ina h alph & bet
    function a_hna=ahna(v)
        a_hna=.128./exp((v+41)./18);
    end
    function b_hna=bhna(v)
        b_hna=4./(1+exp(-(v+18)./5));
    end
    %InaP
    function m_inf=minf(v)
        m_inf=1./(1+exp(-(v+50)./5));
    end
    %I_A
    function a_minf=aminf(v)
       a_minf=1./(1+exp(-(v-17.5)./14));
    end
    function a_tau=atau(v)
        a_tau=25*exp((v+45)./13.3)./(1+exp((v+45)./10));
    end
    function a_hinf=ahinf(v)
        a_hinf=1./(1+exp((v+41.7)./6));
    end
    function ha_tau=hatau(v)
        ha_tau=55.5*exp((v+70)./5.1)./(1+exp((v+70)./5));
    end
    %I_KS (tau=10)
    function ks_minf=ksminf(v)
        ks_minf=1./(1+exp(-(v+34)./6.5));
    end
    function ks_hinf=kshinf(v)
        ks_hinf=1./(1+exp((v+68)./6.6));
    end
    function ks_htau=kshtau(v)
        ks_htau=200+330./(1+exp(-(v+71.6)./6.85));
    end
    %I_CaL
    function cal_a=cala(v)
        cal_a=7.5./(1+exp(-(v-13)./7));
    end
    function cal_b=calb(v)
        cal_b=1.65./(1+exp((v-14)./4));
    end
    function cal_ha=calha(v)
        cal_ha=.0068./(1+exp((v+30)./12));
    end
    function cal_hb=calhb(v)
        cal_hb=.06./(1+exp(-v./11));
    end
    %I_KCA, beta_m=.05
    function kca_a=kcaa(v,CA)
        kca_a=-500*exp((v-65)./27).*(.015-CA)./(1-exp(-(CA-.015)./.0013));
    end
    %I_DR fast delayed rectifier; parms from getKssParm.m, stored in parm_fstDR.mat
    function xx=kss_dr(v)
        xx=.43315*(1+tanh(-(v+13.925)./13.0215))+.1337;
    end

    function xx=nss_dr(v)
        xx=((v+100)/150).^8.5849./(0.5747^8.5849+((v+100)/150).^8.5849);
    end

    function xx=taun_dr(v)
        xx=1./(.27654./exp((v+29.9998)/66.3783)+2.89./(1+exp(-(v-19.0524)/12.8786)));
    end

% -- subfunctions for GC --
%Ina m alph & bet
    function a_na=ana_g(v)
        a_na=.4*(v+25)./(1-exp(-(v+25)./7.2));
    end
    function b_na=bna_g(v)
        b_na=-.124*(v+25)./(1-exp((v+25)./7.2));
    end
    function x=tauM(v)
        x=max(1./(ana_g(v)+bna_g(v)),0.02);
    end
    %Ina h alph & bet
    function a_hna=ahna_g(v)
        a_hna=.03*(v+40)./(1-exp(-(v+40)./1.5));
    end
    function b_hna=bhna_g(v)
        b_hna=-.01*(v+40)./(1-exp((v+40)./1.5));
    end
    function x=hInf(v)
       x=1./(1+exp((v+45)./4));
    end
    function x=tauH(v)
        x=max(1./(ahna_g(v)+bhna_g(v)),0.5);
    end
    %I_A
    function a_minf=aminf_g(v)
       a_minf=1./(1+exp(-(v-7.6)./14));
    end
    function a_hinf=ahinf_g(v)
        a_hinf=1./(1+exp((v+67.4)./6));
    end
    function ha_tau=hatau_g(v)
        ha_tau=138.8*exp((v+70)./5.1)./(1+exp((v+70)./5));
    end
    %I_M
    function x=minfMusc(v)
        x=1./(1+exp(-(v+35)./5));
    end
    function x=tauMusc(v)
        x=1000./(3.3*exp((v+35)./40)+exp(-(v+35)./20));
    end
    % H current
    function x=hcurminf(v)
        x=1./(1+exp((v+80)./10));
    end
    function x=hcurmtau(v)
        x=1176.5*exp((v+65)./23.5)./(1+exp((v+65)./11.8));
    end
    %I_CaPN
    function x=capminf(v)
        x=1./(1+exp(-(v+10)./4));
    end
    function x=capmtau(v)
        x=.4+.7./(exp(-(v+5)./15)+exp((v+5)./15));
    end
    function x=caphinf(v)
        x=1./(1+exp((v+25)./2));
    end
    function x=caphtau(v)
        x=300+100./(exp(-(v+40)./9.5)+exp((v+40)./9.5));
    end
    %I_DR fast delayed rectifier
    function xx=mss_dr(v)
        xx=1./(1+exp(-(v-21)./10));
    end
    function xx=taum_dr(v)
        xx=285.7*exp((v+50)./36.4)./(1+exp((v+50)./18.2));
    end
    %I_CaT
    function x=cathinf(v)
        x=1./(1+exp((v+70)./4)); 
    end
    function x=cattauh(v) 
        x=10+40./(exp(-(v+50)./15)+exp((v+50)./15));
    end
    %I_CAN
    function x=canminf(v)
        x=1./(1+exp(-(v+43)./5.2));
    end
    function x=cantau(v)
        x=1.6+2.7./(exp(-(v+55)./15)+exp((v+55)./15));
    end
    
% -- subfunctions for PGC (mostly same as GC) --
    %I_CaT for PGC
    function x=catminfP(v) %change for rev PGC!
        x=1./(1+exp(-(v+59)./5.5)); %(v+59) for PGC
    end
    function x=cattaumP(v) %change for rev PGC!
        x=1.5+3.5./(exp(-(v+45)./15)+exp((v+30)./15)); %(v+45) for PGC
    end

end %end of main function