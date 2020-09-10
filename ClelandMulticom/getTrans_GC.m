function [fr_MC,fr_GC,v_mcAvg,v_gcAvg]=getTrans_GC(tEnd,Ibg,wghts)
% [fr_MC,fr_GC,v_mcAvg,v_gcAvg]=getTrans_GC(tEnd,Ibg,wghts);
% only saving counts, no PGC cells; getting MC Transfer function coupled to GC
% tEnd=time in ms; Ibg=[ MC ; GC ];
% Ibg/ev in muAmps/cm^2 (smaller), dc current and random AMPA bg for all cells
% wghts=[wM2G;wM2Gn;wG2M]
% wghts=[100;100;100];
%n SUBnetwork GC/MC (Li & Cleland); Ibg:(MC)
% !!! t is in ms  !!! return spike times (spTimes) 
% meant for SHORT times to show spikes AND voltage (not 200,000ms or crazyish)
% using Crank-Nicholson/ImplicTrapezoidalRule and RK4 time-stepping

rng('shuffle'); %random seed

Nrlz=10; %number of realizations

numGC=1; %assuming only 1 GC cell

dt=0.05; %assuming equally spaced
Lt=round(tEnd/dt)+1;

F=9.64853329e4; %Coulombs per Mole
T=273.15+35; %Kelvin
RTovF=8314.472*T/F; %R*1000 so in mV
thick=1; %thickeness of membrane shell, 1micron for Mitral Cell
thick_g=0.2; %thickeness of membrane shell, 0.2micron for GC/PGC

% PARAMS 
C=1.2; %micro-F/cm^2 for MC & PGC
Cg=2; %micro-F/cm^2 for GC
Rm=30; %kilo-Ohm*cm^2 for MC & GC
gL=1/Rm; %leak conduct (mS/cm^2)
tau_m=C*Rm; %membrane time constant (ms)
diam=20; %microns
diam_g=8; %microns
Ra=70; %Ohm*cm
taug_m=Cg*Rm; %membrane time constant (ms)

%MC: Tuft(3=2+1), PDend(4=2+2), Soma(5=2+3), Dend(6=2+4)
%GC: Spine(7=6+1), Dend(8=6+2), Soma(9=6+3)
Len_m=[4;74;5;100]; 
Len_g=[1;30;2];
dx = .25; %1 / (#pts per micron)
Spcv=(0:dx:sum(Len_m)-dx)'; %true length (microns)

Ntt=length(Spcv); %total # compartments in MC
Nttg=length( (0:dx:sum(Len_g)-dx)' );
Ncm=Len_m./dx; %# compartments IF dx divides all entries of Len
Ncmg=Len_g./dx;
Ncmsum=cumsum(Ncm);
Ncmsumg=cumsum(Ncmg);
if(Ncmsum(end)~=Ntt || Ncmsumg(end)~=Nttg)
    disp('index mismatch!');
    return
end

%record all midpoints
idRec=round(Ncm./2); %mid-point
idRec(2:end)=cumsum(Ncm(1:end-1))+idRec(2:end); %true index
idRec_g=round(Ncmg./2); %mid-point
idRec_g(2:end)=cumsum(Ncmg(1:end-1))+idRec_g(2:end); %true index
%only record at soma for each cell
idRec=idRec(3);
idRec_g=idRec_g(3);
lenRec=length(idRec);
lenRec_g=length(idRec_g);
% specify where put gating vars
idGv=(1:Ntt)';
nmGv=length(idGv); %# gating variables
idSegTyp=zeros(nmGv,1); %identity of segment type: Tuft(1), PDend(2), Soma(3), Dend(4)
idSegTyp(1:Ncmsum(1))=1; %putting them everywhere
for j=2:4
    idSegTyp(Ncmsum(j-1)+1:Ncmsum(j))=j;
end
idGv_g=(1:Nttg)';
nmGv_g=length(idGv_g); %# gating variables
idSegTyp_g=zeros(nmGv_g,1); %identity of segment type: Spine(1), Dend(2), Soma(3)
idSegTyp_g(1:Ncmsumg(1))=1; %putting them everywhere
for j=2:3
    idSegTyp_g(Ncmsumg(j-1)+1:Ncmsumg(j))=j;
end 
%...specify injected current here..
IdcM=sparse(zeros(nmGv,1));
IdcM(1:Ncm(1))=Ibg(1);
IdcG=sparse(zeros(nmGv_g,1));
IdcG(1:Ncmg(1))=Ibg(2);

%maximal conduct
% -- for MC --
gNa_1=[20;20;50;25]./gL; %(mS/cm^2) dimLess, Tuft, PDend, Soma, Dend
gNaP_1=[0.1;0.1;0.2;0.15]./gL; 
gDR_1=[10;10;10;20]./gL; 
gA_1=[0;0;1;0]./gL; %only in soma
gKS_1=[18;18;20;20]./gL;
gCaL_1=[0.2;0.2;0.4;.05]./gL;
gKCa_1=[0;0;1;0]./gL; %only in soma
% -- for GC --
gNa_2=[20;20;50]./gL; %(mS/cm^2) dimLess, Spine, Dend, Soma
gA_2=[60;60;20]./gL;
gM_2=[0;0;.5]./gL; %only in soma
gDR_2=[5;5;20]./gL; 
gCaPN_g_2=[.2;.2;0]./gL; %only in dendrite
gCaT_g_2=[0.1;0.1;0]./gL;
gCAN_2=[1;1;0]./gL; %only in dendrite
gKCa_g_2=[.5;.5;0]./gL; %only in soma
%-- augment to lineup with idGv & segment identity (idSegTyp) ---
% -- for MC --
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
% -- repeat for GC --
gNa_g=zeros(nmGv_g,1);
gDR_g=zeros(nmGv_g,1);
gA_g=zeros(nmGv_g,1);
gM_g=zeros(nmGv_g,1);
gCaPN_g=zeros(nmGv_g,1);
gCaT_g=zeros(nmGv_g,1);
gCAN=zeros(nmGv_g,1);
gKCa_g=zeros(nmGv_g,1);
for j=1:3
    gNa_g(idSegTyp_g==j)=gNa_2(j);
    gDR_g(idSegTyp_g==j)=gDR_2(j);
    gA_g(idSegTyp_g==j)=gA_2(j);
    gM_g(idSegTyp_g==j)=gM_2(j);
    gCaPN_g(idSegTyp_g==j)=gCaPN_g_2(j);
    gCaT_g(idSegTyp_g==j)=gCaT_g_2(j);
    gCAN(idSegTyp_g==j)=gCAN_2(j);
    gKCa_g(idSegTyp_g==j)=gKCa_g_2(j);
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

lambd=500*.5*sqrt(diam*Rm*1000/Ra)/100; %space-length in microns
lambd_g=500*.5*sqrt(diam_g*Rm*1000/Ra)/100; %space-length in microns
%matrices for diffusion part
MD=-2*diag(ones(Ntt,1))+diag(ones(Ntt-1,1),1)+diag(ones(Ntt-1,1),-1);
MD(1,1)=-1; MD(Ntt,Ntt)=-1; %no-flux BC
Al = sparse( eye(Ntt) - .5*dt/tau_m*(lambd/dx)^2*MD );
Ar = sparse( eye(Ntt) + .5*dt/tau_m*(lambd/dx)^2*MD );
bv = zeros(Ntt,1);
MD=-2*diag(ones(Nttg,1))+diag(ones(Nttg-1,1),1)+diag(ones(Nttg-1,1),-1);
MD(1,1)=-1; MD(Nttg,Nttg)=-1; %no-flux BC
Al_g = sparse( eye(Nttg) - .5*dt/taug_m*(lambd_g/dx)^2*MD );
Ar_g = sparse( eye(Nttg) + .5*dt/taug_m*(lambd_g/dx)^2*MD );
bv_g = zeros(Nttg,1);
% -- output -- 
if(lenRec==1)
    fr_MC=0;
else
    fr_MC=zeros(lenRec,1);
end
if(lenRec_g==1)
    fr_GC=0;
else
    fr_GC=zeros(numGC,lenRec_g); 
end
v_mcAvg=0; %average MC voltage
v_gcAvg=0; %average PGC voltage
vAll_mc=-70*ones(Ntt,1); %all voltages
vAll_gc=-75*ones(Nttg*numGC,1); %all voltages

minTspk=3; %min time between spikes (ms)
indLst=zeros(lenRec,1); %index of time of last spike
indLst_g=zeros(lenRec_g,1); %index of time of last spike
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

caCon_g=0.05*ones(nmGv_g,1); %not in soma, micrMol/l?
Eca_g=zeros(nmGv_g,1); %dyn variable
gVr=zeros(nmGv_g,13); %gating variables 
giR1=zeros(nmGv_g,13); %aux for RK4
giR2=zeros(nmGv_g,13); %aux for RK4
giR3=zeros(nmGv_g,13); %aux for RK4
giR4=zeros(nmGv_g,13); %aux for RK4
caCr1g=zeros(nmGv_g,1);
caCr2g=zeros(nmGv_g,1);
caCr3g=zeros(nmGv_g,1);
caCr4g=zeros(nmGv_g,1);
rhsV1g=zeros(nmGv_g,1);
rhsV2g=zeros(nmGv_g,1);
rhsV3g=zeros(nmGv_g,1);
rhsV4g=zeros(nmGv_g,1);

%synapses
sidP2M=( idGv <= Ncm(1) ); %start of MC Compart (Tuft) NEED THIS for BACKGROUND AMPA
synG2M=0;  %#inputs to MC 
    %sidG2M=( idGv >= Ncmsum(3)+.98*Ncm(4) ); %end of MC Compart (lateral Dend); move down another GC in another Glom
    sidG2M=( idGv >= Ncmsum(3)+.9*Ncm(4) )&(idGv<=Ncmsum(3)+.92*Ncm(4)); %near the end (same # pts 9 as before)
    numG2M=sum(sidG2M);
synG2M=zeros(numG2M,1); %update size synG2M, # GC inputs to MC
synM2G=zeros(Ncmg(1),1); %#inputs to GC (AMPA)
synM2Gn=zeros(Ncmg(1),1); %#inputs to GC (NMDA)
    sidM2G=( idGv_g <= Ncmg(1) ); %start of GC (spine)
    numM2G=sum(sidM2G);
%weights for synapses
wM2G=wghts(1)/length(synM2G);
wM2Gn=wghts(2)/length(synM2G);
wG2M=wghts(3)/length(synG2M);
%background ampa synapses
bgAmp_m=zeros(Ncm(1),1); 
bgAmp_g=zeros(Ncmg(1),1);
wM_bg=2;
wG_bg=2;

% temp vars for presyn volt
tmpPreVg=zeros(numG2M,1);
tmpPreVmg=zeros(numM2G,1);
    
%only used for for Presv Interp 
    idPre_M2G=(1:numM2G)';
    idPre_G2M=idGv(sidG2M); %end of MC
%make interp biophys, stayin bounds of presyn
    idP_G2M=(1 : (numM2G-1)/(numG2M-1) : numM2G)';
    tmpId=find(sidG2M);
    idP_M2Gend=(tmpId(end)-numG2M+1 : (numG2M-1)/(numM2G-1) : tmpId(end))'; %volt at end of MC
%params for background noisy inputs
lam=1/25; %kHz
lamG=1/25; %kHz
bgJmp=1; %size of Background Jump

%run initial loop to get rid of transients
ppInM=(rand(Ncm(1),round(1000/dt)) < dt*lam);
ppInG=(rand(Ncmg(1),round(1000/dt)) < dt*lamG);
for j=1:round(1000/dt)
    bgAmp_m=bgAmp_m + dt/5*(-bgAmp_m) + bgJmp*ppInM(:,j);
    bgAmp_g=bgAmp_g + dt/5*(-bgAmp_g) + bgJmp*ppInG(:,j);
    %get presynaptic voltages
    tmpPreVg=interp1(idPre_M2G,vAll_gc(sidM2G),idP_G2M,'nearest');
    tmpPreVmg=interp1(idPre_G2M,vAll_mc(sidG2M),idP_M2Gend,'nearest'); %volt at end of MC
    %synapses
    synG2M=synG2M+dt*(1./(1.25*(1+exp(-((tmpPreVg+40)/2)))).*(1-synG2M)-1/18*synG2M);
    synM2G=synM2G+dt*(1./(1+exp(-(tmpPreVmg)/0.2)).*(1-synM2G)-1/5.5*synM2G);
    synM2Gn=synM2Gn+dt*(1/52./(1+exp(-(tmpPreVmg)/0.2)).*(1-synM2Gn)-1/343*synM2Gn);
    % --- step 1 ---
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
    rhsV1g=(vAll_gc-El)+gNa_g.*gVr(:,1).^3.*gVr(:,2).*(vAll_gc-Ena)+gA_g.*gVr(:,3).*gVr(:,4).*(vAll_gc-Ek)+...
        gM_g.*gVr(:,5).*(vAll_gc-Ek)+...
        gCaPN_g.*gVr(:,7).^2.*gVr(:,8).*(vAll_gc-Eca_g)+gKCa_g.*gVr(:,9).*(vAll_gc-Ek)+...
        gDR_g.*gVr(:,10).*(vAll_gc-Ek)+gCaT_g.*(gVr(:,11)).^2.*gVr(:,12).*(vAll_gc-Eca_g)+...
        gCAN.*caCon_g./(200+caCon_g).*gVr(:,13).*(vAll_gc-Ecan);
    rhsV1g=dt/taug_m*(IdcG-rhsV1g);
    giR1(:,1)=dt*2.1*(ana_g(vAll_gc)./(ana_g(vAll_gc)+bna_g(vAll_gc))-gVr(:,1))./tauM(vAll_gc);
    giR1(:,2)=dt*2.1*(hInf(vAll_gc)-gVr(:,2))./tauH(vAll_gc);
    giR1(:,3)=dt*3.3./atau(vAll_gc).*(aminf_g(vAll_gc)-gVr(:,3)); %I_A
    giR1(:,4)=dt*3.3./hatau_g(vAll_gc).*(ahinf_g(vAll_gc)-gVr(:,4)); %I_A
    giR1(:,5)=dt*(minfMusc(vAll_gc)-gVr(:,5))./tauMusc(vAll_gc); %I_M
    giR1(:,6)=dt*2.1./hcurmtau(vAll_gc).*(hcurminf(vAll_gc)-gVr(:,6));
    giR1(:,7)=dt*(capminf(vAll_gc)-gVr(:,7))./capmtau(vAll_gc); %I_CaPN
    giR1(:,8)=dt*(caphinf(vAll_gc)-gVr(:,8))./caphtau(vAll_gc);
    giR1(:,9)=dt*(kcaa(vAll_gc,caCon_g).*(1-gVr(:,9)) - 0.05*gVr(:,9));
    giR1(:,10)=dt*3.3*(mss_dr(vAll_gc)-gVr(:,10))./taum_dr(vAll_gc);    %I_DR
    giR1(:,11)=dt*(catminfG(vAll_gc)-gVr(:,11))./cattaumG(vAll_gc); %I_CaT
    giR1(:,12)=dt*(cathinf(vAll_gc)-gVr(:,12))./cattauh(vAll_gc); 
    giR1(:,13)=dt*(canminf(vAll_gc)-gVr(:,13))./cantau(vAll_gc); %I_CAN
    %update calcium concentration; use both Ca-currents
    Ica=gCaPN_g.*gVr(:,7).^2.*gVr(:,8).*(vAll_gc-Eca_g)+gCaT_g.*(gVr(:,11)).^2.*gVr(:,12).*(vAll_gc-Eca_g); 
    caCr1g=dt*(-Ica./(2*F*thick_g) + (.05-caCon_g)/10 );
    % --- step 2 ---
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
    rhsV2g=(vAll_gc+.5*rhsV1g-El)+gNa_g.*(gVr(:,1)+.5*giR1(:,1)).^3.*(gVr(:,2)+.5*giR1(:,2)).*(vAll_gc+.5*rhsV1g-Ena)+...
        gA_g.*(gVr(:,3)+.5*giR1(:,3)).*(gVr(:,4)+.5*giR1(:,4)).*(vAll_gc+.5*rhsV1g-Ek)+...
        gM_g.*(gVr(:,5)+.5*giR1(:,5)).*(vAll_gc+.5*rhsV1g-Ek)+...
        gCaPN_g.*(gVr(:,7)+.5*giR1(:,7)).^2.*(gVr(:,8)+.5*giR1(:,8)).*(vAll_gc+.5*rhsV1g-Eca_g)+gKCa_g.*(gVr(:,9)+...
        .5*giR1(:,9)).*(vAll_gc+.5*rhsV1g-Ek)+gDR_g.*(gVr(:,10)+.5*giR1(:,10)).*(vAll_gc+.5*rhsV1g-Ek)+...
        gCaT_g.*(gVr(:,11)+.5*giR1(:,11)).^2.*(gVr(:,12)+.5*giR1(:,12)).*(vAll_gc+.5*rhsV1g-Eca_g)+...
        gCAN.*(caCon_g+.5*caCr1g)./(200+(caCon_g+.5*caCr1g)).*(gVr(:,13)+.5*giR1(:,13)).*(vAll_gc+.5*rhsV1g-Ecan);
    rhsV2g=dt/taug_m*(IdcG-rhsV2g);
    giR2(:,1)=dt*2.1*(ana_g(vAll_gc+.5*rhsV1g)./(ana_g(vAll_gc+.5*rhsV1g)+bna_g(vAll_gc+.5*rhsV1g)) -(gVr(:,1)+.5*giR1(:,1)) )./tauM(vAll_gc+.5*rhsV1g);
    giR2(:,2)=dt*2.1*(hInf(vAll_gc+.5*rhsV1g)-(gVr(:,2)+.5*giR1(:,2)))./tauH(vAll_gc+.5*rhsV1g);
    giR2(:,3)=dt*3.3./atau(vAll_gc+.5*rhsV1g).*(aminf_g(vAll_gc+.5*rhsV1g)-(gVr(:,3)+.5*giR1(:,3))); %I_A
    giR2(:,4)=dt*3.3./hatau_g(vAll_gc+.5*rhsV1g).*(ahinf_g(vAll_gc+.5*rhsV1g)-(gVr(:,4)+.5*giR1(:,4))); %I_A
    giR2(:,5)=dt*(minfMusc(vAll_gc+.5*rhsV1g)-(gVr(:,5)+.5*giR1(:,5)))./tauMusc(vAll_gc+.5*rhsV1g); %I_M
    giR2(:,6)=dt*2.1./hcurmtau(vAll_gc+.5*rhsV1g).*(hcurminf(vAll_gc+.5*rhsV1g)-(gVr(:,6)+.5*giR1(:,6)));
    giR2(:,7)=dt*(capminf(vAll_gc+.5*rhsV1g)-(gVr(:,7)+.5*giR1(:,7)))./capmtau(vAll_gc+.5*rhsV1g); %I_CaPN
    giR2(:,8)=dt*(caphinf(vAll_gc+.5*rhsV1g)-(gVr(:,8)+.5*giR1(:,8)))./caphtau(vAll_gc+.5*rhsV1g);
    giR2(:,9)=dt*(kcaa(vAll_gc+.5*rhsV1g,caCon_g+.5*caCr1g).*(1-(gVr(:,9)+.5*giR1(:,9))) - 0.05*(gVr(:,9)+.5*giR1(:,9)));
    giR2(:,10)=dt*3.3*(mss_dr(vAll_gc+.5*rhsV1g)-(gVr(:,10)+.5*giR1(:,10)))./taum_dr(vAll_gc+.5*rhsV1g);    %I_DR
    giR2(:,11)=dt*(catminfG(vAll_gc+.5*rhsV1g)-(gVr(:,11)+.5*giR1(:,11)))./cattaumG(vAll_gc+.5*rhsV1g); %I_CaT
    giR2(:,12)=dt*(cathinf(vAll_gc+.5*rhsV1g)-(gVr(:,12)+.5*giR1(:,12)))./cattauh(vAll_gc+.5*rhsV1g); 
    giR2(:,13)=dt*(canminf(vAll_gc+.5*rhsV1g)-(gVr(:,13)+.5*giR1(:,13)))./cantau(vAll_gc+.5*rhsV1g); %I_CAN
    Ica=gCaPN_g.*(gVr(:,7)+.5*giR1(:,7)).^2.*(gVr(:,8)+.5*giR1(:,8)).*(vAll_gc+.5*rhsV1g-Eca_g)+...
        gCaT_g.*(gVr(:,11)+.5*giR1(:,11)).^2.*(gVr(:,12)+.5*giR1(:,12)).*(vAll_gc+.5*rhsV1g-Eca_g);
    caCr2g=dt*(-Ica./(2*F*thick_g) + (.05-(caCon_g+.5*caCr1g))/10 );
    % --- step 3 ---
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
    rhsV3g=(vAll_gc+.5*rhsV2g-El)+gNa_g.*(gVr(:,1)+.5*giR2(:,1)).^3.*(gVr(:,2)+.5*giR2(:,2)).*(vAll_gc+.5*rhsV2g-Ena)+...
        gA_g.*(gVr(:,3)+.5*giR2(:,3)).*(gVr(:,4)+.5*giR2(:,4)).*(vAll_gc+.5*rhsV2g-Ek)+...
        gM_g.*(gVr(:,5)+.5*giR2(:,5)).*(vAll_gc+.5*rhsV2g-Ek)+...
        gCaPN_g.*(gVr(:,7)+.5*giR2(:,7)).^2.*(gVr(:,8)+.5*giR2(:,8)).*(vAll_gc+.5*rhsV2g-Eca_g)+gKCa_g.*(gVr(:,9)+...
        .5*giR2(:,9)).*(vAll_gc+.5*rhsV2g-Ek)+gDR_g.*(gVr(:,10)+.5*giR2(:,10)).*(vAll_gc+.5*rhsV2g-Ek)+...
        gCaT_g.*(gVr(:,11)+.5*giR2(:,11)).^2.*(gVr(:,12)+.5*giR2(:,12)).*(vAll_gc+.5*rhsV2g-Eca_g)+...
        gCAN.*(caCon_g+.5*caCr2g)./(200+(caCon_g+.5*caCr2g)).*(gVr(:,13)+.5*giR2(:,13)).*(vAll_gc+.5*rhsV2g-Ecan);
    rhsV3g=dt/taug_m*(IdcG-rhsV3g);
    giR3(:,1)=dt*2.1*(ana_g(vAll_gc+.5*rhsV2g)./(ana_g(vAll_gc+.5*rhsV2g)+bna_g(vAll_gc+.5*rhsV2g)) -(gVr(:,1)+.5*giR2(:,1)) )./tauM(vAll_gc+.5*rhsV2g);
    giR3(:,2)=dt*2.1*(hInf(vAll_gc+.5*rhsV2g)-(gVr(:,2)+.5*giR2(:,2)))./tauH(vAll_gc+.5*rhsV2g);
    giR3(:,3)=dt*3.3./atau(vAll_gc+.5*rhsV2g).*(aminf_g(vAll_gc+.5*rhsV2g)-(gVr(:,3)+.5*giR2(:,3))); %I_A
    giR3(:,4)=dt*3.3./hatau_g(vAll_gc+.5*rhsV2g).*(ahinf_g(vAll_gc+.5*rhsV2g)-(gVr(:,4)+.5*giR2(:,4))); %I_A
    giR3(:,5)=dt*(minfMusc(vAll_gc+.5*rhsV2g)-(gVr(:,5)+.5*giR2(:,5)))./tauMusc(vAll_gc+.5*rhsV2g); %I_M
    giR3(:,6)=dt*2.1./hcurmtau(vAll_gc+.5*rhsV2g).*(hcurminf(vAll_gc+.5*rhsV2g)-(gVr(:,6)+.5*giR2(:,6)));
    giR3(:,7)=dt*(capminf(vAll_gc+.5*rhsV2g)-(gVr(:,7)+.5*giR2(:,7)))./capmtau(vAll_gc+.5*rhsV2g); %I_CaPN
    giR3(:,8)=dt*(caphinf(vAll_gc+.5*rhsV2g)-(gVr(:,8)+.5*giR2(:,8)))./caphtau(vAll_gc+.5*rhsV2g);
    giR3(:,9)=dt*(kcaa(vAll_gc+.5*rhsV2g,caCon_g+.5*caCr2g).*(1-(gVr(:,9)+.5*giR2(:,9))) - 0.05*(gVr(:,9)+.5*giR2(:,9)));
    giR3(:,10)=dt*3.3*(mss_dr(vAll_gc+.5*rhsV2g)-(gVr(:,10)+.5*giR2(:,10)))./taum_dr(vAll_gc+.5*rhsV2g);    %I_DR
    giR3(:,11)=dt*(catminfG(vAll_gc+.5*rhsV2g)-(gVr(:,11)+.5*giR2(:,11)))./cattaumG(vAll_gc+.5*rhsV2g); %I_CaT
    giR3(:,12)=dt*(cathinf(vAll_gc+.5*rhsV2g)-(gVr(:,12)+.5*giR2(:,12)))./cattauh(vAll_gc+.5*rhsV2g); 
    giR3(:,13)=dt*(canminf(vAll_gc+.5*rhsV2g)-(gVr(:,13)+.5*giR2(:,13)))./cantau(vAll_gc+.5*rhsV2g); %I_CAN
    Ica=gCaPN_g.*(gVr(:,7)+.5*giR2(:,7)).^2.*(gVr(:,8)+.5*giR2(:,8)).*(vAll_gc+.5*rhsV2g-Eca_g)+...
        gCaT_g.*(gVr(:,11)+.5*giR2(:,11)).^2.*(gVr(:,12)+.5*giR2(:,12)).*(vAll_gc+.5*rhsV2g-Eca_g);
    caCr3g=dt*(-Ica./(2*F*thick_g) + (.05-(caCon_g+.5*caCr2g))/10 );
    % --- step 4 ---
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
    rhsV4g=(vAll_gc+rhsV3g-El)+gNa_g.*(gVr(:,1)+giR3(:,1)).^3.*(gVr(:,2)+giR3(:,2)).*(vAll_gc+rhsV3g-Ena)+...
        gA_g.*(gVr(:,3)+giR3(:,3)).*(gVr(:,4)+giR3(:,4)).*(vAll_gc+rhsV3g-Ek)+...
        gM_g.*(gVr(:,5)+giR3(:,5)).*(vAll_gc+rhsV3g-Ek)+...
        gCaPN_g.*(gVr(:,7)+giR3(:,7)).^2.*(gVr(:,8)+giR3(:,8)).*(vAll_gc+rhsV3g-Eca_g)+gKCa_g.*(gVr(:,9)+...
        giR3(:,9)).*(vAll_gc+rhsV3g-Ek)+gDR_g.*(gVr(:,10)+giR3(:,10)).*(vAll_gc+rhsV3g-Ek)+...
        gCaT_g.*(gVr(:,11)+giR3(:,11)).^2.*(gVr(:,12)+giR3(:,12)).*(vAll_gc+rhsV3g-Eca_g)+...
        gCAN.*(caCon_g+caCr3g)./(200+(caCon_g+caCr3g)).*(gVr(:,13)+giR3(:,13)).*(vAll_gc+rhsV3g-Ecan);
    rhsV4g=dt/taug_m*(IdcG-rhsV4g);
    giR4(:,1)=dt*2.1*(ana_g(vAll_gc+rhsV3g)./(ana_g(vAll_gc+rhsV3g)+bna_g(vAll_gc+rhsV3g)) -(gVr(:,1)+giR3(:,1)) )./tauM(vAll_gc+rhsV3g);
    giR4(:,2)=dt*2.1*(hInf(vAll_gc+rhsV3g)-(gVr(:,2)+giR3(:,2)))./tauH(vAll_gc+rhsV3g);
    giR4(:,3)=dt*3.3./atau(vAll_gc+rhsV3g).*(aminf_g(vAll_gc+rhsV3g)-(gVr(:,3)+giR3(:,3))); %I_A
    giR4(:,4)=dt*3.3./hatau_g(vAll_gc+rhsV3g).*(ahinf_g(vAll_gc+rhsV3g)-(gVr(:,4)+giR3(:,4))); %I_A
    giR4(:,5)=dt*(minfMusc(vAll_gc+rhsV3g)-(gVr(:,5)+giR3(:,5)))./tauMusc(vAll_gc+rhsV3g); %I_M
    giR4(:,6)=dt*2.1./hcurmtau(vAll_gc+rhsV3g).*(hcurminf(vAll_gc+rhsV3g)-(gVr(:,6)+giR3(:,6)));
    giR4(:,7)=dt*(capminf(vAll_gc+rhsV3g)-(gVr(:,7)+giR3(:,7)))./capmtau(vAll_gc+rhsV3g); %I_CaPN
    giR4(:,8)=dt*(caphinf(vAll_gc+rhsV3g)-(gVr(:,8)+giR3(:,8)))./caphtau(vAll_gc+rhsV3g);
    giR4(:,9)=dt*(kcaa(vAll_gc+rhsV3g,caCon_g+caCr3g).*(1-(gVr(:,9)+giR3(:,9))) - 0.05*(gVr(:,9)+giR3(:,9)));
    giR4(:,10)=dt*3.3*(mss_dr(vAll_gc+rhsV3g)-(gVr(:,10)+giR3(:,10)))./taum_dr(vAll_gc+rhsV3g);    %I_DR
    giR4(:,11)=dt*(catminfG(vAll_gc+rhsV3g)-(gVr(:,11)+giR3(:,11)))./cattaumG(vAll_gc+rhsV3g); %I_CaT
    giR4(:,12)=dt*(cathinf(vAll_gc+rhsV3g)-(gVr(:,12)+giR3(:,12)))./cattauh(vAll_gc+rhsV3g); 
    giR4(:,13)=dt*(canminf(vAll_gc+rhsV3g)-(gVr(:,13)+giR3(:,13)))./cantau(vAll_gc+rhsV3g); %I_CAN
    Ica=gCaPN_g.*(gVr(:,7)+giR3(:,7)).^2.*(gVr(:,8)+giR3(:,8)).*(vAll_gc+rhsV3g-Eca_g)+...
        gCaT_g.*(gVr(:,11)+giR3(:,11)).^2.*(gVr(:,12)+giR3(:,12)).*(vAll_gc+rhsV3g-Eca_g);
    caCr4g=dt*(-Ica./(2*F*thick_g) + (.05-(caCon_g+caCr3g))/10 );
    % FINAL Step... update voltage; add implicit trap (Crank-nicholson) diffusion
    bv(idGv)=1/6*(rhsV1+2*rhsV2+2*rhsV3+rhsV4);
    bv(sidP2M)=bv(sidP2M)+dt/C*(-wM_bg*bgAmp_m.*(vAll_mc(sidP2M)-Excs));%synap coupling
    bv(sidG2M)=bv(sidG2M)+dt/C*(-wG2M.*synG2M.*(vAll_mc(sidG2M)-Einh));    
    vAll_mc = Al\( Ar*vAll_mc + bv );
    bv_g(idGv_g)=1/6*(rhsV1g+2*rhsV2g+2*rhsV3g+rhsV4g);
    bv_g(sidM2G)=bv_g(sidM2G)+dt/C*(-(wM2G.*synM2G+wM2Gn.*synM2Gn+wG_bg*bgAmp_g).*(vAll_gc(sidM2G)-Excs)); %synap coupling
    vAll_gc = Al_g\( Ar_g*vAll_gc + bv_g );
    %update gating variables with RK4
    xVr=xVr+1/6*(xiR1+2*xiR2+2*xiR3+xiR4);
    gVr=gVr+1/6*(giR1+2*giR2+2*giR3+giR4);
    %update calcium concentration
    caCon=caCon+1/6*(caCr1+2*caCr2+2*caCr3+caCr4);
    caCon_g=caCon_g+1/6*(caCr1g+2*caCr2g+2*caCr3g+caCr4g);
    %update calcium reversal potential, Nernst eqn
    Eca=RTovF/2*log(10./caCon); 
    Eca_g=RTovF/2*log(10./caCon_g);
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
    if(lenRec_g==1 && numGC==1)
        if(vAll_gc(idRec_g)>vltThres && (dt*(j-indLst_g)>minTspk))
            indLst_g=j;
        end
    else
        indCurSpk=find((vAll_gc(idRec_g)>vltThres).*(dt*(j-indLst_g)>minTspk));  %in GC
        for jCurSpk=1:length(indCurSpk)
            indLst_g(indCurSpk(jCurSpk))=j; %update index of last spike (could have done this with a vector)..
        end
    end
end

for k=1:Nrlz
% --- reset all before new time loop ---
ppInM=(rand(Ncm(1),Lt) < dt*lam);
ppInG=(rand(Ncmg(1),Lt) < dt*lamG);
if(k==1)
    indLst=indLst-round(1000/dt);     %index of time of last spike
    indLst_g=indLst_g-round(1000/dt); %index of time of last spike
else
    indLst=indLst-Lt;     %index of time of last spike
    indLst_g=indLst_g-Lt; %index of time of last spike
end
%--- main time-loop --
for j=1:Lt
    bgAmp_m=bgAmp_m + dt/5*(-bgAmp_m) + bgJmp*ppInM(:,j);
    bgAmp_g=bgAmp_g + dt/5*(-bgAmp_g) + bgJmp*ppInG(:,j);
    
    %get presynaptic voltages
    tmpPreVg=interp1(idPre_M2G,vAll_gc(sidM2G),idP_G2M,'nearest');
    tmpPreVmg=interp1(idPre_G2M,vAll_mc(sidG2M),idP_M2Gend,'nearest'); %volt at end of MC
    %synapses
    synG2M=synG2M+dt*(1./(1.25*(1+exp(-((tmpPreVg+40)/2)))).*(1-synG2M)-1/18*synG2M);
    synM2G=synM2G+dt*(1./(1+exp(-(tmpPreVmg)/0.2)).*(1-synM2G)-1/5.5*synM2G);
    synM2Gn=synM2Gn+dt*(1/52./(1+exp(-(tmpPreVmg)/0.2)).*(1-synM2Gn)-1/343*synM2Gn);
    
    % --- step 1 ---
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
    
    rhsV1g=(vAll_gc-El)+gNa_g.*gVr(:,1).^3.*gVr(:,2).*(vAll_gc-Ena)+gA_g.*gVr(:,3).*gVr(:,4).*(vAll_gc-Ek)+...
        gM_g.*gVr(:,5).*(vAll_gc-Ek)+...
        gCaPN_g.*gVr(:,7).^2.*gVr(:,8).*(vAll_gc-Eca_g)+gKCa_g.*gVr(:,9).*(vAll_gc-Ek)+...
        gDR_g.*gVr(:,10).*(vAll_gc-Ek)+gCaT_g.*(gVr(:,11)).^2.*gVr(:,12).*(vAll_gc-Eca_g)+...
        gCAN.*caCon_g./(200+caCon_g).*gVr(:,13).*(vAll_gc-Ecan);
    rhsV1g=dt/taug_m*(IdcG-rhsV1g);
    giR1(:,1)=dt*2.1*(ana_g(vAll_gc)./(ana_g(vAll_gc)+bna_g(vAll_gc))-gVr(:,1))./tauM(vAll_gc);
    giR1(:,2)=dt*2.1*(hInf(vAll_gc)-gVr(:,2))./tauH(vAll_gc);
    giR1(:,3)=dt*3.3./atau(vAll_gc).*(aminf_g(vAll_gc)-gVr(:,3)); %I_A
    giR1(:,4)=dt*3.3./hatau_g(vAll_gc).*(ahinf_g(vAll_gc)-gVr(:,4)); %I_A
    giR1(:,5)=dt*(minfMusc(vAll_gc)-gVr(:,5))./tauMusc(vAll_gc); %I_M
    giR1(:,6)=dt*2.1./hcurmtau(vAll_gc).*(hcurminf(vAll_gc)-gVr(:,6));
    giR1(:,7)=dt*(capminf(vAll_gc)-gVr(:,7))./capmtau(vAll_gc); %I_CaPN
    giR1(:,8)=dt*(caphinf(vAll_gc)-gVr(:,8))./caphtau(vAll_gc);
    giR1(:,9)=dt*(kcaa(vAll_gc,caCon_g).*(1-gVr(:,9)) - 0.05*gVr(:,9));
    giR1(:,10)=dt*3.3*(mss_dr(vAll_gc)-gVr(:,10))./taum_dr(vAll_gc);    %I_DR
    giR1(:,11)=dt*(catminfG(vAll_gc)-gVr(:,11))./cattaumG(vAll_gc); %I_CaT
    giR1(:,12)=dt*(cathinf(vAll_gc)-gVr(:,12))./cattauh(vAll_gc); 
    giR1(:,13)=dt*(canminf(vAll_gc)-gVr(:,13))./cantau(vAll_gc); %I_CAN
    %update calcium concentration; use both Ca-currents
    Ica=gCaPN_g.*gVr(:,7).^2.*gVr(:,8).*(vAll_gc-Eca_g)+gCaT_g.*(gVr(:,11)).^2.*gVr(:,12).*(vAll_gc-Eca_g); 
    caCr1g=dt*(-Ica./(2*F*thick_g) + (.05-caCon_g)/10 );
    
    % --- step 2 ---
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
    
    rhsV2g=(vAll_gc+.5*rhsV1g-El)+gNa_g.*(gVr(:,1)+.5*giR1(:,1)).^3.*(gVr(:,2)+.5*giR1(:,2)).*(vAll_gc+.5*rhsV1g-Ena)+...
        gA_g.*(gVr(:,3)+.5*giR1(:,3)).*(gVr(:,4)+.5*giR1(:,4)).*(vAll_gc+.5*rhsV1g-Ek)+...
        gM_g.*(gVr(:,5)+.5*giR1(:,5)).*(vAll_gc+.5*rhsV1g-Ek)+...
        gCaPN_g.*(gVr(:,7)+.5*giR1(:,7)).^2.*(gVr(:,8)+.5*giR1(:,8)).*(vAll_gc+.5*rhsV1g-Eca_g)+gKCa_g.*(gVr(:,9)+...
        .5*giR1(:,9)).*(vAll_gc+.5*rhsV1g-Ek)+gDR_g.*(gVr(:,10)+.5*giR1(:,10)).*(vAll_gc+.5*rhsV1g-Ek)+...
        gCaT_g.*(gVr(:,11)+.5*giR1(:,11)).^2.*(gVr(:,12)+.5*giR1(:,12)).*(vAll_gc+.5*rhsV1g-Eca_g)+...
        gCAN.*(caCon_g+.5*caCr1g)./(200+(caCon_g+.5*caCr1g)).*(gVr(:,13)+.5*giR1(:,13)).*(vAll_gc+.5*rhsV1g-Ecan);
    rhsV2g=dt/taug_m*(IdcG-rhsV2g);
    giR2(:,1)=dt*2.1*(ana_g(vAll_gc+.5*rhsV1g)./(ana_g(vAll_gc+.5*rhsV1g)+bna_g(vAll_gc+.5*rhsV1g)) -(gVr(:,1)+.5*giR1(:,1)) )./tauM(vAll_gc+.5*rhsV1g);
    giR2(:,2)=dt*2.1*(hInf(vAll_gc+.5*rhsV1g)-(gVr(:,2)+.5*giR1(:,2)))./tauH(vAll_gc+.5*rhsV1g);
    giR2(:,3)=dt*3.3./atau(vAll_gc+.5*rhsV1g).*(aminf_g(vAll_gc+.5*rhsV1g)-(gVr(:,3)+.5*giR1(:,3))); %I_A
    giR2(:,4)=dt*3.3./hatau_g(vAll_gc+.5*rhsV1g).*(ahinf_g(vAll_gc+.5*rhsV1g)-(gVr(:,4)+.5*giR1(:,4))); %I_A
    giR2(:,5)=dt*(minfMusc(vAll_gc+.5*rhsV1g)-(gVr(:,5)+.5*giR1(:,5)))./tauMusc(vAll_gc+.5*rhsV1g); %I_M
    giR2(:,6)=dt*2.1./hcurmtau(vAll_gc+.5*rhsV1g).*(hcurminf(vAll_gc+.5*rhsV1g)-(gVr(:,6)+.5*giR1(:,6)));
    giR2(:,7)=dt*(capminf(vAll_gc+.5*rhsV1g)-(gVr(:,7)+.5*giR1(:,7)))./capmtau(vAll_gc+.5*rhsV1g); %I_CaPN
    giR2(:,8)=dt*(caphinf(vAll_gc+.5*rhsV1g)-(gVr(:,8)+.5*giR1(:,8)))./caphtau(vAll_gc+.5*rhsV1g);
    giR2(:,9)=dt*(kcaa(vAll_gc+.5*rhsV1g,caCon_g+.5*caCr1g).*(1-(gVr(:,9)+.5*giR1(:,9))) - 0.05*(gVr(:,9)+.5*giR1(:,9)));
    giR2(:,10)=dt*3.3*(mss_dr(vAll_gc+.5*rhsV1g)-(gVr(:,10)+.5*giR1(:,10)))./taum_dr(vAll_gc+.5*rhsV1g);    %I_DR
    giR2(:,11)=dt*(catminfG(vAll_gc+.5*rhsV1g)-(gVr(:,11)+.5*giR1(:,11)))./cattaumG(vAll_gc+.5*rhsV1g); %I_CaT
    giR2(:,12)=dt*(cathinf(vAll_gc+.5*rhsV1g)-(gVr(:,12)+.5*giR1(:,12)))./cattauh(vAll_gc+.5*rhsV1g); 
    giR2(:,13)=dt*(canminf(vAll_gc+.5*rhsV1g)-(gVr(:,13)+.5*giR1(:,13)))./cantau(vAll_gc+.5*rhsV1g); %I_CAN
    %update calcium concentration
    Ica=gCaPN_g.*(gVr(:,7)+.5*giR1(:,7)).^2.*(gVr(:,8)+.5*giR1(:,8)).*(vAll_gc+.5*rhsV1g-Eca_g)+...
        gCaT_g.*(gVr(:,11)+.5*giR1(:,11)).^2.*(gVr(:,12)+.5*giR1(:,12)).*(vAll_gc+.5*rhsV1g-Eca_g);
    caCr2g=dt*(-Ica./(2*F*thick_g) + (.05-(caCon_g+.5*caCr1g))/10 );
    
    % --- step 3 ---
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
    
    rhsV3g=(vAll_gc+.5*rhsV2g-El)+gNa_g.*(gVr(:,1)+.5*giR2(:,1)).^3.*(gVr(:,2)+.5*giR2(:,2)).*(vAll_gc+.5*rhsV2g-Ena)+...
        gA_g.*(gVr(:,3)+.5*giR2(:,3)).*(gVr(:,4)+.5*giR2(:,4)).*(vAll_gc+.5*rhsV2g-Ek)+...
        gM_g.*(gVr(:,5)+.5*giR2(:,5)).*(vAll_gc+.5*rhsV2g-Ek)+...
        gCaPN_g.*(gVr(:,7)+.5*giR2(:,7)).^2.*(gVr(:,8)+.5*giR2(:,8)).*(vAll_gc+.5*rhsV2g-Eca_g)+gKCa_g.*(gVr(:,9)+...
        .5*giR2(:,9)).*(vAll_gc+.5*rhsV2g-Ek)+gDR_g.*(gVr(:,10)+.5*giR2(:,10)).*(vAll_gc+.5*rhsV2g-Ek)+...
        gCaT_g.*(gVr(:,11)+.5*giR2(:,11)).^2.*(gVr(:,12)+.5*giR2(:,12)).*(vAll_gc+.5*rhsV2g-Eca_g)+...
        gCAN.*(caCon_g+.5*caCr2g)./(200+(caCon_g+.5*caCr2g)).*(gVr(:,13)+.5*giR2(:,13)).*(vAll_gc+.5*rhsV2g-Ecan);
    rhsV3g=dt/taug_m*(IdcG-rhsV3g);
    giR3(:,1)=dt*2.1*(ana_g(vAll_gc+.5*rhsV2g)./(ana_g(vAll_gc+.5*rhsV2g)+bna_g(vAll_gc+.5*rhsV2g)) -(gVr(:,1)+.5*giR2(:,1)) )./tauM(vAll_gc+.5*rhsV2g);
    giR3(:,2)=dt*2.1*(hInf(vAll_gc+.5*rhsV2g)-(gVr(:,2)+.5*giR2(:,2)))./tauH(vAll_gc+.5*rhsV2g);
    giR3(:,3)=dt*3.3./atau(vAll_gc+.5*rhsV2g).*(aminf_g(vAll_gc+.5*rhsV2g)-(gVr(:,3)+.5*giR2(:,3))); %I_A
    giR3(:,4)=dt*3.3./hatau_g(vAll_gc+.5*rhsV2g).*(ahinf_g(vAll_gc+.5*rhsV2g)-(gVr(:,4)+.5*giR2(:,4))); %I_A
    giR3(:,5)=dt*(minfMusc(vAll_gc+.5*rhsV2g)-(gVr(:,5)+.5*giR2(:,5)))./tauMusc(vAll_gc+.5*rhsV2g); %I_M
    giR3(:,6)=dt*2.1./hcurmtau(vAll_gc+.5*rhsV2g).*(hcurminf(vAll_gc+.5*rhsV2g)-(gVr(:,6)+.5*giR2(:,6)));
    giR3(:,7)=dt*(capminf(vAll_gc+.5*rhsV2g)-(gVr(:,7)+.5*giR2(:,7)))./capmtau(vAll_gc+.5*rhsV2g); %I_CaPN
    giR3(:,8)=dt*(caphinf(vAll_gc+.5*rhsV2g)-(gVr(:,8)+.5*giR2(:,8)))./caphtau(vAll_gc+.5*rhsV2g);
    giR3(:,9)=dt*(kcaa(vAll_gc+.5*rhsV2g,caCon_g+.5*caCr2g).*(1-(gVr(:,9)+.5*giR2(:,9))) - 0.05*(gVr(:,9)+.5*giR2(:,9)));
    giR3(:,10)=dt*3.3*(mss_dr(vAll_gc+.5*rhsV2g)-(gVr(:,10)+.5*giR2(:,10)))./taum_dr(vAll_gc+.5*rhsV2g);    %I_DR
    giR3(:,11)=dt*(catminfG(vAll_gc+.5*rhsV2g)-(gVr(:,11)+.5*giR2(:,11)))./cattaumG(vAll_gc+.5*rhsV2g); %I_CaT
    giR3(:,12)=dt*(cathinf(vAll_gc+.5*rhsV2g)-(gVr(:,12)+.5*giR2(:,12)))./cattauh(vAll_gc+.5*rhsV2g); 
    giR3(:,13)=dt*(canminf(vAll_gc+.5*rhsV2g)-(gVr(:,13)+.5*giR2(:,13)))./cantau(vAll_gc+.5*rhsV2g); %I_CAN
    %update calcium concentration
    Ica=gCaPN_g.*(gVr(:,7)+.5*giR2(:,7)).^2.*(gVr(:,8)+.5*giR2(:,8)).*(vAll_gc+.5*rhsV2g-Eca_g)+...
        gCaT_g.*(gVr(:,11)+.5*giR2(:,11)).^2.*(gVr(:,12)+.5*giR2(:,12)).*(vAll_gc+.5*rhsV2g-Eca_g);
    caCr3g=dt*(-Ica./(2*F*thick_g) + (.05-(caCon_g+.5*caCr2g))/10 );
    
    % --- step 4 ---
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
    
    rhsV4g=(vAll_gc+rhsV3g-El)+gNa_g.*(gVr(:,1)+giR3(:,1)).^3.*(gVr(:,2)+giR3(:,2)).*(vAll_gc+rhsV3g-Ena)+...
        gA_g.*(gVr(:,3)+giR3(:,3)).*(gVr(:,4)+giR3(:,4)).*(vAll_gc+rhsV3g-Ek)+...
        gM_g.*(gVr(:,5)+giR3(:,5)).*(vAll_gc+rhsV3g-Ek)+...
        gCaPN_g.*(gVr(:,7)+giR3(:,7)).^2.*(gVr(:,8)+giR3(:,8)).*(vAll_gc+rhsV3g-Eca_g)+gKCa_g.*(gVr(:,9)+...
        giR3(:,9)).*(vAll_gc+rhsV3g-Ek)+gDR_g.*(gVr(:,10)+giR3(:,10)).*(vAll_gc+rhsV3g-Ek)+...
        gCaT_g.*(gVr(:,11)+giR3(:,11)).^2.*(gVr(:,12)+giR3(:,12)).*(vAll_gc+rhsV3g-Eca_g)+...
        gCAN.*(caCon_g+caCr3g)./(200+(caCon_g+caCr3g)).*(gVr(:,13)+giR3(:,13)).*(vAll_gc+rhsV3g-Ecan);
    rhsV4g=dt/taug_m*(IdcG-rhsV4g);
    giR4(:,1)=dt*2.1*(ana_g(vAll_gc+rhsV3g)./(ana_g(vAll_gc+rhsV3g)+bna_g(vAll_gc+rhsV3g)) -(gVr(:,1)+giR3(:,1)) )./tauM(vAll_gc+rhsV3g);
    giR4(:,2)=dt*2.1*(hInf(vAll_gc+rhsV3g)-(gVr(:,2)+giR3(:,2)))./tauH(vAll_gc+rhsV3g);
    giR4(:,3)=dt*3.3./atau(vAll_gc+rhsV3g).*(aminf_g(vAll_gc+rhsV3g)-(gVr(:,3)+giR3(:,3))); %I_A
    giR4(:,4)=dt*3.3./hatau_g(vAll_gc+rhsV3g).*(ahinf_g(vAll_gc+rhsV3g)-(gVr(:,4)+giR3(:,4))); %I_A
    giR4(:,5)=dt*(minfMusc(vAll_gc+rhsV3g)-(gVr(:,5)+giR3(:,5)))./tauMusc(vAll_gc+rhsV3g); %I_M
    giR4(:,6)=dt*2.1./hcurmtau(vAll_gc+rhsV3g).*(hcurminf(vAll_gc+rhsV3g)-(gVr(:,6)+giR3(:,6)));
    giR4(:,7)=dt*(capminf(vAll_gc+rhsV3g)-(gVr(:,7)+giR3(:,7)))./capmtau(vAll_gc+rhsV3g); %I_CaPN
    giR4(:,8)=dt*(caphinf(vAll_gc+rhsV3g)-(gVr(:,8)+giR3(:,8)))./caphtau(vAll_gc+rhsV3g);
    giR4(:,9)=dt*(kcaa(vAll_gc+rhsV3g,caCon_g+caCr3g).*(1-(gVr(:,9)+giR3(:,9))) - 0.05*(gVr(:,9)+giR3(:,9)));
    giR4(:,10)=dt*3.3*(mss_dr(vAll_gc+rhsV3g)-(gVr(:,10)+giR3(:,10)))./taum_dr(vAll_gc+rhsV3g);    %I_DR
    giR4(:,11)=dt*(catminfG(vAll_gc+rhsV3g)-(gVr(:,11)+giR3(:,11)))./cattaumG(vAll_gc+rhsV3g); %I_CaT
    giR4(:,12)=dt*(cathinf(vAll_gc+rhsV3g)-(gVr(:,12)+giR3(:,12)))./cattauh(vAll_gc+rhsV3g); 
    giR4(:,13)=dt*(canminf(vAll_gc+rhsV3g)-(gVr(:,13)+giR3(:,13)))./cantau(vAll_gc+rhsV3g); %I_CAN
    %update calcium concentration
    Ica=gCaPN_g.*(gVr(:,7)+giR3(:,7)).^2.*(gVr(:,8)+giR3(:,8)).*(vAll_gc+rhsV3g-Eca_g)+...
        gCaT_g.*(gVr(:,11)+giR3(:,11)).^2.*(gVr(:,12)+giR3(:,12)).*(vAll_gc+rhsV3g-Eca_g);
    caCr4g=dt*(-Ica./(2*F*thick_g) + (.05-(caCon_g+caCr3g))/10 );
    
    % FINAL Step... update voltage; add implicit trap (Crank-nicholson) diffusion
    bv(idGv)=1/6*(rhsV1+2*rhsV2+2*rhsV3+rhsV4);
    bv(sidP2M)=bv(sidP2M)+dt/C*(-wM_bg*bgAmp_m.*(vAll_mc(sidP2M)-Excs)); %synap coupling
    bv(sidG2M)=bv(sidG2M)+dt/C*(-wG2M.*synG2M.*(vAll_mc(sidG2M)-Einh));    
    vAll_mc = Al\( Ar*vAll_mc + bv );
    
    bv_g(idGv_g)=1/6*(rhsV1g+2*rhsV2g+2*rhsV3g+rhsV4g);
    bv_g(sidM2G)=bv_g(sidM2G)+dt/C*(-(wM2G.*synM2G+wM2Gn.*synM2Gn+wG_bg*bgAmp_g).*(vAll_gc(sidM2G)-Excs)); %synap coupling
    vAll_gc = Al_g\( Ar_g*vAll_gc + bv_g );
    %update gating variables with RK4
    xVr=xVr+1/6*(xiR1+2*xiR2+2*xiR3+xiR4);
    gVr=gVr+1/6*(giR1+2*giR2+2*giR3+giR4);
    
    %update calcium concentration
    %Ica=gCaL*xVr(:,7).*xVr(:,8).*(V(:,j)-Eca);
    caCon=caCon+1/6*(caCr1+2*caCr2+2*caCr3+caCr4);
    caCon_g=caCon_g+1/6*(caCr1g+2*caCr2g+2*caCr3g+caCr4g);
    %update calcium reversal potential, Nernst eqn
    Eca=RTovF/2*log(10./caCon); 
    Eca_g=RTovF/2*log(10./caCon_g);
    
    %diff depending on # recording sites (faster w/ 1?)
    if(lenRec==1) %MC spike?
        if(vAll_mc(idRec)>vltThres && (dt*(j-indLst)>minTspk))
            fr_MC=fr_MC+1; %counts
            indLst=j;
        end
    else
        indCurSpk=find((vAll_mc(idRec)>vltThres).*(dt*(j-indLst)>minTspk));  %in MC
        for jCurSpk=1:length(indCurSpk)
            fr_MC(indCurSpk(jCurSpk))=fr_MC(indCurSpk(jCurSpk))+1; %counts
            indLst(indCurSpk(jCurSpk))=j; %update index of last spike (could have done this with a vector)..
        end
    end
    if(lenRec_g==1 && numGC==1)
        if(vAll_gc(idRec_g)>vltThres && (dt*(j-indLst_g)>minTspk))
            fr_GC=fr_GC+1;
            indLst_g=j;
        end
    else
        indCurSpk=find((vAll_gc(idRec_g)>vltThres).*(dt*(j-indLst_g)>minTspk));  %in GC
        for jCurSpk=1:length(indCurSpk)
            gcId=ceil(indCurSpk(jCurSpk)/lenRec_g); %which GC cell
            compId=mod(indCurSpk(jCurSpk),lenRec_g)+1; %which compart in gcId cell
            fr_GC(gcId,compId)=fr_GC(gcId,compId)+1; %counts only
            indLst_g(indCurSpk(jCurSpk))=j; %update index of last spike (could have done this with a vector)..
        end
    end 
    %store average voltage
   v_mcAvg=v_mcAvg+vAll_mc(idRec)/(Lt*Nrlz); %running sum
   v_gcAvg=v_gcAvg+vAll_gc(idRec_g)/(Lt*Nrlz); %running sum
  
end %time-loop j=1:Lt
end %realization loop k=1:Nrlz

%scale firing rates, in Hz
fr_MC=fr_MC/(tEnd*Nrlz)*1000;
fr_GC=fr_GC/(tEnd*Nrlz)*1000;

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
    function x=catminfG(v) %change for rev PGC!
        x=1./(1+exp(-(v+44)./5.5)); %(v+59) for PGC
    end
    function x=cattaumG(v) %change for rev PGC!
        x=1.5+3.5./(exp(-(v+30)./15)+exp((v+30)./15)); %(v+45) for PGC
    end
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
   
end %end of main function