%script to run SS approx method

which_net=1; %1=vary sig (n), 2=smaller sig (sc), 3=larger sig (cs)

sig_vec1=0.5*[.75; 1; 1.5; 1.5; 1.5; 1; .75]; %spontaneous noise
sig_vec2=sig_vec1;
sig_vec2([1 2 6 7])=1.5*sig_vec2([1 2 6 7]); %evoked noise

Nc=7; %7 cell
tauVec=ones(Nc,1);
Gm=zeros(Nc,Nc);


CinMat=diag(ones(Nc,1)); %blocks of weakly correlated noise
CinMat(1,2)=0.1; CinMat(6,7)=0.1; CinMat(3,4)=0.05; CinMat(3,5)=0.05; CinMat(4,5)=0.05;
    CinMat(2,6)=0.05;
CinMat=CinMat+tril(CinMat',-1); %make symmetric (lower triangular)

cellTyp=[2; 0; 1; 1; 1; 0; 2]; %0=Excit (MC), 1=Inhib (GC), 2=Inhib (PGC)
C=1.2; %micro-F/cm^2 for MC & PGC
Cg=2; %micro-F/cm^2 for GC
Rm=30; %kilo-Ohm*cm^2 for MC & GC
Rm_p=20; %kilo-Ohm*cm^2 for PGC
tau_m=C*Rm; %membrane time constant (ms)
taug_m=Cg*Rm; %membrane time constant (ms)
taup_m=C*Rm_p; %membrane time constant (ms)
for j=1:Nc
    switch cellTyp(j)
        case 0 %MC
            tauVec(j)=tau_m/1000; %in sec now
        case 1 %GC
            tauVec(j)=taug_m/1000; %in sec now
        case 2 %PGC
            tauVec(j)=taup_m/1000; %in sec now
    end
end

load dGs_n %gets Gs_good matrix

gP2M=-5;
gM2P=5;
gG2M=Gs_good(10,1);  % !!! set 3 values here !!!, can alter
gGcin=Gs_good(10,2);
gM2G=Gs_good(10,3);
gG2G=0;%not in Cleland, MarellaErmentrout, etc.
Gm=[0   gM2P 0     0 0 0 0; 
   gP2M 0    gG2M gGcin 0 0 0;
   0  gM2G 0   gG2G 0 0 0;
   0  gM2G   gG2G  0 gG2G gM2G 0;
   0 0 0       gG2G  0      gM2G 0;
   0 0 0      gGcin gG2M 0  gP2M;
   0 0 0 0 0           gM2P 0];
Gm(4,:)=.5*Gm(4,:); %common GC getting 2 copies of E/I inputs

muMat=zeros(Nc,2); %1st col Spont; 2nd col is EVOKED
muMat([1 7],1)=6*ones(2,1); %PGC background
muMat([3 4 5],1)=2*ones(3,1); %GC
muMat([2 6],1)=14*ones(2,1); %MC
%odor inputs; SS
evokMu=10;
muMat([2 6],2)=muMat([2 6],1)+evokMu;
muMat([1 7],2)=muMat([1 7],1)+.5*evokMu;
% GC get odor too? unmodeled
muMat(3:5,2)=muMat(3:5,1)+.5*evokMu;

Pstruct=struct('Gm',Gm,'cellTyp',cellTyp); %for anMeanOnly.m to re-scale mu's

%set Struct of Params for an_MethdOrd1
Pstr=struct('tau_vec',tauVec,'CinMat',CinMat,'Gm',Gm,'cellTyp',cellTyp);

%get re-scaled mean(X), for spontaneous
[convged1,mnXSp]=anMeanOnly(Nc,Pstruct,muMat(:,1));
[convged2,mnXCvSp]=anMeanOnly(Nc,Pstruct,muMat(:,1)-mean(muMat(:,1)));
%get re-scaled mean(X), for Evoked
[convged3,mnXEv]=anMeanOnly(Nc,Pstruct,muMat(:,2));
[convged4,mnXCvEv]=anMeanOnly(Nc,Pstruct,muMat(:,2)-mean(muMat(:,2)));
if(convged1==0 || convged2==0 || convged3==0 || convged4==0)
    disp('issue with rescaled mean(s)!')
    return;
end

%spontaneous
[convged1,Corr_val1,cov_Fa,~,cov_Xa,~]=an_MethOrd1(Nc,Pstr,muMat(:,1),mnXCvSp,sig_vec1);
[convged2,Corr_val2,~,mn_Fa,~,mn_Xa]=an_MethOrd1(Nc,Pstr,muMat(:,1),mnXSp,sig_vec1);

%evoked
[convged3,Corr_val3,cov_FaE,~,cov_XaE,~]=an_MethOrd1(Nc,Pstr,muMat(:,2),mnXCvEv,sig_vec2);
[convged4,Corr_val4,~,mn_FaE,~,mn_XaE]=an_MethOrd1(Nc,Pstr,muMat(:,2),mnXEv,sig_vec2);
