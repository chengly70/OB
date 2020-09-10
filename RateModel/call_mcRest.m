function []=call_mcRest(indR)
% Function to run fr model saving the rest of 
% stats (PGC, GC, cov's), for varying sigma
% calling mc_Rest.m 

N_relz=100000;

Nc=7; %1=left PGC, 2=left MC, 3=left GC, 4=common GC, 5=right GC, 6=right MC, 7=right PGC
tauVec=zeros(Nc,1); %time constants, %set below

%time vector
dt=0.001;
tmv=(-2 : dt : 2)'; %in seconds
idSt=round( (.05-tmv(1))/dt )+1; %changed from 0
%idZ=round( (0-tmv(1))/dt )+1; %index of 0

muMat=zeros(Nc,length(tmv));
muMat([1 7],:)=6*ones(2,length(tmv)); %PGC background
muMat([3 4 5],:)=2*ones(3,length(tmv)); %GC
muMat([2 6],:)=14*ones(2,length(tmv)); %MC
%odor inputs
evokMu=10+380*(tmv(idSt:end)-.05).*exp(-(tmv(idSt:end)-.05)./.085);
muMat([2 6],idSt:end)=muMat([2 6],idSt-1)+repmat(evokMu',2,1);
muMat(6,idSt:end)=muMat(2,idSt:end); %same MC input
muMat([1 7],idSt:end)=muMat([1 7],idSt-1)+.5*repmat(evokMu',2,1);
% GC get odor too? unmodeled
muMat(3:5,idSt:end)=muMat(3:5,idSt-1)+.5*repmat(evokMu',3,1);

sig_vec1=0.5*[.75; 1; 1.5; 1.5; 1.5; 1; .75]; %spontaneous noise
sig_vec2=sig_vec1;
sig_vec2([1 2 6 7])=1.5*sig_vec2([1 2 6 7]); %evoked noise
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

%set Struct of Params
Pstr_Par=struct('tau_vec',tauVec,'CinMat',CinMat,'cellTyp',cellTyp);

nameStr='nRst';

tic
    % vary how Gm changes with indR in function
    mc_Rest(Nc,Pstr_Par,tmv,muMat,sig_vec1,sig_vec2,N_relz,nameStr,indR);
toc

end