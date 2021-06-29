function []=cal_mcAwake(indRun)
% Function to run fr model like covSubs_[] 
% using mc_RateAwk.m (just MC stats, no other cells)
% UNLIKE prior files, can calc all 3 (and each has 100 total, easier to code)
% vary indRun=1:300 in .submit files

N_relz=100000;

Nc=7; %1=left PGC, 2=left MC, 3=left GC, 4=common GC, 5=right GC, 6=right MC, 7=right PGC
tauVec=zeros(Nc,1); %time constants, %set below

%time vector
dt=0.001;
tmv=(-0.5 : dt : 0.5)'; %in seconds
idZ=round( (0-tmv(1))/dt )+1; %index of 0

muMat=zeros(Nc,length(tmv));
muMat([1 7],:)=5*ones(2,length(tmv)); %PGC background, GC-2X,MC-2.5X cf. anesth (all rest same)
muMat([3 4 5],:)=4*ones(3,length(tmv)); %GC
muMat([2 6],:)=25*ones(2,length(tmv)); %MC
%odor inputs
evokMu=-2+150*tmv(idZ:end).*exp(-tmv(idZ:end)./.1); %change 2 to -2 in front, 270 to 150 height-exp
muMat([2 6],idZ:end)=muMat([2 6],idZ-1)+repmat(evokMu',2,1);
muMat(6,idZ:end)=muMat(2,idZ:end); %same MC input
muMat([1 7],idZ:end)=muMat([1 7],idZ-1)+.6*repmat(evokMu',2,1);
% GC get odor too? unmodeled
muMat(3:5,idZ:end)=muMat(3:5,idZ-1)+.6*repmat(evokMu',3,1);

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

whchSet=floor((indRun-1)/100)+1; %1-100,_sc, 101-200,_n, 201-300,_cs
indRun=mod(indRun,100); 
if(indRun==0)
    indRun=100;
end

tic
    % vary how Gm changes with indR in function
    mc_RateAwk(Nc,Pstr_Par,tmv,muMat,N_relz,indRun,whchSet);
toc

end