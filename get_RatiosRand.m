% script to get random samples of ratios

% INPUTS: ratcVr,ratCov; ratios from Monte Carlo simulations of firing rate model
ratcVr=cell(3,1); %to match cov # elem
ratCov=cell(3,1);

% time vectors
dt=0.001;
tmv=(-2 : dt : 2)'; %in seconds
Lt=length(tmv); 
Twin=0.1; %large windows
twinv=(-2 : Twin : 2-Twin)';
tspe=-1; %ending time for spontaneous
indS=round((tspe-tmv(1))/dt)+1; %index start spont
tsek=1;
indEs=round((tsek-tmv(1))/dt)+1; %index start evoked
allInd=[1:indS indEs:Lt]; %combine spon & evoked
    
% example code to get ratio from 1 simulation file
for which_net=1:3 
    %loading same data, just to demonstrate
    load nRst242.mat %sims from 1 run
    varTw=interp1(twinv',vrAl_tw(2,:),tmv','pchip'); %interp between Twin to get dt spacing
    covTw=interp1(twinv',covAl_tw(1,:),tmv','pchip');

    ratcVr{which_net,1}=[ratcVr{which_net,1} vrAlt(2,allInd)./varTw(allInd)];
    ratCov{which_net,1}=[ratCov{which_net,1} covAlt(1,allInd)./covTw(allInd)];
end

%combine from all 3 noise regimes, unlike coupling strengths kept separate
varcRat_all=[ratcVr{1,1} ratcVr{2,1} ratcVr{3,1}]'; 
covRati_all=[ratCov{1,1} ratCov{2,1} ratCov{3,1}]';
nSmp_Rat=500; %total # samples; Ratios to test

rndId=ceil(numel(covRati_all)*rand(nSmp_Rat,1)); %uniform random sample

%using same rndId index keeps (Var,Cov) structure
vrRat_smpl=varcRat_all(rndId); 
cvRat_smpl=covRati_all(rndId);

save dRats_smp vrRat_smpl cvRat_smpl nSmp_Rat