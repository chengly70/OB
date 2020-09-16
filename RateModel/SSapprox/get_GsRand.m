% script to get random samples of coupling strengths

% --- Input: Gs_good ----
Gs_good = []; %coupling strength from simulations, 10 x 3
load dGs_n %contains Gs_good

% generate random sample with same stats as Gs_[n/sc/cs]Cm
nSmp_Gs = 1000 - size(Gs_good,1);  %size of random sample is 1000, include orig Gs_good

%get means, set random sample to have same means
mn_Gs = mean(Gs_good);
covM_Gs=cov(Gs_good); %get cov, set random sample to have same cov

R_n=chol(covM_Gs);
Gsmp_n=repmat(mn_Gs,nSmp_Gs,1) + .7*randn(nSmp_Gs,3)*R_n; %nSmp_Gs x 3
Gsmp_n=[Gs_good; Gsmp_n]; %include original

nSmp_g=nSmp_Gs+size(Gs_good,1);

save dGs_rndsmp Gsmp_n covM_Gs mn_Gs R_* nSmp_g