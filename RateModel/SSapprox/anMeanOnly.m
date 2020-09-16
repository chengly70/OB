function [convged,mnXan,mean_all]=anMeanOnly(Nc,Pstruct,mu_vec)
%[convged,mnXan,mean_all]=anMeanOnly(Nc,Pstruct,mu_vec)
% like an_Method.m BUT only solving for mean_X/F assuming sig=0

Num_it=200; %(total) MAX # of iterations
Tolnc=1e-4; %tolerance level for convergence

%synaptic parameters from LiCleland (1/sec)
alph_gaba=1000/1.25;
bet_gaba=1000/18;
alph_ampa=1000/1;
bet_ampa=1000/5.5;
alph_nmda=1000/52;
bet_nmda=1000/343;

% Other parameters
Gm          = Pstruct.Gm;
cellTyp     = Pstruct.cellTyp; %0=Excit (MC), 1=Inhib (GC), 2=Inhib (PGC)

% load precomputed FI (transfer) curves
load('dFI_coupldGc.mat','frateMC','Ibg_v','frateGC');
input_F=Ibg_v./100; %scale down, extra intrins currents
FI_mg=frateMC;
FI_gc=frateGC;
load('dFI_coupldPgc.mat','frateMC'); %assuming SAME Ibg_v
%input_pg=Ibg_v; %assuming same as dFI_coupldGc
FI_mp=frateMC;

barR=max([FI_mg;FI_mp;FI_gc]); %max firing rate
%padd transfer/FI curves with 0s & large FR
tmp_f1=(-10:1:input_F(1)-1)';
input_F=[tmp_f1; input_F];
FI_mg=[zeros(size(tmp_f1)); FI_mg];
FI_gc=[zeros(size(tmp_f1)); FI_gc];
FI_mp=[zeros(size(tmp_f1)); FI_mp];
tmp_f1=(input_F(end)+2 : 2 : input_F(end)+10)';
input_F=[input_F; tmp_f1];
FI_mg=[FI_mg;FI_mg(end)*ones(size(tmp_f1))];
FI_gc=[FI_gc;FI_gc(end)*ones(size(tmp_f1))];
FI_mp=[FI_mp;FI_mp(end)*ones(size(tmp_f1))];

%EASIER/FASTEF when Ibg_v are the same (assuming it is)
F_Mat=zeros(length(input_F),Nc);
for j=1:Nc
    switch cellTyp(j)
        case 0 %MC cell
            F_Mat(:,j)=.5*FI_mp+.5*FI_mg;
        case 1 %GC cell
            F_Mat(:,j)=FI_gc;
        case 2 %PGC cell
            F_Mat(:,j)=FI_mp*exp(-.3);
    end
end

% --- 
synM=zeros(Nc,1);
frt=zeros(Nc,1);

%OUTPUTS; iterative solver, start at uncoupled values
%simulation does not converge by default
convged=(1==0);

%save all of the means at each step
mean_all=mu_vec;
mean_tmp=zeros(Nc,1);
%only save previous var/cov to check convergence

for j=2:Num_it
    %instance of F(sig(j-1)*y+mu(j-1)), len_nrm x Nc matrix for 1D-integ
    for k=1:Nc
        frt(k)=interp1(input_F,F_Mat(:,k),mean_all(k,j-1),'pchip');
        %update synapses
        if(cellTyp(k)==0)%excit
            synM(k)=alph_ampa/barR*frt(k)/(bet_ampa+alph_ampa/barR*frt(k)) + ...
                alph_nmda/barR*frt(k)/(bet_nmda+alph_nmda/barR*frt(k));
        else %inhib; PGC or GC
            synM(k)=alph_gaba/barR*frt(k)/(bet_gaba+alph_gaba/barR*frt(k));
        end
    end
    
    mean_tmp=mu_vec+Gm*synM;
   
    %set tmp variables and updates
    mean_all=[mean_all mean_tmp]; %dynamic increase; can overwrite if want
    
    %check if iteration converges
    if(j>4 && norm(mean_all(:,j-1)-mean_all(:,j))<Tolnc )
        convged=~convged;
        break; %all quantities have converged
    end 
    
end

mnXan=mean_all(:,end);