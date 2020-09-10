%script to get transfer of MC coupled to PGC

dI=50;
Ibg_v=(100 : dI : 2900)';
len_I=length(Ibg_v);
%-outputs--
frateMC=zeros(len_I,1);
fratePGC=zeros(len_I,1);
voltMC=zeros(len_I,1);
voltPGC=zeros(len_I,1);

wghts=100*ones(3,1);

tEnd=5000; %5 secs of biol time * 10 realz = 50sec for each parm

for j=1:len_I
    Ibg=[0.4*Ibg_v(j);Ibg_v(j)];
    [fr_MC,fr_PGC,v_mcAvg,v_pgcAvg]=getTrans_PGC(tEnd,Ibg,wghts);
    
   frateMC(j)=fr_MC;
   fratePGC(j)=fr_PGC;
   voltMC(j)=v_mcAvg;
   voltPGC(j)=v_pgcAvg;
   %save results
   save('dFI_coupldPgc','frateMC','fratePGC','voltMC','voltPGC','Ibg_v')
end