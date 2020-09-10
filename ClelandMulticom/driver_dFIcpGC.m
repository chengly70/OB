%script to get transfer of MC coupled to GC

dI=50;
Ibg_v=(100 : dI : 2900)';
Ibg_GC=9000; %fixed
len_I=length(Ibg_v);
%-outputs--
frateMC=zeros(len_I,1);
frateGC=zeros(len_I,1);
voltMC=zeros(len_I,1);
voltGC=zeros(len_I,1);

wghts=100*ones(3,1);

tEnd=5000; %5 secs of biol time * 10 realz = 50sec for each parm

for j=1:len_I
    Ibg=[Ibg_v(j);Ibg_GC];
    [fr_MC,fr_GC,v_mcAvg,v_gcAvg]=getTrans_GC(tEnd,Ibg,wghts);
    
   frateMC(j)=fr_MC;
   frateGC(j)=fr_GC;
   voltMC(j)=v_mcAvg;
   voltGC(j)=v_gcAvg;
   
   %save results
   save('dFI_coupldGc','frateMC','frateGC','voltMC','voltGC','Ibg_v')
end