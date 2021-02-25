Matlab implementation of Li & Cleland (2013, 2017) OB smaller sub-networks.

Li, G. and T. A. Cleland (2013). A two-layer biophysical model of cholinergic neuromodulation in olfactory bulb. Journal of Neuroscience 33(7), 3037–3058.
Li, G. and T. A. Cleland (2017). A coupled-oscillator model of olfactory bulb gamma oscillations. PLoS computational biology 13(11), e1005760.

To obtain curves in manuscript "Odor-evoked Increases in Spiking Variability in the Olfactory Bulb", set 
input current of GC to 9000 and PGC to 0.4*(current for MC), and set weights to: wghts=wg*ones(3,1); where wg=60,80,100,120,140 
 
driver_dFIcpGC.m -- script that creates and saves dFI_coupldGc.mat, calls function getTrans_GC.m 
driver_dFIcpPGC.m -- script that creates and saves dFI_coupldPgc.mat, calls function getTrans_PGC.m

getTrans_GC.m -- function that implements multicompartment model with coupling between a single GC and single MC cell, reciprocally coupled. 
MC provides both AMPA and NMDA input to GC, GC provides GABA_A input; both cells subject to background AMPA inputs 
(low-pass filtered syn with Poisson Process kicks) and constant input (I in FI curve).
Read comments for inputs/outputs, further details. 

getTrans_PGC.m -- function that implements multicompartment model with coupling between a single PGC and single MC cell. Read comments for inputs/outputs, further details
MC provides both AMPA and NMDA input to PGC, PGC provides GABA_A input; both cells subject to background AMPA inputs 
(low-pass filtered syn with Poisson Process kicks) and constant input (I in FI curve).
Read comments for inputs/outputs, further details. 
NOTE: PGC and GC have different morphologies, ionic currents, etc., all highlighted in Li&Cleland
