clearvars;      
close all;

%% Run simulations
[bs,delta,psi,mu,taut,nvf,nstf,taun,taua,mun,nrst,asp,runTm,tmSt,svPt,edpc,rpc]=...
    deal(6,10,1,0.5,10,0,50,1,1,0,1,0.5,5,0.01,10,0.15,1);

ns_simDynamics(bs,delta,psi,mu,taut,nvf,...
    nstf,taun,taua,mun,nrst,asp,runTm,tmSt,svPt,edpc,rpc);
