function [gBest,gbestvalue,FEcount,etime,MaxFEs,Convergence_curve,iter] = maineHGS(Lmax1,level,Pxy,MaxFEs,Population_size,objectiveFunction,Iteration)
d = level;
NumberOfThresholdValues=2*d;
thdim=NumberOfThresholdValues;
Lmax=Lmax1;
XVmax = Lmax*ones(1,NumberOfThresholdValues);
XVmin = zeros(1,NumberOfThresholdValues);
% Population_size=(thdim/2)*10;
% MaxFEs=(thdim/2)*10000;
% objectiveFunction='renyi2d';
tic
[Fbest,Lbest,FE,MaxFEs,Convergence_curve,iter]=HGS(objectiveFunction,Population_size,thdim,XVmin,XVmax,MaxFEs,Pxy,Iteration);
% [Fbest,Lbest,FE,iteration]=eKGSA(objectiveFunction,Population_size,thdim,XVmin,XVmax,max_it,ElitistCheck,flag,Rpower,h,Pxy);
toc
gBest=Lbest;
gbestvalue=Fbest;
FEcount=FE;
etime = toc;
end