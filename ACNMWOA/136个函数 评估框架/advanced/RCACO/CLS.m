function [ BestPosition,BestFitness,fes ] = CLS(K,Max_iter,fes,lb,ub,BestPosition,BestFitness,fobj)
%CHAOS 此处显示有关此函数的摘要
%   此处显示详细说明
        
	% m = 2500;
	
% setCan = 1-power(abs((FEs-1)/FEs),m);

setCan = (Max_iter-fes+1)/Max_iter;
    x = rand();
    while(~(x~=0.25 && x~=0.5 && x~=0.75 && x~=1))
        x=rand();
    end
    ch(1) = x;
for ant=1:K
    ch(ant+1)=4*ch(ant)*(1-ch(ant));
    CH(ant,:) = lb+ch(ant)*(ub-lb);    
    V = (1-setCan)*BestPosition+setCan*CH(ant);
    %% 边界控制
    Flag4ub=V>ub';
    Flag4lb=V<lb';
    V=(V.*(~(Flag4ub+Flag4lb)))+ub'.*Flag4ub+lb'.*Flag4lb;
    %% 边界控制结束
    fes = fes + 1;
    FitnessV=fobj(V);%计算适应度值
    if (FitnessV<BestFitness)
        BestFitness = FitnessV;
        BestPosition = V;
        break;
    end
 end
     

end

