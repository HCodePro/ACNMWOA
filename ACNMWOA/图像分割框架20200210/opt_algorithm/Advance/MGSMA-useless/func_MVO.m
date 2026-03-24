function [Universes]=func_MVO(Inflation_rates,Universes,Best_universe,N,FEs,MaxFEs,ub,lb)
    WEP_Max=1;
    WEP_Min=0.2;
    %Eq. (3.3) in the paper
    WEP=WEP_Min+FEs*((WEP_Max-WEP_Min)/MaxFEs);
    
    %Travelling Distance Rate (Formula): Eq. (3.4) in the paper
    TDR=1-((FEs)^(1/6)/(MaxFEs)^(1/6));
    
    [sorted_Inflation_rates,sorted_indexes]=sort(Inflation_rates);
    
    for newindex=1:N
        Sorted_universes(newindex,:)=Universes(sorted_indexes(newindex),:);
    end
    
    %Normaized inflation rates (NI in Eq. (3.1) in the paper)
    normalized_sorted_Inflation_rates=normr(sorted_Inflation_rates);
    
    Universes(1,:)= Sorted_universes(1,:);
    
    %Update the Position of universes
    for i=2:size(Universes,1)%Starting from 2 since the firt one is the elite
        Back_hole_index=i;
        for j=1:size(Universes,2)
            r1=rand();
            if r1<normalized_sorted_Inflation_rates(i)
                White_hole_index=RouletteWheelSelection(-sorted_Inflation_rates);% for maximization problem -sorted_Inflation_rates should be written as sorted_Inflation_rates
                if White_hole_index==-1
                    White_hole_index=1;
                end
                %Eq. (3.1) in the paper
                Universes(Back_hole_index,j)=Sorted_universes(White_hole_index,j);
            end
            
            if (size(lb,2)==1)
                %Eq. (3.2) in the paper if the boundaries are all the same
                r2=rand();
                if r2<WEP
                    r3=rand();
                    if r3<0.5
                        Universes(i,j)=Best_universe(1,j)+TDR*((ub-lb)*rand+lb);
                    end
                    if r3>0.5
                        Universes(i,j)=Best_universe(1,j)-TDR*((ub-lb)*rand+lb);
                    end
                end
            end
            
            if (size(lb,2)~=1)
                %Eq. (3.2) in the paper if the upper and lower bounds are
                %different for each variables
                r2=rand();
                if r2<WEP
                    r3=rand();
                    if r3<0.5
                        Universes(i,j)=Best_universe(1,j)+TDR*((ub(j)-lb(j))*rand+lb(j));
                    end
                    if r3>0.5
                        Universes(i,j)=Best_universe(1,j)-TDR*((ub(j)-lb(j))*rand+lb(j));
                    end
                end
            end
            
        end
    end

end