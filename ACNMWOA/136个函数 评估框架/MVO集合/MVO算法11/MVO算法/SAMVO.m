function [Best_universe,Convergence_curve]=saMVO(N,MaxFEs,lb,ub,dim,fobj)

%Two variables for saving the position and inflation rate (fitness) of the best universe
Best_universe=zeros(1,dim);
Best_universe_Inflation_rate=inf;

%Initialize the positions of universes
Universes=initialization(N,dim,ub,lb);

%Minimum and maximum of Wormhole Existence Probability (min and max in
% Eq.(3.3) in the paper
WEP_Max=1;
WEP_Min=0.2;
FEs=0;
Convergence_curve=[];

%Iteration(time) counter
Time=1;

%Main loop
while FEs<MaxFEs
    
    %Eq. (3.3) in the paper
    WEP=WEP_Min+FEs*((WEP_Max-WEP_Min)/MaxFEs);
    
    %Travelling Distance Rate (Formula): Eq. (3.4) in the paper
    TDR=1-((FEs)^(1/6)/(MaxFEs)^(1/6));
    
    %Inflation rates (I) (fitness values)
    Inflation_rates=zeros(1,size(Universes,1));
    
    T=(min(Inflation_rates)-max(Inflation_rates))/log(0.8);%??�??温度
    r=0.95;%??温�????
    c=0;
    while c<N
        for i=1:size(Universes,1)
            
            %Boundary checking (to bring back the universes inside search
            % space if they go beyoud the boundaries
            Flag4ub=Universes(i,:)>ub;
            Flag4lb=Universes(i,:)<lb;
            %TempUniverse=(Universes(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            TempUniverse=0.5.*(Universes(i,:)-(ub-lb)/(FEs*FEs)+(ub-lb)*2*rand/(FEs*FEs));
            %Calculate the inflation rate (fitness) of universes
            Inflation_rates(1,i)=fobj(Universes(i,:));
            FEs=FEs+1;
            TempUniverse_rate=fobj(TempUniverse);
            FEs=FEs+1;

            dE=Inflation_rates(1,i)-TempUniverse_rate;
            if dE>0
                Universes(i,:)=TempUniverse;
                Inflation_rates(1,i)=TempUniverse_rate;
                if TempUniverse_rate<Best_universe_Inflation_rate
                    Best_universe_Inflation_rate=TempUniverse_rate;
                    Best_universe=TempUniverse;
                end
            else
                if exp(-(dE)/T)>rand
                    Universes(i,:)=TempUniverse;
                    Inflation_rates(1,i)=TempUniverse_rate;
                end
            end
            T=T*r;
      c=c+1;
        end
    end
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
    
    Convergence_curve(Time)=Best_universe_Inflation_rate;
    
    %Print the best universe details after every 50 iterations
%     if mod(Time,50)==0
%         display(['At iteration ', num2str(Time), ' the best universes fitness is ', num2str(Best_universe_Inflation_rate)]);
%     end
    Time=Time+1;
end
