function [Fbest,Lbest,FE,MaxFEs,Convergence_curve,iter]=EvolutionLGCMFO(fname,N,thdim,lb,ub,MaxFEs,Pxy,Iteration)
Fbest=-inf;
dim=thdim/2;
FEs=0;
%Initialize the positions of moths
Moth_pos=random_initialization(N,dim,ub,lb);
Convergence_curve=[];

it=1;

% Main loop
while FEs < MaxFEs || it<=Iteration
    
    % Number of flames Eq. (3.14) in the paper
    Flame_no=round(N-FEs*((N-1)/MaxFEs));
    
    for i=1:size(Moth_pos,1)
        
        % Check if moths go out of the search spaceand bring it back
        Flag4ub=Moth_pos(i,:)>ub;
        Flag4lb=Moth_pos(i,:)<lb;
        Moth_pos(i,:)=(Moth_pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
        
        % Calculate the fitness of moths
		if FEs < MaxFEs
			[Moth_fitness(1,i),FEs,Moth_pos(i,:)] =sigle_evaluation(Moth_pos(i,:),dim,thdim,fname,Pxy,FEs);
		else
			break;
		end
        
    end
       
    if it==1
        % Sort the first population of moths
        [fitness_sorted I]=sort(Moth_fitness,'descend');
        sorted_population=Moth_pos(I,:);
        
        % Update the flames
        best_flames=sorted_population;
        best_flame_fitness=fitness_sorted;
    else
        
        % Sort the moths
        double_population=[previous_population;best_flames];
        double_fitness=[previous_fitness best_flame_fitness];
        
        [double_fitness_sorted I]=sort(double_fitness,'descend');
        double_sorted_population=double_population(I,:);
        
        fitness_sorted=double_fitness_sorted(1:N);
        sorted_population=double_sorted_population(1:N,:);
        
        % Update the flames
        best_flames=sorted_population;
        best_flame_fitness=fitness_sorted;
    end
    
    % Update the position best flame obtained so far
    Best_flame_score=fitness_sorted(1);
    Best_flame_pos=sorted_population(1,:);
      
    previous_population=Moth_pos;
    previous_fitness=Moth_fitness;
    
    % a linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a=-1+FEs*((-1)/MaxFEs);
    k=1;
    for i=1:size(Moth_pos,1)
        
        for j=1:size(Moth_pos,2)
            if i<=Flame_no % Update the position of the moth with respect to its corresponsing flame
                
                % D in Eq. (3.13)
                %sorted_population火焰

                distance_to_flame=abs(sorted_population(i,j)-Moth_pos(i,j));
                b=1;
                t=(a-1)*rand+1;
                
                % Eq. (3.12)
                %飞蛾随火焰更新位置
                Moth_pos_s(i,j)=distance_to_flame*exp(b.*t).*cos(t.*2*pi)+sorted_population(i,j);
            end
            
            if i>Flame_no % Upaate the position of the moth with respct to one flame
                
                % Eq. (3.13)
                distance_to_flame=abs(sorted_population(i,j)-Moth_pos(i,j));
                b=1;
                t=(a-1)*rand+1;
                
                % Eq. (3.12)
                Moth_pos_s(i,j)=distance_to_flame*exp(b.*t).*cos(t.*2*pi)+sorted_population(Flame_no,j);
            end          
        end
        
        Moth_pos_m_levy1=Moth_pos_s(i,:)*(1+k*Levy1(1));
		Flag4ub=Moth_pos_m_levy1>ub;
        Flag4lb=Moth_pos_m_levy1<lb;
        Moth_pos_m_levy1=(Moth_pos_m_levy1.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
		
        Moth_pos_m_gaus=Moth_pos_s(i,:)*(1+k*randn(1));
		Flag4ub=Moth_pos_m_gaus>ub;
        Flag4lb=Moth_pos_m_gaus<lb;
        Moth_pos_m_gaus=(Moth_pos_m_gaus.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
		
        r= cauchyrnd(0, 1, 1);t=randperm(1); cauchy=r(t(1),:);gaus=randn(1,1);
        Moth_pos_m_cauchy=Moth_pos_s(i,:)*(1+k*cauchy);
        Flag4ub=Moth_pos_m_cauchy>ub;
        Flag4lb=Moth_pos_m_cauchy<lb;
        Moth_pos_m_cauchy=(Moth_pos_m_cauchy.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
		
        Flag4ub=Moth_pos_s(i,:)>ub;
        Flag4lb=Moth_pos_s(i,:)<lb;
        Moth_pos_s(i,:)=(Moth_pos_s(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb; 
		
		if FEs < MaxFEs
			[Moth_fitness_m_levy1,FEs,Moth_pos_m_levy1] =sigle_evaluation(Moth_pos_m_levy1,dim,thdim,fname,Pxy,FEs);
		else
			break;
		end
		 
		if FEs < MaxFEs
			[Moth_fitness_m_gaus,FEs,Moth_pos_m_gaus] =sigle_evaluation(Moth_pos_m_gaus,dim,thdim,fname,Pxy,FEs);
		else
			break;
		end
		
		if FEs < MaxFEs
			[Moth_fitness_m_cauchy,FEs,Moth_pos_m_cauchy] =sigle_evaluation(Moth_pos_m_cauchy,dim,thdim,fname,Pxy,FEs);
		else
			break;
		end
		
		if FEs < MaxFEs
			[Moth_fitness_s,FEs,Moth_pos_s(i,:)] =sigle_evaluation(Moth_pos_s(i,:),dim,thdim,fname,Pxy,FEs);
		else
			break;
		end
		
        if FEs < MaxFEs
			[Moth_fitness_f,FEs,Moth_pos(i,:)] =sigle_evaluation(Moth_pos(i,:),dim,thdim,fname,Pxy,FEs);
		else
			break;
		end
        Moth_fitness_comb=[Moth_fitness_m_levy1,Moth_fitness_m_gaus,Moth_fitness_m_cauchy,Moth_fitness_s,Moth_fitness_f];
        [~,m]=max(Moth_fitness_comb);
        if m==1   %%听说switch case效率更高啊！
            Moth_pos(i,:)=Moth_pos_m_levy1;
        elseif m==2
            Moth_pos(i,:)=Moth_pos_m_gaus; 
        elseif m==3
            Moth_pos(i,:)=Moth_pos_m_cauchy; 
        elseif m==4
            Moth_pos(i,:)=Moth_pos_s(i,:); 
        else
            Moth_pos(i,:)=Moth_pos(i,:); 
        end
    end 
    Convergence_curve(it)=Best_flame_score;
    if Fbest<Best_flame_score
        FE=FEs;
        iter=it;
    end
    Fbest=Best_flame_score;
    Lbest=Best_flame_pos;	
    it=it+1; 
end
end

function o=Levy1(d)
% d=1;
beta=3/2;
%Eq. (3.10)
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);

% Eq. (3.9)
o=step;

end



