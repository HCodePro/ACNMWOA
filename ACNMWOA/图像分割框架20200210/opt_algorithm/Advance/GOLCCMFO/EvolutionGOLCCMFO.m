function [Fbest,Lbest,FE,MaxFEs,Convergence_curve,iter]=EvolutionGOLCCMFO(fname,N,thdim,lb,ub,MaxFEs,Pxy,Iteration)
Fbest=-inf;
dim=thdim/2;
FEs=0;
type=0;
Cr=0.7;
%Initialize the positions of moths
Moth_pos=random_initialization(N,dim,ub,lb);
for i=1:size(Moth_pos,1)

    % Check if moths go out of the search spaceand bring it back
    Flag4ub=Moth_pos(i,:)>ub;
    Flag4lb=Moth_pos(i,:)<lb;
    Moth_pos(i,:)=(Moth_pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;        
    % Calculate the fitness of moths
    if FEs<MaxFEs
        [Moth_fitness(1,i),FEs,Moth_pos(i,:)] =sigle_evaluation(Moth_pos(i,:),dim,thdim,fname,Pxy,FEs);
    else
        break;
    end
end
[ Moth_pos, FEs, Moth_fitness ] = GOL( Moth_pos,Moth_fitness,dim,ub,lb,thdim,fname,Pxy,FEs,type);             
it=1;
Convergence_curve=[];
type=1;
% Main loop
while  FEs < MaxFEs || it<=Iteration
% while Iteration<Max_iteration+1
    
    % Number of flames Eq. (3.14) in the paper
    Flame_no=round(N-FEs*((N-1)/MaxFEs));
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
    
    for i=1:size(Moth_pos,1)
        
        for j=1:size(Moth_pos,2)
            if i<=Flame_no % Update the position of the moth with respect to its corresponsing flame
                
                % D in Eq. (3.13)
                distance_to_flame=abs(sorted_population(i,j)-Moth_pos(i,j));
                b=1;
                t=(a-1)*rand+1;
                
                % Eq. (3.12)
                Moth_pos(i,j)=distance_to_flame*exp(b.*t).*cos(t.*2*pi)+sorted_population(i,j);
            end
            
            if i>Flame_no % Upaate the position of the moth with respct to one flame
                
                % Eq. (3.13)
                distance_to_flame=abs(sorted_population(i,j)-Moth_pos(i,j));
                b=1;
                t=(a-1)*rand+1;
                
                % Eq. (3.12)
                Moth_pos(i,j)=distance_to_flame*exp(b.*t).*cos(t.*2*pi)+sorted_population(Flame_no,j);
            end
            
        end
        % Check if moths go out of the search spaceand bring it back
        Flag4ub=Moth_pos(i,:)>ub;
        Flag4lb=Moth_pos(i,:)<lb;
        Moth_pos(i,:)=(Moth_pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        if FEs<MaxFEs
            [Moth_fitness(1,i),FEs,Moth_pos(i,:)] =sigle_evaluation(Moth_pos(i,:),dim,thdim,fname,Pxy,FEs);
        else
            break;
        end
       %% Cross
        Bvc = randperm(2*dim); %%此处2*dim替换dim
        Mvc(i,:) = Moth_pos(i,:);
        Mvc(i,:) = (Mvc(i,:)-lb)./(ub-lb);

        p2 = 0.6;  %p2取0.2到0.8之间
        for j=1:(2*dim)/2 %%此处2*dim替换dim
            p = unifrnd(0,1);  %生成服从均匀分布的0-1的随机数
            if  p<p2
                no1= Bvc(2*j-1);
                no2 = Bvc(2*j);
                r = unifrnd(0,1);  %生成服从均匀分布的0-1的随机数
                Mvc(i,no1)=r*Mvc(i,no1)+(1-r)*Mvc(i,no2);
            end
        end
        Mvc(i,:) = Mvc(i,:).*(ub-lb)+lb;
       %% Boundary repair
        Flag4ub=Mvc(i,:)>ub;
        Flag4lb=Mvc(i,:)<lb;
        Mvc(i,:)=(Mvc(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        if FEs<MaxFEs
            [MvcFitness(1,i),FEs,Mvc(i,:)] =sigle_evaluation(Mvc(i,:),dim,thdim,fname,Pxy,FEs);
        else
            break;
        end
    end
    ComFitness1=[Moth_fitness,MvcFitness];
    ComPositions1=[Moth_pos;Mvc];
    [ComFitness1,index]=sort(ComFitness1,'descend');
    for i=1:size(Moth_pos,1)
        Moth_pos(i,:)=ComPositions1(index(i),:);
        Moth_fitness(1,i)=ComFitness1(1,i);
    end
    %% Criss
    Bhc = randperm(N);
    MhcFitness = zeros(1,N);
    for i=1:(N/2)
        no1= Bhc(2*i-1);
        no2 = Bhc(2*i);
        for j=1:2*dim %%此处2*dim替换dim
            r1 = unifrnd(0,1);  %生成服从均匀分布的0-1的随机数
            r2 = unifrnd(0,1);  %生成服从均匀分布的0-1的随机数
            c1 = (rand(1)*2)-1; %生成服从均匀分布的-1到1的随机数
            c2 = (rand(1)*2)-1;
            Mhc(no1,j)=r1*Moth_pos(no1,j)+(1-r1)*Moth_pos(no2,j)+c1*(Moth_pos(no1,j)-Moth_pos(no2,j));
            Mhc(no2,j)=r2*Moth_pos(no2,j)+(1-r2)*Moth_pos(no1,j)+c2*(Moth_pos(no2,j)-Moth_pos(no1,j));
        end
    end
    %% Boundary repair
    for i = 1:N
        Tp=Mhc(i,:)>ub;
        Tm=Mhc(i,:)<lb;
        Mhc(i,:)=(Mhc(i,:).*(~(Tp+Tm)))+ub.*Tp+lb.*Tm;
        if FEs<MaxFEs
            [MhcFitness(1,i),FEs,Mhc(i,:)] =sigle_evaluation(Mhc(i,:),dim,thdim,fname,Pxy,FEs);
        else
            break;
        end
    end
    ComFitness2=[Moth_fitness,MhcFitness];
    ComPositions2=[Moth_pos;Mhc];
    [ComFitness2,index]=sort(ComFitness2,'descend');
    for i =1:size(Moth_pos,1)
        Moth_pos(i,:)=ComPositions2(index(i),:);
        Moth_fitness(1,i)=ComFitness2(1,i);
    end
	if rand() <= Cr
		[ Moth_pos, FEs, Moth_fitness ] = GOL( Moth_pos,Moth_fitness,dim,ub,lb,thdim,fname,Pxy,FEs,type);  
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





