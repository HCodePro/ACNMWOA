
% Cai Z, Gu J, Luo J, Zhang Q, Chen H, Pan Z, Li Y, Li C. Evolving an Optimal Kernel Extreme Learning Machine by Using an Enhanced Grey Wolf Optimization Strategy. Expert Systems with Applications, 2019.

	


function [Fbest,Lbest,FE,MaxFEs,Convergence_curve,iter]=IterationIGWO(fname,SearchAgents_no,thdim,lb,ub,MaxFEs,Pxy,Iteration)
Fbest=-inf;
dim=thdim/2;

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,2*dim);
Alpha_score=-inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,2*dim);
Beta_score=-inf; %change this to -inf for maximization problems

beta_num=10;
omega_num=15;
it=1;
FEs = 0;% Loop counter

%Initialize the positions of search agents
Positions=random_initialization(SearchAgents_no,dim,ub,lb);
Convergence_curve=[];
% Convergence_curve=zeros(1,Max_iter);

% Main loop
while FEs < MaxFEs || it<=Iteration
    for i=1:size(Positions,1) 
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
		% Calculate objective function for each search agent
		[fitness(i),FEs,Positions(i,:)] =sigle_evaluation(Positions(i,:),dim,thdim,fname,Pxy,FEs);
		% Update Alpha
		if fitness(i)>Alpha_score
			Alpha_score=fitness(i); % Update alpha
			Alpha_pos=Positions(i,:);
			best_index=i;
		end
		
		if fitness(i)<Alpha_score && fitness(i)>Beta_score
			Beta_score=fitness(i); % Update beta
			Beta_pos=Positions(i,:);
		end
    end
    
    [~,index]=sort(fitness,'descend');
    Positions=Positions(index,:);
    % Update the Position of search agents
    a=2-it*((2)/Iteration); % a decreases linearly fron 2 to 0

    r=sqrt(sum((Alpha_pos-Beta_pos).^2));
    sa2=repmat(Alpha_pos,[beta_num,1])+((2*r).*rand([beta_num,2*dim])-r);
    sa2_finess=zeros(beta_num,1);   
    for k=1:beta_num
        Flag4ub=sa2(k,:)>ub;
        Flag4lb=sa2(k,:)<lb;
        sa2(k,:)=(sa2(k,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb; 
		% Calculate objective function for each search agent
		[sa2_finess(k),FEs,sa2(k,:)] =sigle_evaluation(sa2(k,:),dim,thdim,fname,Pxy,FEs);
    end
    [sa2_finess,maxIndex]=max(sa2_finess);
    if sa2_finess(1)>Alpha_score%Alpha 设置为局部更优解
        %sprintf('Alpha_score=%s,sa2_finess=%s',num2str(Alpha_score),num2str(sa2_finess(1)));
        Alpha_score=sa2_finess;
        Alpha_pos=sa2(maxIndex,:);
    end
    for i=1:SearchAgents_no-omega_num
        for j=1:size(Positions,2)
            r1=rand();
            r2=rand(); 
            A=2*a*r1-a; % Equation (3.3)
            C=2*r2; % Equation (3.4)

            D=abs(C*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            Positions(i,j)=Alpha_pos(j)-A*D; % Equation (3.5)-part 3  
%                 Positions(i,j)=Alpha_pos(j)+(Alpha_pos(j)-Positions(i,j));
        end
    end
    %update Deta
    Omega=random_initialization(omega_num,dim,ub,lb);
    Positions([SearchAgents_no-omega_num+1:SearchAgents_no],:)=Omega;
    Convergence_curve(it)=Alpha_score;  
	if Fbest<Alpha_score
        FE=FEs;
        iter=it;
    end
    Fbest=Alpha_score;
    Lbest=Alpha_pos;	
    it=it+1; 
end

end