% function [Destination_fitness,bestPositions,Convergence_curve]=HGS(N,Max_iter,lb,ub,dim,fobj)
% function [bestPositions,Convergence_curve]=HGS(N,MaxFEs,lb,ub,dim,fobj)
function [Fbest,Lbest,FE,MaxFEs,Convergence_curve,iter]=HGS(fname,N,thdim,lb,ub,MaxFEs,Pxy,Iteration)
if Iteration~=0
    MaxFEs=Iteration;
end
dim=thdim/2;
Fbest=-inf;

bestPositions=zeros(1,2*dim);
tempPosition=zeros(N,2*dim);

Destination_fitness=-inf;%change this to -inf for maximization problems
Worstest_fitness=inf;
AllFitness = -inf*ones(N,1);%record the fitness of all positions
VC1 = ones(N,1);%record the variation control of all positions

weight3 = ones(N,2*dim);%hungry weight of each position
weight4 = ones(N,2*dim);%hungry weight of each position

%Initialize the set of random solutions
% X=initialization(N,dim,ub,lb);
X=random_initialization(N,dim,ub,lb);
Convergence_curve=[];
it=1; %Number of iterations
FEs=0;

hungry = zeros(1,size(X,1));%record the hungry of all positions
count=0;

% Main loop
while  FEs <= MaxFEs
    if Iteration~=0
        FEs=it;
    end
    
    VC2 = 0.03; %The variable of variation control 

    sumHungry = 0;%record the sum of each hungry 
    
    %sort the fitness
    for i=1:size(X,1)
        % Check if solutions go outside the search space and bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
%         AllFitness(i) = fobj(X(i,:));
        if FEs < MaxFEs
            [AllFitness(i),FEs,X(i,:)] =sigle_evaluation(X(i,:),dim,thdim,fname,Pxy,FEs,Iteration);
        else
            break;
        end
        FEs=FEs+1;
    end
    [AllFitnessSorted,IndexSorted] = sort(AllFitness,'descend');
    bestFitness = AllFitnessSorted(1);
    worstFitness = AllFitnessSorted(size(X,1));
    
    %update the best fitness value and best position
    if bestFitness > Destination_fitness
        bestPositions=X(IndexSorted(1),:);
        Destination_fitness = bestFitness;
        count=0;
    end
    
    if worstFitness < Worstest_fitness
        Worstest_fitness = worstFitness;
    end
    
    for i = 1:size(X,1)
         %calculate the variation control of all positions
         VC1(i) = sech(abs(AllFitness(i)-Destination_fitness));    
         %calculate the hungry of each position
        if Destination_fitness == AllFitness(i)
            hungry(1,i) = 0;
            count = count+1;
            tempPosition(count,:)=X(i,:);
        else
            temprand = rand();
            c = (AllFitness(i)-Destination_fitness)/(Worstest_fitness-Destination_fitness)*temprand*2*(ub-lb);
            if c<100
                b=100*(1+temprand);
            else
                b=c;
            end   
            hungry(1,i) = hungry(1,i)+ max(b); 
            sumHungry = sumHungry + hungry(1,i);
        end
    end 
    
    %calculate the hungry weight of each position
    for i=1:size(X,1)
        for j=2:size(X,2)
                weight3(i,j) = (1-exp(-abs(hungry(1,i)-sumHungry)))*rand()*2;
                if rand()<VC2
                    weight4(i,j) = hungry(1,i)*size(X,1)/sumHungry*rand();
                else
                    weight4(i,j) = 1;
                end
        end
        
    end
    
    
    % Update the Position of search agents
    shrink=2*(1-FEs/MaxFEs); % a decreases linearly fron 2 to 0
    for i=1:size(X,1)
        if rand<VC2
            X(i,:) = X(i,j)*(1+randn(1));
        else
             if count>=1
                 A = randi([1,count]);
                 for j=1:size(X,2)
                     r = rand();
                     vb = 2*shrink*r-shrink;%[-a,a]
                     % Moving based on the bestPosition
                     % The transformation range is controlled by weight3,bestPositions and X
                     if r>VC1(i)
                         X(i,j) = weight4(i,j)*tempPosition(A,j)+vb*weight3(i,j)*abs(tempPosition(A,j)-X(i,j));
                     else
                         X(i,j) = weight4(i,j)*tempPosition(A,j)-vb*weight3(i,j)*abs(tempPosition(A,j)-X(i,j));
                     end
                 end
             end
        end
    end
    
    Convergence_curve(it)=Destination_fitness;
    if Fbest<Destination_fitness
        FE=FEs;
        iter=it;
    end
    Fbest=Destination_fitness;
    Lbest=bestPositions;
    it=it+1;
end
% toc
end



