function [best_pos,Convergence_curve]=MSFOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
%%Zhang, Y., et al. (2016). "A novel multi-scale cooperative mutation Fruit Fly Optimization Algorithm." Knowledge-Based Systems 114: 24-35.
tic
%% 果蝇的相关参数
sizepop = SearchAgents_no;  %种群规模
FEs=0;
it=1;
w0 = 1;
alpha = 0.95;
M = 5;%%Scale number
W = max(ub)-min(lb);
%%Convergence_curve=zeros(1,Max_iter);
%% 初始果蝇群体位置
foa_option.sx = 1;
foa_option.sy = 1;
sx = foa_option.sx;

X_axis = sx*rands(1,dim);
X = zeros(sizepop, dim);
D = ones(sizepop, dim);
S = zeros(sizepop, dim);
Smell = zeros(sizepop, 1); 
delta = ones(M,1);
fit = zeros(M,1);
sub_group = zeros(sizepop/M,1,M);


%% 果蝇寻优开始
%利用嗅觉寻找食物
w = w0*alpha;
for i = 1:sizepop
    for j = 1:dim
        X(i,j) = X_axis(1,j) + w*(min(lb) + (max(ub) - min(lb)) * rand());
        D(i,j) = X(i,j);
        S(i,j) = D(i,j);
    end 
    FEs = FEs+1;
    Smell(i) = fobj(S(i,1:dim));
end
[bestSmell, bestindex] = min(Smell);
X_axis = X(bestindex,:);
best_pos = S(bestindex,:);
bestCVaccuarcy = bestSmell;

%% 果蠅迭代尋優
times = 0;
 while FEs<MaxFEs
    w = w0*alpha^FEs;
    for i = 1:sizepop
        for j = 1:dim
            X(i,j) = X_axis(1,j) + w*(min(lb) + (max(ub) - min(lb)) * rand());
            D(i,j) = X(i,j);
            S(i,j) = D(i,j);
        end 
        FEs = FEs+1;
        Smell(i) = fobj(S(i,1:dim));
    end

	[bestSmell, bestindex] = min(Smell);
    X_axis = X(bestindex,:);
    best_pos=S(bestindex,:);
    if bestSmell == bestCVaccuarcy
        times = times+1;
    end
    %%开始MSCM策略
    if times>=dim/2
        sort(Smell);
        for k=1:M   %%根据适应度分组
            if k==1
                for p = 1:M
                    sub_group(:,1,k) = Smell(p);
                end
            else
                for p = ((k-1)*M)+1:k*M 
                    sub_group(:,1,k) = Smell(p);
                end
            end
        end
        for k = 1:M   %%计算各子群适应度平均值
            fit(k) = sum(sub_group(:,1,k))/(sizepop/M);
        end
        fit_max = max(fit);
        fit_min = min(fit);
        for k = 1:M   %%计算各子群的标准差
            delta(k) = delta(k)*exp((M*fit(k)-sum(fit))/(fit_max-fit_min));
        end
        for k = 1:M    %%标准化
            if delta(k)>W/4
                delta(k) = abs((W/4)-delta(k));
            end
        end
        %%更新新种群的位置
        if FEs<MaxFEs
        FEs = FEs+6;
        f_min = fobj(X_axis(1,:)+rand(1)*delta(1,1));
        for i = 1:M
            if fobj(X_axis(1,:)+rand(1)*delta(i,1))<f_min
                f_min = fobj(X_axis(1,:)+rand(1)*delta(i,1));
            end
        end
        if (fobj(X_axis(1,:)+rand(1)*delta(M,1))<fobj(X_axis(1,:)+rand(1)*W))&&(fobj(X_axis(1,:)+rand(1)*delta(M,1))==f_min)
            X_axis(1,:) = X_axis(1,:)+rand(1)*delta(M,1);
        else
            X_axis(1,:) = X_axis(1,:)+rand(1)*W;
        end
        end
      
    end
	if bestSmell < bestCVaccuarcy
        bestCVaccuarcy = bestSmell;
        best_pos=S(bestindex,:);
%         X_axis = X(bestindex,:);
	end
     Convergence_curve(it)=(bestCVaccuarcy);
     it=it+1;
 end
toc
end

