function [best_pos,Convergence_curve]=CEFOA(SearchAgents_no,MaxFEs,lb,ub,dim,fobj)
tic
%%Han, X., et al. (2018). "Novel fruit fly optimization algorithm with trend search and co-evolution." Knowledge-Based Systems 141: 1-17.
%% 果蝇的相关参数
sizepop = SearchAgents_no;  %种群规模
FEs=0;
it=1;

%%Convergence_curve=zeros(1,Max_iter);
%% 初始果蝇群体位置
foa_option.sx = 1;
foa_option.sy = 1;
foa_option.ax = 20;
foa_option.ay = 20;
foa_option.bx = 10;
foa_option.by = 10;
sx = foa_option.sx;
sy = foa_option.sy;
ax = foa_option.ax;
ay = foa_option.ay;
bx = foa_option.bx;
by = foa_option.by;
X_axis = sx*rands(1,dim);
Y_axis = sy*rands(1,dim);
X = zeros(sizepop, dim);
Y = zeros(sizepop, dim);
D = ones(sizepop, dim);
S = zeros(sizepop, dim);
Smell = zeros(sizepop, 1); 

f_ijx = zeros(sizepop,dim);
f_ijy = zeros(sizepop,dim);
c_ijx = zeros(sizepop,dim);
c_ijy = zeros(sizepop,dim);
k = zeros(sizepop,1);

%% 果蝇寻优开始
%利用嗅觉寻找食物
for i = 1:sizepop
    tempx = ax*rand() - bx;
    tempy = ay*rand() - by;
    for j = 1:dim
        X(i,j) = X_axis(1,j) + tempx;
        Y(i,j) = Y_axis(1,j) + tempy;
        D(i,j) = (X(i,j).^2 + Y(i,j).^2).^0.5;
        S(i,j) = 1./D(i,j);
    end
	Flag4ub=S(i,:)>ub;
    Flag4lb=S(i,:)<lb;
    S(i,:)=(S(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;  
    FEs=FEs+1;
    Smell(i) = fobj(S(i,1:dim));
end

[bestSmell, bestindex] = min(Smell);
X_axis = X(bestindex,:);
Y_axis = Y(bestindex,:);
best_pos=S(bestindex,:);
bestCVaccuarcy = bestSmell;

%% 果蠅迭代尋優
    c_ix = zeros(sizepop,1);   
    c_iy = zeros(sizepop,1);
    r_ix = ones(sizepop,1);
    r_iy = ones(sizepop,1);
 while FEs<MaxFEs
    for i = 1:sizepop              %%sub_algorithm 1
        for j = 1:dim
            f_ijx(i,j) = (Smell(i)-bestCVaccuarcy)./(X(i,j)-X_axis(1,j));
            f_ijy(i,j) = (Smell(i)-bestCVaccuarcy)./(Y(i,j)-Y_axis(1,j));
            if f_ijx(i,j)>0
                c_ijx(i,j) = 1;
                c_ix(i) = c_ix(i) + 1;
            else
                c_ijx(i,j) = 0;
            end
            if f_ijy(i,j)>0
                c_ijy(i,j) = 1;
                c_iy(i) = c_iy(i) + 1;
            else
                c_ijy(i,j) = 0; 
            end
        end
        k(i) = dim*log((Smell(i)/sum(Smell))/(FEs/MaxFEs));
        if c_ix(i)>k(i)
            c_ix(i) = -c_ix(i);
        end
        if c_iy(i)>k(i)
            c_iy(i) = -c_iy(i);
        end
        r_ix(i) = r_ix(i)*exp(c_ix(i));
        r_iy(i) = r_iy(i)*exp(c_iy(i));
        X_axis(1,:) = X_axis(1,:)+(ub-lb)*r_ix(i)*randn();
        Y_axis(1,:) = Y_axis(1,:)+(ub-lb)*r_iy(i)*randn();
    end
    for i = 1:sizepop    %%compute initial smell concentration with(X_axis,Y_axis);
    tempx = ax*rand() - bx;
    tempy = ay*rand() - by;
    for j = 1:dim
        X(i,j) = X_axis(1,j) + tempx;
        Y(i,j) = Y_axis(1,j) + tempy;
        D(i,j) = (X(i,j).^2 + Y(i,j).^2).^0.5;
        S(i,j) = 1./D(i,j);
    end  
    FEs = FEs+1;
    Smell(i) = fobj(S(i,1:dim));
    end
    

	[bestSmell, bestindex] = min(Smell);    %%Sub_algorithm 2
    [worstSmell,worstindex] = max(Smell);
    X_axis = X(bestindex,:);
    Y_axis = Y(bestindex,:);
    X_axis_worst = X(worstindex,:);
    Y_axis_worst = Y(worstindex,:);
    dis = (sum(X_axis.^2-X_axis_worst.^2)+sum(Y_axis.^2-Y_axis_worst.^2)).^0.5;
    valt = bestSmell-worstSmell;
    beta = 1/(exp(-(dis/valt)));
    if beta > 0.5
        X_axis = X_axis_worst;
        Y_axis = Y_axis_worst;
    end
    for i = 1:sizepop
        tempx = ax*rand() - bx;
        tempy = ay*rand() - by;
        for j = 1:dim
            X(i,j) = X_axis(1,j)+tempx;
            Y(i,j) = Y_axis(1,j)+tempy;
            D(i,j) = (X(i,j).^2 + Y(i,j).^2).^0.5;
            S(i,j) = 1./D(i,j);
        end
        FEs = FEs+1;
        Smell(i) = fobj(S(i,:));
    end
    [bestSmell, bestindex] = min(Smell);
    if bestSmell<bestCVaccuarcy
        bestCVaccuarcy = bestSmell;
        X_axis = X(bestindex,:);
        Y_axis = Y(bestindex,:);
        best_pos=S(bestindex,:);
    end
    
     Convergence_curve(it)=(bestCVaccuarcy);
     it=it+1;
 end
toc
end

