function [h]= renyi2d(xI,Pxy)    
    [M,N]=size(xI);
    alpha=0.5;
    ind = Pxy == 0;
    ind = ind .* eps;
    Pxy = Pxy + ind;
    clear ind
    %% 计算s,t
    for i=1:N/2
        eval([['s',num2str(i)],'=round(xI(M,i));']);
    end
    for i=N/2+1:N
        eval([['t',num2str(i-N/2)],'=round(xI(M,i));']);
    end
    %% 计算P
    P0=sum(sum(Pxy(1:s1,1:t1)));    
    for l=1:N/2-1       
        eval([['P',num2str(l)],'=0;']);
        i=(eval(['s',num2str(l)])+1):1:eval(['s',num2str(l+1)]);
        j=(eval(['t',num2str(l)])+1):1:eval(['t',num2str(l+1)]);
        eval([['P',num2str(l)],'=sum(sum(Pxy(i,j)));']);     
    end
    eval([['P',num2str(N/2)],'=0;']);
    i=(eval(['s',num2str(N/2)])+1):1:256;
    j=(eval(['t',num2str(N/2)])+1):1:256;
    eval([['P',num2str(N/2)],'=sum(sum(Pxy(i,j)));']);  
    %% 计算H
    H0=0;
    H0=sum(sum((Pxy(1:s1,1:t1)./P0).^alpha));
    H0=(1/(1-alpha))*log(H0+eps);
    for l=1:N/2-1      
        eval([['H',num2str(l)],'=0;']);
        i=(eval(['s',num2str(l)])+1):eval(['s',num2str(l+1)]);
        j=(eval(['t',num2str(l)])+1):eval(['t',num2str(l+1)]);
        res=eval(['P',num2str(l)]);
        eval([['H',num2str(l)],'=sum(sum((Pxy(i,j)./res).^alpha));']); 
        res=eval(['H',num2str(l)]);
        eval([['H',num2str(l)],'=(1/(1-alpha))*log(res+eps);']); 
    end
    eval([['H',num2str(N/2)],'=0;']);
    i=(eval(['s',num2str(N/2)])+1):256;
    j=(eval(['t',num2str(N/2)])+1):256;
    res=eval(['P',num2str(N/2)]);
    eval([['H',num2str(N/2)],'=sum(sum((Pxy(i,j)./res).^alpha));']); 
    res=eval(['H',num2str(N/2)]);
    eval([['H',num2str(N/2)],'=(1/(1-alpha))*log(res+eps);']);  
    %% 计算h
    h=0;
    for i=1:N/2+1
        h=h+eval(['H',num2str(i-1)]);
    end
   if ((isinf(h))|| (h==0))
        h=0;
    end
end