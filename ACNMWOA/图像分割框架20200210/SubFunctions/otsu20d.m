function [h]= otsu2d(xI,Pxy)
    ind = Pxy == 0;
    ind = ind .* eps;
    Pxy = Pxy + ind;
    clear ind
    [M,N]=size(xI);
    for i=1:N/2
        eval([['s',num2str(i)],'=round(xI(M,i));']);
    end
    for i=N/2+1:N
        eval([['t',num2str(i-N/2)],'=round(xI(M,i));']);
    end
    
    %% ¼ÆËãw
    w0=sum(sum(Pxy(1: s1,1:t1)));
    for l=1:N/2-1 
        eval([['w',num2str(l)],'=0;']);
        i=(eval(['s',num2str(l)])+1):1:eval(['s',num2str(l+1)]);
        j=(eval(['t',num2str(l)])+1):1:eval(['t',num2str(l+1)]);
        eval([['w',num2str(l)],'=sum(sum(Pxy(i,j)));']);  
    end    
    eval([['w',num2str(N/2)],'=0;']);
    i=(eval(['s',num2str(N/2)])+1):256;
    j=(eval(['t',num2str(N/2)])+1):256;
    eval([['w',num2str(N/2)],'=sum(sum(Pxy(i,j)));']);

    %% ¼ÆËãux,uy
    ux = 0;
    uy = 0;
    i = 1:256;
    j = 1:256;
    ux=sum(sum(i * Pxy(i,j)));
    uy=sum(sum(j .* Pxy(i,j)));
    %% ¼ÆËãuix,uiy
    u0x=0;
    u0y=0;
    if w0~=0
        i = 1:s1;
        j = 1:t1;
        u0x=sum(sum(i * Pxy(i,j)))/w0;
        u0y=sum(sum(j .* Pxy(i,j)))/w0;
    end
    for l=1:N/2-1       
        eval([['u',num2str(l),'x'],'=0;']);
        eval([['u',num2str(l),'y'],'=0;']);
        if eval(['w',num2str(l)])~=0
            i = (eval(['s',num2str(l)])+1):eval(['s',num2str(l+1)]);
            j = (eval(['t',num2str(l)])+1):eval(['t',num2str(l+1)]);
            res=sum(sum(i * Pxy(i,j)))/eval(['w',num2str(l)]);
            eval([['u',num2str(l),'x'],'=res;']);
            res1=sum(sum(j .* Pxy(i,j)))/eval(['w',num2str(l)]);
            eval([['u',num2str(l),'y'],'=res1;']);
        end
    end
    eval([['u',num2str(N/2),'x'],'=0;']);
    eval([['u',num2str(N/2),'y'],'=0;']);
    if eval(['w',num2str(N/2)])~=0
        i = (eval(['s',num2str(N/2)])+1):256;
        j = (eval(['t',num2str(N/2)])+1):256;
        res=sum(sum(i * Pxy(i,j)))/eval(['w',num2str(N/2)]);
        eval([['u',num2str(N/2),'x'],'=res;']);
        res1=sum(sum(j .* Pxy(i,j)))/eval(['w',num2str(N/2)]);
        eval([['u',num2str(N/2),'y'],'=res1;']);
    end
    h=0;
    for i=1:N/2+1
        h=h+eval(['w',num2str(i-1)])*((eval(['u',num2str(i-1),'x'])-ux)^2+(eval(['u',num2str(i-1),'y'])-uy)^2);
    end

    if ((isinf(h))|| (h==0))
        h=0;
    end
end
