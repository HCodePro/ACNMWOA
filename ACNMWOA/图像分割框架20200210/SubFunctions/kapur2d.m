function [h]= kapur2d(xI,Pxy)     
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
    
    %% º∆À„P
    P0=sum(sum(Pxy(1:s1,1:t1)));
    for l=1:N/2-1
        eval([['P',num2str(l)],'=0;']);
        eval([['w',num2str(l)],'=0;']);
        i=(eval(['s',num2str(l)])+1):1:eval(['s',num2str(l+1)]);
        j=(eval(['t',num2str(l)])+1):1:eval(['t',num2str(l+1)]);
        eval([['P',num2str(l)],'=sum(sum(Pxy(i,j)));']);
    end    
    eval([['P',num2str(N/2)],'=0;']);
    i=(eval(['s',num2str(N/2)])+1)+1:1:256;
    j=(eval(['t',num2str(N/2)])+1)+1:1:256;
    eval([['P',num2str(N/2)],'=sum(sum(Pxy(i,j)));']);
    
    %% º∆À„H
    H0=0;
    i=1:s1;
    j=1:t1;
    if P0>0
        H0=-sum(sum((Pxy(i,j)/P0).*log(Pxy(i,j)/P0)));
    end
    
    for l=1:N/2-1       
        eval([['H',num2str(l)],'=0;']);
        i=(eval(['s',num2str(l)])+1):eval(['s',num2str(l+1)]);
        j=(eval(['t',num2str(l)])+1):eval(['t',num2str(l+1)]);
        P=eval(['P',num2str(l)]);
        if P>0
            eval([['H',num2str(l)],'=-sum(sum((Pxy(i,j)/P).*log(Pxy(i,j)/P)));']);
        end
    end
    eval([['H',num2str(N/2)],'=0;']);
    i=(eval(['s',num2str(N/2)])+1):256;
    j=(eval(['t',num2str(N/2)])+1):256;
    P=eval(['P',num2str(N/2)]);
    if P>0
        eval([['H',num2str(N/2)],'=-sum(sum((Pxy(i,j)/P).*log(Pxy(i,j)/P)));']);
    end
    
    %% º∆À„h
    h=0;
    for i=1:N/2+1
        h=h+eval(['H',num2str(i-1)]);
    end    
   if ((isinf(h))|| (h==0))
        h=0;
    end
    
end

    

