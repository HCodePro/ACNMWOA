function [h]= otsu2d(xI,Pxy)
%     clear;
%     clc;
%     load('Pxy.mat');
%     load('x1.mat');
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
    %% 计算背景区的概率
    w0=0;
    i=1:s1;
    j=1:t1;
    w0=sum(sum(Pxy(i,j)));
    
    
    
%     for i=1:s1
%         for j=1:t1
%             w0=w0+Pxy(i,j);
%         end
%     end
    %% 计算目标区概率
    for l=1:N/2-1
        eval([['w',num2str(l)],'=0;']);
        i=eval(['s',num2str(l)])+1:eval(['s',num2str(l+1)]);
        j=eval(['t',num2str(l)])+1:eval(['t',num2str(l+1)]);
        eval([['w',num2str(l)],'=sum(sum(Pxy(i,j)));']);
%         res=0;
%         for i=eval(['s',num2str(l)])+1:eval(['s',num2str(l+1)])
%             for j=eval(['t',num2str(l)])+1:eval(['t',num2str(l+1)])
%                 res=res+Pxy(i,j);
%             end
%         end
%         eval([['w',num2str(l)],'=res;']);
    end
    
    eval([['w',num2str(N/2)],'=0;']);
    i=eval(['s',num2str(N/2)])+1:256;
    j=eval(['t',num2str(N/2)])+1:256;
    eval([['w',num2str(N/2)],'=sum(sum(Pxy(i,j)));']);
    
%     res=0;
%     for i=eval(['s',num2str(N/2)])+1:256
%         for j=eval(['t',num2str(N/2)])+1:256
%             res=res+Pxy(i,j);
%         end
%     end
%     eval([['w',num2str(N/2)],'=res;']);
    
    %% 计算u
    ux=0;
    uy=0;
    i=1:256;
    j=1:256;
    ux=sum(sum(i*Pxy(i,j)));
    uy=sum(sum(j.*Pxy(i,j)));
%     for i=1:256
%         for j=1:256
%             ux=ux+i*Pxy(i,j);
%             uy=uy+j*Pxy(i,j);
%         end
%     end
    
    u0x=0;
    u0y=0;
    if w0~=0
        i=1:s1;
        j=1:t1;
        u0x=sum(sum(i*Pxy(i,j)))/w0;
        u0y=sum(sum(j.*Pxy(i,j)))/w0;
        
%         for i=1:s1
%             for j=1:t1
%                 u0x=u0x+i*Pxy(i,j)/w0;
%                 u0y=u0y+j*Pxy(i,j)/w0;
%             end
%         end
    end
    
    
    for l=1:N/2-1
        eval([['u',num2str(l),'x'],'=0;']);
        eval([['u',num2str(l),'y'],'=0;']);
        if eval(['w',num2str(l)])~=0
            i=eval(['s',num2str(l)])+1:eval(['s',num2str(l+1)]);
            j=eval(['t',num2str(l)])+1:eval(['t',num2str(l+1)]);
            resx=sum(sum(i*Pxy(i,j)))/eval(['w',num2str(l)]);
            resy=sum(sum(j.*Pxy(i,j)))/eval(['w',num2str(l)]);
            eval([['u',num2str(l),'x'],'=resx;']);
            eval([['u',num2str(l),'y'],'=resy;']);          
%             resx=0;
%             resy=0;
%             for i=eval(['s',num2str(l)])+1:eval(['s',num2str(l+1)])
%                 for j=eval(['t',num2str(l)])+1:eval(['t',num2str(l+1)])
%                     resx=resx+i*Pxy(i,j);
%                     resy=resy+j*Pxy(i,j);
%                 end
%             end
%             resx=resx/eval(['w',num2str(l)]);
%             resy=resy/eval(['w',num2str(l)]);
%             eval([['u',num2str(l),'x'],'=resx;']);
%             eval([['u',num2str(l),'y'],'=resy;']);
        end
    end
    
    eval([['u',num2str(N/2),'x'],'=0;']);
    eval([['u',num2str(N/2),'y'],'=0;']);
    if eval(['w',num2str(N/2)])~=0
        i=eval(['s',num2str(N/2)])+1:256;
        j=eval(['t',num2str(N/2)])+1:256;
        resx=sum(sum(i*Pxy(i,j)))/eval(['w',num2str(N/2)]);
        resy=sum(sum(j.*Pxy(i,j)))/eval(['w',num2str(N/2)]);
        eval([['u',num2str(N/2),'x'],'=resx;']);
        eval([['u',num2str(N/2),'y'],'=resy;']);    
%         resx=0;
%         resy=0;
%         for i=eval(['s',num2str(N/2)])+1:256
%             for j=eval(['t',num2str(N/2)])+1:256
%                 resx=resx+i*Pxy(i,j);
%                 resy=resy+j*Pxy(i,j);
%             end
%         end
%         resx=resx/eval(['w',num2str(N/2)]);
%         resy=resy/eval(['w',num2str(N/2)]);
%         eval([['u',num2str(N/2),'x'],'=resx;']);
%         eval([['u',num2str(N/2),'y'],'=resy;']);
    end
    
    %% 计算h
    h=0;
    for i=1:N/2+1
        h=h+eval(['w',num2str(i-1)])*((eval(['u',num2str(i-1),'x'])-ux)^2+(eval(['u',num2str(i-1),'y'])-uy)^2);
    end
    if ((isinf(h))|| (h==0))
        h=0;
    end