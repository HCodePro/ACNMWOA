function [X] = random_initialization(N,dim,ub,lb)
%RANDOM_INITIALIZATION 此处显示有关此函数的摘要
%   此处显示详细说明
    dim1=dim;
    X1=initialization(dim1,N,ub,lb); 
    for si=1:size(X1,1)
       X1(si,:)=sort(X1(si,:)); 
    end
    dim2=dim;
    X2=initialization(dim2,N,ub,lb); 
    for si=1:size(X2,1)
       X2(si,:)=sort(X2(si,:)); 
    end
    X=[X1 X2];
end

