% function [outputArg1,outputArg2] = Friedman(inputArg1,inputArg2)
%FRIEDMAN 此处显示有关此函数的摘要
%   此处显示详细说明
function Friedman1(xlsfilename,method,runs,dir_name,image_filename)
% 
% clear all;
% close all;
% clc;
% xlsfilename='test_seg.xlsx';
% 
% method={'WOA','CLS_Replace_ACO_R'};
% runs=2;

%% 读取所有评估数据
% for i=1:size(method,2)
%     for j=1:runs
%         sheetname=[method{i},'_',num2str(j),'run_eval'];
%         [~,~,sheet] = xlsread(xlsfilename,sheetname);
% %         eval([sheetname,'=sheet;']);
%         eval_val(i,j,:,:)=sheet;
%     end
% end
[eval_val,sheet] = readData(method,xlsfilename,image_filename,runs);


[method_num,run_num,img_num,evl_num]=size(eval_val);
img_num=img_num-1;
evl_num=evl_num-1;


%% 提取每一种优化方法利用所有图片评估runs次后的所有值
for k=1:evl_num
    for i=1:method_num
        for j=1:run_num
            tempdata(i,j,:)=eval_val(i,j,:,k+1);
    %         data(i+i*(j-1),:)=eval_val(i,j,:,k+1);
        end
    end

    meanranks = [];
    Data = [];
    for i=1:img_num
        for j=1:method_num
            ttData(:,j)=cell2mat(tempdata(j,:,i+1))';   
        end
%         x = data1;
        Data = [Data; ttData];   
        
        
%         [p(i),table,stat(i)] = friedman(x);             %Frideman秩方差分析
%         meanranks = [meanranks;stat(i).meanranks];
    end
    
    ans{1,1} = 'method_name';
    % ans{2,1} = 'p';
    ans{2,1} = 'mean_level';
    ans{3,1} = 'rank';

    x =  Data;
    [p,table,stat] = friedman(x);
    mean_level = stat.meanranks;
    
    
    sheetname1=cell2mat(sheet(1,k+1));
    
%     if sheetname1=='PRI'||strcmp(sheetname1,'SSIM')||sheetname1=='FSIM'||sheetname1=='PSNR'||sheetname1=='NCC'
    if strcmp(sheetname1,'PRI')||strcmp(sheetname1,'SSIM')||strcmp(sheetname1,'FSIM')||strcmp(sheetname1,'PSNR')||strcmp(sheetname1,'NCC')
        [~,index] = sort(mean_level,'descend');
    else
        [~,index] = sort(mean_level);
    end
    
    
    index2{k+1,1}=sheetname1;

    for c=1:size(index,2)
        index2{1,c+1}=cell2mat(method(c));
%         index2{l,c+1}=index1(l,c);
        index2{k+1,c+1}=index(1,c);
    end
    
%     [~,index] = sort(mean_level);
    
    
    
    
    for i = 1:method_num
        ans{1,i+1} = cell2mat(method(i)); 
    %     ans{2,i+1} = p(i);
        ans{2,i+1} = mean_level(i);
        if i~=1 && mean_level(index(i)) == mean_level(index(i-1))
            ans{3 ,index(i)+1} = ans{3 ,index(i-1)+1};
        else
            ans{3, index(i) + 1} = i;
        end

    end
    ans{1,method_num+2} = 'p';
    ans{2,method_num+2} = p;    
    
    
    sheetname2=[sheetname1,'_Fried'];
    xls_filename=[dir_name,'\',cell2mat(sheet(1,k+1)),'_statistics.xlsx'];
    xlswrite(xls_filename,ans,sheetname2);

    close all;

end
% index2{1,1}='eval_name';
% xlswrite(xlsfilename,index2,'eval_res_fried1');

