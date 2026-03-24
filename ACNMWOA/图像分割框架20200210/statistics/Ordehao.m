% function  Orderhao(filename)
function Ordehao(eval_names,dir_name)


%对不同优化方法对同一个函数进行均值排序
% clear all;
% close all;
% clc;
% filename='test_seg.xlsx';
% dir_name=pwd;
% eval_names={'BDE','PRI','VOI','GCE','SSIM','FSIM','RMSE','PSNR','NCC','AD','MD','NAE'};
											

for p=1:size(eval_names,2)
    xls_filename=[dir_name,'\',eval_names{1,p},'_statistics.xlsx'];
    
    sheetname1=eval_names{1,p};
    
    overall_name=['overall',eval_names{1,p}];
    
    
%   此处显示详细说明
    [~,~,rawdata] =  xlsread(xls_filename,overall_name);
    num = size(rawdata);
    flag = 0;
    cur = 1;
    for i = 2:num(1)
        if flag == 0
            fun(cur).name = cell2mat(rawdata(i,2));
            cur = cur + 1;
        else
            break;
        end
        if i+1<=num(1) && ~strcmp(rawdata(i,1),rawdata(i+1,1))
            flag = 1;
        end
    end
    funcnum = [num(1) - 1]/(cur - 1);%%函数的数目
    %%函数的下标获取为 （k - 1）*（cur - 1） +2 —— （k - 1）*（cur - 1) + cur
    for i = 1:funcnum
        k = 1;
        for j = (i - 1) * (cur - 1) + 2 : (i - 1)*(cur - 1) + cur
            f(k) = cell2mat(rawdata(j,5));      %取均值
            k = k + 1;
        end
        
        

    
        if strcmp(sheetname1,'PRI')||strcmp(sheetname1,'SSIM')||strcmp(sheetname1,'FSIM')||strcmp(sheetname1,'PSNR')||strcmp(sheetname1,'NCC')
            [~,index] = sort(f,'descend');
        else
            [~,index] = sort(f);
        end

        
        
%         [~,index] = sort(f);
        level = 1;
        for j = 1:cur - 1
            fun(index(j)).func(i) = level ;
            level = level + 1;
            if j>1 && f(index(j)) == f(index(j-1))
                fun(index(j)).func(i) = fun(index(j-1)).func(i);
            end
        end
    end
    for i = 1 : cur - 1
        fun(i).mean = mean(fun(i).func);
        temp(i) = fun(i).mean;
    end
    
    
    [~,index] = sort(temp);
    result{1,1} = 'imgName';
    for i = 1:funcnum
        result{i+1,1} = cell2mat(rawdata((i - 1)*(cur - 1) +2,1));
    end
    result{funcnum+2, 1} = 'mean_level';
    result{funcnum+3, 1} = '平均结果';
    for i = 1: cur-1
        result{1,i+1} = fun(i).name;
        for j = 1:funcnum
            result{j+1,i+1} = fun(i).func(j);
        end
        result{funcnum+2,i+1} = (fun(i).mean);
        result{funcnum+3,index(i)+1} =i;
        if i>1 && temp(index(i)) == temp(index(i-1))
            result{funcnum+3,index(i)+1} = result{funcnum+3,index(i - 1)+1};
        end
    end
    write_sheetname1=['result',eval_names{1,p}];
    xlswrite(xls_filename,result,write_sheetname1);
    write_sheetname2=['result & pValue',eval_names{1,p}];
    % 合并pValue与result
    algorithm_num = cur - 1;
    xlswrite(xls_filename,result(:,1),write_sheetname2,'A1');
    xlswrite(xls_filename,result(:,2),write_sheetname2,'B1');
    for i = 1 : algorithm_num - 1
        
        number = 'B'+(i*3-2) - 'A';
        numberOfA = floor(number / 26);
        numberOfVar = mod(number, 26);

        frontName = [];
        if numberOfA ~= 0
            for j = 1 : numberOfA
                frontName = [frontName, 'A'];
            end
        end
        behindName = char('A' + numberOfVar);
        
        xlswrite(xls_filename,result(:,i+2),write_sheetname2,[[frontName, behindName],'1']);
    end
    
end
% end
% close all;
% clear all;
% clc
% for i=1:50000
%    for j=1:10
%        rand_value(i,j)=rand;
%    end
% end
% xlswrite('rand_value',rand_value,'result_val');

