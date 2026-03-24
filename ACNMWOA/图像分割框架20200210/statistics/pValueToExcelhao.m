function pValueToExcelhao(eval_names,method,xlsfilename,runs,dir_name,image_filename)
%将第一个优化方法与其余四个优化函数进行威尔科克森符号秩检验（signrank），并将第一个优化方法与其他的均值和signrank比较
%   此处显示详细说明
% clear all;
% close all;
% clc;
% xlsfilename='test_seg.xlsx';
% runs=10;
%     
% dir_name=pwd;
% method={'CC_Replace_ACO_R','WOA','DE'};
% eval_names={'BDE','PRI','VOI','GCE','SSIM','FSIM','RMSE','PSNR','NCC','AD','MD','NAE'};
% 




readFilename=xlsfilename;
%% 读取所有评估数据
% for i=1:size(method,2)
%     for j=1:runs
%         sheetname=[method{i},'_',num2str(j),'run_eval'];
%         [~,~,sheet] = xlsread(readFilename,sheetname);
% %         eval([sheetname,'=sheet;']);
%         eval_val(i,j,:,:)=sheet;
%     end
% end
[eval_val,sheet] = readData(method,xlsfilename,image_filename,runs);


[method_num,run_num,img_num,evl_num]=size(eval_val);
img_num=img_num-1;
evl_num=evl_num-1;



for k=1:size(eval_names,2)
    
    writeFilename=[dir_name,'\',eval_names{1,k},'_statistics.xlsx'];

%     sheetname1=eval_names{1,k};
    
    overall_name=['overall',eval_names{1,k}];

    [~,~,rawdata] = xlsread(writeFilename,overall_name);
    num = size(rawdata);
    flag = 0;
    cur = 1;
    for i = 2:num(1)
        if flag == 0
            cur = cur + 1; 
        else
            break;
        end
        if i+1<=num(1) && ~strcmp(rawdata(i,1),rawdata(i+1,1))
            flag = 1;
        end
    end

    functioon_num = (num(1) - 1)/(cur - 1);         
    algorithm_num = cur - 1;         
    fold = runs;                  
    write_sheetname1=['pValue',eval_names{1,k}];
    sheet_name=write_sheetname1;
    tag=['B3:B',num2str(1+algorithm_num)]; 
    [~,label] = xlsread(writeFilename,overall_name,tag);
    label=label';
    xlswrite(writeFilename, label, sheet_name, 'B1')          %产生pValue表中的第一行

    for i=1:size(image_filename,2)
%         Function_name=['F',num2str(i)]; 
        Function_name={image_filename{1,i}};       
        xlswrite(writeFilename,Function_name , sheet_name, ['A',num2str(i+1)]);   %pValue表中第一列
    end


    for i=1:method_num
        for j=1:run_num
            tempdata(i,j,:)=eval_val(i,j,:,k+1);
    %         data(i+i*(j-1),:)=eval_val(i,j,:,k+1);
        end
    end
    




    % for i=1:functioon_num 
    %     temp = xlsread(readFilename,['F',num2str(i)],'','basic');
    %     data(:,i) = temp(:,end);
    % end
    Data1=[];
    for i=1:img_num
        for j=1:method_num
            ttData(:,j)=cell2mat(tempdata(j,:,i+1))';
        end
    %         x = data1;
        Data1 = [Data1, ttData]; 
    end
    data=[];
    for i=1:img_num
        Data2=[];
        for j=1:method_num
%             Data2=[Data2;Data1(:,(j-1)*img_num+1)];        
            Data2=[Data2;Data1(:,method_num*(i-1)+j)];
        end
        data(:,i) = Data2;
    end




    pvalue=zeros(functioon_num,(algorithm_num-1));
    for i=1:functioon_num
        for j=1:(algorithm_num-1)
            pvalue(i,j)=signrank(data(1:fold,i),data(j*fold+1:j*fold+fold,i));  %威尔科克森符号秩检验
        end
    end
    xlswrite(writeFilename, pvalue, sheet_name, 'B2')              %产生F1中的所有测试数据，从C1开始
    write_sheetname2=['result & pValue',eval_names{1,k}];
    % 合并pValue与result
    for i = 1 : algorithm_num - 1

        number = 'C'+(i*3-2) - 'A';
        numberOfA = floor(number / 26);
        numberOfVar = mod(number, 26);

        frontName = [];
        if numberOfA ~= 0
            for j = 1 : numberOfA
                frontName = [frontName, 'A'];
            end
        end
        behindName = char('A' + numberOfVar);


        xlswrite(writeFilename, pvalue(:,i), write_sheetname2,[[frontName, behindName],'2']);
        xlswrite(writeFilename, {'pvalue'}, write_sheetname2,[[frontName, behindName],'1']);
    end
    % 统计强弱
    xlswrite(writeFilename, {'CountOfBetter'}, write_sheetname2,['A',num2str(functioon_num + 4)]);
    xlswrite(writeFilename, {'CountOfWorse'}, write_sheetname2,['A',num2str(functioon_num + 5)]);
    xlswrite(writeFilename, {'CountOfEqual'}, write_sheetname2,['A',num2str(functioon_num + 6)]);
    tmpAlpha = ' ';
    tmpResult = repmat(tmpAlpha, functioon_num, algorithm_num - 1);
    tmpdata = xlsread(writeFilename, write_sheetname2);
    for i = 1 : algorithm_num - 1
        cntOfBetter = 0;
        cntOfWorse = 0;
        for j = 1 : functioon_num
            if tmpdata(j, i*3) < 0.05                   %威尔科克森检验统计
                 if tmpdata(j, 1) < tmpdata(j, i*3-1)        %均值比较
                    tmpResult(j, i) = '+';
                    cntOfBetter = cntOfBetter + 1;
                else
                    tmpResult(j, i) = '-';
                    cntOfWorse = cntOfWorse + 1;
                 end

            end
        end
        cntOfEqual = functioon_num - cntOfBetter - cntOfWorse;

        number = 'C'+(i*3-1) - 'A';
        numberOfA = floor(number / 26);
        numberOfVar = mod(number, 26);

        frontName = [];
        if numberOfA ~= 0
            for j = 1 : numberOfA
                frontName = [frontName, 'A'];
            end
        end
        behindName = char('A' + numberOfVar);


        xlswrite(writeFilename, tmpResult(:, i), write_sheetname2,[[frontName, behindName],'2']);
        xlswrite(writeFilename, cntOfBetter, write_sheetname2,[[frontName, behindName],num2str(functioon_num + 4)]);
        xlswrite(writeFilename, cntOfWorse, write_sheetname2,[[frontName, behindName],num2str(functioon_num + 5)]);
        xlswrite(writeFilename, cntOfEqual, write_sheetname2,[[frontName, behindName],num2str(functioon_num + 6)]);
    end

end