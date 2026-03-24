function [eval_names]=statistics(xlsfilename,method,runs,dir_name,image_filename) 
% clear;
% close all;
% clc;
% xlsfilename='test_seg.xlsx';
% method={'CC_Replace_ACO_R','WOA','DE'};
% runs=10;

%% 读取所有评估数据
% for i=1:size(method,2)
%     for j=1:runs
%         sheetname=[method{i},'_',num2str(j),'run_eval'];
%         [~,~,sheet] = xlsread(xlsfilename,sheetname);
%         %% 过滤图片
% %         it=1;
%         t=1;
%         for m=2:size(sheet(:,1),1)
%             flag=0;
%             for n=1:size(image_filename,2)
%                 if strcmp(sheet{m,1},image_filename{n})
%                     flag=1;
%                 end
%             end
% %             if flag==1
% %                 temp{it,1}=sheet{m,1};
% %                 it=it+1;
% %             end
%             if flag==1
%                 selected(t)=m;
%                 t=t+1;
%             end
%             flag=0;
%         end
%         selected=[1,selected];
% %         sheet1=sheet;
%         sheet=sheet(selected,:);
%         selected=[];
%         %%
%         
% %         eval([sheetname,'=sheet;']);
%         eval_val(i,j,:,:)=sheet;
%     end
% end
[eval_val,sheet] = readData(method,xlsfilename,image_filename,runs);

[method_num,run_num,img_num,evl_num]=size(eval_val);
img_num=img_num-1;
evl_num=evl_num-1;
eval_names={};
for k=1:evl_num
    
    eval_names{k}=sheet{1,k+1};
    startLineNum=2;
    %% 提取每一种优化方法利用所有图片评估runs次后的所有值
    for i=1:method_num
        for j=1:run_num
            tempdata(i,j,:)=eval_val(i,j,:,k+1);
    %         data(i+i*(j-1),:)=eval_val(i,j,:,k+1);
        end
    end
    

    for imgNum=1:img_num
        for i=1:method_num
            data(i,:)=cell2mat(tempdata(i,:,imgNum+1));
    %         data(i,:)=tempdata(i,:,imgNum+1);
        end
        methodNum=size(method,2);
        medName_labels2=method';
        statistic_values=zeros(methodNum,4);
        for it=1:methodNum
            statistic_values(it,:)=[max(data(it,:)),min(data(it,:)),mean(data(it,:)),std(data(it,:))];
        end
        xls_filename=[dir_name,'\',cell2mat(sheet(1,k+1)),'_statistics.xlsx'];
        
        
        xlsname=[];
        xlsname=[xlsname,'overall',sheet{1,k+1}];
        imNum_label=repmat({sheet{imgNum+1,1}},methodNum,1);
        xlswrite(xls_filename, imNum_label, xlsname, ['A',num2str(startLineNum)]);
        xlswrite(xls_filename, medName_labels2, xlsname, ['B',num2str(startLineNum)]);
        xlswrite(xls_filename, statistic_values, xlsname, ['C',num2str(startLineNum)]);
        startLineNum=startLineNum+methodNum;
    end
    title={'img','method','max','min','mean','std'};
    
    
    xlswrite(xls_filename,title, xlsname, 'A1');
end