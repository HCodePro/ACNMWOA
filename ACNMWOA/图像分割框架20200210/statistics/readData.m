function [eval_val,sheet] = readData(method,xlsfilename,image_filename,runs)
    for i=1:size(method,2)
        for j=1:runs
            sheetname=[method{i},'_',num2str(j),'run_eval'];
            [~,~,sheet] = xlsread(xlsfilename,sheetname);
            %% ¹ýÂËÍ¼Æ¬
    %         it=1;
            t=1;
            for m=2:size(sheet(:,1),1)
                flag=0;
                for n=1:size(image_filename,2)
                    if strcmp(sheet{m,1},image_filename{n})
                        flag=1;
                    end
                end
    %             if flag==1
    %                 temp{it,1}=sheet{m,1};
    %                 it=it+1;
    %             end
                if flag==1
                    selected(t)=m;
                    t=t+1;
                end
                flag=0;
            end
            selected=[1,selected];
    %         sheet1=sheet;
            sheet=sheet(selected,:);
            selected=[];
            %%
    %         eval([sheetname,'=sheet;']);
            eval_val(i,j,:,:)=sheet;
        end
    end
end

