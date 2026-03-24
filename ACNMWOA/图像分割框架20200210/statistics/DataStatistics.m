function DataStatistics(xlsxname,method,runs,dir_name,image_filename)
%% 数据统计
[eval_names]=statistics(xlsxname,method,runs,dir_name,image_filename); 
Ordehao(eval_names,dir_name);
pValueToExcelhao(eval_names,method,xlsxname,runs,dir_name,image_filename);
Friedman(xlsxname,method,runs,dir_name,image_filename);
%Friedman1(xlsxname,method,runs,dir_name,image_filename);
end

