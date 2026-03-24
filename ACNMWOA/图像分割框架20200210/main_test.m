%% 对比原始类的算法运行
clear;
close all;
clc;
addpath([pwd, '\SubFunctions']);
addpath([pwd, '\statistics']);
addpath(genpath(pwd))
addpath(genpath(strcat(pwd, '\opt_algorithm\Advance')))
%% 自定义参数
runs = 30; %跑的次数
% runs = 10;
%% 一般运行分割阈值：10、15、20、25
partitionLevels = 30; %与优化算法维度成正比thdim=(partitionLevels-1)*2；取值2,3,4,5,6,10,15,20,25,30
% Population_size=((round(partitionLevels-1)/2))*10;      %种群大小
Population_size = 30;
% 图片大小
%% 肾病理
% ImageSize.width = 400;
% ImageSize.height = 240;
%% 乳腺癌
% ImageSize.width = 512;
% ImageSize.height = 512;
%% 遥感图像
% ImageSize.width = 512;
% ImageSize.height = 512;
%% COVID19
ImageSize.width = 512;
ImageSize.height = 400;

%% 评估模式和迭代模式
%%%%%%%%%%%%%%%%%%%%二选一  非选中项必须为‘0’%%%%%%%%%%%%%%%%%%
% MaxFEs=(partitionLevels-1)*1000;
%方案一：
MaxFEs = 0;
Iteration = 50;

%方案二：
% Iteration = 0;
% MaxFEs = 2000;

%% 评估方法
%可选方法：'BDE','PRI','VOI','GCE','SSIM','FSIM','RMSE','MSE','PSNR','NCC','AD','SC','MD','NAE'
%  evaluate_name={'BDE','PRI','VOI','GCE','SSIM','FSIM','RMSE','MSE','PSNR','NCC','AD','SC','MD','NAE'};
% evaluate_name={'BDE','PRI','VOI','GCE'};
evaluate_name = {'PSNR', 'SSIM', 'FSIM'};

%% 图片列表
%   image_filename={'100080.jpg','147091.jpg','170057.jpg','216081.jpg','241048.jpg',...
%                  '253027.jpg','260058.jpg','35008.jpg','35010.jpg','35049.jpg',...
%                  '38082.jpg','76053.jpg'};
% image_filename={'100080.jpg','101085.jpg','101087.jpg','106024.jpg','108005.jpg',...
%                 '113016.jpg','119082.jpg','12003.jpg','124084.jpg','138078.jpg',...
%                 '14037.jpg','147091.jpg','157055.jpg','161062.jpg','170057.jpg',...
%                 '175032.jpg','175043.jpg','19021.jpg','208001.jpg','216081.jpg',...
%                 '219090.jpg','220075.jpg','223061.jpg','231015.jpg','232038.jpg',...
%                 '24063.jpg','241004.jpg','241048.jpg','253027.jpg','260058.jpg',...
% 				  '291000.jpg','296059.jpg','304034.jpg','3063.jpg','3096.jpg',...
%                 '314016.jpg','317080.jpg','35008.jpg','35010.jpg','35049.jpg',...
%                 '35070.jpg','351093.jpg','37073.jpg','38082.jpg','38092.jpg',...
%                 '45096.jpg','61060.jpg','65010.jpg','66053.jpg','69020.jpg',...
%                 '76053.jpg','86068.jpg'};
% image_filename={'175043.jpg','219090.jpg','220075.jpg','231015.jpg','291000.jpg',...
%                 '304034.jpg','35049.jpg','37073.jpg','45096.jpg','69020.jpg'};
%% 肾病理
% image_filename = {'10A.jpg', '111A.jpg', '11A.jpg', '11C.jpg', '128C.jpg', '130B.jpg', '132B.jpg', ...
%                       '133D.jpg', '141B.jpg', '143B.jpg', '151D.jpg', '157C.jpg', '163C.jpg', '168D.jpg', '172A.jpg', ...
%                       '172B.jpg', '173D.jpg', '193A.jpg', '1A.jpg', '200A.jpg', '27A.jpg', '29D.jpg', '2A.jpg', '30C.jpg', ...
%                       '39C.jpg', '3A.jpg', '42B.jpg', '43C.jpg', '4A.jpg', '51B.jpg', '56D.jpg', '57D.jpg', '5A.jpg', '60B.jpg', ...
%                       '62B.jpg', '64C.jpg', '65B.jpg', '6A.jpg', '6C.jpg', '71C.jpg', '73B.jpg', '7A.jpg', '7D.jpg', '85B.jpg', ...
%                       '85D.jpg', '8A.jpg', '8C.jpg', '94D.jpg', '96B.jpg', '9A.jpg'};
%% 乳腺癌
% image_filename = {'rxa_1.JPG', 'rxa_2.JPG ', 'rxa_3.jpg ', 'rxa_4.JPG ', 'rxa_5.JPG ', ...
%                       'rxa_6.JPG ', 'rxa_7.JPG ', 'rxa_8.jpg ', 'rxa_9.JPG ', 'rxa_10.JPG', ...
%                       'rxa_11.jpg', 'rxa_12.jpg', 'rxa_13.JPG', 'rxa_14.jpg', 'rxa_15.JPG', ...
%                       'rxa_16.JPG', 'rxa_17.JPG', 'rxa_18.jpg', 'rxa_19.JPG', 'rxa_20.JPG', ...
%                       'rxa_21.JPG', 'rxa_22.jpg', 'rxa_23.JPG', 'rxa_24.JPG', 'rxa_25.JPG', ...
%                       'rxa_26.JPG', 'rxa_27.JPG', 'rxa_28.JPG', 'rxa_29.JPG', 'rxa_30.JPG', ...
%                       'rxa_31.JPG', 'rxa_32.JPG', 'rxa_33.jpg', 'rxa_34.jpg', 'rxa_35.JPG', ...
%                       'rxa_36.JPG', 'rxa_37.JPG', 'rxa_38.JPG', 'rxa_39.jpg', 'rxa_40.jpg', ...
%                       'rxa_41.JPG', 'rxa_42.JPG', 'rxa_43.jpg', 'rxa_44.JPG', 'rxa_45.jpg', ...
%                       'rxa_46.JPG', 'rxa_47.jpg', 'rxa_48.JPG', 'rxa_49.jpg', 'rxa_50.JPG', ...
%                       'rxa_51.jpg', 'rxa_52.jpg', 'rxa_53.jpg', 'rxa_54.jpg', 'rxa_55.JPG', ...
%                       'rxa_56.JPG', 'rxa_57.jpg', 'rxa_58.JPG', 'rxa_59.JPG', 'rxa_60.JPG'};
%% 遥感图像
% image_filename = {'rs_1.tiff ', 'rs_2.tiff ', 'rs_3.tiff ', 'rs_4.tiff ', 'rs_5.tiff ', ...
%                       'rs_6.tiff ', 'rs_7.tiff ', 'rs_8.tiff ', 'rs_9.tiff ', 'rs_10.tiff', ...
%                       'rs_11.tiff', 'rs_12.tiff', 'rs_13.tiff', 'rs_14.tiff', 'rs_15.tiff', ...
%                       'rs_16.tiff', 'rs_17.tiff', 'rs_18.tiff', 'rs_19.tiff', 'rs_20.tiff', ...
%                       'rs_21.tiff', 'rs_22.tiff', 'rs_23.tiff', 'rs_24.tiff', 'rs_25.tiff', ...
%                       'rs_26.tiff', 'rs_27.tiff', 'rs_28.tiff', 'rs_29.tiff', 'rs_30.tiff', ...
%                       'rs_31.tiff', 'rs_32.tiff', 'rs_33.tiff', 'rs_34.tiff', 'rs_35.tiff', ...
%                       'rs_36.tiff', 'rs_37.tiff', 'rs_38.tiff'};
%% COVID19
image_filename = {'01.jpg', '02.jpg', '03.jpg', '04.jpg', '05.jpg', '06.jpg', '07.jpg', '08.jpg', '09.jpg', '10.jpg', ...
                      '11.jpg', '12.jpg', '13.jpg', '14.jpg', '15.jpg', '16.jpg', '17.jpg', '18.jpg', '19.jpg', '20.jpg', ...
                      '21.jpg', '22.jpg', '23.jpg', '24.jpg', '25.jpg', '26.jpg', '27.jpg', '28.jpg', '29.jpg', '30.jpg', ...
                      '31.jpg', '32.jpg', '33.jpg', '34.jpg', '35.jpg', '36.jpg', '37.jpg', '38.jpg', '39.jpg', '40.jpg', ...
                      '41.jpg', '42.jpg', '43.jpg', '44.jpg', '45.jpg', '46.jpg', '47.jpg', '48.jpg', '49.jpg', '50.jpg', ...
                      '51.jpg', '52.jpg', '53.jpg', '54.jpg', '55.jpg'};

% image_filename = {'01_BC_G1_9057_40x_2.JPG'};
%% 对比算法列表
% method={'BLRIME','RIME', 'DE', 'IGWO', 'IWOA'};
%% ACNMWOA
% method = {'ACNMWOA', 'WOA', 'SSA', 'RIME', 'DE', 'PSO', 'IGWO', 'IWOA', 'CLSGMFO', 'LGCMFO'};
%% CCRIME
% method = {'CCRIME', 'RIME', 'SSA', 'HHO', 'WOA', 'SMA', 'LGCMFO', 'CLPSO', 'CLSGMFO', 'DECLS', 'WDE'};
%% ACDE
% method = {'ACDE', 'DE', 'SSA', 'HHO', 'WOA', 'SMA', 'LGCMFO', 'SCADE', 'CLSGMFO', 'DECLS', 'WDE'};
%% ACRIME
method = {'ACRIME', 'RIME', 'SSA', 'HHO', 'DE', 'SMA', 'LGCMFO', 'WDE', 'CLSGMFO', 'm_SCA', 'DECLS'};
%% 分割方法列表
% obj_name={'otsu2d','renyi2d','kapur2d'};
% obj_name={'kapur2d'};
obj_name = {'renyi2d'};

%% 准备工作
% 构建图片结构体
image_files.name = [];

for i = 1:size(image_filename, 2)
    image_files(i).name = cell2mat(image_filename(i));
end

% 构建路径
timestr = datestr(now, 'mm-dd-HH-MM');

if size(method, 2) == 1
    obj_path = [pwd, '\exp_result\HCode-', method{1}, '_Level', num2str(partitionLevels), '_result-set_' timestr];
else
    obj_path = [pwd, '\exp_result\HCode-', method{1}, '_', method{2}, '_Level', num2str(partitionLevels), '_result-set_' timestr];
end

%% 分割方法遍历
for obj = 1:size(obj_name, 2)
    %% 初始化
    objectiveFunction = obj_name{obj};
    dir_name = [obj_path, '\' objectiveFunction];
    xlsxname = [dir_name, '\FES-', timestr, '.xlsx'];
    NumofRecord = 40;
    measure1_name = {'BDE', 'PRI', 'VOI', 'GCE'};
    measure2_name = {'SSIM', 'FSIM', 'RMSE', 'MSE', 'PSNR', 'NCC', 'AD', 'SC', 'MD', 'NAE'};
    measure1_name = intersect(evaluate_name, measure1_name);
    measure2_name = intersect(evaluate_name, measure2_name);
    nLines = size(method, 2);
    basic_linestyles = cellstr(char('-', ':', '-.', '--'));
    basic_Markers = cellstr(char('o', 'x', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p', 'h', '.'));
    MarkerEdgeColors = hsv(nLines);
    linestyles = repmat(basic_linestyles, ceil(nLines / numel(basic_linestyles)), 1);
    Markers = repmat(basic_Markers, ceil(nLines / numel(basic_Markers)), 1);

    %% 开始主循环
    for i = 1:length(method)
        methodname = method{i};
        addpath([pwd '\opt_algorithm\' methodname]);
        pathName = [dir_name, '\', method{i}, '_result-set'];

        if ~exist('pathName', 'dir')
            mkdir(pathName);
        end

        disp(['*****************  ', 'Folds of ', method{i}, ' in ', objectiveFunction, '  *****************']);

        for im_len = 1:length(image_filename)
            Find_res(i, im_len).fitness = -inf;
            Find_res(i, im_len).imgname = {};
            Find_res(i, im_len).thresholds = [];
            Find_res(i, im_len).methodName = method{i};
            Find_res(i, im_len).run = 0;
            Find_res(i, im_len).FEs = 0;
            Find_res(i, im_len).Iter = 0;
            Find_res(i, im_len).AvgTime = 0;
        end

        parfor fold = 1:runs
            disp(['Fold', num2str(fold)]);
            [fd_methodata, fd_xlsimagename, fd_cg_curves, fd_FEcount, fd_image_info] = callmain(methodname, pathName, partitionLevels, MaxFEs, Population_size, image_files, NumofRecord, objectiveFunction, ImageSize, Iteration);
            fold_methodata(fold, :, :) = fd_methodata;
            fold_xlsimagename(fold, :) = fd_xlsimagename;
            fold_cg_curves(fold, :, :) = fd_cg_curves;
            fold_FEcount(fold, :) = fd_FEcount;
            fold_image_info(fold, :, :) = fd_image_info;
        end

        for j = 1:runs
            image_info(:, :) = fold_image_info(j, :, :);

            for im_sign = 1:length(image_info)
                % For mat
                sampleLabels = image_info(im_sign).mat_sampleLabels;
                superpixelLabels = image_info(im_sign).mat_superpixelLabels;
                refImage = image_info(im_sign).mat_refImage;
                dataFilename1 = image_info(im_sign).mat_dataFilename1;
                dataFilename2 = image_info(im_sign).mat_dataFilename2;
                save (dataFilename1, 'sampleLabels', 'superpixelLabels');
                save (dataFilename2, 'refImage');
            end

            sheetname = strcat(methodname, '_', num2str(j), 'run');
            methodata(:, :) = fold_methodata(j, :, :);
            xlsimagename = fold_xlsimagename(j, :);
            cg_curves(:, :) = fold_cg_curves(j, :, :);
            FEcount = fold_FEcount(j, :);

            all_curves(i, j, :, :) = cg_curves; % method  runs image NumofRecord
            xlsmethodata = {};

            for l = 1:size(methodata, 2) - 5
                xlsmethodata{1, l} = ['gBest', num2str(l)];
            end

            xlsmethodata{1, size(methodata, 2) - 4} = 'gbestvalue';
            xlsmethodata{1, size(methodata, 2) - 3} = 'FEscount';
            xlsmethodata{1, size(methodata, 2) - 2} = 'etime';
            xlsmethodata{1, size(methodata, 2) - 1} = 'MaxFEs';
            xlsmethodata{1, size(methodata, 2)} = 'Iteration';
            xlsimagename = xlsimagename';
            xlsmethodata = [xlsmethodata; methodata];
            xlswrite(xlsxname, xlsimagename, sheetname, 'A2');
            xlswrite(xlsxname, xlsmethodata, sheetname, 'B1');

            for im_sign = 1:length(image_filename)
                Find_res(i, im_sign).AvgTime = Find_res(i, im_sign).AvgTime + methodata{im_sign, size(methodata, 2) - 2};

                if methodata{im_sign, size(methodata, 2) - 4} > Find_res(i, im_sign).fitness
                    Find_res(i, im_sign).fitness = methodata{im_sign, size(methodata, 2) - 4};
                    %                     Find_res(i,im_sign).thresholds=xlsmethodata{im_sign+1,partitionLevels+1:size(methodata,2)-5};
                    Find_res(i, im_sign).thresholds = cell2mat(methodata(im_sign, partitionLevels + 1:size(methodata, 2) - 5));
                    Find_res(i, im_sign).imgname{1} = image_filename{im_sign};
                    Find_res(i, im_sign).FEs = methodata{im_sign, size(methodata, 2) - 3};
                    Find_res(i, im_sign).Iter = methodata{im_sign, size(methodata, 2)};
                    Find_res(i, im_sign).run = j;
                end

            end

            % For measuring the BDE, VoI, PRE, GCE
            if isempty(measure1_name) == 0
                seg_measure1(sheetname, xlsxname, xlsimagename, pathName, image_files, measure1_name, ImageSize);
            end

            % For measuring the other parameters, e.g. SSIM, RMSE, etc.
            if isempty(measure2_name) == 0
                seg_measure2(sheetname, xlsxname, xlsimagename, pathName, image_files, measure1_name, measure2_name);
            end

            close all;
        end

        %% 输出图像
        for im_sign = 1:length(image_filename)
            image_info(:, :) = fold_image_info(Find_res(i, im_sign).run, :, :);
            % 绘制灰度图
            figure;
            imshow(image_info(im_sign).gray);
            save_pic(image_info(im_sign).gray_filename);

            % 保存一维直方图
            figure;
            imhist(image_info(im_sign).hist);
            save_pic(image_info(im_sign).hist_filename);

            % 保存非均值滤波图像
            figure;
            imshow(image_info(im_sign).NLM);
            save_pic(image_info(im_sign).NLM_filename);

            % 保存二维直方图
            figure;
            mesh(image_info(im_sign).TwoD);
            save_pic(image_info(im_sign).TwoD_filename);

            % For displaying the segmented image
            figure;
            imshow(image_info(im_sign).seg_gray);
            save_pic(image_info(im_sign).seg_gray_filename);

            % For displaying the segmented image
            figure;
            imshow(image_info(im_sign).seg_color);
            save_pic(image_info(im_sign).seg_color_filename);

            % For seg image
            figure;
            save_pic_seg(image_info(im_sign).seg, image_info(im_sign).seg_filename, image_info(im_sign).seg_gBest, image_info(im_sign).seg_level);
            close all;
        end

    end

    %% 存阈值
    for im_sign = 1:length(image_filename)

        for method_sign = 1:length(method)
            xlsThData{method_sign + 1, 1} = Find_res(method_sign, im_sign).methodName;
            xlsThData{1, partitionLevels + 2} = 'Fitness';
            xlsThData{1, partitionLevels + 3} = 'FEs';
            xlsThData{1, partitionLevels + 4} = 'Iter';
            xlsThData{1, partitionLevels + 5} = 'AvgTime';

            for th_sign = 1:partitionLevels
                xlsThData{1, th_sign + 1} = ['Thresh', num2str(th_sign)];
                xlsThData{method_sign + 1, th_sign + 1} = Find_res(method_sign, im_sign).thresholds(th_sign);
            end

            xlsThData{method_sign + 1, partitionLevels + 2} = Find_res(method_sign, im_sign).fitness;
            xlsThData{method_sign + 1, partitionLevels + 3} = Find_res(method_sign, im_sign).FEs;
            xlsThData{method_sign + 1, partitionLevels + 4} = Find_res(method_sign, im_sign).Iter;
            xlsThData{method_sign + 1, partitionLevels + 5} = Find_res(method_sign, im_sign).AvgTime / runs;
        end

        Thresh_data = [dir_name, '\', 'Thresholds_data.xlsx'];
        Thresh_sheet = image_filename{im_sign};
        xlswrite(Thresh_data, xlsThData, Thresh_sheet);
    end

    %% 迭代曲线
    for l = 1:size(image_filename, 2)
        clf
        set(gcf, 'Position', [200, 200, 700, 400]);
        imgname = image_filename{1, l};

        for it = 1:size(method, 2)
            xls_curves(:, :) = all_curves(it, :, l, :); % method  runs image NumofRecord
            iter_data = [dir_name, '\', 'Iteration_data.xlsx'];
            folds_labels = [1:runs]';
            xls_method = repmat({method{it}}, [runs, 1]);
            xlsnum = runs * (it - 1) + 1;
            xlswrite(iter_data, xls_method, imgname, ['A', num2str(xlsnum)]);
            xlswrite(iter_data, folds_labels, imgname, ['B', num2str(xlsnum)]);
            xlswrite(iter_data, xls_curves, imgname, ['C', num2str(xlsnum)])
            yy(it, :) = mean(all_curves(it, 1:runs, l, :)); % method  runs image NumofRecord

            if MaxFEs ~= 0
                xx = [1:NumofRecord] * (MaxFEs / NumofRecord);
            else
                xx = [1:NumofRecord] * (Iteration / NumofRecord);
            end

        end

        for it = 1:size(method, 2)
            semilogy(xx, yy(it, :), [linestyles{it} Markers{it}], 'LineWidth', 1.5, 'Color', MarkerEdgeColors(it, :));
            hold on;
        end

        %去除文件名后缀
        flag = [];
        sign = 0;

        while ~strcmp(flag, '.')
            flag = imgname(end - sign);
            sign = sign + 1;
        end

        imgname = imgname(1:end - sign);
        title(imgname);
        set(gcf, 'color', 'white')
        set(gca, 'YScale', 'log', 'YLimMode', 'auto')

        if MaxFEs ~= 0
            xlabel('FEs');
        else
            xlabel('Iter');
        end

        ylabel('Best Value');
        legend(method, 'Location', 'Best');
        legend('boxoff')
        filename = [dir_name, '\', imgname, '_carve'];
        save_pic(filename);
    end

    close all;
    %     %% 数据统计
    %     [eval_names]=statistics(xlsxname,method,runs,dir_name);
    %     Ordehao(eval_names,dir_name);
    %     pValueToExcelhao(eval_names,method,xlsxname,runs,dir_name,image_filename);
    %     Friedman(xlsxname,method,runs,dir_name);
    %     Friedman1(xlsxname,method,runs,dir_name);
end

%% 数据统计
parfor obj = 1:size(obj_name, 2)
    tic;
    disp(['Fold', num2str(obj)]);
    objectiveFunction = obj_name{obj};
    dir_name = [obj_path, '\' objectiveFunction];
    xlsxname = [dir_name, '\FES-', timestr, '.xlsx'];
    DataStatistics(xlsxname, method, runs, dir_name, image_filename)
    toc;
end
