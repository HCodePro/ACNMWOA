clear all
close all
clc
tic
% add function
% addpath([pwd,'/CCMWOA']);
%  addpath([pwd,'/Origin']);
%  addpath([pwd,'/statics']);
%  addpath([pwd,'/subFunction']);
%  addpath([pwd,'/cec05']);
addpath([pwd, '/HCode_Algorithm']);
addpath(genpath(strcat(pwd, '/HCode_Algorithm')));
addpath(genpath(pwd));
%%parm setting
% option setting
SearchAgents_no = 30; % Number of search agents
NumofRecord = 40;
dim = 30; % 30 50 100
MaxFEs = dim * 5000 * 2;

%% basic algoritem
%% advanced algorithm
%% ACNMWOA
algorithmName = {'ACNMWOA', 'WDE', 'LSHADE', 'SADE', 'SHADE', 'LWOA', 'IWOA', 'CLACO', 'LSHADE_cnEpSi', 'CLPSO'};

Flod = 30;

%% perpare xlsfile
timestr = datestr(now, 'mm-dd-HH-MM');
dirname = ['exp_result/HCode-Advance-', algorithmName{1}, '-', algorithmName{2}, '-', timestr, '-Dim', num2str(dim)];
mkdir(dirname);
filename = [dirname, '/Function23-', timestr];

te = {'F', 'Algrithm', 'max', 'min', 'mean', 'std'};
xlsfilename = [filename, '.xlsx'];
xlswrite(xlsfilename, te, 'overall')
startLineNum = 2; % the startLineNum of overall sheet to write data

algrithmNum = size(algorithmName, 2);
% generate the space of linestyles, MarkerEdgeColors,Markers
nLines = algrithmNum;
basic_linestyles = cellstr(char('-', ':', '-.', '--'));
basic_Markers = cellstr(char('o', 'x', '+', '*', 's', 'd', 'v', '^', '<', '>', 'p', 'h', '.'));
MarkerEdgeColors = hsv(nLines);
linestyles = repmat(basic_linestyles, ceil(nLines / numel(basic_linestyles)), 1);
Markers = repmat(basic_Markers, ceil(nLines / numel(basic_Markers)), 1);
%% =======uncomment for specify function===============
%% CEC14:24-53;  CEC17:107-136; CEC19:137-146; CEC20:147-156
Function_name_all = {'F107', 'F109', 'F110', 'F111', 'F112', 'F113', 'F114', 'F115', 'F116', 'F117', 'F118', 'F119', 'F120', 'F121', 'F122', 'F123', 'F124', 'F125', 'F126', 'F127', 'F128', 'F129', 'F130', 'F131', 'F132', 'F133', 'F134', 'F135', 'F136'};

for funcNum = 1:size(Function_name_all, 2)
    dim = 30;
    Function_name = Function_name_all{funcNum};
    [lb, ub, dim, fhd] = Get_Functions(Function_name, dim);
    disp(['----------------', Function_name, '--------------------']);
    Function_name = ['F', num2str(funcNum)];
    % benchmark function
    cg_curves = zeros(algrithmNum, Flod, NumofRecord);

    parfor cflod = 1:Flod
        display(['flod', num2str(cflod)]);

        for cnum = 1:algrithmNum
            alg_fhd = str2func(algorithmName{cnum});
            [~, cg_curve] = alg_fhd(SearchAgents_no, MaxFEs, lb, ub, dim, fhd);
            cg_curves(cnum, cflod, :) = MySampling(cg_curve, NumofRecord);
        end

    end

    %% write data to excel sheet
    algName_labels = cell(algrithmNum * Flod, 1);
    all_curves = zeros(algrithmNum * Flod, NumofRecord);
    folds_labels = repmat(int32(1:Flod)', [algrithmNum, 1]);

    for it = 1:algrithmNum
        algName_labels((it - 1) * Flod + 1:(it - 1) * Flod + Flod) = algorithmName(it);
        all_curves((it - 1) * Flod + 1:(it - 1) * Flod + Flod, :) = cg_curves(it, :, :);
    end

    xlswrite(xlsfilename, algName_labels, Function_name, 'A1')
    xlswrite(xlsfilename, folds_labels, Function_name, 'B1')
    xlswrite(xlsfilename, all_curves, Function_name, 'C1')
    %write overall
    algName_labels2 = algorithmName';
    statistic_values = zeros(algrithmNum, 4);

    for it = 1:algrithmNum
        statistic_values(it, :) = [max(cg_curves(it, :, end)), min(cg_curves(it, :, end)), mean(cg_curves(it, :, end)), std(cg_curves(it, :, end))];
    end

    funcNum_label = repmat({Function_name}, algrithmNum, 1);
    xlswrite(xlsfilename, funcNum_label, 'overall', ['A', num2str(startLineNum)])
    xlswrite(xlsfilename, algName_labels2, 'overall', ['B', num2str(startLineNum)])
    xlswrite(xlsfilename, statistic_values, 'overall', ['C', num2str(startLineNum)])
    startLineNum = startLineNum + algrithmNum;
    %% plot cg_curveline
    clf
    set(gcf, 'Position', [0, 0, 1000, 600])

    for it = 1:algrithmNum
        yy(it, :) = mean(all_curves((it - 1) * Flod + 1:(it - 1) * Flod + Flod, :));
    end

    xx = [1:NumofRecord] * (MaxFEs / NumofRecord);

    for it = 1:algrithmNum
        semilogy(xx, yy(it, :), [linestyles{it} Markers{it}], 'LineWidth', 1.5, 'Color', MarkerEdgeColors(it, :));
        hold on;
    end

    hold off;
    set(gcf, 'color', 'white')
    set(gca, 'YScale', 'log', 'YLimMode', 'auto')
    title(Function_name)
    xlabel('FEs');
    ylabel('Best Score');
    algorithmName1 = strrep(algorithmName, '_', '\_');
    legend(algorithmName1);
    %     legend(algorithmName,'Best');
    legend('boxoff')

    a = findobj(gcf); % get the handles associated with the current figure
    allaxes = findall(a, 'Type', 'axes');
    % alllines=findall(a,'Type','line');
    alltext = findall(a, 'Type', 'text');
    set(allaxes, 'FontName', 'Times', 'LineWidth', 1, 'FontSize', 14, 'FontWeight', 'bold');
    % set(alllines,'Linewidth',1);
    set(alltext, 'FontName', 'Times', 'FontSize', 14, 'FontWeight', 'bold')
    %
    set(gcf, 'PaperUnits', 'inches');
    % set(gcf, 'PaperUnits', 'centimeters');
    krare = 3.5;
    % x_width=krare*1.618 ;
    x_width = krare * 5/3;
    %  x_width=3*1;
    y_width = krare * 4/3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    % set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 4])
    set(gca, ...
        'Box', 'on', ...
        'TickDir', 'in', ...
        'TickLength', [.02 .02], ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'YGrid', 'on', ...
        'XGrid', 'off', ...
        'XColor', [.3 .3 .3], ...
        'YColor', [.3 .3 .3], ...
        'LineWidth', 1);
    axis tight
    grid on
    box on
    saveas(gcf, [filename, '-', Function_name, '-carve'], 'fig')
    filename1 = [filename, '-', Function_name, '-carve'];
    print(filename1, '-dtiff', '-r300'); %<-Save as PNG with 300 DPI
end

%     Orderhao(xlsfilename);
%     pValueToExcelKJ1(xlsfilename,Flod);
%     FridTest3( xlsfilename ,Flod)
%     FridTest4hao( xlsfilename ,Flod)
% Order(xlsfilename);
% pValueToExcelKJ(xlsfilename,Flod);
% FridTest3( xlsfilename ,Flod)
% FridTest4( xlsfilename ,Flod)
Orderhao(xlsfilename);
pValueToExcelhao(xlsfilename, Flod);
FridTest3(xlsfilename, Flod)
FridTest4(xlsfilename, Flod)

toc
close all
