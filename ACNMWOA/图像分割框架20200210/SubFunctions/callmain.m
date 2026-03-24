function [methodata, xlsimagename, cg_curves, FEcount, image_info] = callmain(methodname, pathName, partitionLevels, MaxFEs, Population_size, image_files, NumofRecord, objectiveFunction, ImageSize, Iteration)
    NORMALIZE_IMAGE = true;
    %     NumofRecord=40;
    method = methodname;
    methodata = [];
    filePath = 'Berkeley-dataset/';
    %     directoryFiles = dir(filePath);
    directoryFiles = image_files;
    it = 1;

    for i = 1:length(directoryFiles)
        xlscurrentFilename = directoryFiles(i).name;
        %         if (length(xlscurrentFilename)>4) && strcmp(xlscurrentFilename(end-3:end),'.jpg')

        if (length(xlscurrentFilename) > 4)
            xlsimagename{it} = xlscurrentFilename;
            it = it + 1;
        end

    end

    im_sign = 1;

    for fileIndex = 1:length(directoryFiles)
        currentFilename = directoryFiles(fileIndex).name;
        %         disp(currentFilename);
        %         if (length(currentFilename)>4) && strcmp(currentFilename(end-3:end),'.jpg')
        if (length(currentFilename) > 4)
            % It is an image file

            %去除文件名后缀
            flag = [];
            sign = 0;

            while ~strcmp(flag, '.')
                flag = currentFilename(end - sign);
                sign = sign + 1;
            end

            imageName = currentFilename(1:end - sign);
            dataFilename1 = [pathName '\' imageName '.mat'];

            dataFilename2_path = [pwd '\BW_dataset'];
            %             exist dataFilename2_path dir
            if exist(dataFilename2_path, 'dir') == 0
                mkdir(dataFilename2_path);
            end

            dataFilename2 = [dataFilename2_path '\' imageName '.mat'];
            imageFullPath = [filePath currentFilename];
            image = imread(imageFullPath);

            if ~isempty(whos('image'))
                [imageX, imageY, imageDepth] = size(image);

                if NORMALIZE_IMAGE
                    % Resize the image to 320*240
                    %                     imageSize = 320;
                    %                     if imageX>=imageY
                    %                         resizedImageY = round(imageY/imageX*imageSize);
                    %                         image = imresize(image,[imageSize,resizedImageY]);
                    %                     else
                    %                         resizedImageX = round(imageX/imageY*imageSize);
                    %                         image = imresize(image,[resizedImageX,imageSize]);
                    %                     end

                    %                     imageSize = 320;
                    width = ImageSize.width;
                    height = ImageSize.height;

                    if imageX >= imageY
                        %                         resizedImageY = round(imageY/imageX*imageSize);
                        image = imresize(image, [width, height]);
                    else
                        %                         resizedImageX = round(imageX/imageY*imageSize);
                        image = imresize(image, [height, width]);
                    end

                end

                %% Convert to GreyScale
                if (imageDepth == 3)
                    image = rgb2gray(image);
                end

                % 保存灰度图

                %                 figure;
                %                 imshow(image);
                filename = [[pathName '\'], [currentFilename(1:end - sign) '_gray']];
                %                 save_pic(filename);

                image_info(im_sign).gray = image;
                image_info(im_sign).gray_filename = filename;

                % 保存一维直方图
                %                 figure;
                %                 imhist(image);
                filename = [[pathName '\'], [currentFilename(1:end - sign) '_hist']];
                %                 save_pic(filename);

                image_info(im_sign).hist = image;
                image_info(im_sign).hist_filename = filename;

                %% Compute the non-local means (NLM) of an image
                refImage = image;
                I = image;
                number_of_levels = partitionLevels + 1;
                level = partitionLevels;
                [m, n] = size(I);
                a = I;
                a0 = im2double(a);
                t = 7;
                f = 2;
                a3 = nlmeans(a0, t, f);
                % 保存非均值滤波图像
                %                 figure;
                %                 imshow(a3);
                filename = [[pathName '\'], [currentFilename(1:end - sign) '_NLM']];
                %                 save_pic(filename);
                image_info(im_sign).NLM = a3;
                image_info(im_sign).NLM_filename = filename;

                for i = 1:m

                    for j = 1:n
                        a4(i, j) = a3(i, j);
                    end

                end

                a4 = (a4 - min(a4(:))) / (max(a4(:)) - min(a4(:)));
                a4 = im2uint8(a4);
                [m, n] = size(a);
                a0 = (a);
                fxy = zeros(256, 256);

                for i = 1:m

                    for j = 1:n
                        c = a0(i, j);
                        d = (a4(i, j));
                        fxy(c + 1, d + 1) = fxy(c + 1, d + 1) + 1;
                    end

                end

                Pxy = fxy / m / n;
                %% Display the 2D-histogram
                % 保存二维直方图
                %                 figure;
                %                 mesh(Pxy);
                filename = [[pathName '\'], [currentFilename(1:end - sign) '_2D']];
                %                 save_pic(filename);
                image_info(im_sign).TwoD = Pxy;
                image_info(im_sign).TwoD_filename = filename;

                Lmax1 = 254;

                if size(I, 3) == 1 %grayscale image
                    mainename = ['maine' method];
                    mainename = str2func(mainename(1, :));
                    [gBest, gbestvalue, FEcount, etime, MaxFEs, Convergence_curve, iter] = mainename(Lmax1, level, Pxy, MaxFEs, Population_size, objectiveFunction, Iteration);
                    %% return optimal intensity
                    intensity = round(gBest);
                    %% return fitness value
                    fitness = gbestvalue;
                end

                %% Extracting the threshold values
                Thresholds = intensity((number_of_levels):end);
                %% Sorting the thresholds
                srtThr = sort(Thresholds);
                %% Checking the case of similar threshold values
                v = find(diff(srtThr) == 0);

                while (size(v, 2) > 0)

                    if size(I, 3) == 1 %grayscale image
                        [gBest, gbestvalue, FEcount, etime, MaxFEs, Convergence_curve, iter] = mainename(Lmax1, level, Pxy, MaxFEs, Population_size, objectiveFunction, Iteration);
                        %% return optimal intensity
                        intensity = round(gBest);
                        %% return fitness value
                        fitness = gbestvalue;
                    end

                    %% Extracting the threshold values
                    Thresholds = intensity((number_of_levels):end);
                    %% Sorting the thresholds
                    srtThr = sort(Thresholds);
                    %% Checking the case of similar threshold values
                    v = find(diff(srtThr) == 0);
                end

                %% 原程序缺少此文件
                X = imageGRAY(I, Thresholds);
                %% For displaying the segmented image
                %                  figure;
                %                  imshow(X);
                filename = [[pathName '\'], [currentFilename(1:end - sign) '_seg_gray']];
                %                  save_pic(filename);
                image_info(im_sign).seg_gray = X;
                image_info(im_sign).seg_gray_filename = filename;

                X1 = X;
                X1 = double(X1);
                Thresholds = sort(Thresholds);
                X = imquantize(I, Thresholds);
                Y = label2rgb(X);
                %% For displaying the segmented image
                %                  figure;
                %                  imshow(Y);
                % 保存图片
                filename = [[pathName '\'], [currentFilename(1:end - sign) '_seg_color']];
                %                 save_pic(filename);

                image_info(im_sign).seg_color = Y;
                image_info(im_sign).seg_color_filename = filename;

                filename = [[pathName '\'], [currentFilename(1:end - sign) '_seg']];
                image_info(im_sign).seg = image;
                image_info(im_sign).seg_filename = filename;
                image_info(im_sign).seg_level = level;
                image_info(im_sign).seg_gBest = gBest;

                %                 figure;

                %                 set (gcf,'position',[500,500,800,200] );
                %                 plot(imhist(image),'LineWidth', 1);
                %                 set(gca,'XTick',0:100:300)
                %                 set(gca,'YTick',0:round(max(imhist(image))/3):max(imhist(image)))
                %                 axis([0,300,0,max(imhist(image))]);
                %                 hold on;
                %                 for c=1:level-1
                %                     plot([gBest(level+c-1) gBest(level+c-1)], get(gca, 'YLim'), '-r', 'LineWidth', 1);
                %                 end
                %                 hold off;
                %                 filename=[[pathName '\'],[currentFilename(1:end-sign) '_seg']];
                %                 a=findobj(gcf); % get the handles associated with the current figure
                %                 allaxes=findall(a,'Type','axes');
                %                 % alllines=findall(a,'Type','line');
                %                 alltext=findall(a,'Type','text');
                %                 set(allaxes,'FontName','Times','LineWidth',1,'FontSize',14,'FontWeight','bold');
                %                 % set(alllines,'Linewidth',1);
                %                 set(alltext,'FontName','Times','FontSize',14,'FontWeight','bold')
                %                 %
                %                 krare=3.5;
                %                 set(gca, ...
                %                     'Box'         , 'on'     , ...
                %                     'TickDir'     , 'in'     , ...
                %                     'TickLength'  , [.02 .02] , ...
                %                     'XMinorTick'  , 'on'      , ...
                %                     'YMinorTick'  , 'on'      , ...
                %                     'YGrid'       , 'off'      , ...
                %                     'XGrid'       , 'off'      , ...
                %                     'XColor'      , [.3 .3 .3], ...
                %                     'YColor'      , [.3 .3 .3], ...
                %                     'LineWidth'   , 1         );
                %                 axis tight
                %             %     grid on
                %                 box on
                %                 saveas(gcf, filename,'fig')
                %                 print(filename,'-dtiff', '-r300'); %<-Save as PNG with 300 DPI

                %% Post-Processing steps
                for l = 1:2 * level
                    outdata{1, l} = gBest(l);
                end

                %                 outdata={gBest(1),gBest(2),gBest(3),gBest(4),gBest(5),gBest(6),gBest(7),gBest(8), gbestvalue, FEcount,etime, MaxFEs};
                outdata{1, 2 * level + 1} = gbestvalue;
                outdata{1, 2 * level + 2} = FEcount;
                outdata{1, 2 * level + 3} = etime;
                outdata{1, 2 * level + 4} = MaxFEs;
                outdata{1, 2 * level + 5} = iter;
                methodata = [methodata; outdata];
                %                 xlsimagename=xlsimagename1;
                cg_curves(fileIndex, :) = MySampling(Convergence_curve, NumofRecord);
                sampleLabels = X;
                superpixelLabels = X1;

                image_info(im_sign).mat_sampleLabels = sampleLabels;
                image_info(im_sign).mat_superpixelLabels = superpixelLabels;
                image_info(im_sign).mat_refImage = refImage;
                image_info(im_sign).mat_dataFilename1 = dataFilename1;
                image_info(im_sign).mat_dataFilename2 = dataFilename2;
                %                 save (dataFilename1, 'sampleLabels', 'superpixelLabels');
                %                 save (dataFilename2, 'refImage');
                im_sign = im_sign + 1;
                clearvars -except im_sign image_info FEcount Iteration ImageSize objectiveFunction cg_curves NumofRecord refImage dataFilename2 Population_size MaxFEs partitionLevels xlsimagename method X1 X J dataFilename1 filePath directoryFiles pathName fileIndex NORMALIZE_IMAGE methodata;
            end

        end

    end

end
