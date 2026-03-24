%===================================================================================
% MATLAB code for multi-level image thresholding segmentation using 2DNLMeKGSA.
% Author: Himanshu Mittal (himanshu.mittal224@gmail.com), 
%           Mukesh Saraswat (saraswatmukesh@gmail.com)
% Modified this file for the non-commercial purpose only.
% (https://people.eecs.berkeley.edu/~yang/software/lossy_segmentation/)
%
% Developed in MATLAB R2015a
%
% Reference: "An optimum multi-level image thresholding segmentation using
%            non-local means 2D histogram and exponential Kbest gravitational 
%            search algorithm." Engineering Applications of Artificial 
%            Intelligence, Volume 71, Pages 226-235, Elsevier, 2018. 
%            https://doi.org/10.1016/j.engappai.2018.03.001
%
% File purpose: Measuring the BDE, VoI, PRE, and GCE values.
%===================================================================================
 
function seg_measure1(sheetname,xlsxname,xlsimagename,pathName,image_files,measure1_name,ImageSize)

testImageList.name=[];
for i=1:size(image_files,2)
    image_name=image_files(i).name;
    
    flag=[];
    sign=0;
    while ~strcmp(flag,'.') 
        flag=image_name(end-sign); 
        sign=sign+1;
    end
    
    
    testImageList(i).name=[image_name(1,1:end-sign),'.mat'];
end


benchPath = 'benchset/';
strc=strcat(pathName,'/');
testPath = strc;
cropBenchImage = true;

% testImageList = dir(testPath);
if isempty(testImageList)
    error('Cannot find the image directories.');
end

% resizeScale = 320;
width=ImageSize.width;
height=ImageSize.height;

averageBoundaryError = 0;
averageRI = 0;
averageVOI = 0;
averageGCE = 0;

imageCount = 0;
for imageIndex=1:length(testImageList)
    if testImageList(imageIndex).name(1)=='.'
        continue;
    end
    
    testFilename = [testPath testImageList(imageIndex).name];
    if ~strcmp(lower(testFilename(end-3:end)),'.mat')
        % Not a valid data file
        continue;
    end
    
    benchFilename = [benchPath testImageList(imageIndex).name];
    if isempty(dir(benchFilename))
        % The bench data for this image is not available
        warning(['Cannot find the bench file for ' testImageList(imageIndex).name]);
    end
    
    disp(['Processing ' testImageList(imageIndex).name]);
    superpixelLabels = [];
    load(benchFilename);
    load(testFilename);
    imageCount = imageCount + 1;
    
    % Comparison script
    totalBoundaryError = 0;
    sumRI = 0;
    sumVOI = 0;
    sumGCE = 0;

    [imageX, imageY] = size(sampleLabels);
    [benchX, benchY] = size(imageLabelCell{1});
    
    for benchIndex=1:length(imageLabelCell)
        benchLabels = imageLabelCell{benchIndex};
        
        % Rescale benchimage
%         if benchX>benchY
%             benchY=benchY*resizeScale/benchX;
%             benchX=resizeScale;
%         else
%             benchX=benchX*resizeScale/benchY;
%             benchY=resizeScale;
%         end
        if benchX>=benchY
%                         resizedImageY = round(imageY/imageX*imageSize);
            benchLabels = imresize(benchLabels,[width,height],'nearest');
        else
%                         resizedImageX = round(imageX/imageY*imageSize);
            benchLabels = imresize(benchLabels,[height,width],'nearest');

        end
        
        
        
        
%         benchLabels=imresize(benchLabels, [benchX, benchY],'nearest');
        
        cropBoundarySize = (size(benchLabels,1)-size(sampleLabels,1))/2;
            
        % update the four error measures:        
        totalBoundaryError = totalBoundaryError + compare_image_boundary_error(benchLabels, sampleLabels);        
        
        [curRI,curGCE,curVOI] = compare_segmentations(sampleLabels,benchLabels);       
        sumRI = sumRI + curRI;
        sumVOI = sumVOI + curVOI;
        sumGCE = sumGCE + curGCE;        
    end
    BE(imageCount)=totalBoundaryError;
    RI(imageCount)=sumRI;
    VOI(imageCount)=sumVOI;
    GCE(imageCount)=sumGCE;
    % update the averages... note that sumRI / length(imageLabelCell) is
    % equivalent to the PRI.
    averageBoundaryError = averageBoundaryError + totalBoundaryError / length(imageLabelCell);
    averageRI = averageRI + sumRI / length(imageLabelCell);
    averageVOI = averageVOI + sumVOI / length(imageLabelCell);
    averageGCE = averageGCE + sumGCE / length(imageLabelCell);

%     disp(['Current err:  Boundary  RI  VOI  GCE:']);
%     disp([num2str(averageBoundaryError/imageCount) '  ' num2str(averageRI/imageCount) ...
%          '  ' num2str(averageVOI/imageCount) '  ' num2str(averageGCE/imageCount)]);
end
for i=1:size(BE,2)
    BDE1{i}=BE(i);
    PRI1{i}=RI(i);
    VOI1{i}=VOI(i);
    GCE1{i}=GCE(i);
end

A=[];
for l=1:size(measure1_name,2)
    A=[A;eval([measure1_name{l},'1'])];   
end

A=A';
% xlsA={'BDE','PRI','VOI','GCE'};
xlsA=measure1_name;
xlsA=[xlsA;A];
% st1=strcat(sheetname_eval);
st1=strcat(sheetname,'_eval');
xlswrite(xlsxname,xlsimagename,st1,'A2');
xlswrite(xlsxname,xlsA,st1,'B1');
% xlswrite(xlsxname,xlsA,st1);

AvgBDE1 = averageBoundaryError / imageCount;
AvgPRI1 = averageRI / imageCount;
AvgVOI1 = averageVOI / imageCount;
AvgGCE1 = averageGCE / imageCount;


B=[];
xlsB={};
for l=1:size(measure1_name,2)
    B{1,l}=eval(['Avg',measure1_name{l},'1']);
    xlsB{1,l}=['Avg',measure1_name{l}];
end



% B={AvgBDE1,AvgPRI1,AvgVOI1,AvgGCE1};

% xlsB={'AvgBDE','AvgPRI','AvgVOI','AvgGCE'};
xlsB=[xlsB;B];
st2=strcat(sheetname,'_eval_avg');
xlswrite(xlsxname,xlsB,st2);
end