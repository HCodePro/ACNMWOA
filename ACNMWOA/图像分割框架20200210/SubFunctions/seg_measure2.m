function seg_measure2(sheetname,xlsxname,xlsimagename,pathName,image_files,measure1_name,measure2_name)

directoryFiles.name=[];
for i=1:size(image_files,2)
    image_name=image_files(i).name;
    
    flag=[];
    sign=0;
    while ~strcmp(flag,'.') 
        flag=image_name(end-sign); 
        sign=sign+1;
    end
    
    
    directoryFiles(i).name=[image_name(1,1:end-sign),'.mat'];
end



filePath = 'BW_dataset/';
% directoryFiles = dir(filePath);
% strc=strcat(methodname,'_result-set');
% pathName = strc;
counter=1;
for fileIndex=1:length(directoryFiles)
    currentFilename = directoryFiles(fileIndex).name;
    disp(currentFilename);
    
    if (length(currentFilename)>4) && strcmp(currentFilename(end-3:end),'.mat')
        % It is an image file        
        imageName = currentFilename(1:end-4);
        dataFilename = [filePath imageName '.mat'];
        dataFilename1 = [pathName '/' imageName '.mat'];
        load(dataFilename);
        ref=refImage;
        
        
        
        
        
        load(dataFilename1);
        A=uint8(superpixelLabels);
        
        %% Calculate the SSIM
        [ssimval, ssimmap] = ssim(A,ref);
        SSIM(counter)=ssimval;

        %% Calculate the FSIM
        [fsim, fsimc] = FeatureSIM(ref, A);
        FSIM(counter)=fsim;

        %Mean Square Error 
        MSE(counter) = MeanSquareError(ref, A);

        %Root Mean Square Error 
        RMSE(counter) = sqrt(MeanSquareError(ref, A));

        %Peak Signal to Noise Ratio 
        PSNR(counter) = PeakSignaltoNoiseRatio(ref, A);

        %Normalized Cross-Correlation 
        NK(counter) = NormalizedCrossCorrelation(ref, A);

        %Average Difference 
        AD(counter) = AverageDifference(ref, A);

        %Structural Content 
        SC(counter) = StructuralContent(ref, A);

        %Maximum Difference 
        MD(counter) = MaximumDifference(ref, A);

        %Normalized Absolute Error
        NAE(counter) = NormalizedAbsoluteError(ref, A);

        counter=counter+1;
        
    end
end

for i=1:size(SSIM,2)    
    SSIM1{i}=SSIM(i);
    FSIM1{i}=FSIM(i);
    RMSE1{i}=RMSE(i);
    MSE1{i}=MSE(i);
    PSNR1{i}=PSNR(i);
    NCC1{i}=NK(i);
    AD1{i}=AD(i);
    SC1{i}=SC(i);
    MD1{i}=MD(i);
    NAE1{i}=NAE(i);
end


A=[];
for l=1:size(measure2_name,2)
    A=[A;eval([measure2_name{l},'1'])];   
end
A=A';


% A=[SSIM1;FSIM1;RMSE1;PSNR1;NCC1;AD1;MD1;NAE1]';
% xlsA={'SSIM','FSIM','RMSE','PSNR','NCC','AD','MD','NAE'};

xlsA=measure2_name;
xlsA=[xlsA;A];
% st1=strcat(sheetname,'_II');
st1=strcat(sheetname,'_eval');

xlswrite(xlsxname,xlsimagename,st1,'A2');
% xlswrite(xlsxname,xlsimagename,st1,'A2');

loaction='A'+length(measure1_name)+1;
loaction=char(loaction);


xlswrite(xlsxname,xlsA,st1,[loaction,'1']);

% xlswrite(xlsxname,xlsA,st1);



% A=[SSIM;FSIM;MSE;RMSE;PSNR;NK;AD;SC;MD;NAE];
% st1=strcat(sheetname,'_II');
% xlswrite(xlsxname,A',st1);

        %% Average
        AvgSSIM=mean(SSIM)        ;
        AvgFSIM=mean(FSIM);
        AvgMSE=mean(MSE);
        AvgRMSE=mean(RMSE);
        AvgPSNR=mean(PSNR);
        AvgNCC=mean(NK);
        AvgAD=mean(AD);
        AvgSC=mean(SC);
        AvgMD=mean(MD);
        AvgNAE=mean(NAE);
    
% 'SSIM','FSIM','RMSE','MSE','PSNR','NCC','AD','SC','MD','NAE'        
        
B=[];
xlsB={};
for l=1:size(measure2_name,2)
    B{1,l}=eval(['Avg',measure2_name{l}]);
    xlsB{1,l}=['Avg',measure2_name{l}];
end
        
        
        
% B={AvgSSIM,AvgFSIM,AvgRMSE,AvgPSNR,AvgNCC,AvgAD,AvgMD,AvgNAE};
% xlsB={'AvgSSIM','AvgFSIM','AvgRMSE','AvgPSNR','AvgNCC','AvgAD','AvgMD','AvgNAE'};
xlsB=[xlsB;B];
% st2=strcat(sheetname,'_avg_II');
st2=strcat(sheetname,'_eval_avg');

loaction='A'+length(measure1_name);
loaction=char(loaction);

xlswrite(xlsxname,xlsB,st2,[loaction,'1']);

% st2=strcat(sheetname,'_avg_II');
% xlswrite(xlsxname,B,st2);
end