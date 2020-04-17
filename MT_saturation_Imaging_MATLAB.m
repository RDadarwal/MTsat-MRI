%---------------------------------------------------------------------------------------------------------------------------------
%%                    Magnetization Transfer saturation Imaging
%                            Author: Rakshit Dadarwal; update: 17/04/2020

% Reference: MT saturation (MT sat) parameter estimation and model fitting is based on the method
%described by Helms et al. 2008 (https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.21732).
%---------------------------------------------------------------------------------------------------------------------------------
clc;clear;

% subject id
sub=297;
% MRI scanner dicom storage directory
dcm_dir='/mnt/scanner/siemens/human/sub_0';
sub_dcm=strcat(dcm_dir,num2str(sub),'/dicom');
% measurement numbers
pd_meas=34; t1_meas=28; mt_meas=31;
% save directory
save_dir=strcat('/home/user/sub_0',num2str(sub),'/');

if exist(sub_dcm) == 7
    sprintf('start magnetization transfer saturation model fitting...')
    cd(sub_dcm)
    % select dicom files
    [filenamePD, pathnamePD] = uigetfile('multiselect', 'on', strcat('*-00',num2str(pd_meas),'-*.*'), 'select PD-weighted images');
    [filenameMT, pathnameMT] = uigetfile('multiselect', 'on', strcat('*-00',num2str(mt_meas),'-*.*'),'select MT-weighted images');
    [filenameT1, pathnameT1] = uigetfile('multiselect', 'on',strcat('*-00',num2str(t1_meas),'-*.*'), 'select T1-weighted images');
    
    % read dicom header information
    InfoFirstPD = dicominfo(fullfile(pathnamePD,filenamePD{1}));
    W = InfoFirstPD.Height; H = InfoFirstPD.Width;
    TRpd = InfoFirstPD.RepetitionTime; TEpd = InfoFirstPD.EchoTime;
    FApd = InfoFirstPD.FlipAngle* pi / 180; Navpd=InfoFirstPD.NumberOfAverages;
    
    InfoFirstMT = dicominfo(fullfile(pathnameMT,filenameMT{1}));
    TRmt = InfoFirstMT.RepetitionTime; TEmt = InfoFirstMT.EchoTime;
    FAmt = InfoFirstMT.FlipAngle* pi / 180; Navmt=InfoFirstMT.NumberOfAverages;
    
    InfoFirstT1 = dicominfo(fullfile(pathnameT1,filenameT1{1}));
    TRt1 = InfoFirstT1.RepetitionTime; TEt1 = InfoFirstT1.EchoTime;
    FAt1 = InfoFirstT1.FlipAngle* pi / 180; Navt1=InfoFirstT1.NumberOfAverages;
    
    % change directory
    cd(save_dir)
    
    % read dicom data
    for n=1:length(filenamePD)
        PD(:,:,n) = double(dicomread(fullfile(pathnamePD,filenamePD{n})));
    end
    
    for n=1:length(filenameMT)
        MT(:,:,n) = double(dicomread(fullfile(pathnameMT,filenameMT{n})));
    end
    
    for n=1:length(filenameT1)
        T1(:,:,n) = double(dicomread(fullfile(pathnameT1,filenameT1{n})));
    end
    
    % create brain mask
    mask = PD / max(PD(:));
    thr = graythresh(mask);
    mask(mask<thr/2) = 0;
    mask(mask>0) = 1;
    
    % apply brain mask on T1-w, MT-w and PD-w images
    T1 = T1.*mask; PD = PD.*mask; MT = MT.*mask;
    
    %% Calculate MTR, apparent T1 and, MT saturation maps
    
    % calcualte apparent T1
    T1_app = 2 .* ((PD./FApd - T1./FAt1)./(T1.*FAt1./TRt1 - PD.* FApd./TRpd));
    
    % calcualte apparent A
    A_app  = T1 .* PD ./ (PD * FApd - T1 * FAt1) * (FApd/FAt1 - (FAt1/FApd));
    
    % calcualte apparent Magnetization Transfer saturation
    MT_sat = (A_app * FAmt ./ MT -1) ./T1_app * TRmt - FAmt^2/2;
    
    % calculate Magnetization Transfer Ratio (MTR)
    MTR = (PD - MT)./PD * 100;
    
    %% save results
    niftiwrite(abs(MT_sat), 'MT_sat','Compressed',true); niftiwrite(abs(MTR), 'MTR','Compressed',true);
    niftiwrite(abs(T1_app), 'T1_app','Compressed',true); niftiwrite(abs(A_app), 'A_app','Compressed',true);
end
