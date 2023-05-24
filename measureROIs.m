%% Measure ROIs %%

% A semi-automated Z-stack analysis pipeline to:
%   1) Save Z-stacks as individual .tifs (for each channel)
%   2) Manually set the threshold and generate binary mask for each Z-plane for each channel
%   3) Select and save ROIs from masks (to be applied to all channels that will be analyzed) 
%   4) Measure and save results of user-selected parameters for each ROI (for each Z-plane for each channel) 
%   5) Calculate the corrected total cell fluorescence (CTCF) for each channel and perform ratiometric analysis

% NOTES:
%   - Requires installation of MIJ and ImageJ from MATLAB File Exchange: 
%       https://www.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlabImageJ
%   - Requires download of Bio-Formats Package to ImageJ .jars folder:
%       https://www.openmicroscopy.org/bio-formats/downloads/
%   - Requires download of the most updated NetCDF java source files (.zip) to ImageJ .jars folder:
%       https://downloads.unidata.ucar.edu/netcdf/
%   - ImageJ setup prior to running code: 
%       - spatially calibrate images for correct "Area" measurement for CTCF
%       - calibrate gray values of images for correct "Mean" measurement for CTCF
%       - set results in ROI Manager to include "Integrated Density" measurement for CTCF
%       - directory to custom macro "measureROIs.ijm" to do the measurements using ROI Manager
%   - Run each section separately
%   - Wand tool for cell ROI selection, Oval tool for background ROI selection

% Written by Kayla Fernando (kayla.fernando@duke.edu) (5/17/23)

%% Set paths and use MIJ to upload image to ImageJ

close all 
clear all 
clc

mouse = 'mouse';
sample = 1; % 1 = EBC side; 2 = non-EBC side
NumberOfZPoints = 21; % from image metadata 
analysisPath = 'Y:\\'; % directory to image analysis folder
applicationPath = 'Y:\\'; % directory to image processing applications
filePath = [analysisPath mouse '\' mouse '_' num2str(sample)]; % used later for CTCF calculation

vers = ['R' version('-release')];
javaaddpath(['C:\Program Files\MATLAB\' vers '\java\mij.jar']);
javaaddpath(['C:\Program Files\MATLAB\' vers '\java\ij.jar']);
MIJ.start(java.lang.String(applicationPath)) % loads macros
disp('Open an .ims file')
    MIJ.run('Bio-Formats Importer') % Data Browser; Group files with similar names; Dimensions; Series 3

% Make stacks folder to store the Z-stacks for each channel
mkdir([analysisPath mouse '\' mouse '_' num2str(sample) '\stacks']);
newFolder = ([analysisPath mouse '\' mouse '_' num2str(sample) '\stacks']);
cd(newFolder)
disp('Save as TIFF in image analysis folder and ignore the big scary exception that will pop up')
    MIJ.run('Bio-Formats Exporter') % Name appropriately; Tagged Image File Format; Write each channel to a separate file 

%% Each channel will be saved as a grayscale image stack, pseudo-color each channel accordingly

MIJ.run("RGB Color");
MIJ.run("Make Composite");
MIJ.run("Split Channels");

%% Separate each channel into individual .tifs and make a copy that will be used for making masks

newFolder = ([analysisPath mouse '\' mouse '_' num2str(sample) '\stacks']);
cd(newFolder)
MIJ.run('Stack to Images')

for k = 1:NumberOfZPoints
    MIJ.run('Save')
    MIJ.run('Duplicate...')
    MIJ.run('Save')
    MIJ.run('Close')
    MIJ.run('Close')
end 
              
%% Create binary mask and manually select ROIs at each Z-plane that will be applied to all channels

% Make rois folder to store the ROI coordinates for each Z-plane 
mkdir([analysisPath mouse '\' mouse '_' num2str(sample) '\rois']);
newFolder = ([analysisPath mouse '\' mouse '_' num2str(sample) '\rois']);
cd(newFolder)
        
for k = 1:NumberOfZPoints
    disp(['Open *copy* of GFP Z-plane number ' num2str(k)])
        MIJ.run('Open...')
    disp('Manually adjust brightness and contrast and Apply')   
        MIJ.run('Brightness/Contrast...') % okay to update LUT
    pause
    disp(['Open *copy* of RFP Z-plane number ' num2str(k)])
        MIJ.run('Open...')
    disp('Manually adjust brightness and contrast and Apply')   
        MIJ.run('Brightness/Contrast...') % okay to update LUT
    pause
    MIJ.run('Merge Channels...')
    MIJ.run('RGB Color')
    MIJ.run('8-bit')
    MIJ.run('Gaussian Blur...') % sigma = 1.25
    disp('Manually adjust threshold and Apply')
        MIJ.run('Threshold...')
    pause
    disp(['Watershed, then save as mask' num2str(k)])
        MIJ.run('Convert to Mask')
        MIJ.run('Watershed')
        MIJ.run('Tiff...')
    pause
    disp('Select and save cell and background ROIs from this mask')
    disp('Clear values before continuing')
        MIJ.run('ROI Manager...')
    pause
end

%% Measure values of all selected ROIs at each Z-plane for the current channel

% Analyze the current active channel in ImageJ and set folder naming conventions
channel = 'GFP';

% Make results folder for current channel
mkdir([analysisPath mouse '\' mouse '_' num2str(sample) '\' channel 'results']);
newFolder = ([analysisPath mouse '\' mouse '_' num2str(sample) '\' channel 'results']);
cd(newFolder)
              
for k = 1:NumberOfZPoints
    disp(['Open *original* Z-plane number ' num2str(k)])
        MIJ.run('Open...') 
    disp(['Open RoiSet' num2str(k)])
        MIJ.run('Open...') 
    disp(['Save results of Z-plane number ' num2str(k) ...
        ', then clear Results, clear ROI Manager, and close current Z-plane']) % Select last ROI; Measure; Save As
    pause
end

% Repeat for all desired channels

%% Calculate CTCF and do ratiometric analysis

[ratio,gfpCTCF,rfpCTCF] = CTCF(filePath,NumberOfZPoints);

% Make histograms of RFP/GFP ratio for each Z-plane
for k = 1:length(ratio)
    figure; 
    histogram(ratio{k},'BinWidth',0.5); 
    hold on; 
    title(['Z-plane ' num2str(k)]);
    xlabel('CTCF RFP/CTCF GFP ratio');
    ylabel('# of ROIs');
    xline(1); 
    hold off;
end

% Close ImageJ
MIJ.closeAllWindows
MIJ.exit
