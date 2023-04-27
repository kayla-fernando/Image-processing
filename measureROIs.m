%% Measure ROIs %%

% A semi-automated Z-stack analysis pipeline to:
%   1) Save Z-stacks as individual .tifs (for each channel)
%   2) Select and save ROIs (to be applied to all channels that will be analyzed) 
%   3) Measure and save results of user-selected parameters for each ROI (for each Z-plane for each channel) 
%   4) Calculate the corrected total cell fluorescence (CTCF) for each channel and perform ratiometric analysis

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

% Written by Kayla Fernando (kayla.fernando@duke.edu) (4/26/23)

%% Set paths and use MIJ to upload image to ImageJ

close all; clear all; clc

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
    MIJ.run('Bio-Formats Importer') 
    MIJ.run('Split Channels')

% Make stacks folder to store the Z-stacks for each channel
mkdir([analysisPath mouse '\' mouse '_' num2str(sample) '\stacks']);
disp('Save channels in image analysis folder')
    MIJ.run('Save')

%% Manually select ROIs at each Z-plane that will be applied to all channels

% Make rois folder to store the ROI coordinates for each Z-plane 
mkdir([analysisPath mouse '\' mouse '_' num2str(sample) '\rois']);
newFolder = ([analysisPath mouse '\' mouse '_' num2str(sample) '\rois']);
cd(newFolder)
disp('Select and save ROIs at each Z-plane')
    MIJ.run('ROI Manager...') % check Show All

% Select ROIs at this Z-plane, then save the RoiSet to the rois folder
% Clear values and move to the next Z-plane

%% Measure values of all selected ROIs at each Z-plane for the current channel

% Analyze the current active channel in ImageJ and set folder naming conventions
channel = 'GFP';

% Separate current channel into individual .tifs
newFolder = ([analysisPath mouse '\' mouse '_' num2str(sample) '\stacks']);
cd(newFolder)
MIJ.run('Stack to Images')

% Save each image
for k = 1:NumberOfZPoints
    MIJ.run('Save')
    MIJ.run('Close')
end 

% Make results folder for current channel
mkdir([analysisPath mouse '\' mouse '_' num2str(sample) '\' channel 'results']);
newFolder = ([analysisPath mouse '\' mouse '_' num2str(sample) '\' channel 'results']);
cd(newFolder)
for k = 1:NumberOfZPoints
    disp(['Open Z-plane number ' num2str(k)])
        MIJ.run('Open...') 
    disp(['Open RoiSet' num2str(k)])
        MIJ.run('Open...') 
    disp(['Save results of Z-plane number ' num2str(k) ...
        ', then clear Results, clear ROI Manager, and close current Z-plane'])
        MIJ.run('measureROIs') % custom macro
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
