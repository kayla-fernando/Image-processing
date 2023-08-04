close all; clear all; clc

mouse = 'KF124';
slide = 1; % slide number
slice = 6; % slice number
sample = 1; % 1 = EBC side; 2 = non-EBC side
condition = 'stained'; % stained or unstained 
analysisPath = 'Y:\\home\kayla\Histology analysis\'; % directory to image analysis folder
filepath = [analysisPath mouse ' ' condition '\' mouse '_' num2str(slide) '_slice' num2str(slice) '_' num2str(sample)]; % used later for CTCF calculation

ZPlaneToUse = 5;

gfpCTCF = cell(1,numel(ZPlaneToUse));
rfpCTCF = cell(1,numel(ZPlaneToUse));
ratio = cell(1,numel(ZPlaneToUse));

% For the selected Z-plane
for k = ZPlaneToUse
    % Get the GFP results
    gfpResults = readmatrix([filepath '\GFP results non PCs.csv']); 
    gfpBackground = readmatrix([filepath '\GFP results background.csv']);
        gfpNames = detectImportOptions([filepath '\GFP results non PCs.csv']).VariableNames;
    % Get the RFP results
    rfpResults = readmatrix([filepath '\RFP results non PCs.csv']); % get the RFP results
    rfpBackground = readmatrix([filepath '\RFP results background.csv']);
        rfpNames = detectImportOptions([filepath '\RFP results non PCs.csv']).VariableNames;
    % How many ROIs are at this Z-plane
    NumberofROIs = size(gfpResults,1);
    % Go through each ROI
    for n = 1:NumberofROIs
        % GFP: (integrated density - (area of selected cell * mean fluorescence of background readings)) / area
        gfpCTCF{k}(n) = (gfpResults(n,find(strcmp(gfpNames,'IntDen'))) - ...
            (gfpResults(n,find(strcmp(gfpNames,'Area'))) * mean(gfpBackground(end-2:end,find(strcmp(gfpNames,'Mean')))))) / ...
            gfpResults(n,find(strcmp(gfpNames,'Area')));
        % RFP: (integrated density - (area of selected cell * mean fluorescence of background readings)) / area
        rfpCTCF{k}(n) = (rfpResults(n,find(strcmp(rfpNames,'IntDen'))) - ...
            (rfpResults(n,find(strcmp(rfpNames,'Area'))) * mean(rfpBackground(end-2:end,find(strcmp(rfpNames,'Mean')))))) / ...
            rfpResults(n,find(strcmp(rfpNames,'Area')));
        % Ratiometric analysis with CTCF
        ratio{k}(n) = rfpCTCF{k}(n)/gfpCTCF{k}(n);
    end
end

% Make histograms of RFP/GFP ratio for each Z-plane
for k = 1:length(ratio)
    figure; 
    histogram(ratio{k},'BinWidth',0.1,'FaceColor','k'); 
    hold on; 
    title(['Z-plane ' num2str(k)]);
    xlabel('CTCF_{RFP} / CTCF_{GFP}');
    ylabel('# ROIs');
    xline(1); xlim([0 2]);
    hold off;
end

ratio{k}';

%% Input data for EBC PCs, non-EBC PCs, & non-PCs histograms

EBC_PCs = [];
nonEBC_PCs = [];
non_PCs = [];

%% Make EBC PCs, non-EBC PCs, & non-PCs histograms

figure; histogram(EBC_PCs,'BinWidth',0.1,'FaceColor','k'); hold on; title('EBC PCs');
xlabel('CTCF_{RFP} / CTCF_{GFP}'); ylabel('# ROIs'); xline(1); xlim([0 2]); hold off;

figure; histogram(nonEBC_PCs,'BinWidth',0.1,'FaceColor','k'); hold on; title('Non-EBC PCs');
xlabel('CTCF_{RFP} / CTCF_{GFP}'); ylabel('# ROIs'); xline(1); xlim([0 2]); hold off;

figure; histogram(non_PCs,'BinWidth',0.1,'FaceColor','k'); hold on; title('Non-PCs');
xlabel('CTCF_{RFP} / CTCF_{GFP}'); ylabel('# ROIs'); xline(1); xlim([0 2]); hold off;
