%% ===================== Script Overview =====================
% File: iEEG_data_validation.m
%
% This script is designed to analyze intracranial electroencephalography
% (iEEG) data by computing dynamic network biomarker (DNB) and causal
% network markers (CNMs). The analysis focuses on identifying critical
% transitions during epileptic seizures through network dynamics. The
% output includes comprehensive visualizations of DNB trends, CNM-GC
% (Granger causality), CNM-TE (transfer entropy) patterns, and the variance
% of network dynamics.
%
% Input Data: 
%   - Path: ..iEEG dataset/ID../Sz..mat
%   - Variable: EEG (N*M matrix, N=data length, M=electrode channels)
%
% Outputs:
%   - Plot: Diagram for iEEG datasets in Main text and Supporting
%   Information
%
% Dependencies:
%   - MATLAB R2023a on author's computer (or any suitable versions)
%   - Parallel Computing Toolbox
%   - Custom functions (Granger_causality, Transfer_entropy)
%
% Usage:
%   1. Ensure data files is downloaded
%   2. Modify ID and Sz to analyze different patient/seizure cases
%   3. Run this script directly (click "Run" or type script name in the command window)
%
% Author: Zezhou Wang <zezhou_wang@foxmail.com> | Last Updated: 2025-05-26
% Project: https://github.com/wzzzzzyb/CNMs
% =========================================================

%% Step 0: Download the iEEG Dataset
% Read the "iEEG Dataset Download Instuction.pdf" to prepare for the
% necessary data

%% Step 1: Initialize Environment and Load iEEG Data
clear; close all; clc;  % Clear workspace, figures, and command window
Path = pwd;             % Specify the path of the m file (default to the current working directory)

ID = 1; Sz = 4;         % Select your interested dataset
% In Figure 4, the cases are a) ID1Sz4, b) ID7Sz2, c) ID5Sz2, d) ID16Sz1,
% you can use the following statement, for each experiment:
% ID = 1; Sz = 4;
% ID = 7; Sz = 2;
% ID = 5; Sz = 2;
% ID = 16; Sz = 1;

load([Path '\iEEG dataset\ID' num2str(ID) '\Sz' num2str(Sz)]);  % load the mat
disp(['Processing: ID' num2str(ID) '_Sz'  num2str(Sz) '...']);
%% Step 2: Data Preprocessing
f = 512;            % Write frequency (see source website for details)
[N M] = size(EEG);  % Data row and column size, N represents data length, and M represents the electrode channels
time_window = f;    % Set the data per second as a time window
T = floor(N/time_window);  % The number of time windows after dividing the data, and each time window represents 1 second

% Convert EEG data to T * timew_window * M size
EEG_tw = NaN(T,time_window,M);
for i = 1:M
    for j = 1:T
        EEG_tw(j,:,i) = EEG((j-1)*time_window+1:j*time_window, i);
    end
end

% In this iEEG short-term dataset, data from 3 minutes before the onset of the 
% seizure and 3 minutes after the end of the seizure are included.
seizure_begin = 180*512/time_window;    % Time point of seizure begin
seizure_end = T-180*512/time_window;    % Time point of seizure end

% Calculate the covariance between nodes in the i-th time window and store it as a 3-order tensor as a whole
aux = permute(EEG_tw, [2 3 1]);         % Transpose EEG_tw to facilitate covariance calculation
covariance_tensor = NaN(M,M,T-1);       % Covariance tensor
for i = 1:T-1
    covariance_tensor(:,:,i) = cov(aux(:,:,i));
end

%% Step 3: Determine the Dominant Group (DG) and Non-dominant Group (NDG)
% We use the K-means clustering algorithm here to cluster the average variance 
% of nodes during the onset of the seizure. For a comparison of the performance
% of other clustering algorithms, please refer to Supplemental Material.

% Calculate the mean variance during seizure
mean_var = NaN(M,1);
for i = 1:M
    mean_var(i,1) = mean(covariance_tensor(i,i,seizure_begin:seizure_end));
end

% Obtain DG (high variance) and NDG (low variance) via K-means clustering algorithm
[idx C] = kmeans(mean_var,2);
A = 1:M;
if C(1) > C(2)
    dominant_group = A(idx == 1)';
    non_dominant_group = A(idx == 2)';
else
    dominant_group = A(idx == 2)';
    non_dominant_group = A(idx == 1)';
end

% Provide the position combinations between DG and DG, as well as between DG and NDG
s1 = length(dominant_group)*(length(dominant_group)-1)/2;   % Size of combinations between DG-DG
s2 = length(dominant_group)*length(non_dominant_group);     % Size of combinations between DG-NDG

DG_DG_pos  = NaN(s1,2);     % Position correspondence between DG-DG
k = 1;
for i = 1:length(dominant_group)
    for j = i+1:length(dominant_group)
        DG_DG_pos(k,1) = dominant_group(i);
        DG_DG_pos(k,2) = dominant_group(j);
        k = k + 1;
    end
end

DG_NDG_pos = NaN(s2,2);     % Position correspondence between DG-NDG
k = 1;
for i = 1:length(dominant_group)
    for j = 1:length(non_dominant_group)
        DG_NDG_pos(k,1) = dominant_group(i);
        DG_NDG_pos(k,2) = non_dominant_group(j);
        k = k + 1;
    end
end
%% Step 4: Dynamical Network Marker (DNB)
% DNB(Net) := SD_d .* PCC_d ./ PCC_o
SD_d = zeros(T-1,1); % Calculate the SD_d
for i = 1:length(dominant_group)
    SD_d = SD_d + reshape(covariance_tensor(dominant_group(i), dominant_group(i),:),T-1,1);
end
SD_d = SD_d / length(dominant_group);

PCC_d = zeros(T-1,1); % Calculate the PCC_d

for i = 1:s1
    PCC_d = PCC_d + abs(reshape(covariance_tensor(DG_DG_pos(i,1), DG_DG_pos(i,2),:),T-1,1));
end
PCC_d = PCC_d / s1;

PCC_o = zeros(T-1,1); % Calculate the PCC_o
for i = 1:s2
    PCC_o = PCC_o + abs(reshape(covariance_tensor(DG_NDG_pos(i,1), DG_NDG_pos(i,2),:),T-1,1));
end
PCC_d = PCC_d / s2;
DNB = SD_d .* PCC_d ./ PCC_o;

%% Step 5: Causal Network Marker - Granger Causality and Transfer Entropy (CNM-GC and CNM-TE)
% Open Parallel Computing Toolbox to speed up, please refer to your own computer 
% configuration determine the number of thread pools (default maximum to 12).
tic; 
parpool('Processes', min(12, feature('numcores')));   

GC = NaN(length(DG_NDG_pos),T-1);
TE = NaN(length(DG_NDG_pos),T-1);
parfor i = 1:length(DG_NDG_pos)
    aux1 = NaN(T-1,1); % Auxiliary variables introduced to execute parfor
    aux2 = NaN(T-1,1);
    for j = 1:T-1
        % X: DG data, Y: NDG data, X_t, Y_t correspond to the time window
        % data of the i-th pair at time point j, and X_T, Y_T are their
        % data for time point j+1.
        X_t = EEG_tw(j,  :,DG_NDG_pos(i,1))';
        X_T = EEG_tw(j+1,:,DG_NDG_pos(i,1))';
        Y_t = EEG_tw(j,  :,DG_NDG_pos(i,2))';
        Y_T = EEG_tw(j+1,:,DG_NDG_pos(i,2))';

        % Calculate GC and TE by Granger_causality and Transfer_entropy
        % function
        [aux1(j), ~] = Granger_causality(X_t, X_T, Y_t, Y_T);
        [aux2(j), ~] = Transfer_entropy(X_t, X_T, Y_t, Y_T, 50);
    end
    GC(i,:) = aux1;     % GC(X -> Y)
    TE(i,:) = aux2;     % TE(X -> Y)
end
CNM_GC = 1 ./ mean(GC(:,:), 1)';
CNM_TE = 1 ./ mean(TE(:,:), 1)';

delete(gcp('nocreate')); toc;    % Close the Parallel Computing Toolbox
%% Step 6: Visualization
% Output subplot in Supplemental Material, including CNM-GC, CNM-TE, DNB,
% node variance, DG node variance, and NDG node variance.
short_cycle = 5;        % Short cycle moving average length (blue line)
long_cycle = 12;        % Long cycle moving average length (red line)

% Center the image and set the length and width
scnsize  = get(0,'ScreenSize'); 
len = 800;
wid = 600;
set(gcf,'Position', [scnsize(1,3)/2-len/2, scnsize(1,4)/2-wid/2,len, wid]);

sgtitle(['ID' num2str(ID) 'Sz' num2str(Sz)]) % General Title
% ===================== GC diagram =====================
subplot(3,2,1)
plot(1:T-1,CNM_GC,'Color',[0.7,0.7,0.7])
hold on
plot(1:T-1,movmean(CNM_GC,short_cycle),'Color',[0.2,0.3,0.49],'LineWidth',1)
plot(1:T-1,movmean(CNM_GC,long_cycle),'Color',[0.85,0.33,0.1],'LineWidth',2)
% Seizure begin annotation (blue pointed corner, consistent below)
scatter(seizure_begin, max(CNM_GC,[],'all'),100, 'vr', 'filled','MarkerFaceColor',[0.2,0.3,0.49]);
% Seizure end annotation (red pointed corner, consistent below)
scatter(seizure_end, max(CNM_GC,[],'all'),100, 'vr', 'filled','MarkerFaceColor',[0.85,0.33,0.1]);
box off
ax = gca;
ax.FontSize = 9;
ax.LineWidth = 1;
xlim([1 T-1]);
xlabel('t(s)');
ylabel('CNM-GC');
hold off

% ===================== TE diagram =====================
subplot(3,2,3)
plot(1:T-1,CNM_TE,'Color',[0.7,0.7,0.7])
hold on
plot(1:T-1,movmean(CNM_TE,short_cycle),'Color',[0.2,0.3,0.49],'LineWidth',1)
plot(1:T-1,movmean(CNM_TE,long_cycle),'Color',[0.85,0.33,0.1],'LineWidth',2)
scatter(seizure_begin, max(CNM_TE,[],'all'),100, 'vr', 'filled','MarkerFaceColor',[0.2,0.3,0.49]);
scatter(seizure_end, max(CNM_TE,[],'all'),100, 'vr', 'filled','MarkerFaceColor',[0.85,0.33,0.1]);
box off
ax = gca;
ax.FontSize = 9;
ax.LineWidth = 1;
xlim([1 T-1]);
xlabel('t(s)');
ylabel('CNM-TE');
hold off

% ===================== DNB diagram =====================
subplot(3,2,5)
plot(1:T-1,DNB,'Color',[0.7,0.7,0.7])
hold on
plot(1:T-1,movmean(DNB,short_cycle),'Color',[0.2,0.3,0.49],'LineWidth',1)
plot(1:T-1,movmean(DNB,long_cycle),'Color',[0.85,0.33,0.1],'LineWidth',2)
scatter(seizure_begin, max(DNB,[],'all'),100, 'vr', 'filled','MarkerFaceColor',[0.2,0.3,0.49]);
scatter(seizure_end, max(DNB,[],'all'),100, 'vr', 'filled','MarkerFaceColor',[0.85,0.33,0.1]);
box off
ax = gca;
ax.FontSize = 9;
ax.LineWidth = 1;
xlim([1 T-1]);
xlabel('t(s)');
ylabel('DNB');
hold off

% ===================== Variance diagram =====================
subplot(3,2,2)
hold on;
max_element = -inf;
for i = 1:M
    plot(1:T-1,reshape(covariance_tensor(i,i,:),T-1,1))
    if max_element < max(covariance_tensor(i,i,:))
        max_element = max(covariance_tensor(i,i,:));
    end
end
scatter(seizure_begin, max_element,100, 'vr', 'filled','MarkerFaceColor',[0.2,0.3,0.49]);
scatter(seizure_end, max_element,100, 'vr', 'filled','MarkerFaceColor',[0.85,0.33,0.1]);
box off
ax = gca;
ax.FontSize = 9;
ax.LineWidth = 1;
xlim([1 T-1]);
xlabel('t(s)');
ylabel('Variance of Each Node');
hold off

% ===================== DG Variance diagram =====================
subplot(3,2,4)
hold on;
for i = dominant_group'
    plot(1:T-1,reshape(covariance_tensor(i,i,:),T-1,1))
end
scatter(seizure_begin, max_element,100, 'vr', 'filled','MarkerFaceColor',[0.2,0.3,0.49]);
scatter(seizure_end, max_element,100, 'vr', 'filled','MarkerFaceColor',[0.85,0.33,0.1]);
box off
ax = gca;
ax.FontSize = 9;
ax.LineWidth = 1;
xlim([1 T-1]);
xlabel('t(s)');
ylabel('Variance of Nodes in DG');
hold off

% ===================== NDG Variance diagram =====================
subplot(3,2,6)
hold on;
for i = non_dominant_group'
    plot(1:T-1,reshape(covariance_tensor(i,i,:),T-1,1))
end
scatter(seizure_begin, max_element,100, 'vr', 'filled','MarkerFaceColor',[0.2,0.3,0.49]);
scatter(seizure_end, max_element,100, 'vr', 'filled','MarkerFaceColor',[0.85,0.33,0.1]);
box off
ax = gca;
ax.FontSize = 9;
ax.LineWidth = 1;
xlim([1 T-1]);
xlabel('t(s)');
ylabel('Variance of Nodes in NDG');
hold off
% ==============================================================