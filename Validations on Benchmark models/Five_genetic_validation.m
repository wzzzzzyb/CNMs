%% ===================== Script Overview =====================
% File: Five_genetic_validation.m
%
% Description: This script analyzes the causal
% network markers (CNMs) for five-gene genetic network, a typical benchmark
% model for detecting early-warning signal ahead of critical transition.
%
% Inputs: 
%   - Variable, Dynamics, and parameters
%
% Outputs:
%   - Plot: Diagram for CNMs of five-gene genetic network in Main text and supplemental material
%
% Dependencies:
%   - MATLAB R2023a on author computer (or any suitable versions)
%
% Usage:
%   Run this script directly (click "Run" or type script name in the command window)
% =========================================================

%% Step 0: Identify the dynamics
clear;close all;clc
dt = 0.008;
steps = 1e3;
par_num = 2*1e3; % Parallel trajectories
P = [0.01:0.01:0.1,0.2:0.1:1];
p_num = length(P);

t = 0:dt:dt*(steps-1); % Time
D = 1e-5; % Noise strength

stab = [1,0,1,3,2]; % Attractor

%% Step 1: Langevin simulation
S = zeros(steps,5,par_num,p_num);
S(1,:,:,:) = repmat(stab,[1,1,par_num,p_num]) + 0.01 * rand(1,5,par_num,p_num); % Initial value

for pnum = 1:p_num
    for num = 1:par_num
        for i = 1:steps-1
            S(i+1,:,num,pnum) = S(i,:,num,pnum)...
                +force(S(i,:,num,pnum),P(pnum),0)*dt...
                +normrnd(0,sqrt(2*D*dt),1,5); % Simulation
        end
    end
end

barS = S - repmat(stab,[steps,1,par_num,p_num]); % Centralization
toc

%% Step 2: Verify the dominant group and the non-dominant group

dominant_group = [1,2];
non_dominant_group = [3,4,5];

% Position correspondence between DG-NDG
DG_NDG_pos = NaN(length(dominant_group)*length(non_dominant_group),2);
k = 1;
for i = 1:length(dominant_group)
    for j = 1:length(non_dominant_group)
        DG_NDG_pos(k,1) = dominant_group(i);
        DG_NDG_pos(k,2) = non_dominant_group(j);
        k = k + 1;
    end
end

%% Step 3: Causal Network Marker - Granger Causality and Transfer Entropy (CNM-GC and CNM-TE)
% Open Parallel Computing Toolbox to speed up, please refer to your own computer 
% configuration determine the number of thread pools (default to 12).
tic;
% parpool('Processes', 12);

CNM_GC = NaN(p_num,1);
CNM_TE = NaN(p_num,1);

for pnum = 1:p_num
    GC = NaN(steps-1,length(DG_NDG_pos));
    TE = NaN(steps-1,length(DG_NDG_pos));
    for i = 1:length(DG_NDG_pos)
        aux1 = NaN(steps-1,1); % Auxiliary variables introduced to execute parfor
        aux2 = NaN(steps-1,1);
        for j = 1:steps-1
            % X: DG data, Y: NDG data, X_t, Y_t correspond to the time window
            % data of the i-th pair at time point j, and X_T, Y_T are their
            % data for time point j+1.
            X_t = reshape(barS(j,DG_NDG_pos(i,1),:,pnum),[par_num,1]);
            X_T = reshape(barS(j+1,DG_NDG_pos(i,1),:,pnum),[par_num,1]);
            Y_t = reshape(barS(j,DG_NDG_pos(i,2),:,pnum),[par_num,1]);
            Y_T = reshape(barS(j+1,DG_NDG_pos(i,2),:,pnum),[par_num,1]);
    
            % Calculate GC and TE by Granger_causality and Transfer_entropy function

            % GC
            [~,~,residual] = regress(Y_T,Y_t);          % Fitting of H0 model
            var_xi0_x2y = residual' * residual / length(residual);
            [~,~,residual] = regress(Y_T,[Y_t,X_t]);    % Fitting of H1 model
            var_xi1_x2y = residual' * residual / length(residual);
            aux1(j) = log(var_xi0_x2y / var_xi1_x2y);

            % TE
            Tra = [X_T X_t Y_T Y_t];  % matrix for joint distribution estimation
            pieces = 50;

            % Discretize data into bins
            Tra_max = max(Tra,[],1);
            Tra_min = min(Tra,[],1);
            delta = (Tra_max-Tra_min)/pieces;  % Bin width calculation (50 pieces)

            % Convert to integer indices for accumarray
            Tra_acc = floor((Tra-repmat(Tra_min,[length(X_t),1]))./repmat(delta,[length(X_t),1]))+1;
            Tra_acc(Tra_acc == pieces+1) = pieces;  % Handle edge cases
            % 1D count, p(x_t), p(y_t)
            dist1 = zeros(pieces,2);  
            dist1(:,1) = accumarray(Tra_acc(:,2),1);
            dist1(:,2) = accumarray(Tra_acc(:,4),1);
            dist1 = dist1./sum(dist1);

            % 2D count, p(x_T,x_t), p(y_T,y_t), p(x_t,y_t)
            dist2 = zeros(pieces,pieces,3);
            dist2(:,:,1) = accumarray(Tra_acc(:,[1,2]),1);
            dist2(:,:,2) = accumarray(Tra_acc(:,[3,4]),1);
            dist2(:,:,3) = accumarray(Tra_acc(:,[2,4]),1);
            dist2 = dist2./sum(dist2, [1 2]);

            % 3D count, p(x_T,x_t,y_t), p(y_T,x_t,y_t)
            dist3 = zeros(pieces,pieces,pieces,2);
            dist3(:,:,:,1) = accumarray(Tra_acc(:,[1,2,4]),1);
            dist3(:,:,:,2) = accumarray(Tra_acc(:,[2,3,4]),1);
            dist3 = dist3./sum(dist3, [1 2 3]);

            % Entropy terms
            total_f = 0;
            for k1 = 1:pieces
                for k2 = 1:pieces
                    if dist2(k1,k2,1)~=0 && dist1(k2,1)~=0
                       total_f = total_f - dist2(k1,k2,1)*log2(dist2(k1,k2,1)/dist1(k2,1));
                    end
                end
            end
            total_s = 0;
            for k1 = 1:pieces
                for k2 = 1:pieces
                    for k3 = 1:pieces
                        if dist3(k1,k2,k3,1)~=0 && dist2(k2,k3,3)~=0
                           total_s = total_s - dist3(k1,k2,k3,1)*log2(dist3(k1,k2,k3,1)/dist2(k2,k3,3));
                        end
                    end
                end
            end
            aux2(j) = total_f - total_s;
        end
        GC(:,i) = aux1;    % GC(X -> Y)
        TE(:,i) = aux2;     % TE(X -> Y)
    end

    mean_GC = mean(GC(steps/2:end,:),1);
    mean_TE = mean(TE(steps/2:end,:),1);
    CNM_GC(pnum) = 1 ./ mean(mean_GC);
    CNM_TE(pnum) = 1 ./ mean(mean_TE);
end

toc;

%% Step 4: Figure plot
figure()

plot(P,CNM_GC,'o-','Color',[160/255,201/255,235/255],'LineWidth',2);

xlabel('Parameter P','FontSize',18)
ylabel('Causality Biomarker','FontSize',18)
xlim([0,1.05])
legend('CNM-GC','FontSize',16,'Box','off','Location','NorthWest')

% 调整轴
box off
ax = gca;
ax.FontSize = 16;
ax.LineWidth = 1;
ax.XDir = 'reverse';

figure()

plot(P,CNM_TE,'o-','Color',[249/255,218/255,173/255],'LineWidth',2);

% 调整标签
xlabel('Parameter P','FontSize',18)
ylabel('Causality Biomarker','FontSize',18)
xlim([0,1.05])
legend('CNM-TE','FontSize',16,'Box','off','Location','NorthWest')

% 调整轴
box off
ax = gca;
ax.FontSize = 16;
ax.LineWidth = 1;
ax.XDir = 'reverse';
