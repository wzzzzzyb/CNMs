%% ===================== Script Overview =====================
% File: Five_genetic_validation.m
%
% Description: This script analyzes the causal network markers (CNMs) for
% five-gene genetic network, a typical benchmark model for detecting
% early-warning signal ahead of critical transition.
%
% Inputs: 
%   - No additional input required
%
% Outputs:
%   - Plot: Diagram for CNMs of five-gene genetic network
%
% Dependencies:
%   - MATLAB R2023a on author computer (or any suitable versions)
%   - Parallel Computing Toolbox
%   - Custom functions (force.m)
%
% Usage:
%   Run this script directly (click "Run" or type script name in the command window)
%
% Project: https://github.com/wzzzzzyb/CNMs | Last Updated: 2025-05-29
% =========================================================
%% Step 0: Setting the dynamics
clear;close all;clc
dt = 0.008;                     % Time steps
steps = 1000;                   % Simulation steps
par_num = 2000;                 % Parallel trajectories
D = 1e-5;                       % Noise strength
P = [0.01:0.01:0.1,0.2:0.1:1];  % Selected parameters

p_num = length(P);
t = 0:dt:dt*(steps-1);          % Time
stab = [1,0,1,3,2];             % Attractor

%% Step 1: Langevin simulation
% It costs nearly 1-2mins. Please note that each result has a certain degree
% of randomness.
tic
S = zeros(steps,5,par_num,p_num);
S(1,:,:,:) = repmat(stab,[1,1,par_num,p_num]) + 0.01 * rand(1,5,par_num,p_num); % Initial value
for pnum = 1:p_num
    fprintf('Simulating %d/%d\n', pnum, p_num);
    for num = 1:par_num
        for i = 1:steps-1
            S(i+1,:,num,pnum) = S(i,:,num,pnum)...
                + force(S(i,:,num,pnum),P(pnum),0)*dt...
                + normrnd(0,sqrt(2*D*dt),1,5);      % Simulation
        end
    end
end
barS = S - repmat(stab,[steps,1,par_num,p_num]);    % Centralization
toc
fprintf('Simulations Complete!');

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
tic; parpool('Processes', min(12, feature('numcores')));

CNM_GC = NaN(p_num,1);
CNM_TE = NaN(p_num,1);
parfor pnum = 1:p_num
    GC = NaN(steps-1,length(DG_NDG_pos));
    TE = NaN(steps-1,length(DG_NDG_pos));
    for i = 1:length(DG_NDG_pos)
        Aux1 = NaN(steps-1,1); % Auxiliary variables introduced to execute parfor
        Aux2 = NaN(steps-1,1);
        for j = steps/2:steps-1
            % X: DG data, Y: NDG data, X_t, Y_t correspond to the time window
            % data of the i-th pair at time point j, and Y_T are the data
            % for time point j+1.
            X_t = reshape(barS(j,DG_NDG_pos(i,1),:,pnum),[par_num,1]);
            Y_t = reshape(barS(j,DG_NDG_pos(i,2),:,pnum),[par_num,1]);
            Y_T = reshape(barS(j+1,DG_NDG_pos(i,2),:,pnum),[par_num,1]);
    
            % GC
            [~,~,residual] = regress(Y_T,Y_t);          % Fitting of H0 model
            var_xi0_x2y = residual' * residual / length(residual);
            [~,~,residual] = regress(Y_T,[Y_t,X_t]);    % Fitting of H1 model
            var_xi1_x2y = residual' * residual / length(residual);
            Aux1(j) = log(var_xi0_x2y / var_xi1_x2y);

            % TE
            pieces = 50;
            Tra = [X_t Y_T Y_t];
            scale = max(max(Tra,[],'all'),-min(Tra,[],'all'));
            Tra_max = scale * ones(1,3);
            Tra_min = -scale * ones(1,3);
            delta = (Tra_max-Tra_min)/pieces;  % Bin width calculation
            
            % Convert to integer indices for accumarray
            Tra_acc = floor((Tra-repmat(Tra_min,[length(X_t),1]))./repmat(delta,[length(X_t),1]))+1;
            Tra_acc(Tra_acc == pieces+1) = pieces;  % Handle edge cases
            
            % 1D marginal distributions
            dist1 = zeros(pieces,2);
            aux1 = accumarray(Tra_acc(:,1),1); % p(x_t)
            dist1(1:size(aux1,1),1) = aux1;  
            aux1 = accumarray(Tra_acc(:,3),1); % p(y_t)
            dist1(1:size(aux1,1),2) = aux1;  
            dist1 = dist1./sum(dist1);  % Normalize
        
            % 2D joint distributions
            dist2 = zeros(pieces,pieces,2);
            aux2 = accumarray(Tra_acc(:,[2,3]),1);  % p(y_T,y_t)
            dist2(1:size(aux2,1),1:size(aux2,2),1) = aux2;
            aux2 = accumarray(Tra_acc(:,[1,3]),1);  % p(x_t,y_t)
            dist2(1:size(aux2,1),1:size(aux2,2),2) = aux2;
            dist2 = dist2./sum(dist2, [1 2]);  % Normalize
            
            % 3D joint distributions
            dist3 = zeros(pieces,pieces,pieces);
            aux3 = accumarray(Tra_acc(:,[1,2,3]),1);  % p(x_t,y_T,y_t)
            dist3(1:size(aux3,1),1:size(aux3,2),1:size(aux3,3)) = aux3;
            dist3 = dist3./sum(dist3, [1 2 3]);  % Normalize
                  
            % Conditional entropy terms
            total_f = 0;  % H(Y_T|Y_t)
            for k1 = 1:pieces
                for k2 = 1:pieces
                    if dist2(k1,k2,1)~=0 && dist1(k2,2)~=0
                       total_f = total_f - dist2(k1,k2,1)*log2(dist2(k1,k2,1)/dist1(k2,2));
                    end
                end
            end
            
            % Joint entropy terms
            total_s = 0;  % H(Y_T|Y_t,X_t)
            for k1 = 1:pieces
                for k2 = 1:pieces
                    for k3 = 1:pieces
                        if dist3(k1,k2,k3)~=0 && dist2(k1,k3,2)~=0
                           total_s = total_s - dist3(k1,k2,k3)*log2(dist3(k1,k2,k3)/dist2(k1,k3,2));
                        end
                    end
                end
            end
            % Final TE Calculation
            Aux2(j) = total_f - total_s;  % TE_{X→Y} = H(Y_T|Y_t) - H(Y_T|Y_t,X_t)
        end
        GC(:,i) = Aux1;     % GC(X -> Y)
        TE(:,i) = Aux2;     % TE(X -> Y)
    end
    mean_GC = mean(GC(steps/2:end,:),1);
    mean_TE = mean(TE(steps/2:end,:),1);
    CNM_GC(pnum) = 1 ./ mean(mean_GC);
    CNM_TE(pnum) = 1 ./ mean(mean_TE);
end
delete(gcp('nocreate')); toc;    % Close the Parallel Computing Toolbox

%% Step 4: Figure plot
figure()    % GC
plot(P,CNM_GC,'o-','Color',[160/255,201/255,235/255],'LineWidth',2);
xlabel('Parameter P','FontSize',18)
ylabel('Causality Biomarker','FontSize',18)
xlim([0,1.05])
legend('CNM-GC','FontSize',16,'Box','off','Location','NorthWest')
box off
ax = gca;
ax.FontSize = 16;
ax.LineWidth = 1;
ax.XDir = 'reverse';

figure()    % TE
plot(P,CNM_TE,'o-','Color',[249/255,218/255,173/255],'LineWidth',2);
xlabel('Parameter P','FontSize',18)
ylabel('Causality Biomarker','FontSize',18)
xlim([0,1.05])
legend('CNM-TE','FontSize',16,'Box','off','Location','NorthWest')
box off
ax = gca;
ax.FontSize = 16;
ax.LineWidth = 1;
ax.XDir = 'reverse';