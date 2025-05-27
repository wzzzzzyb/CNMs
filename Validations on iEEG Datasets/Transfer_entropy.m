function [TE_X2Y, TE_Y2X] = Transfer_entropy(X_t, X_T, Y_t, Y_T, pieces)
% ===================== Function Overview =====================
% Transfer_entropy - Calculate bidirectional transfer entropy (TE) between variables X and Y
%  
% Detailed Description:
%   This function computes the transfer entropy (TE) between two variables X and Y in both 
%   directions (X→Y and Y→X). Transfer entropy quantifies the amount of information transferred 
%   from one process to another, capturing nonlinear causal relationships. The implementation 
%   follows the discrete binning approach to estimate probability distributions.
%
% Algorithm:
%   1. Data Discretization: 
%      - Discretize X_t, X_T, Y_t, Y_T into "pieces" bins each.
%      - Generate matrix Tra = [X_T, X_t, Y_T, Y_t].
%   2. Probability Estimation:
%      - Compute 1D marginal probabilities (p(x_t), p(y_t)).
%      - Compute 2D joint probabilities (p(x_T,x_t), p(y_T,y_t), p(x_t,y_t)).
%      - Compute 3D joint probabilities (p(x_T,x_t,y_t), p(x_t,y_T,y_t)).
%   3. Entropy Calculation:
%      - TE_{X→Y} = Σ p(y_T,y_t,x_t) * log2(p(y_T|y_t,x_t)/p(y_T|y_t))
%      - TE_{Y→X} = Σ p(x_T,x_t,y_t) * log2(p(x_T|x_t,y_t)/p(x_T|x_t))
%
% Input Arguments:  
%   X_t      - [N×1 vector] Observations of X at time t.
%   X_T      - [N×1 vector] Observations of X at time T (T > t).
%   Y_t      - [N×1 vector] Observations of Y at time t.
%   Y_T      - [N×1 vector] Observations of Y at time T (T > t).
%   pieces   - [Integer] Number of bins for discretizing data (typically 30-50).
%
% Output Arguments:  
%   TE_X2Y   - [Non-negative scalar] Transfer entropy from X to Y. Larger values indicate stronger information flow.
%   TE_Y2X   - [Non-negative scalar] Transfer entropy from Y to X.
%
% Example (Y transfers information to X):
%   X_t = randn(512,1);          % Random historical data for X
%   Y_t = randn(512,1);          % Random historical data for Y
%   X_T = X_t + Y_t + 0.1*randn(512,1);  % X_T depends on Y_t
%   Y_T = Y_t + 0.1*randn(512,1);
%   [TE_X2Y, TE_Y2X] = Transfer_entropy(X_t, X_T, Y_t, Y_T, 30);
%   % Expected result: TE_Y2X >> TE_X2Y
%
% Notes:
%   - Requires sufficient data length for stable probability estimation.
%   - Larger "pieces" values improve resolution but require more data.
%   - Negative TE values may arise from numerical errors; treat values < 0 as 0.
%
% Project: https://github.com/wzzzzzyb/CNMs
% =============================================================
    %% Data Preparation
    Tra = [X_T X_t Y_T Y_t];  % matrix for joint distribution estimation
    
    % Discretize data into bins
    Tra_max = max(Tra,[],1);
    Tra_min = min(Tra,[],1);
    delta = (Tra_max-Tra_min)/pieces;  % Bin width calculation
    
    % Convert to integer indices for accumarray
    Tra_acc = floor((Tra-repmat(Tra_min,[length(X_t),1]))./repmat(delta,[length(X_t),1]))+1;
    Tra_acc(Tra_acc == pieces+1) = pieces;  % Handle edge cases
    
    %% Probability Estimation
    % 1D marginal distributions
    dist1 = zeros(pieces,2);  
    dist1(:,1) = accumarray(Tra_acc(:,2),1);  % p(x_t)
    dist1(:,2) = accumarray(Tra_acc(:,4),1);  % p(y_t)
    dist1 = dist1./sum(dist1);  % Normalize
    
    % 2D joint distributions
    dist2 = zeros(pieces,pieces,3);
    dist2(:,:,1) = accumarray(Tra_acc(:,[1,2]),1);  % p(x_T,x_t)
    dist2(:,:,2) = accumarray(Tra_acc(:,[3,4]),1);  % p(y_T,y_t)
    dist2(:,:,3) = accumarray(Tra_acc(:,[2,4]),1);  % p(x_t,y_t)
    dist2 = dist2./sum(dist2, [1 2]);  % Normalize
    
    % 3D joint distributions
    dist3 = zeros(pieces,pieces,pieces,2);
    dist3(:,:,:,1) = accumarray(Tra_acc(:,[1,2,4]),1);  % p(x_T,x_t,y_t)
    dist3(:,:,:,2) = accumarray(Tra_acc(:,[2,3,4]),1);  % p(x_t,y_T,y_t)
    dist3 = dist3./sum(dist3, [1 2 3]);  % Normalize
    
    %% Entropy Calculation
    % Conditional entropy terms
    total_f_1 = 0;  % H(Y_T|Y_t)
    total_f_2 = 0;  % H(X_T|X_t)
    for k1 = 1:pieces
        for k2 = 1:pieces
            if dist2(k1,k2,2)~=0 && dist1(k2,2)~=0
               total_f_1 = total_f_1 - dist2(k1,k2,2)*log2(dist2(k1,k2,2)/dist1(k2,2));
            end
            if dist2(k1,k2,1)~=0 && dist1(k2,1)~=0
               total_f_2 = total_f_2 - dist2(k1,k2,1)*log2(dist2(k1,k2,1)/dist1(k2,1));
            end
        end
    end
    
    % Joint entropy terms
    total_s_1 = 0;  % H(Y_T|Y_t,X_t)
    total_s_2 = 0;  % H(X_T|X_t,Y_t)
    for k1 = 1:pieces
        for k2 = 1:pieces
            for k3 = 1:pieces
                if dist3(k1,k2,k3,2)~=0 && dist2(k1,k3,3)~=0
                   total_s_1 = total_s_1 - dist3(k1,k2,k3,2)*log2(dist3(k1,k2,k3,2)/dist2(k1,k3,3));
                end
                if dist3(k1,k2,k3,1)~=0 && dist2(k2,k3,3)~=0
                   total_s_2 = total_s_2 - dist3(k1,k2,k3,1)*log2(dist3(k1,k2,k3,1)/dist2(k2,k3,3));
                end
            end
        end
    end
    
    %% Final TE Calculation
    TE_X2Y = total_f_1 - total_s_1;  % TE_{X→Y} = H(Y_T|Y_t) - H(Y_T|Y_t,X_t)
    TE_Y2X = total_f_2 - total_s_2;  % TE_{Y→X} = H(X_T|X_t) - H(X_T|X_t,Y_t)
end