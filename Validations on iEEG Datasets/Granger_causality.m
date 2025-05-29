function [GC_X2Y, GC_Y2X] = Granger_causality(X_t, X_T, Y_t, Y_T)
% ===================== Function Overview =====================
% Granger_causality - Calculate bidirectional Granger causality (GC) between variables X and Y
%
% Detailed Description:
%   This function computes the Granger causality (GC) between two variables X and Y in both
%   directions (X→Y and Y→X). Granger causality assesses whether the historical values of one
%   variable (e.g., X) improve the prediction of another variable (e.g., Y) beyond what can be
%   predicted by the historical values of Y alone. The GC value is derived from the ratio of
%   prediction errors between restricted (H0) and unrestricted (H1) linear regression models.
%
% Algorithm:
%   1. For X→Y:
%     - H0 Model: Predict Y_T using only Y_t (linear regression: Y_T ~ a*Y_t + ξ₀).
%     - H1 Model: Predict Y_T using both Y_t and X_t (multivariate regression: Y_T ~ b*Y_t + c*X_t + ξ₁).
%     - Compute Granger causality as log(E[ξ₀²]/E[ξ₁²]), where E[ξ₀²] and E[ξ₁²] are the residual 
%     variances of H0 and H1 models (E[ξ]≈0, so E[ξ²]≈Var(ξ)), respectively.
%   2. For Y→X: Repeat the above steps with X and Y swapped.
%
% Input Arguments:
%   X_t - [N×1 vector] Observations of X at time t.
%   X_T - [N×1 vector] Observations of X at time T (T > t).
%   Y_t - [N×1 vector] Observations of Y at time t.
%   Y_T - [N×1 vector] Observations of Y at time T (T > t).
%
% Output Arguments:
%   GC_X2Y - [Non-negative scalar] Granger causality from X to Y. Larger values indicate stronger causal influence.
%   GC_Y2X - [Non-negative scalar] Granger causality from Y to X.
%
% Example (Y Granger-causes X):
%   X_t = rand(512,1); % Random historical data for X
%   Y_t = rand(512,1); % Random historical data for Y
%   X_T = Y_t + 0.1*rand(512,1); % X_T depends on Y_t (Y causes X)
%   Y_T = rand(512,1); % Y_T is random noise
%   [GC_X2Y, GC_Y2X] = Granger_causality(X_t, X_T, Y_t, Y_T)
%   % Expected result: GC_Y2X >> GC_X2Y
%
% Notes:
%   - Assumes linear relationships and uses 1st-order autoregressive models (only t and T time points).
%   - The intercept term of the data is approximately non-existent, as paper presents a linear 
%   expansion around zero steady state.
%   - Larger sample sizes (N) improve reliability.
%
% Project: https://github.com/wzzzzzyb/CNMs | Last Updated: 2025-05-29
% =============================================================
    % ------------------------ X->Y ------------------------
    [~,~,residual] = regress(Y_T,Y_t);          % Fitting of H0 model
    var_xi0_x2y = residual' * residual / length(residual);
    [~,~,residual] = regress(Y_T,[Y_t,X_t]);    % Fitting of H1 model
    var_xi1_x2y = residual' * residual / length(residual);
    GC_X2Y = log(var_xi0_x2y / var_xi1_x2y);    % log(E[ξ₀²]/E[ξ₁²])
    % ------------------------ Y->X ------------------------
    [~,~,residual] = regress(X_T,X_t);          % Fitting of H0 model
    var_xi0_y2x = residual' * residual / length(residual);
    [~,~,residual] = regress(X_T,[X_t,Y_t]);    % Fitting of H1 model
    var_xi1_y2x = residual' * residual / length(residual);
    GC_Y2X = log(var_xi0_y2x / var_xi1_y2x);    % log(E[ξ₀²]/E[ξ₁²])
    % -------------------------------------------------------
end