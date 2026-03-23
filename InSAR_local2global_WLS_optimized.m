function [V_global, V_global_q] = InSAR_local2global_WLS_optimized( ...
    Lon, Lat, V_local, V_local_q, VUmed_int, VUmed_std_int, ndegree)
% ==============================================================
% InSAR_local2global_WLS_optimized
% --------------------------------------------------------------
% Purpose:
% Transforms local InSAR-derived vertical land motion (VLM) data
% into the global GNSS reference frame (IGS20) using a weighted
% polynomial regression model with full uncertainty propagation.
%
% Description:
% This function estimates and applies a correction surface (Δ-field)
% between InSAR local velocities and GNSS-interpolated velocities
% (VUmed_int). The correction is modeled as a polynomial surface
% fitted by Weighted Least Squares (WLS), where each observation
% weight reflects both InSAR and GNSS uncertainty.
%
% The method ensures a statistically consistent transformation from
% relative InSAR velocities (local reference frame) to absolute global
% velocities (GNSS reference frame), with rigorous propagation of
% uncertainty at each pixel, including the covariance between the
% InSAR input and the fitted correction.
%
% Transformation Model:
% Δ_i = VUmed_int_i − V_local_i
% w_i = 1 / (σ_InSAR_i² + σ_GNSS_i²)
% β̂ = (AᵀWA)⁻¹AᵀWΔ
% correction_i = B_iβ̂
% V_global_i = V_local_i + correction_i
%
% Uncertainty Propagation (COMPLETE):
% Cov(β̂) = (AᵀWA)⁻¹ · MSE
% Var(correction_i) = B_i · Cov(β̂) · B_iᵀ
% Cov(V_local_i, correction_i) = -σ²_InSAR_i · w_i · B_i · (AᵀWA)⁻¹ · A_iᵀ
% Var(V_global_i) = σ²_InSAR_i + Var(correction_i) + 2·Cov(V_local_i, correction_i)
%
% The covariance term accounts for the statistical dependency between
% the InSAR input (V_local) and the correction estimate, which arises
% because V_local appears in both the fitting residuals and the final
% transformation. This correlation is negative and reduces the total
% uncertainty for pixels used in the fitting process.
%
% Inputs:
% Lon, Lat           - InSAR pixel coordinates [deg] (N×1 vectors)
% V_local            - InSAR vertical velocities [mm/yr] (N×1 vector)
% V_local_q          - InSAR velocity standard deviation [mm/yr] (N×1)
% VUmed_int          - GNSS-interpolated velocities at pixel locations [mm/yr] (N×1)
% VUmed_std_int      - Interpolated GNSS standard deviation [mm/yr] (N×1)
% ndegree            - Polynomial degree (1, 2, or 3)
%                      1: planar surface
%                      2: quadratic surface
%                      3: cubic surface
%
% Outputs:
% V_global           - Global-frame InSAR vertical velocities [mm/yr] (N×1)
% V_global_q         - Propagated uncertainty per pixel [mm/yr] (N×1)
%
% Notes:
% • The weighting accounts for both GNSS and InSAR uncertainty.
% • Supports polynomial degrees 1 (planar), 2 (quadratic), and 3 (cubic).
% • Outliers in InSAR velocities are removed using the IQR method.
% • Fully vectorized; no loops over pixels for optimal performance.
% • Memory-efficient computation of covariance terms using diagonal extraction.
% • Covariance propagation ensures pixel-level uncertainty realism.
% • Accounts for correlation between InSAR input and correction estimate.
%
% Example:
% [V_global, V_global_q] = InSAR_local2global_weighted_uncertainty_ed3( ...
%     Lon, Lat, V_local, V_local_q, VUmed_int, VUmed_std_int, 2);
%
% Reference:
% Reshadati, M. (2025). "Transforming Local InSAR-Derived VLM to a
% Global Reference Frame Using Weighted GNSS-Integrated Uncertainty
% with Full Covariance Propagation."
%
% Author: Mahmoud Reshadati
% Date: 2025-10-30
% Updated: 2025-11-10 (Added complete covariance term in uncertainty propagation)
% Version: ed3
% ==============================================================

fprintf('Running Local-to-Global InSAR conversion (ed3) with full uncertainty propagation...\n\n');

%% Step 1: Outlier removal (InSAR)
Q1 = quantile(V_local, 0.25);
Q3 = quantile(V_local, 0.75);
IQR = Q3 - Q1;
lower_bound = Q1 - 1.5 * IQR;
upper_bound = Q3 + 1.5 * IQR;
outlier_idx = (V_local < lower_bound | V_local > upper_bound);

filtered_V_local = V_local;
filtered_V_local(outlier_idx) = nan;

fprintf('Outlier removal: %d pixels flagged as outliers\n', sum(outlier_idx));

%% Step 2: Compute Δ = GNSS − InSAR
df = VUmed_int - filtered_V_local;
valid_idx = ~isnan(df);

df_clean = df(valid_idx);
Lon_clean = Lon(valid_idx);
Lat_clean = Lat(valid_idx);
V_local_q_clean = V_local_q(valid_idx);
VUmed_std_clean = VUmed_std_int(valid_idx);

fprintf('Valid pixels for fitting: %d / %d\n', sum(valid_idx), numel(V_local));

%% Step 3: Combined variance and weights
sigma2_combined = (V_local_q_clean.^2) + (VUmed_std_clean.^2);
w = 1 ./ sigma2_combined;

%% Step 4: Design matrices for polynomial surface
switch ndegree
    case 1
        A = [ones(length(Lon_clean),1), Lon_clean, Lat_clean];
        B = [ones(length(Lon),1), Lon, Lat];
    case 2
        A = [ones(length(Lon_clean),1), Lon_clean, Lat_clean, ...
             Lon_clean.^2, Lat_clean.^2, Lon_clean.*Lat_clean];
        B = [ones(length(Lon),1), Lon, Lat, ...
             Lon.^2, Lat.^2, Lon.*Lat];
    case 3
        A = [ones(length(Lon_clean),1), Lon_clean, Lat_clean, ...
             Lon_clean.^2, Lat_clean.^2, Lon_clean.^3, Lat_clean.^3, ...
             Lon_clean.*Lat_clean.^2, Lon_clean.^2.*Lat_clean];
        B = [ones(length(Lon),1), Lon, Lat, ...
             Lon.^2, Lat.^2, Lon.^3, Lat.^3, ...
             Lon.*Lat.^2, Lon.^2.*Lat];
    otherwise
        error('ndegree must be 1, 2, or 3');
end

n_params = size(A, 2);
fprintf('Polynomial degree: %d (%d parameters)\n', ndegree, n_params);

%% Step 5: Weighted least squares solution
% Efficient WLS without full diag(W)
Aw = A .* sqrt(w);
bw = df_clean .* sqrt(w);
X = (Aw' * Aw) \ (Aw' * bw);

% Residuals and weighted MSE
residuals = df_clean - A * X;
MSE = sum(w .* residuals.^2) / (length(df_clean) - n_params);

% Covariance of coefficients
AtWA_inv = pinv(Aw' * Aw);
Cov_beta = AtWA_inv * MSE;

fprintf('Weighted MSE: %.4f (mm/yr)²\n', MSE);
fprintf('Weighted RMSE: %.4f mm/yr\n', sqrt(MSE));

%% Step 6: Apply correction and propagate uncertainty (VECTORIZED)
shift = B * X;
V_global = V_local + shift;

% --- 6.1 Variance from polynomial model uncertainty ---
Var_corr = sum((B * Cov_beta) .* B, 2);  % [N_total × 1]

% --- 6.2 Covariance correction term (VECTORIZED, memory-efficient) ---
% Extract design matrix rows for pixels used in fitting
B_valid = B(valid_idx, :);  % [N_valid × p]
A_valid = A;                 % [N_valid × p] (already cleaned)

% Compute diagonal of B_valid * AtWA_inv * A_valid' efficiently:
% Using the identity: diag(B*M*A') = sum((B*M) .* A, 2)
diag_BA = sum((B_valid * AtWA_inv) .* A_valid, 2);  % [N_valid × 1]

% Covariance term: -2 * σ²_InSAR * w * diag_BA
cov_terms = -2 * (V_local_q_clean.^2) .* w .* diag_BA;  % [N_valid × 1]

% --- 6.3 Combine all variance components ---
Var_global = V_local_q.^2 + Var_corr;  % Initialize with InSAR + model variance
Var_global(valid_idx) = Var_global(valid_idx) + cov_terms;  % Add covariance correction

% --- 6.4 Final uncertainty (ensure non-negative) ---
V_global_q = sqrt(max(Var_global, 0));

%% Diagnostics
fprintf('\n');
fprintf('========== Transformation Summary ==========\n');
fprintf('✅ Transformation complete: %d InSAR pixels converted to global frame\n', numel(V_local));
fprintf('Pixels used in fitting: %d\n', sum(valid_idx));
fprintf('Pixels with covariance correction: %d\n', sum(valid_idx));

% Calculate average uncertainty change from covariance term
avg_uncertainty_before = mean(sqrt(V_local_q(valid_idx).^2 + Var_corr(valid_idx)));
avg_uncertainty_after = mean(V_global_q(valid_idx));
avg_reduction = avg_uncertainty_before - avg_uncertainty_after;

fprintf('\n--- Uncertainty Statistics (fitted pixels) ---\n');
fprintf('Average σ_InSAR:                %.4f mm/yr\n', mean(V_local_q(valid_idx)));
fprintf('Average σ_model:                %.4f mm/yr\n', mean(sqrt(Var_corr(valid_idx))));
fprintf('Average σ_global (w/o cov):     %.4f mm/yr\n', avg_uncertainty_before);
fprintf('Average σ_global (with cov):    %.4f mm/yr\n', avg_uncertainty_after);
fprintf('Average uncertainty reduction:  %.4f mm/yr (%.1f%%)\n', ...
    avg_reduction, 100*avg_reduction/avg_uncertainty_before);

fprintf('\n--- Velocity Statistics ---\n');
fprintf('V_local range:   [%.2f, %.2f] mm/yr\n', min(V_local), max(V_local));
fprintf('V_global range:  [%.2f, %.2f] mm/yr\n', min(V_global), max(V_global));
fprintf('Correction range: [%.2f, %.2f] mm/yr\n', min(shift), max(shift));

fprintf('\n--- Sample Output (first 8 pixels) ---\n');
fprintf('Columns: [GNSS, InSAR_local, InSAR_global, σ_local, σ_global]\n');
disp(array2table([VUmed_int(1:8), V_local(1:8), V_global(1:8), V_local_q(1:8), V_global_q(1:8)], ...
    'VariableNames', {'GNSS', 'InSAR_local', 'InSAR_global', 'sigma_local', 'sigma_global'}));

fprintf('=============================================\n');

end