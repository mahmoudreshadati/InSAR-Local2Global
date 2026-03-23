function [VLM_GPS_interp, VLM_GPS_interp_std] = InSAR_OverSampling_IDW_optimized( ...
    LATI, LONI, VLM_GPS, VLM_GPS_std, lat_insar, lon_insar, power, AssignedRAM)
% ==============================================================
%  InSAR_OverSampling_IDW_optimized
%  --------------------------------------------------------------
%  Purpose:
%    Performs spatial interpolation of GNSS vertical land motion (VLM)
%    data onto InSAR pixel coordinates using inverse-distance weighting
%    (IDW) combined with GNSS uncertainty weighting.
%
%  Description:
%    • Each GNSS station contributes to nearby InSAR pixels with a weight
%      proportional to (1 / distance^power) × (1 / σ²), where σ is the
%      GNSS velocity uncertainty.
%    • The method operates in local tangent-plane coordinates for accurate
%      geometric distances (km scale).
%    • To handle large InSAR datasets (millions of pixels), the process
%      runs in memory-managed chunks based on the available RAM.
%
%  Mathematical Formulation:
%    For each pixel j:
%        wᵢⱼ = ( (1 / dᵢⱼ^p) / σᵢ² )
%        V_interpⱼ = Σ( wᵢⱼ × Vᵢ ) / Σ( wᵢⱼ )
%        σ_interpⱼ = √[ Σ( (wᵢⱼ / Σwᵢⱼ)² × σᵢ² ) ]
%
%  Inputs:
%    LATI, LONI        - GNSS station coordinates (vectors, degrees)
%    VLM_GPS           - GNSS vertical velocities (mm/yr)
%    VLM_GPS_std       - GNSS velocity uncertainties (mm/yr)
%    lat_insar, lon_insar - InSAR pixel coordinates (vectors, degrees)
%    power (optional)  - IDW power parameter (default = 2)
%    AssignedRAM       - RAM limit for chunking (GB, default = 8)
%
%  Outputs:
%    VLM_GPS_interp     - Interpolated GNSS VLM at InSAR pixels (mm/yr)
%    VLM_GPS_interp_std - Propagated uncertainties (mm/yr)
%
%  Notes:
%    • Distances are computed in kilometers using latlon2local (ENU).
%    • Handles NaN or zero uncertainties by assigning conservative values.
%    • Chunk processing avoids out-of-memory errors even for >10⁶ pixels.
%    • Fully vectorized inside each chunk; no pixel-wise loops.
%
%  Example:
%    [V_interp, V_std] = InSAR_Claude_WeightedVLMInterp_Chunk( ...
%        lat_gnss, lon_gnss, Vu_gnss, Su_gnss, lat_insar, lon_insar, 2, 32);
%
%  Author: Mahmoud Reshadati
%  Date:   2025-10-30
% ==============================================================

fprintf('----------------\nOversampling GNSS stations at InSAR pixel locations:\n\n')

%% Input arguments check
if nargin < 7
    power = 2; % distance metric power
end
if nargin < 8
    AssignedRAM = 8; % 8GB RAM default
end

%% Memory estimation and chunk size
% Detailed explanation inside comments — estimates how many pixels can
% be processed per chunk given AssignedRAM and matrix count assumptions.
N_gps = numel(LATI);
N_insar = numel(lat_insar);
bytes_per_double = 8;
num_matrices = 5;
mem_overhead_factor = 1.3;

memNeededPerMatrix = bytes_per_double * N_gps * N_insar;
memNeeded = memNeededPerMatrix * num_matrices * mem_overhead_factor;
memNeededGB = memNeeded / 1e9;
fprintf('Estimated peak RAM for full interpolation: %.2f GB\n', memNeededGB);

memPerChunkPerMatrix = bytes_per_double * N_gps;
chunk_size = floor((AssignedRAM * 1e9) / (memPerChunkPerMatrix * num_matrices * mem_overhead_factor));
if chunk_size < 1
    error('AssignedRAM too low for even one chunk. Increase AssignedRAM or reduce data size.');
end
Nchunks = ceil(N_insar / chunk_size);
fprintf('\nStarting interpolation in %d chunks (~%.0f pixels each)\n', Nchunks, chunk_size);
memChunkPerMatrixGB = (N_gps * chunk_size * bytes_per_double) / 1e9;
memChunkGB = memChunkPerMatrixGB * num_matrices * mem_overhead_factor;
fprintf('Each chunk will use ~%.2f GB RAM (including overhead)\n', memChunkGB);

%% Coordinate preparation and input sanitization
LATI = LATI(:);
LONI = LONI(:);
VLM_GPS = VLM_GPS(:);
VLM_GPS_std = VLM_GPS_std(:);
lat_insar = lat_insar(:);
lon_insar = lon_insar(:);

valid_std = VLM_GPS_std > 0 & isfinite(VLM_GPS_std);
if any(~valid_std)
    max_valid_std = max(VLM_GPS_std(valid_std));
    if isempty(max_valid_std) || ~isfinite(max_valid_std)
        max_valid_std = 1.0;
    end
    VLM_GPS_std(~valid_std) = max_valid_std * 10;
    fprintf('Warning: %d GNSS stations had invalid std; replaced with %.2f mm/yr\n', ...
            sum(~valid_std), max_valid_std * 10);
end

w_gps = 1 ./ (VLM_GPS_std.^2);

%% Convert to local Cartesian coordinates (ENU)
lat0 = mean([LATI; lat_insar], 'omitnan');
lon0 = mean([LONI; lon_insar], 'omitnan');
alt0 = 0;
[x_gps, y_gps, ~] = latlon2local(LATI, LONI, zeros(size(LATI)), [lat0, lon0, alt0]);
[x_insar, y_insar, ~] = latlon2local(lat_insar, lon_insar, zeros(size(lat_insar)), [lat0, lon0, alt0]);
x_gps = x_gps / 1000; y_gps = y_gps / 1000;
x_insar = x_insar' / 1000; y_insar = y_insar' / 1000;

%% Interpolation (chunked)
VLM_GPS_interp = NaN(N_insar, 1);
VLM_GPS_interp_std = NaN(N_insar, 1);

for c = 1:Nchunks
    t0 = tic;
    idx_start = (c-1)*chunk_size + 1;
    idx_end = min(c*chunk_size, N_insar);
    idx_chunk = idx_start:idx_end;

    x_chunk = x_insar(idx_chunk);
    y_chunk = y_insar(idx_chunk);

    d = sqrt((x_gps - x_chunk).^2 + (y_gps - y_chunk).^2);
    d(d == 0) = 1e-6;

    W = (w_gps ./ (d.^power));
    Wsum = sum(W, 1);
    valid = Wsum > 0 & isfinite(Wsum);

    if any(~valid)
        fprintf('  Warning: %d pixels in chunk %d have invalid weights\n', sum(~valid), c);
    end

    W_normalized = W ./ Wsum;
    W_normalized(:, ~valid) = NaN;

    VLM_GPS_interp(idx_chunk) = (VLM_GPS' * W_normalized)';

    W_normalized_sq = W_normalized.^2;
    VLM_GPS_std_sq = VLM_GPS_std.^2;
    variance_interp = sum(W_normalized_sq .* VLM_GPS_std_sq, 1);
    VLM_GPS_interp_std(idx_chunk) = sqrt(variance_interp)';

    fprintf('  Chunk %d/%d done in %.1f s\n', c, Nchunks, toc(t0));
end

fprintf('✅ Interpolation complete: %d GNSS → %d InSAR pixels (distances in km)\n', ...
        N_gps, N_insar);
fprintf('Valid interpolations: %d (%.1f%%)\n', sum(isfinite(VLM_GPS_interp)), ...
        100*sum(isfinite(VLM_GPS_interp))/N_insar);
fprintf('------------------------------------------------------\n');

end
