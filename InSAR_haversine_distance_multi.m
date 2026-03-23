function D = InSAR_haversine_distance_multi(lat_ref, lon_ref, lat_gps, lon_gps)
% ==============================================================
%  InSAR_haversine_distance_multi
%  --------------------------------------------------------------
%  Purpose:
%    Computes great-circle distances (in meters) between multiple
%    reference points (e.g., regional centers or InSAR pixels)
%    and multiple GNSS station coordinates using the haversine formula.
%
%  Description:
%    This function returns a full [Nref × Ngps] distance matrix,
%    where each entry D(i,j) represents the great-circle distance
%    between reference point i and GNSS station j.
%
%  Equation (Haversine formula):
%    a = sin²(Δφ/2) + cos(φ₁)·cos(φ₂)·sin²(Δλ/2)
%    c = 2·atan2(√a, √(1−a))
%    d = R·c
%
%  Inputs:
%    lat_ref, lon_ref : [Nref × 1] Reference point coordinates (degrees)
%    lat_gps, lon_gps : [Ngps × 1] GNSS station coordinates (degrees)
%
%  Output:
%    D : [Nref × Ngps] matrix of great-circle distances (meters)
%
%  Notes:
%    • The computation assumes a spherical Earth (mean radius R = 6371 km).
%    • Fully vectorized implementation (no explicit loops).
%    • Used in GNSS spatial filtering and region selection.
%
%  Example:
%    % Compute distances from two region centers to three GNSS stations
%    lat_ref = [34.05; 40.7];
%    lon_ref = [-118.25; -74.0];
%    lat_gps = [33.9; 41.0; 39.8];
%    lon_gps = [-118.5; -73.9; -75.1];
%    D = InSAR_haversine_distance_multi(lat_ref, lon_ref, lat_gps, lon_gps);
%
%  Author: Mahmoud Reshadati
%  Date:   2025-10-30
% ==============================================================

    % ---- Constants ----
    R = 6371000; % mean Earth radius [m]

    % ---- Convert to radians ----
    lat_ref = deg2rad(lat_ref(:));
    lon_ref = deg2rad(lon_ref(:));
    lat_gps = deg2rad(lat_gps(:));
    lon_gps = deg2rad(lon_gps(:));

    % ---- Compute latitude/longitude differences ----
    % Using broadcasting (outer subtraction)
    dLat = lat_gps' - lat_ref; % [Nref × Ngps]
    dLon = lon_gps' - lon_ref; % [Nref × Ngps]

    % ---- Haversine formula ----
    a = sin(dLat/2).^2 + cos(lat_ref) .* cos(lat_gps') .* sin(dLon/2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));

    % ---- Distance in meters ----
    D = R * c;
end
