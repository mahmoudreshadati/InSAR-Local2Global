function fnc_GNSS_near_Region(MidasMatFile, outputFilename, outputAddress, ...
                    minYears_period, maxYears_period, StartYear, ...
                    lat_ref, lon_ref, max_distance_meter)
% ==============================================================
%  fnc_GNSS_near_Region
%  --------------------------------------------------------------
%  Purpose:
%    Filters GNSS stations from the global MIDAS IGS20 dataset
%    to extract only those that are:
%       - within a specified temporal range (duration and start year)
%       - within a given radial distance of a region of interest
%
%  This function is used to prepare regional GNSS subsets for
%  subsequent InSAR local-to-global transformation analyses.
%
%  Inputs:
%    MidasMatFile       - Path to the global MIDAS IGS20 velocity file (.mat)
%    outputFilename     - Name of the output regional GNSS file (.mat)
%    outputAddress      - Directory where the regional file will be saved
%    minYears_period    - Minimum GNSS data duration (years)
%    maxYears_period    - Maximum GNSS data duration (years)
%    StartYear          - Minimum acceptable installation/start year
%    lat_ref, lon_ref   - Reference coordinates of the region center
%    max_distance_meter - Maximum search radius (m) from region center
%
%  Outputs (saved in MAT file):
%    MidasTable_Region      - Filtered table of regional GNSS stations
%    station_GNSS_region    - Station names (string array)
%    lat_GNSS_region        - Latitudes of regional GNSS stations
%    lon_GNSS_region        - Longitudes of regional GNSS stations
%    height_GNSS_region     - Heights of GNSS stations (m)
%    Ve_GNSS_region, Vn_GNSS_region, Vu_GNSS_region  - East/North/Up velocities (m/yr)
%    Se_GNSS_region, Sn_GNSS_region, Su_GNSS_region  - Corresponding uncertainties (m/yr)
%
%  Example:
%    fnc_GNSS_near_Region('GlobalMaps/MIDAS_IGS20_Velocities.mat', ...
%        'NewYork_Local2Global_MidasStations_IGS20.mat', 'GlobalMaps', ...
%        10, 25, 2000, 40.7, -74.0, 5e5);
%
%  Author: Mahmoud Reshadati
%  Date:   2025-10-30
% ==============================================================

fprintf('Running fnc_GNSS_near_Region for output file: %s\n', outputFilename);
fprintf('Filtering GNSS stations near the area of interest...\n\n');

%% =============================================================
%  Step 1: Data Loading and Validation
% =============================================================

if ~isfile(MidasMatFile)
    error('Input file "%s" not found. Run download_midas_IGS20.m first.', MidasMatFile);
end

% The MIDAS dataset (TableData) contains GNSS velocities in m/yr.
load(MidasMatFile, 'TableData');
T = TableData;   % Create a local copy for convenience

lat_gps = T.lat;
lon_gps = T.lon;

%% =============================================================
%  Step 2: Apply Temporal Filters
% =============================================================
% Retain stations based on observation duration and start year

idx_duration = (T.duration >= minYears_period) & (T.duration <= maxYears_period);
idx_start    = (T.t_start >= StartYear);
idx_inrange_temporal = idx_duration & idx_start;

%% =============================================================
%  Step 3: Apply Spatial (Distance) Filter
% =============================================================
% Compute great-circle distance (in meters) from region center
% using haversine formula
D = InSAR_haversine_distance_multi(lat_ref, lon_ref, lat_gps, lon_gps);

% Keep all stations within max_distance_meter
idx_inrange_radius = any(D <= max_distance_meter, 1);

%% =============================================================
%  Step 4: Combine Filters and Extract Data
% =============================================================
Valid_idx = idx_inrange_temporal(:) & idx_inrange_radius(:);

MidasTable_Region = T(Valid_idx, :);

% Extract relevant columns into named variables
station_GNSS_region = string(MidasTable_Region.station);
lat_GNSS_region     = MidasTable_Region.lat;
lon_GNSS_region     = MidasTable_Region.lon;
height_GNSS_region  = MidasTable_Region.height;

Ve_GNSS_region = MidasTable_Region.Ve;
Vn_GNSS_region = MidasTable_Region.Vn;
Vu_GNSS_region = MidasTable_Region.Vu;

Se_GNSS_region = MidasTable_Region.Se;
Sn_GNSS_region = MidasTable_Region.Sn;
Su_GNSS_region = MidasTable_Region.Su;

% adding parameters into a struct for saving
RegionalGNSSParams.minYears_period = minYears_period;
RegionalGNSSParams.maxYears_period = maxYears_period;
RegionalGNSSParams.StartYear = StartYear;
RegionalGNSSParams.max_distance_meter = max_distance_meter;
RegionalGNSSParams.lat_ref = lat_ref;
RegionalGNSSParams.lon_ref = lon_ref;
%% =============================================================
%  Step 5: Save Filtered Dataset
% =============================================================
OutputFile = fullfile(outputAddress, outputFilename);

save(OutputFile, ...
    'MidasTable_Region', 'station_GNSS_region', ...
    'lat_GNSS_region', 'lon_GNSS_region', 'height_GNSS_region', ...
    'Ve_GNSS_region', 'Vn_GNSS_region', 'Vu_GNSS_region', ...
    'Se_GNSS_region', 'Sn_GNSS_region', 'Su_GNSS_region',...
    'RegionalGNSSParams');


%% =============================================================
%  Step 6: Reporting Summary
% =============================================================
fprintf('✅ Filtered dataset saved successfully:\n -> %s\n', OutputFile);
fprintf('Total stations in MIDAS file: %d\n', height(T));
fprintf('Stations after region/time filtering: %d\n', height(MidasTable_Region));
fprintf('All velocity components are in meters/year.\n');
fprintf('----------------------------------------------------------\n\n');

end
