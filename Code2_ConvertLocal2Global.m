% ==============================================================
%      InSAR Local-to-Global VLM Conversion Workflow
% ==============================================================
%
%  Overview:
%  ----------
%  This script performs the complete transformation of local InSAR
%  line-of-sight (LOS) velocities into a globally referenced vertical
%  land-motion (VLM) frame using regional GNSS (MIDAS IGS20) data.
%
%  The workflow consists of four key stages:
%
%    1️⃣  **GNSS Oversampling (Interpolation)**
%        - The sparse GNSS vertical velocities (Vu, Su) are spatially
%          interpolated over the dense InSAR pixel grid using an
%          inverse-distance weighting (IDW) approach.
%        - This produces a continuous GNSS reference surface at every
%          InSAR pixel location.
%
%    2️⃣  **Local-to-Global Transformation**
%        - A polynomial correction (of degree `PN_Degree`) is fitted to
%          align the local InSAR VLM field with the oversampled GNSS
%          velocities.
%        - This step removes long-wavelength biases and converts the
%          InSAR VLM field into the IGS20 global reference frame.
%
%    4️⃣  **Result Saving**
%        - The following outputs are written per city:
%             a) Oversampled GNSS fields  (IDW surface)
%             b) Global-referenced InSAR VLM maps
%
%  Input Requirements:
%  -------------------
%    • Regional GNSS subsets from `fnc_GNSS_near_Region` (IGS20) meter/yr
%    • InSAR LOS velocity products {elpx_ll, MVel} (cm/yr)
%
%  Output Units:
%    • All final VLM values and uncertainties are stored in **mm/yr**.
%
%  Author: Mahmoud Reshadati
%  Date:   2025-10-30
% ==============================================================

clc; clear; % close all;

%% =============================================================
%  PART 1: PARAMETERS AND DATA PATHS
% =============================================================

CityName = {'LosAngele','NewYork'};   % Target cities (used for Variable names)

% --------- City definitions and GNSS input files ---------------
local2Global_MidasStations_Files = {  % GNSS (IGS20) input files per city
    'GlobalMaps/LosAngele_Local2Global_MidasStations_IGS20.mat'
    'GlobalMaps/NewYork_Local2Global_MidasStations_IGS20.mat'
    }; 

% --------- InSAR input data directories ------------------------
InSARPath = {                          
    'C:\Users\mahmo\Desktop\Paper-Local2Global with uncertainty\CodeReady2Publish\InSARLosAngeles'
    'C:\Users\mahmo\Desktop\Paper-Local2Global with uncertainty\CodeReady2Publish\InSARNewYork'
    };

% --------- Output directories ---------------------------------
OutputDir = {                          
    'LosAngeles_Outputs'
    'NewYork_Outputs'
    };

% --------- Processing parameters -------------------------------
IDW_power = 1;           % Power parameter for inverse-distance weighting
PN_Degree = 1;           % Polynomial degree for local-to-global correction
available_memory = 8;   % Available RAM for chunked interpolation (GB)

%% =============================================================
%  PART 2: LOCAL-TO-GLOBAL CONVERSION FOR EACH CITY
% =============================================================

for i = 1:numel(CityName)
    fprintf('\nProcessing City: %s\n\n', CityName{i});

    %% -------------------- Load GNSS data ---------------------
    % Load GNSS stations (regional subset, IGS20 reference)
    GlobalGNSSInfo = load(local2Global_MidasStations_Files{i});
    lat_GNSS = GlobalGNSSInfo.lat_GNSS_region;
    lon_GNSS = GlobalGNSSInfo.lon_GNSS_region;
    VLM_GNSS = GlobalGNSSInfo.Vu_GNSS_region;  
    VLMq_GNSS = GlobalGNSSInfo.Su_GNSS_region; 

    % Convert GNSS data from m/yr → mm/yr
    VLM_GNSS  = VLM_GNSS * 1000;
    VLMq_GNSS = VLMq_GNSS * 1000;
    fprintf('GNSS data from IGS20 loaded (units converted to mm/yr).\n');

    %% -------------------- Load InSAR data --------------------
    % InSAR data are stored in cm/yr (LOS direction)
    load([InSARPath{i} '/elpx_ll']);
    load([InSARPath{i} '/MVel']);
    load([InSARPath{i} '/MVelq']);
    load([InSARPath{i} '/uimages']);

    lat_InSAR = elpx_ll(:,2);
    lon_InSAR = elpx_ll(:,1);

    % Convert LOS → vertical using incidence angle
    VLM_InSAR  = MVel  ./ cosd(elpx_ll(:,3));
    VLMq_InSAR = MVelq ./ cosd(elpx_ll(:,3));

    % Convert cm/yr → mm/yr
    VLM_InSAR  = VLM_InSAR  * 10;
    VLMq_InSAR = VLMq_InSAR * 10;
    fprintf('InSAR data loaded and converted to vertical mm/yr.\n\n');

    %% =========================================================
    %  STEP 1: INVERSE-DISTANCE OVERSAMPLING OF GNSS DATA
    % =========================================================

    if ~exist(OutputDir{i}, "dir")
        mkdir(OutputDir{i});
    end

    savename_oversampled = sprintf('%s_GNSS_Oversampled_P%d.mat', CityName{i}, IDW_power);
    saveAddress_oversampledData = fullfile(OutputDir{i}, savename_oversampled);

    if exist(saveAddress_oversampledData, 'file')
        fprintf('Oversampled data already exists:\n -> %s\n', saveAddress_oversampledData);
        fprintf('Process skipped. Delete the file to rerun this step.\n---------\n');
        load(saveAddress_oversampledData);
        savingflag_oversample = 0;
    else
        % Interpolate GNSS VLM onto InSAR pixel locations using IDW weighting
        [VLM_GNSS_Oversampled, VLMq_GNSS_Oversampled] = ...
            InSAR_OverSampling_IDW_optimized(lat_GNSS, lon_GNSS, VLM_GNSS, VLMq_GNSS, ...
            lat_InSAR, lon_InSAR, IDW_power, available_memory);
        
        lat_oversampledGNSS = lat_InSAR;
        lon_oversampledGNSS = lon_InSAR;
        savingflag_oversample = 1;
    end

    %% =========================================================
    %  STEP 2: TRANSFORM LOCAL InSAR → GLOBAL VLM FRAME
    % =========================================================

    saveNameGlobal = sprintf('%s_GlobalInSAR_L%d_P%d.mat', CityName{i}, IDW_power, PN_Degree);
    saveAddress_GlobalData = fullfile(OutputDir{i}, saveNameGlobal);

    if exist(saveAddress_GlobalData, 'file')
        fprintf('Global VLM data already exists:\n -> %s\n', saveAddress_GlobalData);
        fprintf('Process skipped. Delete the file to rerun this step.\n---------\n');
        load(saveAddress_GlobalData);
        savingflag_local2global = 0;
    else
        % Apply local-to-global transformation (polynomial correction)
        [VLM_GlobalInSAR, VLMq_GlobalInSAR] = InSAR_local2global_WLS_optimized( ...
            lon_InSAR, lat_InSAR, VLM_InSAR, VLMq_InSAR, ...
            VLM_GNSS_Oversampled, VLMq_GNSS_Oversampled, PN_Degree);
            
        lat_GlobalInSAR = lat_InSAR;
        lon_GlobalInSAR = lon_InSAR;
        savingflag_local2global = 1;
    end

    %% =========================================================
    %  STEP 3: SAVE RESULTS
    % =========================================================

    % Save oversampled GNSS data
    if savingflag_oversample
        save(saveAddress_oversampledData, ...
            'VLM_GNSS_Oversampled','VLMq_GNSS_Oversampled', ...
            'lat_oversampledGNSS','lon_oversampledGNSS', ...
            'IDW_power');
    end

    % Save globally transformed InSAR data
    if savingflag_local2global
        save(saveAddress_GlobalData, ...
            'VLM_GlobalInSAR','VLMq_GlobalInSAR', ...
            'lon_GlobalInSAR','lat_GlobalInSAR', ...
            'IDW_power','PN_Degree', ...
            'VLM_InSAR','VLMq_InSAR');
    end

end