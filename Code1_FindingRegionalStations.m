% ==============================================================
%        Extract Regional GNSS Stations for Global Conversion
% ==============================================================
%  This script identifies GNSS (MIDAS IGS20) stations located near
%  each study region (city) and prepares their data for subsequent
%  InSAR local-to-global velocity transformation.
%
%  Author: Mahmoud Reshadati
%  Date:   2025-10-30
% ==============================================================

% --------- Initialization ---------
clc; clear; close all;

%% =============================================================
%  PART 1: PARAMETER DEFINITIONS
% =============================================================

% --------- Global GNSS MIDAS input file -----------------------
MidasMatFile = 'GlobalMaps/MIDAS_IGS20_Velocities.mat';   % Contains GNSS velocities in IGS20 frame

% --------- Target regions ------------------------------------
CityName = {'NewYork', 'LosAngele'};      % Target cities (name of files to generate)
MidasOutputPath = 'GlobalMaps';            % Directory for saving regional GNSS subsets

% --------- GNSS station selection criteria -------------------
minYears_period   = [10, 10];   % Minimum operation duration (years) for GNSS inclusion
maxYears_period   = [25, 25];   % Maximum allowable data span (years)
StartYear         = [2000, 2000]; % Exclude GNSS stations installed after this year

% --------- Reference coordinates for each region -------------
% Approximate central coordinates of the area (used for distance filtering)
Region_Reflat     = [40.7 , 34.05];   % Reference latitude for each region
Region_Reflon     = [-74.0 , -118.24]; % Reference longitude for each region

% --------- Distance threshold (search radius) ----------------
% Maximum radial distance (m) around region center for GNSS inclusion
MaxDistanceRef_meter = [0.5e+6, 0.5e+6];  % 500 km radius per city

% --------- Visualization flag --------------------------------
% Set to 1 to visualize GNSS station distribution per city
plot_flag = 1;

%% =============================================================
%  PART 2: REGIONAL GNSS STATION EXTRACTION
% =============================================================

for i = 1:numel(CityName)

    fprintf('\nExtracting GNSS stations for city: %s\n', CityName{i});

    % --------- Output file definition -------------------------
    outputFilename = [CityName{i},'_Local2Global_MidasStations_IGS20.mat'];

    % --------- Run GNSS extraction function -------------------
    % This function filters GNSS stations from the global MIDAS file
    % based on:
    %   - Data longevity (min/max years)
    %   - Start year constraint
    %   - Geographical distance from region center
    %
    % The output includes:
    %   lat_GNSS_region, lon_GNSS_region, Vu_GNSS_region, Su_GNSS_region
    %   saved in a MAT file for later use.

    fnc_GNSS_near_Region(MidasMatFile, outputFilename, MidasOutputPath, ...
        minYears_period(i), maxYears_period(i), StartYear(i), ...
        Region_Reflat(i), Region_Reflon(i), MaxDistanceRef_meter(i));

    % --------- Optional visualization --------------------------
    % Plot the spatial distribution of GNSS stations selected for the region
    if plot_flag ~= 0
        DataRegionalGNSSFile = fullfile(MidasOutputPath, outputFilename);
        DataRegionalGNSS = load(DataRegionalGNSSFile);

        figure('Name', sprintf('GNSS Map - %s', CityName{i}), 'Color', 'w');
        geoscatter(DataRegionalGNSS.lat_GNSS_region, DataRegionalGNSS.lon_GNSS_region, ...
            25, DataRegionalGNSS.Vu_GNSS_region * 1000, 'filled');
        hold on; geoplot(Region_Reflat(i),Region_Reflon(i),...
            'rv','markersize',15,'MarkerFaceColor','w');
        geobasemap satellite;
        colormap(jet); colorbar;
        clim([-5, 5]);
        title(sprintf('GNSS stations for %s (velocity in mm/yr)', CityName{i}));
    end
end