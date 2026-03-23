
% ==============================================================
%           DOWNLOAD AND PARSE MIDAS IGS20 VELOCITY FILE
% ==============================================================
%
%  Overview:
%  ----------
%  This script downloads, parses, and structures the global GNSS
%  velocity dataset from the MIDAS (Multivariate Integrated Data
%  Analysis System) IGS20 solution, maintained by the University of
%  Nevada, Reno (UNR).
%
%  The MIDAS dataset provides long-term GNSS station velocities in
%  the IGS20 reference frame, including vertical (Vu) and horizontal
%  (Ve, Vn) components and their corresponding uncertainties.
%
%  The script performs the following main steps:
%
%    1️⃣  **Download**
%        - Retrieves the most recent MIDAS IGS20 file from:
%          https://geodesy.unr.edu/gps_timeseries/IGS20/midas/midas.IGS.txt
%        - Skips the download if the file already exists locally.
%
%    2️⃣  **Parse**
%        - Reads the space-delimited text file according to the
%          structure described in:
%          https://geodesy.unr.edu/velocities/midas.readme.txt
%        - Extracts key columns: station name, position, time span,
%          velocities (Ve, Vn, Vu), and uncertainties (Se, Sn, Su).
%
%    3️⃣  **Convert and Save**
%        - Stores all fields as a MATLAB table `TableData` and saves
%          it to `MIDAS_IGS20_Velocities.mat` for later use in
%          regional GNSS filtering (`fnc_GNSS_near_Region.m`).
%
%    4️⃣  **Quick Visualization (optional)**
%        - Displays a global scatter map of vertical velocities Vu
%          (converted to mm/yr) for quick sanity checking.
%
%  Output:
%  --------
%    • MIDAS_IGS20_Velocities.mat containing:
%        - TableData: full parsed table
%        - station, lat, lon, height
%        - Ve, Vn, Vu, Se, Sn, Su (in meters/year)
%
%  Reference for Midas data:
%    Blewitt, G., et al. (2016), "GNSS Velocity Field for the United
%    States and Globally: The MIDAS Solution," *JGR Solid Earth*.
%
%  Author: Mahmoud Reshadati
%  Date:   2025-10-30
% ==============================================================


% ---- Step 0: Initialization ----
clear; clc;close all;

MidasFile_URL = 'https://geodesy.unr.edu/gps_timeseries/IGS20/midas/midas.IGS.txt';
OutputAddress = 'GlobalMaps';
OutputMidasTextFile = 'midas.IGS20.txt'; % text filename to be saved
OutputGNSSMatFile = 'MIDAS_IGS20_Velocities.mat'; % mat filename to be saved

fprintf('Starting GNSS Midas Downloader...\n\n');

% ---- Step 1: Download file ----

% Creating Destination Folder
if ~exist(OutputAddress,'dir')
    mkdir(OutputAddress);
end

% Checking Midas Text File Exist:
MidasTextFileAddress = fullfile(OutputAddress,OutputMidasTextFile);

if isfile(MidasTextFileAddress)
    fprintf(['Download Midas text File skipped.' ...
        '\nFile already exists locally:\n  ''%s'' \n'], ...
        MidasTextFileAddress);
    fprintf('To download a new version, delete %s\n\n',OutputMidasTextFile);
else
    fprintf('Downloading MIDAS IGS20 velocity file ...\n');
    try
        websave(MidasTextFileAddress, MidasFile_URL);
    catch ME
        error('Failed to download MIDAS file: %s', ME.message);
    end
    fprintf('Download complete, Saved File: %s\n\n', MidasTextFileAddress);
end

% ---- Step 2: Define variable names and format ----
% According to the README file (25–27 columns)
varNames = { ...
    'station', 'version', 't_start', 't_end', 'duration', ...
    'n_epochs_total', 'n_epochs_good', 'n_vel_samples', ...
    'Ve', 'Vn', 'Vu', ...
    'Se', 'Sn', 'Su', ...
    'E0', 'N0', 'U0', ...
    'frac_outlier_E', 'frac_outlier_N', 'frac_outlier_U', ...
    'std_pair_E', 'std_pair_N', 'std_pair_U', ...
    'n_steps', 'lat', 'lon', 'height' ...
    };

% The file uses space-separated columns
fmt = [ ...
    '%4s %s %f %f %f %f %f %f ' ...    % up to n_vel_samples
    '%f %f %f %f %f %f ' ...           % velocities + uncertainties
    '%f %f %f ' ...                    % offsets
    '%f %f %f ' ...                    % outlier fractions
    '%f %f %f ' ...                    % std velocity pairs
    '%f %f %f %f' ];                   % n_steps, lat, lon, height

% ---- Step 3: Read data ----
fprintf('Reading Midas text file ...\n');
fid = fopen(MidasTextFileAddress, 'r');
C = textscan(fid, fmt, 'CommentStyle', '#', 'Delimiter', ' ', 'MultipleDelimsAsOne', true);
fclose(fid);

% ---- Step 4: Convert to table ----
TableData = table(C{:}, 'VariableNames', varNames);

% ---- Step 5: Create arrays for easy access ----
station = string(TableData.station);
lat = TableData.lat;
lon = TableData.lon;
height = TableData.height;
Ve = TableData.Ve;  % Velocities in m/yr
Vn = TableData.Vn;
Vu = TableData.Vu;
Se = TableData.Se;  % uncertainties in m/yr
Sn = TableData.Sn;
Su = TableData.Su;

% ---- Step 6: Save locally ----
GNSS_matFile_address = fullfile(OutputAddress,OutputGNSSMatFile);
save(GNSS_matFile_address, ...
    'TableData','station','lat','lon','height', ...
    'Ve','Vn','Vu','Se','Sn','Su');

fprintf('Saved GNSS data to MIDAS_IGS20_Velocities.mat\n\n');
fprintf('Total stations loaded: %d\n', size(TableData,1));
fprintf('Unit of the GNSS data: meter/year \n');

% ---- (Optional) quick visualization ----
figure;
scatter(lon, lat, 10, Vu*1000, 'filled'); % Vu in mm/yr
colormap(jet); colorbar;
xlabel('Longitude'); ylabel('Latitude');
title('Global Vertical Land Motion (MIDAS IGS20, mm/yr)');


