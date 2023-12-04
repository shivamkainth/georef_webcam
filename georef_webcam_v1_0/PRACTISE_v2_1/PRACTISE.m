% Haerer, Bernhardt and Schulz (2016)
% "PRACTISE - Photo Rectification And ClassificaTIon SoftwarE (V.2.1)"
%
    % incorporating adapted versions of
    %
    % Aronica, Bates and Horritt, 2002
    % "Assessing the uncertainty in distributed model predictions using
    % observed binary pattern information within GLUE" (Hydrol Process)
    %
    % Corripio, 2004
    % "Snow surface albedo estimation using terrestrial photography"
    % (Int J Rem Sens)
    %
    % Haerer, Bernhardt, Schulz and Corripio (2012)
    % "PRACTISE - Photo Rectification And ClassificaTIon SoftwarE (V.1.0)"
    % (GMD)
    %
    % Salvatori et al., 2011
    % "Snow cover monitoring with images from digital camera systems"
    % (ItJRS)
    %
    % Tolson and Shoemaker, 2007
    % "Dynamically dimensioned search algorithm for computationally
    % efficient watershed model calibration" (WRR)
    %
    % Wang, Robinson and White, 2000
    % "Generating Viewsheds without Using Sightliines" (PE&RS)
%
%   written by
%   Stefan Haerer (LMU Munich)
%   08/2012
%   contact: stefan.haerer@boku.ac.at
%       updated by Stefan Haerer (BOKU Vienna, 22/12/2015)
%       automatically edited by georef_webcam by Sebastian Buchelt (University of Würzburg, 15/06/2020)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Name:       PRACTISE
%   Purpose:    Main file of PRACTISE
%   Comment:    This file starts and controls the complete software tool,
%               be careful with any changes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
tic
%
disp(' ')
disp('Starting PRACTISE - Photo Rectification And ClassificaTIon SoftwarE')
% load input configuration
Input_PRACTISE_m='Input_PRACTISE';
run(Input_PRACTISE_m)
%%%%%%%%%%%%%%%%%%%%%%%%%% check & administration %%%%%%%%%%%%%%%%%%%%%%%%%
% check switches
if ~exist('vs', 'var')
    error('No viewshed switch specified, check parameter file')
end
if ~exist('os', 'var')
    error(['No GCP/DDS optimisation switch specified, check ', ...
     'parameter file'])
end
if ~exist('cs', 'var')
    error('No classification switch specified, check parameter file')
end
if ~exist('rs', 'var')
    error('No remote sensing switch specified, check parameter file')
end
if rs>0
    if ~exist('rsps', 'var')
        error(['No remote sensing probability switch specified, ', ...
         'check parameter file'])
    end
    if ~exist('rsms', 'var')
        error('No remote sensing mask switch specified, check ', ...
         'parameter file')
    end
end
if ~exist('is', 'var')
    error('No image switch specified, check parameter file')
end
% Camera positions (viewpoint and targetpoint) defined?
if ~exist('cam', 'var')
    error(['Camera view- and targetpoint: longitude and latitude ', ...
     'positions are not defined.'])
elseif size(cam, 2)<2
    error('Camera targetpoint: longitude and latitude position is ', ...
     'not defined.')
end
% check if Unix and change "\" to "/" in "fin_..." and "fout_..." in case
list=who('-regexp', 'fin_');
for i=1:length(list)
    assignin('base',char(list(i)),Unix_Slash_PRACTISE(evalin('base',char(list(i)))));
end
clear list
if is==0
    fout_folder=Unix_Slash_PRACTISE(fout_folder);
end
% define and create output folder
if is==1
    dt_stamp=strrep(strrep(strrep(datestr(now), '-', '_'), ':', '_'), ...
     ' ', '_');
    fout_folder=[fin_folder, 'output_', strrep(fin_imfolder, filesep, ...
     '_'), dt_stamp, filesep];
end
% Output folder defined and exists?
if ~exist(fout_folder, 'dir')
    mkdir(fout_folder);
else
    fout_folder_msg1=['The selected output folder does already ', ...
         'exist and existing files might be overwritten. Please ', ...
         'check before confirming to continue.'];
    if IsOctave_PRACTISE()
        fout_folder_ok1=questdlg(fout_folder_msg1, 'Hint', ...
         'Continue', 'Stop', 'Continue');
        if fout_folder_ok1(1:4)=='Stop'
            return
        end
        clear fout_folder_ok1 fout_folder_msg1
    else
%         uisetpref('clearall'); % to reset dialog box properties
        fout_folder_ok1=uigetpref('userhelp', 'ctrlc', 'Hint', ...
        fout_folder_msg1, 'Ok');
        clear fout_folder_ok1
    end
end
% create logfile and save input file (*.m)
ls=1;
if ls==1
    disp(['The used input file ''', Input_PRACTISE_m, '.m'' will be ', ...
     'saved and a logfile will be produced in ''', fout_folder, ''''])
    if ~exist('dt_stamp', 'var')
        dt_stamp=strrep(strrep(strrep(datestr(now), '-', '_'), ...
         ':', '_'), ' ', '_');
    end
    diary([fout_folder, 'logfile_PRACTISE_', dt_stamp, '.txt']);
    copyfile([Input_PRACTISE_m, '.m'], fout_folder);
end
if is==1 || ls==1
    clear dt_stamp
end
clear Input_PRACTISE_m
%
disp('PRACTISE will be run using the following options:')
% display and check chosen switches
%   viewshed
if vs == 0
    disp('- use existing viewshed (Arc/Info ASCII Grid)')
elseif vs == 1
    disp('- generate viewshed')
else
    error('Viewshed switch option not available, check parameter file')
end
%   ground control points (GCPs) & optimisation
if os == 0
    disp('- w/o GCPs & w/o DDS optimisation')
elseif os == 1
    disp('- w GCPs & w/o DDS optimisation')
elseif os == 2
    disp('- w GCPs & w DDS optimisation')
elseif os == 3
    disp('- w GCPs & w DDS optimisation in an interactive mode')
else
    error(['GCP/DDS optimisation switch option not available, check ', ...
     'parameter file'])
end
%   image
if is == 0
    disp('- classify a single image')
elseif is == 1
    if rs > 0 || os > 1
        disp(['- classify a single image with automatically derived ', ...
         'file name'])
    elseif rs == 0 && os <= 1
        disp(['- classify all ', fin_imformat(2:end), '-images in ''', ...
         fin_folder, fin_imfolder, ''''])
    end
    if cs == 0
        disp('using the manual classification mode')
    elseif cs==1
        disp('using the automatic blue band classification')
    elseif cs==2
        disp('using the automatic blue band + pca classification')
    elseif cs==3
        disp(['using an interactive classification mode (default: ', ...
         'automatic blue band classification)'])
    else
        error(['Classification switch option not available, check ', ...
         'parameter file'])
    end
else
    error('Image switch option not available, check parameter file')
end
%   satellite image
if rs == 0
    disp('- w/o satellite image classification')
elseif rs == 1
    disp(['- w satellite image classification (of unzipped raw ', ...
     'Landsat data) using a NDSI threshold '])
elseif rs == 2
    disp(['- w satellite image classification (of preprocessed NDSI ', ...
     'data) using a NDSI threshold '])
else
    error(['Remote sensing switch option not available, check ', ...
     'parameter file'])
end
if rs > 0
    if rsps == 1
        if cs == 2
            disp(['calibrated by probability values of the snow/no ', ...
             'snow photo map'])
        elseif cs < 2
            disp(['calibrated by binary values of the snow/no snow ', ...
             'photo map'])
            disp(['(comment: remote sensing probability switch on ', ...
             'only works, when the pca classification switch is or ', ...
             'might be (=3) selected, remote sensing probability ', ...
             'switch will be deactivated now)'])
            rsps=0;
        elseif cs == 3
            disp(['calibrated by binary or probability values in ', ...
             'snow/no snow photo map dependent on the subsequent ', ...
             'selection in the interactive classification'])
        end
    elseif rsps == 0
        disp('calibrated by binary values of the snow/no snow photo map')
    else
        error(['Remote sensing probability option not available, ', ...
         'check parameter file'])
    end
    if rsms==0
        disp('and w/o using a mask (for excluding e.g. clouds)')
    elseif rsms == 1
        disp('and using a binary mask (for excluding e.g. clouds)')
    elseif rsms == 2
        disp(['and using an existing fmask map (excluding e.g. ', ...
         'cloud, cloud shadow and/or water pixel)'])
    else
        error(['Remote sensing mask option not available, check ', ...
         'parameter file'])
    end
else
    if rsms > 0 || rsps > 0
        disp(['The remote sensing switch is deactivated, hence remote ', ...
         'sensing mask and and probability option will be deactivated'])
        rsms=0;
        rsps=0;
    end
end
%   check additional switch combinations
if vs == 0 && os > 1
    error(['An activated GCP/DDS optimisation does not work in ', ...
     'combination with the use of an existing viewshed, check ', ...
     'parameter file'])
end
% Priority of graphics (problematic)
if IsOctave_PRACTISE()
    disp(['The current ''graphics_toolkit'' is ''', graphics_toolkit, '''.'])
    disp(['Hint: Displaying figures in Octave is not always stable. Change ', ...
     'between the ''available_graphics_toolkits'' if problems occur.'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% load %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Loading and preparing input')
% load
%   DEM
%       read LINK file (if existing)
if strcmp(fin_demW(end-8:end), '.dem.LINK')
    fid=fopen(fin_demW);
    if fid==-1
        error('Linking file to DEM could not be found.')
    end
    fin_demW=fgetl(fid);
    fin_demW=[fin_demW, fgetl(fid)];
    fin_demW=Unix_Slash_PRACTISE(fin_demW);
    fclose(fid);
end
%       read DEM file
disp(['Load DEM ''', fin_demW, ''''])
fid=fopen(fin_demW);
if fid==-1
    error('DEM file could not be found.')
end
for i=1:6 % 6 headerlines
    dummy=fgetl(fid);
    dummypos=regexp(dummy, ' ');
    header_W(i,:)=cellstr(dummy);
    headern_W(i,:)=cellstr(dummy(1:dummypos(1)));
    headerv_W(i,:)=str2num(dummy(dummypos(1):end));
end
clear dummy dummypos i
demWrcz = fscanf(fid, '%f', [headerv_W(1,1), inf])';
demWrcz(demWrcz==headerv_W(6,1))=NaN;
fclose(fid);
clear ans fid
%   Photograph (RGB)
%       single photo
if exist('fin_imfolder', 'var')
    fin_imagepath=[fin_folder, fin_imfolder];
    if exist('fin_image', 'var')
        N_images=1;
    else
        if  is==0
            error('Photo file name is not defined.')
%       multiple photos
        elseif is==1
            if exist('fin_imformat', 'var')
                list=dir([fin_imagepath, '*', fin_imformat]);
                if length(list)<1
                    list=dir([fin_imagepath, '*', '.img.LINK']);
                    if length(list)<1
                        error('No photo or link to a photo could be found.')
                    elseif length(list)>1
                        error('Several linking files to photos or photo folders could be found.')
                    end
                    fin_imagelink=list(1).name;
%           read LINK file (if existing)
                    fid=fopen([fin_imagepath, fin_imagelink]);
                    if fid==-1
                        error('Linking file to image could not be found.')
                    end
                    fin_imagepath=fgetl(fid);
                    fin_imagepath=Unix_Slash_PRACTISE(fin_imagepath);
                    if ~feof(fid)
                        list(1).name=fgetl(fid);
                    else
                        list=dir([fin_imagepath, '*', fin_imformat]);
                    end
                    fclose(fid);
                end
                if length(list)>1 && (rs>0 || os>1)
                    error(['Optimising the camera location and orientation ', ...
                     'using GCPs or classifying satellite images while ', ...
                     'processing more than one image does not work.'])
                end
                for i=1:length(list)
                    fin_images{i,1}=list(i).name;
                end
                clear i list
                N_images=length(fin_images);
%           define the first photograph file name (might be the only one)
                fin_image=char(fin_images(1));
                fin_image=Unix_Slash_PRACTISE(fin_image);
            else
                error('Photo format could not be found.')
            end
        end
    end
%      load first photograph (might be the only one)
%           read photo
    disp(['Load photo from folder ''', fin_imagepath, ''''])
    photo = imread([fin_imagepath, fin_image]);
%               photo size (pixel rows and columns)
    [pix_r pix_c dummy]=size(photo);
    clear dummy
    if is==1
%               automatically read root name of file
        fin_imageDot=regexp(fin_image, '\.');
        f_name=fin_image(1:fin_imageDot(end)-1);
        clear fin_imageDot
        if rs>0
%               automatically read date and time from photo file name
%                (for check of time shift to raw Landsat data)
            if ~strcmpi(f_name(end-2:end), 'utc')
                if exist('fin_imagelink', 'var')
                    fin_imagelink=fin_imagelink(1:end-9);
                    fin_imageDot=regexp(fin_imagelink, '\.');
                    fin_imagelink=fin_imagelink(1:fin_imageDot(end)-1);
                    time_pos=fin_imagelink(end-15:end);
                    if ~strcmpi(fin_imagelink(end-2:end), 'utc')
                        error(['Neither the end of the photo filename ', ...
                         '(time zone) nor the linking photo filename ', ...
                         'are ''utc'' or ''UTC'''])
                    end
                    clear fin_imageDot fin_imagelink
                else
                    error(['The end of the photo filename (time ', ...
                     'zone) must be ''utc'' or ''UTC'''])
                end
            else
                time_pos=f_name(end-15:end);
            end
            time=struct('year', str2num(time_pos(1:4)), ...
                        'month', str2num(time_pos(5:6)), ...
                        'day', str2num(time_pos(7:8)), ...
                        'hour', str2num(time_pos(10:11)), ...
                        'min', str2num(time_pos(12:13)), ...
                        'sec', 0, ...
                        'UTC', 0);
            time_names = fieldnames(time);
            for i = 1:length(time_names)
                if isempty(time.(time_names{i,1}))
                    error(['Date and Time format of photo file ', ...
                     'name is not correct, check if file name is ', ...
                     'in the format ''*yyyymmdd_HHMMutc', ...
                     fin_imformat, ''''])
                end
            end
            clear i time_pos time_names
            disp(['Date and time of the photograph in UTC, check ', ...
             'if correct: '])
            disp(time)
        end
        clear fin_imformat
    end
else
    error('Photo folder could not be found.')
end
%   Viewshed
if vs==0
    if is==1
        list=dir([fin_folder, fin_imfolder, '*', fin_vsformat]);
        if length(list)<1
            list=dir([fin_folder, fin_imfolder, '*', '.view.LINK']);
        end
        if length(list)==1
            fin_viewW=list.name;
            clear list
        else
            error(['No or multiple viewsheds or links to viewsheds ', ...
             'could be found.'])
        end
        clear fin_vsformat
    end
%       read LINK file (if existing)
    if strcmp(fin_viewW(end-9:end), '.view.LINK')
        fid=fopen([fin_folder, fin_imfolder, fin_viewW]);
        if fid==-1
            error('Linking file to viewshed could not be found.')
        end
        fin_viewWpath=fgetl(fid);
        fin_viewWpath=Unix_Slash_PRACTISE(fin_viewWpath);
        fin_viewW=fgetl(fid);
        fin_viewW=Unix_Slash_PRACTISE(fin_viewW);
        fclose(fid);
    else
        fin_viewWpath=[fin_folder, fin_imfolder];
    end
%       read viewshed
    disp(['Load viewshed ''', fin_viewWpath, fin_viewW, ''''])
    fid=fopen([fin_viewWpath, fin_viewW]);
    if fid==-1
        error('Viewshed file could not be found.')
    end
    for i=1:6
        dummy=fgetl(fid);
        if dummy==char(header_W(i,:))
            clear dummy
        else
            error('Viewshed and DEM header not equal, please check.')
        end
    end
    clear i
    viewW = fscanf(fid, '%f', [headerv_W(1,1), inf])';
    viewW(viewW==headerv_W(6,1))=NaN;
    fclose(fid);
    clear ans fid
end
%   Ground control points
if os>0
    if is==1
        list=dir([fin_folder, fin_imfolder, '*', fin_gcpformat]);
        if length(list)<1
            list=dir([fin_folder, fin_imfolder, '*', '.gcp.LINK']);
        end
        if length(list)==1
            fin_gcpW=list.name;
            clear list
        else
            error(['No or multiple GCP files or links to GCP files ', ...
             'could be found.'])
        end
        clear fin_gcpformat
    end
%       read LINK file (if existing)
    if strcmp(fin_gcpW(end-8:end), '.gcp.LINK')
        fid=fopen([fin_folder, fin_imfolder, fin_gcpW]);
        if fid==-1
            error('Linking file to GCP file could not be found.')
        end
        fin_gcpWpath=fgetl(fid);
        fin_gcpWpath=Unix_Slash_PRACTISE(fin_gcpWpath);
        fin_gcpW=fgetl(fid);
        fin_gcpW=Unix_Slash_PRACTISE(fin_gcpW);
        fclose(fid);
    else
        fin_gcpWpath=[fin_folder, fin_imfolder];
    end
%       read GCP file
    disp(['Load GCP file ''', fin_gcpWpath, fin_gcpW, ''''])
    fid=fopen([fin_gcpWpath, fin_gcpW]);
    if fid==-1
        error('GCP file could not be found.')
    end
    cell_gcpW=textscan(fid, '%f %f %f %u %u %s', 'headerlines', 1);
    for i=1:5 % longitude, latitude, altitude, photo rows and photo cols
        gcpW(i,:)=cell_gcpW{1,i}';
    end
    gcpW_name(1,:)=cell_gcpW{1,6}';
    fclose(fid);
    clear ans fid i cell_gcpW
end
%   Satellite image
if rs>0
    fin_satpath=[fin_folder, fin_satfolder];
    if ~exist('[fin_satpath, fin_satfolder_image]', 'dir')
        list=dir([fin_satpath, '*', '.satfolder.LINK']);
        if length(list)~=1
            error(['No satellite image folder and no or multiple links to the destination of the satellite image folder could be found.'])
        end
        %       read LINK file (if existing)
        if strcmp(list.name(end-14:end), '.satfolder.LINK')
            fid=fopen([fin_satpath, list.name]);
            if fid==-1
                error('Linking file to satellite folder could not be found.')
            end
            fin_satpath=fgetl(fid);
            fin_satpath=Unix_Slash_PRACTISE(fin_satpath);
            fclose(fid);
        end
    end
    clear fin_satfolder list
end
if rs==1
%     raw unzipped Landsat & Landsat Look data
    if exist('fin_satpath', 'var') && exist('fin_satfolder_image', 'var')
%       define the satellite image input file names
%           Landsat image metadata (TXT-file)
        clear list
        list=dir([fin_satpath, fin_satfolder_image, '*', '_MTL.txt']);
        if length(list)==1
            fin_satname_meta=list(1).name;
            clear list
        else
            error('No or multiple MTL (metadata) file names of the Landsat image are found.')
        end
%           Landsat Look image (JPG-file) if existing
        if exist('fin_satfolder_look', 'var')
            list=dir([fin_satpath, fin_satfolder_look, ...
             fin_satname_meta(1:end-8), '.jpg']);
            if length(list)==1
                fin_satname_look=list(1).name;
            end
            clear list
        end
        if ~exist('fin_satname_look', 'var')
            disp(['No Landsat Look image file name could be found, ', ...
                'NDSI map will be used as background of snow map figure.'])
            clear list fin_satfolder_look
        end
        if strcmp(fin_satname_meta(1:3), 'LT5') || strcmp(fin_satname_meta(1:3), 'LE7')
%           additional Landsat 7  file names (TXT-file)
%               Earth-Sun distance
            fin_satname_DistEtoS='earth_sun_distance.txt';
            fin_satname_DistEtoS=Unix_Slash_PRACTISE(fin_satname_DistEtoS);
%               Solar spectral irradiances
            if strcmp(fin_satname_meta(1:3), 'LT5')
                fin_satname_RadSol='lt5_tm_chkur_solar_spectral_radiances.txt';
            elseif strcmp(fin_satname_meta(1:3), 'LE7')
                fin_satname_RadSol='le7_etm_plus_chkur_solar_spectral_radiances.txt';
            end
            fin_satname_RadSol=Unix_Slash_PRACTISE(fin_satname_RadSol);
        end
%       load Landsat files
%           load Landsat Look image
        if exist('fin_satname_look', 'var')
            disp(['Load Landsat Look image ''', ...
             fin_satpath, fin_satfolder_look, fin_satname_look, ''''])
            satLook=imread([fin_satpath, fin_satfolder_look, ...
             fin_satname_look]);
            if ~exist('satLook', 'var')
                error('Landsat Look image file could not be found.')
            end
        end
        if strcmp(fin_satname_meta(1:3), 'LT5') || strcmp(fin_satname_meta(1:3), 'LE7')
%           load Earth-Sun Distance (according to Landsat 7 handbook)
            disp(['Load Earth-Sun distance from ''', fin_satname_DistEtoS , ''''])
            fid=fopen(fin_satname_DistEtoS);
            if fid==-1
                error('Earth-Sun distance file of the Landsat 7 could not be found.')
            end
            for i=1:2 % 2 headerlines
                fgetl(fid);
            end
            satDistEtoS=fscanf(fid, '%u\t%f', [2, inf])';
            fclose(fid);
            clear ans fid i
%           load ETM+ Solar Spectral Irradiances (according to Landsat 7 handbook)
            disp(['Load Landsat solar spectral irradiances from ''', ...
             fin_satname_RadSol, ''''])
            fid=fopen(fin_satname_RadSol);
            if fid==-1
                error('Solar spectral irradiances file for the Landsat data could not be found.')
            end
            for i=1:2 % 2 headerlines
                fgetl(fid);
            end
            satRadSol=fscanf(fid, '%u\t%f', [2, inf])';
            fclose(fid);
            clear ans fid i
        end
%           load Landsat metadata
        disp(['Load Landsat metadata from ''', fin_satpath, ...
            fin_satfolder_image, fin_satname_meta, ''''])
        fid=fopen([fin_satpath, fin_satfolder_image, fin_satname_meta]);
        if fid==-1
            error('MTL (metadata) file of the Landsat image could not be found.')
        end
        while ~feof(fid)
            dummystr=fgetl(fid);
            dummynr=regexp(dummystr, '=');
            if ~isempty(dummynr)
                dummynames=strtrim(dummystr(1:dummynr-1));
                satMTL.(dummynames)=strtrim(dummystr(dummynr+1:end));
                if ~isempty(regexp(satMTL.(dummynames), '"'))
                    satMTL.(dummynames)=regexprep(satMTL.(dummynames), '"', '');
                end
            end
        end
        fclose(fid);
        satMETA=struct('datestring', satMTL.DATE_ACQUIRED, ...
         'timestring', satMTL.SCENE_CENTER_TIME, ...
         'x_ul', satMTL.CORNER_UL_PROJECTION_X_PRODUCT, ...
         'y_ul', satMTL.CORNER_UL_PROJECTION_Y_PRODUCT, ...
         'sunelev', satMTL.SUN_ELEVATION, ...
         'nrows', satMTL.REFLECTIVE_LINES, ...
         'ncols', satMTL.REFLECTIVE_SAMPLES, ...
         'cellsize', satMTL.GRID_CELL_SIZE_REFLECTIVE, ...
         'nodata', '0');
        if strcmp(fin_satname_meta(1:3), 'LT5') || strcmp(fin_satname_meta(1:3), 'LE7')
            satMETA.bGfnamestring=satMTL.FILE_NAME_BAND_2;
            satMETA.bGlmax=satMTL.RADIANCE_MAXIMUM_BAND_2;
            satMETA.bGlmin=satMTL.RADIANCE_MINIMUM_BAND_2;
            satMETA.bGqcalmax=satMTL.QUANTIZE_CAL_MAX_BAND_2;
            satMETA.bGqcalmin=satMTL.QUANTIZE_CAL_MIN_BAND_2;
            satMETA.bMIRfnamestring=satMTL.FILE_NAME_BAND_5;
            satMETA.bMIRlmax=satMTL.RADIANCE_MAXIMUM_BAND_5;
            satMETA.bMIRlmin=satMTL.RADIANCE_MINIMUM_BAND_5;
            satMETA.bMIRqcalmax=satMTL.QUANTIZE_CAL_MAX_BAND_5;
            satMETA.bMIRqcalmin=satMTL.QUANTIZE_CAL_MIN_BAND_5;
            satMETA.bNIRfnamestring=satMTL.FILE_NAME_BAND_4;
            satMETA.bNIRlmax=satMTL.RADIANCE_MAXIMUM_BAND_4;
            satMETA.bNIRlmin=satMTL.RADIANCE_MINIMUM_BAND_4;
            satMETA.bNIRqcalmax=satMTL.QUANTIZE_CAL_MAX_BAND_4;
            satMETA.bNIRqcalmin=satMTL.QUANTIZE_CAL_MIN_BAND_4;
        elseif strcmp(fin_satname_meta(1:3), 'LC8')
            satMETA.bGfnamestring=satMTL.FILE_NAME_BAND_3;
            satMETA.bGreflmult=satMTL.REFLECTANCE_MULT_BAND_3;
            satMETA.bGrefladd=satMTL.REFLECTANCE_ADD_BAND_3;
            satMETA.bMIRfnamestring=satMTL.FILE_NAME_BAND_6;
            satMETA.bMIRreflmult=satMTL.REFLECTANCE_MULT_BAND_6;
            satMETA.bMIRrefladd=satMTL.REFLECTANCE_ADD_BAND_6;
            satMETA.bNIRfnamestring=satMTL.FILE_NAME_BAND_5;
            satMETA.bNIRreflmult=satMTL.REFLECTANCE_MULT_BAND_5;
            satMETA.bNIRrefladd=satMTL.REFLECTANCE_ADD_BAND_5;
        end
        dummynames=char(fieldnames(satMETA));
        for i=1:size(dummynames,1)
            if isempty(regexp(dummynames(i,:),'string'))
                satMETA.(strtrim(dummynames(i,:)))= ...
                 str2num(satMETA.(strtrim(dummynames(i,:))));
            end
        end
        clear ans fid i dummynames
        if IsOctave_PRACTISE()
            satMETA_t_str_dummy=satMETA.timestring ...
             (1:regexp(satMETA.timestring(1:end-1),'\.')-1);
            satJDAY=datenum([satMETA.datestring, ' ', ...
             satMETA_t_str_dummy], 'yyyy-mm-dd HH:MM:SS');
            clear satMETA_t_str_dummy
        else
            satJDAY=datenum([satMETA.datestring, ' ', ...
             satMETA.timestring(1:end-1)], 'yyyy-mm-dd HH:MM:SS');
        end
%           load Landsat bands (TIF-files)
%               green band
        satG=imread([fin_satpath, fin_satfolder_image, ...
         satMETA.bGfnamestring]);
%               mid infrared band
        satMIR=imread([fin_satpath, fin_satfolder_image, ...
         satMETA.bMIRfnamestring]);
%               near infrared band
        satNIR=imread([fin_satpath, fin_satfolder_image, ...
            satMETA.bNIRfnamestring]);
%               if Mapping Toolbox is available, e.g.
%         [satG, satG_R, satG_bbox]=geotiffread([fin_folder, ...
%             fin_satfolder, fin_satfolder_image, satMETA.bGfnamestring]);
%         geotiffinfo([fin_satpath, fin_satfolder_image, ...
%             satMETA.bGfnamestring])
%         figure()
%         mapshow(satG, satG_R)
%         close all;
    else
        error('Landsat folder name and/or subfolder names do not exist.')
    end
elseif rs==2
%     NDSI map (Arc/Info ASCII Grid)
%         NDSI map metadata information
    satMETA=struct('x_ul', NaN, 'y_ul', NaN, ...
     'nrows', NaN, 'ncols', NaN, 'cellsize', NaN, 'nodata', NaN, ...
     'datetimestring', satdatetime);
    clear satdatetime
    dummynames=fieldnames(satMETA);
    dummynr=[4, 3, 1, 2, 5, 6];
%         load NDSI map metadata information
    disp(['Load NDSI map ''', fin_satpath, fin_satfolder_image, ...
     fin_satname_ndsi, ''''])
    fid=fopen([fin_satpath, fin_satfolder_image, fin_satname_ndsi]);
    if fid==-1
        error('NDSI map could not be found.')
    end
    for i=1:6 % 6 headerlines
        dummy=fgetl(fid);
        dummypos=regexp(dummy, ' ');
        header_sat(i,:)=cellstr(dummy);
        satMETA.(dummynames{dummynr(i)})=str2num(dummy(dummypos(1):end));
    end
%             shift coordinates from corner to pixel centre
    satMETA.x_ul=satMETA.x_ul+(0.5*satMETA.cellsize);
    satMETA.y_ul=satMETA.y_ul-(0.5*satMETA.cellsize)+(satMETA.nrows*satMETA.cellsize);
    clear i dummy dummypos dummynames dummynr i
%         load NDSI map data
    satNDSI = fscanf(fid, '%f', [satMETA.ncols, inf])';
    satNDSI(satNDSI==satMETA.nodata)=NaN;
    i=find(satNDSI<-1 | satNDSI>1);
    satNDSI(i)=2;
    fclose(fid);
    clear ans fid
%         recording time of satellite image
    satJDAY=datenum(satMETA.datetimestring, 'yyyy-mm-dd HH:MM:SS');
end
%   remote sensing mask
if rsms>0
    if rsms==1
%       using a binary mask (e.g. for clouds, Arc/Info ASCII Grid)
%           read mask
        disp(['Load satellite image mask ''', fin_satpath, ...
         fin_satfolder_mask, fin_satname_mask, ''''])
        fid=fopen([fin_satpath, fin_satfolder_mask, fin_satname_mask]);
        if fid==-1
            error('Remote sensing binary mask file could not be found.')
        end
        for i=1:6
            dummy=fgetl(fid);
            dummypos=regexp(dummy, ' ');
            header_satmask(i,:)=cellstr(dummy);
            headerv_satmask(i,:)=str2num(dummy(dummypos(1):end));
        end
        clear i dummy dummypos
        satmask=fscanf(fid, '%f', [headerv_satmask(1,1), inf])';
        fclose(fid);
        clear ans fid
    elseif rsms==2
%       excluding optionally clouds, cloud shadows and/or water (fmask algorithm)
%           fmask file name
        list=dir([fin_satpath, fin_satfolder_mask, ...
         '*', 'Fmask']);
        if length(list)==1
            fin_satname_mask=list(1).name;
        end
        clear list
%           read fmask hdr file
        fid=fopen([fin_satpath, fin_satfolder_mask, fin_satname_mask, ...
         '.hdr']);
        if fid==-1
            error('Remote sensing binary mask file could not be found.')
        end
        for i=1:12
            dummy=fgetl(fid);
            header_satmask(i,:)=cellstr(dummy);
            if i==3 || i==4
                dummypos=regexp(dummy, '=');
                headerv_satmask(i-2,1)=str2num(dummy(dummypos(1)+1:end));
            end
        end
        dummypos=regexp(dummy, ',');
        for i=3:5
            headerv_satmask(i,1)=str2num(dummy(dummypos(i)+1:dummypos(i+1)-1));
        end
        headerv_satmask(4,1)=headerv_satmask(4,1)-(headerv_satmask(2,1)*headerv_satmask(5,1));
        headerv_satmask(6,1)=255;
        fclose(fid);
        clear i ans fid dummy dummypos
%           read fmask file
        disp(['Load fmask map of satellite image ''', fin_satpath, ...
         fin_satfolder_mask, fin_satname_mask, ''''])
        fid=fopen([fin_satpath, fin_satfolder_mask, fin_satname_mask]);
        if fid==-1
            error('Remote sensing fmask file could not be found.')
        end
        satmask=fread(fid, headerv_satmask(1)*headerv_satmask(2), 'uint8', 0, 'ieee-le');
        satmask=reshape(satmask, [headerv_satmask(1), headerv_satmask(2)])';
        fclose(fid);
        clear ans fid
    end
%       check if the satellite data and its mask data have the same
%       cell size and registration point
    if satMETA.cellsize==headerv_satmask(5,1)
        if mod(satMETA.x_ul, satMETA.cellsize)~= ...
        mod(headerv_satmask(3,1)+0.5*headerv_satmask(5,1), headerv_satmask(5,1))
            error(['The longitude coordinate of the registration points of satellite ', ...
             'image and mask data are not identical'])
        end
        if mod(satMETA.y_ul, satMETA.cellsize)~= ...
        mod(headerv_satmask(4,1)+0.5*headerv_satmask(5,1), headerv_satmask(5,1))
            error(['The latitude coordinate of the registration points of satellite ', ...
             'image and mask data are not identical'])
        end
    else
        error('The cell size of the satellite image and the mask data are not identical')
    end
    if rsms==1
%       pixels that are zero or assigned no data are not masked
        satmask_code=unique(satmask);
        satmask_code=satmask_code(satmask_code~=headerv_satmask(6,1));
        satmask(satmask==headerv_satmask(6,1))=0;
    elseif rsms==2
%       check that the set mask code only includes existing and
%       reasonable fmask code values
        if length(satmask_code)>4 || length(satmask_code)~=length([find(satmask_code==1), ...
        find(satmask_code==2), find(satmask_code==4), find(satmask_code==255)])
            error(['The set mask code is not reasonable when using ', ...
                'fmask data, check parameter file'])
        else
            for i=1:length(satmask_code)
                satmask(satmask==satmask_code(i))=100;
            end
            satmask(satmask~=100)=0;
        end
    end
%       convert from double to logical
%           0 = clear
%           1 = not clear or no observation
    satmask(satmask~=0)=1;
    satmask=logical(satmask);
end
%%%%%%%%%%%%%%%%%%%%%%%%% GCPs & DDS optimisation %%%%%%%%%%%%%%%%%%%%%%%%%
if os==0
% w/o GCPs and w/o DDS optimisation
%   create vector for coordinate transformation
    X0=[cam(1,1), cam(2,1), cam_off(1), cam_rol, cam(1,2), cam(2,2), ...
     cam_off(2), cam_foc, cam_hei, cam_wid];
else
% w GCPs and w or w/o DDS optimisation
%   first assumption: interactive optimisation (if os==2 will be changed to non-interactive later)
    os_intact='y';
    os_startval='y';
    os_quest=1;
    gcpRMSE_count=1;
    gcpP_0_dummy=0;
%   interactive loop
    while os_intact=='y'
        if os_startval=='y'
%       create vector for coordinate transformation or seed in DDS optimisation
            X0=[cam(1,1), cam(2,1), cam_off(1), cam_rol, cam(1,2), ...
             cam(2,2), cam_off(2), cam_foc, cam_hei, cam_wid];
        end
%       w GCPs and w or w/o DDS optimisation
%           original RMSE value and projected coordinates of GCPs
        [gcpRMSE_orig, gcpP]=DDS_RMSE_PRACTISE(X0, demWrcz, ...
         headerv_W, gcpW, pix_c, pix_r);
        if os>1
%       w DDS optimisation
            toc
            disp('Starting DDS optimisation')
%           upper and lower boundaries (UB and LB) using UBD and LBD
%               check if a sign is missing in lower boundaries
            if length(LBD)~=sum(LBD<=0)
                i=find(LBD>0);
                LBD(i)=-LBD(i);
                clear i
            end
%               set UB and LB using UBD and LBD and vaerify that overall
%                boundaries are not increased
            if ~exist('UB_gcp', 'var')
                UB_gcp=X0+UBD;
                LB_gcp=X0+LBD;
                UB_gcp_orig=UB_gcp;
                LB_gcp_orig=LB_gcp;
            else
                i=find(X0+UBD<UB_gcp_orig);
                UB_gcp(i)=X0(i)+UBD(i);
                i=find(X0+UBD>UB_gcp_orig);
                UB_gcp(i)=UB_gcp_orig(i);
                i=find(X0+LBD>LB_gcp_orig);
                LB_gcp(i)=X0(i)+LBD(i);
                i=find(X0+LBD<LB_gcp_orig);
                LB_gcp(i)=LB_gcp_orig(i);
                clear i
            end
%           call DDS optimisation
            [Xopt, Fopt]=DDS_PRACTISE('DDS_RMSE_PRACTISE', X0, LB_gcp, ...
             UB_gcp, DDS_R_os, DDS_MaxEval_os, demWrcz, headerv_W, ...
             gcpW, pix_c, pix_r);
%           save original camera parameters (incl. altitude)
            [cam_orig, cam_rc_orig]=Cam_PRACTISE(cam, headerv_W, ...
             demWrcz, cam_off);
%           Store RMSE value before optimisation
            if gcpP_0_dummy==0
                gcpRMSE_0=gcpRMSE_orig;
                gcpP_0=gcpP;
                X_0=X0;
%               save original camera parameters (incl. altitude)
                [cam_0, cam_rc_0]=Cam_PRACTISE(cam, headerv_W, ...
                 demWrcz, cam_off);
            end
            if ~exist('gcpRMSE_opt', 'var') || Fopt<gcpRMSE_orig
                gcpP_orig=gcpP;
            end
%           optimal RMSE value and projected coordinates of GCPs
            [gcpRMSE_opt, gcpP]=DDS_RMSE_PRACTISE(Xopt, demWrcz, headerv_W, ...
                gcpW, pix_c, pix_r);
            clear Fopt
        end
%       plot photograph and hold on for visual investigation of the optimisation
        close(figure(1));
        if IsOctave_PRACTISE()
            figure1=figure(1,'position',get(0,'screensize')+[100 50 -150 -150]);
        else
            figure1=figure(1);
            set(figure1, 'units', 'normalized', 'outerposition', [0 0 1 1], ...
              'PaperUnits', 'centimeters', 'PaperPosition', [0 0 50 35], ...
              'PaperSize', [50 35], 'PaperType', '<custom>');
        end
        if os==1
            set(figure1, 'NumberTitle', 'off', 'Name', 'GCPs');
        elseif os>1
            set(figure1, 'NumberTitle', 'off', 'Name', 'GCPs and DDS optimisation', 'visible', 'off');
        end
        set(figure1, 'Color', [1 1 1], 'InvertHardcopy', 'off');
        axes1 = axes('FontSize', 15.91, 'TickDir', 'out', ...
          'XAxisLocation', 'top', 'YDir', 'reverse', ...
          'DataAspectRatio', [1 1 1], 'Parent', figure1);
		if ~IsOctave_PRACTISE()
			set(axes1, 'FontWeight', 'demi')
		end
        axis(axes1, [0.5 pix_c+0.5 0.5 pix_r+0.5]);
        xlabel(axes1, 'Pixel column', 'FontSize', 18.185);
        ylabel(axes1, 'Pixel row', 'FontSize', 18.185);
        box(axes1, 'on');
        hold(axes1, 'all');
        %
        image(photo, 'Parent', axes1);
        hold on
        plot(gcpW(4,:), gcpW(5,:), 'gx', 'LineWidth', 2, 'MarkerSize', os_MarkSiz+2)
        if os>1
            if gcpP_0_dummy==0
                fig1h1=plot(gcpP_0(1,:), gcpP_0(2,:), 'o', 'MarkerEdgeColor', [1 0 0], 'MarkerSize', os_MarkSiz);
                if IsOctave_PRACTISE()
                    set (fig1h1, 'MarkerSize', os_MarkSiz+8)
                end
            elseif gcpP_0_dummy==1
                fig1h2=plot(gcpP_0(1,:), gcpP_0(2,:), 'o', 'MarkerEdgeColor', [1 0 0], 'MarkerSize', os_MarkSiz+4);
                if IsOctave_PRACTISE()
                    set (fig1h2, 'MarkerSize', os_MarkSiz+8)
                else
                    plot(gcpP_orig(1,:), gcpP_orig(2,:), 'o', 'MarkerEdgeColor', [1 0 0], 'MarkerSize', os_MarkSiz);
                end
            end
        end
        fig1h=plot(gcpP(1,:), gcpP(2,:), '.', 'MarkerEdgeColor', [1 0 0], 'MarkerSize', os_MarkSiz);
        if os==1
            Leg=legend('GCP positions', 'Calculated GCP positions');
        elseif os>1
            if gcpP_0_dummy==1 && ~IsOctave_PRACTISE()
                Leg=legend('GCP positions', ...
                 'Calculated GCP positions before DDS optimisation', ...
                 'Calculated GCP positions before last successful DDS optimisation', ...
                 'Calculated GCP positions after DDS optimisation');
            else
                Leg=legend('GCP positions', ...
                 'Calculated GCP positions before DDS optimisation', ...
                 'Calculated GCP positions after DDS optimisation');
            end
        end
        if IsOctave_PRACTISE()
            set(Leg, 'location', 'southoutside', 'color', [0.75 0.75 0.75])
            fig1_leg=get(figure1,'children');
            fig1_leg_pos=get(fig1_leg(1), 'position');
            set(fig1_leg(1), 'position', fig1_leg_pos+[0, 0, 0.03, 0]);
        end
        %
        clear axes1 Leg fig1_leg fig1_leg_pos
%       interactive questions
        if os==3
            if os_quest==1
                if gcpRMSE_opt<gcpRMSE_optthres
                    disp(['The optimised RMSE distance between real and calculated GCPs (', ...
                        num2str(gcpRMSE_opt), ') is below the set threshold of ', num2str(gcpRMSE_optthres)])
                    if ~IsOctave_PRACTISE()
%                         uisetpref('clearall'); % to reset the dialog box properties
                        oscs_ok1=uigetpref('userhelp', 'ctrlc', 'Hint', ...
                            ['Switching between different Matlab windows ', ...
                            'while a message box is open is not possible. ', ...
                            'Use the keyboard combination ''Ctrl + C'' when the ', ...
                            'message box is active to allow this option.'], ...
                            'Ok');
                    end
                    set(figure1, 'visible', 'on');
                    button=questdlg('Do you want to run the optimisation again?', ...
                        'GCPs and DDS optimisation', 'Yes using start values', ...
                        'Yes using new optimised values', 'No', 'Yes using start values');
                elseif gcpRMSE_count>=gcpRMSE_countthres
                    disp(['The maximum number of DDS optimisation runs (', num2str(gcpRMSE_countthres), ...
                        ') has been reached without any optimised RMSE distance between real and calculated', ...
                        ' GCPs below the set threshold of ', num2str(gcpRMSE_optthres)])
                    gcpRMSE_opttemp(gcpRMSE_count)=gcpRMSE_opt;
                    gcpRMSE_Xopttemp(gcpRMSE_count,:)=Xopt;
                    gcpRMSE_gcpPtemp(gcpRMSE_count,:)=[gcpP(1,:), gcpP(2,:)];
                    gcpRMSE_opt=min(gcpRMSE_opttemp);
                    disp(['The camera parameters of the DDS optimisation run with the lowest RMSE distance', ...
                        ' between real and calculated GCPs (', num2str(gcpRMSE_opt), ') will be used'])
                    gcpRMSE_optrow=find(gcpRMSE_opttemp==gcpRMSE_opt);
                    Xopt=gcpRMSE_Xopttemp(gcpRMSE_optrow(1),:);
                    gcpP(1,:)=gcpRMSE_gcpPtemp(gcpRMSE_optrow(1),1:end/2);
                    gcpP(2,:)=gcpRMSE_gcpPtemp(gcpRMSE_optrow(1),end/2+1:end);
                    if ~IsOctave_PRACTISE()
                        oscs_ok1=uigetpref('userhelp', 'ctrlc', 'Hint', ...
                            ['Switching between different Matlab windows ', ...
                            'while a message box is open is not possible. ', ...
                            'Use the keyboard combination ''Ctrl + C'' when the ', ...
                            'message box is active to allow this option.'], ...
                            'Ok');
                    end
                    delete(fig1h);
                    plot(gcpP(1,:), gcpP(2,:), '.', 'MarkerEdgeColor', [1 0 0], 'MarkerSize', os_MarkSiz)
                    set(figure1, 'visible', 'on');
                    button=questdlg('Do you want to run the optimisation again?', ...
                        'GCPs and DDS optimisation', 'Yes using start values', ...
                        'Yes using new optimised values', 'No', 'Yes using start values');
                else
                    button='Yes using start values';
                end
            elseif os_quest==2
                set(figure1, 'visible', 'on');
                button=questdlg('Do you want to optimise the values further?', ...
                    'GCPs and DDS optimisation', 'Yes using old optimised values', ...
                    'Yes using new optimised values', 'No', 'Yes using old optimised values');
            end
            switch button
                case ''
                    disp('PRACTISE has been terminated by user.')
                    diary off
                    return
                case 'Yes using start values'
                    os_startval='n';
                    if gcpRMSE_count>=gcpRMSE_countthres
                        gcpRMSE_count=1;
                    else
                        gcpRMSE_opttemp(gcpRMSE_count)=gcpRMSE_opt;
                        gcpRMSE_Xopttemp(gcpRMSE_count,:)=Xopt;
                        gcpRMSE_gcpPtemp(gcpRMSE_count,:)=[gcpP(1,:), gcpP(2,:)];;
                        gcpRMSE_count=gcpRMSE_count+1;
                    end
                case 'Yes using new optimised values'
                    if gcpRMSE_opt<gcpRMSE_orig
                        gcpRMSE_dummy=gcpRMSE_orig;
                        cam_dummy=cam_orig;
                        cam_rc_dummy=cam_rc_orig;
                        X0_dummy=X0;
                    end
                    os_startval='y';
                    os_quest=2;
                    if IsOctave_PRACTISE() && gcpP_0_dummy==0
                        os_msg2=['If the multiplication factors for the upper and/or ', ...
                         'lower boundaries are adapted, any changes in specific boundary ', ...
                         'values (+/-delta...) are ignored.'];
                        os_ok2=helpdlg(os_msg2, 'Hint');
                        clear os_ok2 os_msg2
                    elseif ~IsOctave_PRACTISE()
                        os_ok2=uigetpref('userhelp', 'multfact', 'Hint',  ...
                         ['If the multiplication factors for the upper and/or ', ...
                         'lower boundaries are adapted, any changes in ', ...
                         'specific boundary values (+/-delta...) are ignored.'], ...
                         'Ok');
                    end
                    gcpP_0_dummy=1;
%           forward new parameters to existing variables
                    cam(1,1)=Xopt(1,1);
                    cam(2,1)=Xopt(1,2);
                    cam_off(1)=Xopt(1,3);
                    cam_rol=Xopt(1,4);
                    cam(1,2)=Xopt(1,5);
                    cam(2,2)=Xopt(1,6);
                    cam_off(2)=Xopt(1,7);
                    cam_foc=Xopt(1,8);
                    cam_hei=Xopt(1,9);
                    cam_wid=Xopt(1,10);
	                os_enter={'Multiplication factor for upper boundary deviations'; ...
                        'Multiplication factor for lower boundary deviations'; ...
		   	            '+/-delta cam location (longitude) in m:'; ...
                        '+/-delta cam location (latitude) in m:'; ...
                        '+/-delta cam location offset in m:'; ...
                        '+/-delta cam roll in degree:'; ...
                        '+/-delta cam target (longitude) in m:'; ...
                        '+/-delta cam target (latitude) in m:'; ...
                        '+/-delta cam target offset in m:'; ...
                        '+/-delta cam focal length in m:'; ...
                        '+/-delta cam sensor height in m:'; ...
                        '+/-delta cam sensor width in m:'};
                    os_def={num2str(1); num2str(1); ...
                        [num2str(UBD(1)), ' / ', num2str(LBD(1))]; ...
                        [num2str(UBD(2)), ' / ', num2str(LBD(2))]; ...
                        [num2str(UBD(3)), ' / ', num2str(LBD(3))]; ...
                        [num2str(UBD(4)), ' / ', num2str(LBD(4))]; ...
                        [num2str(UBD(5)), ' / ', num2str(LBD(5))]; ...
                        [num2str(UBD(6)), ' / ', num2str(LBD(6))]; ...
                        [num2str(UBD(7)), ' / ', num2str(LBD(7))]; ...
                        [num2str(UBD(8)), ' / ', num2str(LBD(8))]; ...
                        [num2str(UBD(9)), ' / ', num2str(LBD(9))]; ...
                        [num2str(UBD(10)), ' / ', num2str(LBD(10))]};
                    disp('Optimisation values and upper and lower boundaries (in metres unless otherwise noted)')
                    disp(['Values of cam ...     locat long       ', ...
                        'locat lat    locat off    roll (in deg)       ', ...
                        'targ long        targ lat     targ off    ', ...
                        'focal length    sensor height    sensor width'])
                    disp([sprintf('%s\t%12.1f\t%12.1f\t%9.2f\t%12.2f\t%12.1f\t%12.1f\t%9.2f\t%13.5f\t%14.5f\t%14.5f', ...
                        ['before optimisation'], X0)])
                    disp([sprintf('%s\t%12.1f\t%12.1f\t%9.2f\t%12.2f\t%12.1f\t%12.1f\t%9.2f\t%13.5f\t%14.5f\t%14.5f', ...
                        ['after optimisation'], Xopt)])
                    disp([sprintf('%s\t\t%12.1f\t%12.1f\t%9.2f\t%12.2f\t%12.1f\t%12.1f\t%9.2f\t%13.5f\t%14.5f\t%14.5f', ...
                        ['upper boundary'], UB_gcp)])
                    disp([sprintf('%s\t\t%12.1f\t%12.1f\t%9.2f\t%12.2f\t%12.1f\t%12.1f\t%9.2f\t%13.5f\t%14.5f\t%14.5f', ...
                        ['lower boundary'], LB_gcp)])
                    if IsOctave_PRACTISE()
                        os_answer=inputdlg(os_enter, 'Upper and lower boundary deviations', ...
                            1, os_def);
                    else
                        os_answer=inputdlg(os_enter, 'Upper and lower boundary deviations', ...
                            1, os_def, 'on');
                    end
                    if isempty(os_answer)
                        disp('PRACTISE has been terminated by user.')
                        diary off
                        return
                    else
                        if str2num(os_answer{1})~=1 || str2num(os_answer{2})~=1
                            UBD=UBD*str2num(os_answer{1});
                            LBD=LBD*str2num(os_answer{2});
                        else
                            for i=1:length(UBD)
                                pos_slash=strfind(os_answer{i+2}, '/');
                                if isempty(pos_slash)
                                    error(['The input dialog has not been changed correctly. ', ...
                                        'Please use the following format for the boundary ', ...
                                        'deviations e.g. ''100 / 100'''])
                                end
                                UBD(1,i)=str2num(os_answer{i+2}(1:pos_slash-1));
                                LBD(1,i)=str2num(os_answer{i+2}(pos_slash+1:end));
                                clear pos_slash
                            end
                        end
                    end
                    clear os_answer os_def os_enter
                case 'No'
                    if exist('gcpRMSE_dummy', 'var') && ...
                     gcpRMSE_opt==gcpRMSE_orig
                        gcpRMSE_orig=gcpRMSE_dummy;
                        cam_orig=cam_dummy;
                        cam_rc_orig=cam_rc_dummy;
                        X0=X0_dummy;
                    end
                    os_intact='n';
%           forward new parameters to existing variables
                    cam(1,1)=Xopt(1,1);
                    cam(2,1)=Xopt(1,2);
                    cam_off(1)=Xopt(1,3);
                    cam_rol=Xopt(1,4);
                    cam(1,2)=Xopt(1,5);
                    cam(2,2)=Xopt(1,6);
                    cam_off(2)=Xopt(1,7);
                    cam_foc=Xopt(1,8);
                    cam_hei=Xopt(1,9);
                    cam_wid=Xopt(1,10);
                case 'Yes using old optimised values'
                    os_startval='n';
            end
        else
            os_intact='n';
        end
        clear figure1 fig1h fig1h1 fig1h2
    end
    if os>1 && os_intact=='n'
        disp(['Optimisation values and upper and lower boundaries ', ...
            '(in metres unless otherwise noted)'])
        disp(['Values of cam ...     locat long       ', ...
            'locat lat    locat off    roll (in deg)       ', ...
            'targ long        targ lat     targ off    ', ...
            'focal length    sensor height    sensor width'])
        disp([sprintf('%s\t%12.1f\t%12.1f\t%9.2f\t%12.2f\t%12.1f\t%12.1f\t%9.2f\t%13.5f\t%14.5f\t%14.5f', ...
            ['before optimisation'], X0)])
        disp([sprintf('%s\t%12.1f\t%12.1f\t%9.2f\t%12.2f\t%12.1f\t%12.1f\t%9.2f\t%13.5f\t%14.5f\t%14.5f', ...
            ['after optimisation'], Xopt)])
        disp([sprintf('%s\t\t%12.1f\t%12.1f\t%9.2f\t%12.2f\t%12.1f\t%12.1f\t%9.2f\t%13.5f\t%14.5f\t%14.5f', ...
            ['upper boundary'], UB_gcp)])
        disp([sprintf('%s\t\t%12.1f\t%12.1f\t%9.2f\t%12.2f\t%12.1f\t%12.1f\t%9.2f\t%13.5f\t%14.5f\t%14.5f', ...
            ['lower boundary'], LB_gcp)])
        clear gcpP_0_dummy
    end
    if exist('gcpRMSE_opt', 'var')
        disp(['The RMSE between real and calculated GCPs is ', ...
         num2str(gcpRMSE_opt),'px after the optimisation.'])
    else
        disp(['The RMSE between real and calculated GCPs is ', ...
         num2str(gcpRMSE_orig),'px.'])
    end
end
clear UBD LBD os_intact os_startval os_quest oscs_ok1 os_ok2 proxy_n ...
 button gcpRMSE_count gcpRMSE_dummy gcpRMSE_Xopttemp gcpRMSE_opttemp ...
 gcpRMSE_gcpPtemp gcpRMSE_optrow X0_dummy cam_dummy cam_rc_dummy os_MarkSiz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% viewshed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate camera parameters (raster position and altitude)
[cam, cam_rc]=Cam_PRACTISE(cam, headerv_W, demWrcz, cam_off);
% Specified camera angles viewshed
if vs==1
    toc
%   normalize DEM to viewpoint altitude
    demWvs=(demWrcz/headerv_W(5,1))-cam_rc(3,1);
%       create radial buffer zone around camera location
%           j=col i=row cam_rc(1,1)=col cam_rc(2,1)=row
    [j i]=meshgrid(1:size(demWvs, 2), 1:size(demWvs, 1));
    buffer_radius_rc=buffer_radius/headerv_W(5,1);
    buffer_cam=sqrt((j-cam_rc(1,1)).^2+(i-cam_rc(2,1)).^2)<=buffer_radius_rc;
    demWvs(buffer_cam==1)=min(min(demWvs))-1; % smaller than minimum existing value
    clear i j buffer_radius_rc
    demWvs(cam_rc(2,1),cam_rc(1,1))=0;
%   Cam properties & viewshed raster initialisation
    disp('Start calculation of specified camera angles (+ 0.1rad) viewshed')
%       field of view (use maximum range, if cam_rol~=0 + 0.1rad)
    roldummy=cam_rol*pi/180; % deg to rad
    FOV_vert=atan((0.5*cam_hei*abs(cos(roldummy))+0.5*cam_wid*abs(sin(roldummy))) ...
                    /cam_foc)+0.1;
    FOV_hor=atan((0.5*cam_wid*abs(cos(roldummy))+0.5*cam_hei*abs(sin(roldummy))) ...
                    /cam_foc)+0.1;
    clear roldummy
%       calculate azimuth angles (clockwise from Azi(1,1) to Azi(2,1),
%       with N=0, E=pi/2, S=pi, W=3*pi/2)
%           delta DEM in the viewing direction [DEM pixel] -> Cam target-Cam position
    N0_rc(1,1)=cam_rc(1,2)-cam_rc(1,1); % delta longitude (=delta cols)
    N0_rc(2,1)=cam_rc(2,1)-cam_rc(2,2); % delta latitude (=-delta rows)
    N0_rc(3,1)=cam_rc(3,2)-cam_rc(3,1); % delta altitude
%           calculate angles in viewing direction
    if N0_rc(2,1)>0 && N0_rc(1,1)>0 % delta latitude, delta longitude > 0 = NE
        N0_rc_hor=pi/2-atan(N0_rc(2,1)/N0_rc(1,1));
    elseif N0_rc(2,1)<0 && N0_rc(1,1)>0 % delta latitude < 0 & delta longitude > 0 = SE
        N0_rc_hor=pi/2-atan(N0_rc(2,1)/N0_rc(1,1));
    elseif N0_rc(2,1)<0 && N0_rc(1,1)<0 % delta latitude, delta longitude < 0 = SW
        N0_rc_hor=3*pi/2-atan(N0_rc(2,1)/N0_rc(1,1));
    elseif N0_rc(2,1)>0 && N0_rc(1,1)<0 % delta latitude > 0 & delta longitude < 0 = NW
        N0_rc_hor=3*pi/2-atan(N0_rc(2,1)/N0_rc(1,1));
    elseif N0_rc(2,1)==0 && N0_rc(1,1)==0
        error(['Camera position and target at the same xy ', ...
            'coordinate, observing the sky or the soil in nadir ', ...
            'position does not work yet.'])
    elseif N0_rc(2,1)==0
        if N0_rc(1,1)>0
            N0_rc_hor=pi/2;
        elseif N0_rc(1,1)<0
            N0_rc_hor=3*pi/2;
        end
    elseif N0_rc(1,1)==0
        if N0_rc(2,1)>0
            N0_rc_hor=0;
        elseif N0_rc(2,1)<0
            N0_rc_hor=pi;
        end
    end
%           calculate azimuth boundaries (with N=0, E=pi/2, ...)
    Azi(1,1)=N0_rc_hor-FOV_hor;
    if Azi(1,1)<0
        Azi(1,1)=2*pi+Azi(1,1);
    end
    Azi(2,1)=N0_rc_hor+FOV_hor;
    if Azi(2,1)>2*pi
        Azi(2,1)=Azi(2,1)-2*pi;
    end
%           define viewshed sectors with NE=1, SE=2, SW=3, NW=4
    for i=1:2
        if Azi(i,1)>=0 && Azi(i,1)<pi/2 %N-E
            Azi(i,2)=1;
        elseif Azi(i,1)>=pi/2 && Azi(i,1)<pi %E-S
            Azi(i,2)=2;
        elseif Azi(i,1)>=pi && Azi(i,1)<3*pi/2 %S-W
            Azi(i,2)=3;
        elseif Azi(i,1)>=3*pi/2 && Azi(i,1)<2*pi %W-N
            Azi(i,2)=4;
        end
    end
    if Azi(1,1)>Azi(2,1) % Azimuth angle including N (2*pi to 0)
        azi_sec=[Azi(1,2):4, 1:Azi(2,2)];
    else
        azi_sec=[Azi(1,2):Azi(2,2)];
    end
    azi_sec_name=['Northeast'; 'Southeast'; 'Southwest'; 'Northwest'];
    azi_sec_name=azi_sec_name(azi_sec,:);
    disp('The viewshed for sector(s) ')
    disp(azi_sec_name)
    disp(['will be calculated and no sight barriers within ', ...
        num2str(buffer_radius),  'm distance to the camera location are assumed.'])
    clear(azi_sec_name)
%       calculate vertical angles (from Vert(1,1) to Vert(2,1); upwards to downwards)
%           upwards=+, horizontal=0, downwards=-
    N0_rc_vert=atan(N0_rc(3,1)/norm(N0_rc(1:2,1))); % norm(N0_rc(1:2,1)) is never 0
    Vert(1,1)=N0_rc_vert+FOV_vert;
    Vert(2,1)=N0_rc_vert-FOV_vert;
    if Vert(1,1)>pi/2
        Vert(1,1)=pi/2;
        disp(['Warning: Sky observation angle adjusted to pi/2 as ', ...
         'higher values are not usable. If wanted, use 2 different ', ...
         'viewshed calculations and intersect them!'])
    elseif Vert(2,1)<-pi/2
        Vert(2,1)=-pi/2;
        disp(['Warning: Surface observation angle adjusted to -pi/2 ', ...
         'as lower values are not usable. If wanted, use 2 different ', ...
         'viewshed calculations and intersect them!'])
    end
    clear FOV_hor FOV_vert
%       viewshed raster initialization
    viewW=-1*ones(headerv_W(2,1), headerv_W(1,1));
    viewW(isnan(demWvs))=NaN;
%   initialization of temporary rasters and variables
%       R=maximal normalized altitude of visibile pixel
%       index=number-indexed raster (row 1 & col 1 = 1, row 2 & col 1 = 2, ...)
    R=viewW;
    index=reshape([1:headerv_W(2,1)*headerv_W(1,1)], headerv_W(2,1), []);
    r=cam_rc(2,1);
    c=cam_rc(1,1);
%   viewshed
%       2 process steps for specified camera angles computation:
%           1) algorithm assumes that all pixels in the chosen sectors
%              that lie inbetween the vertical angles are visible
%           2) but azimuth angles are saved in case of vertical visibility and
%              taken out later if they do not fulfill the azimuth conditions
%       calculation of the viewpoint and the 8 neighbouring pixels
    viewW(r-1,c+1)=pi/4; %NE
    viewW(r+1,c+1)=3/4*pi; %SE
    viewW(r+1,c-1)=pi+pi/4; %SW
    viewW(r-1,c-1)=pi+3/4*pi; %NW
    viewW(r-1,c)=0; %N
    viewW(r+1,c)=pi; %S
    viewW(r,c-1)=3*pi/2; %W
    viewW(r,c+1)=pi/2; %E
    for i=-1:1
        for j=-1:1
            if i~=0 || j~=0
                vertdummy=atan(demWvs(r+i,c+j)/sqrt(i^2+j^2));
                if vertdummy>Vert(1,1) || vertdummy<Vert(2,1)
                    viewW(r+i,c+j)=-1;
                end
            else
                if Vert(2,1)>-pi/2
                    viewW(r,c)=-1;
                else
                    viewW(r,c)=N0_rc_hor;
                end
            end
        end
    end
    R(r-1:r+1,c-1:c+1)=demWvs(r-1:r+1,c-1:c+1);
%       calculation of the N, E, S and W directions, seen from the viewpoint
%           using proxy vectors
    for dir=1:4
        if dir==1 % N
%               flipud to S direction (=from viewpoint to the south)
            demdummy=flipud(demWvs(1:r-1,c));
            viewdummy=flipud(viewW(1:r-1,c));
            Rdummy=flipud(R(1:r-1,c));
            indexdummy=flipud(index(1:r-1,c));
            azidummy=0;
        elseif dir==2 % E
%               transpose to S direction
            demdummy=demWvs(r,c+1:headerv_W(1,1))';
            viewdummy=viewW(r,c+1:headerv_W(1,1))';
            Rdummy=R(r,c+1:headerv_W(1,1))';
            indexdummy=index(r,c+1:headerv_W(1,1))';
            azidummy=pi/2;
        elseif dir==3 % S
            demdummy=demWvs(r+1:headerv_W(2,1),c);
            viewdummy=viewW(r+1:headerv_W(2,1),c);
            Rdummy=R(r+1:headerv_W(2,1),c);
            indexdummy=index(r+1:headerv_W(2,1),c);
            azidummy=pi;
        elseif dir==4 % W
%               transpose and flipud to S direction
            demdummy=flipud(demWvs(r,1:c-1)');
            viewdummy=flipud(viewW(r,1:c-1)');
            Rdummy=flipud(R(r,1:c-1)');
            indexdummy=flipud(index(r,1:c-1)');
            azidummy=3*pi/2;
        end
%           S direction or rotated/mirrored other directions to S direction
        for l=2:size(demdummy, 1)
            Z=Rdummy(l-1,1)*l/(l-1);
            if isfinite(demdummy(l,1))
                if demdummy(l,1)>Z
                    vertdummy=atan(demdummy(l)/l);
                    if vertdummy>Vert(1,1) || vertdummy<Vert(2,1)
                        viewdummy(l,1)=-1; % invisible
                    else
                        viewdummy(l,1)=azidummy; % vertically visible
                    end
                    Rdummy(l,1)=demdummy(l,1);
                else
                    viewdummy(l,1)=-1; % invisible
                    Rdummy(l,1)=Z;
                end
            else
                viewdummy(l,1)=NaN; % like invisible, but NaN
                Rdummy(l,1)=Z;
            end
        end
        viewW(indexdummy)=viewdummy;
        R(indexdummy)=Rdummy;
        clear demdummy viewdummy Rdummy indexdummy
    end
    clear l dir azidummy
%       calculation of the chosen sectors
%           using proxy arrays
%               W-NW + NW-N = W-N, S-SW + SW-W = S-W,
%               N-NE + NE-E = N-E, E-SE + SE-S = E-S)
    for l=1:length(azi_sec)
        dir=azi_sec(l);
        if dir==1 % NE-E & N-NE
%               leftrright array (NE flipped to NW)
            demdummy=fliplr(demWvs(1:r,c:headerv_W(1,1)));
            viewdummy=fliplr(viewW(1:r,c:headerv_W(1,1)));
            Rdummy=fliplr(R(1:r,c:headerv_W(1,1)));
            indexdummy=fliplr(index(1:r,c:headerv_W(1,1)));
            secdummy=-pi/2; % angle clockwise from N=0
        elseif dir==2 % E-SE & SE-S
%               upsidedown + leftrright array (SE flipped to NW)
            demdummy=flipud(fliplr(demWvs(r:headerv_W(2,1),c:headerv_W(1,1))));
            viewdummy=flipud(fliplr(viewW(r:headerv_W(2,1),c:headerv_W(1,1))));
            Rdummy=flipud(fliplr(R(r:headerv_W(2,1),c:headerv_W(1,1))));
            indexdummy=flipud(fliplr(index(r:headerv_W(2,1),c:headerv_W(1,1))));
            secdummy=pi/2;
        elseif dir==3 % W-SW & S-SW
%               upsidedown array (SW flipped to NW)
            demdummy=flipud(demWvs(r:headerv_W(2,1),1:c));
            viewdummy=flipud(viewW(r:headerv_W(2,1),1:c));
            Rdummy=flipud(R(r:headerv_W(2,1),1:c));
            indexdummy=flipud(index(r:headerv_W(2,1),1:c));
            secdummy=-3*pi/2;
        elseif dir==4 % W-NW & NW-N
            demdummy=demWvs(1:r,1:c);
            viewdummy=viewW(1:r,1:c);
            Rdummy=R(1:r,1:c);
            indexdummy=index(1:r,1:c);
            secdummy=3*pi/2;
        end
        rr=size(demdummy, 1);
        cc=size(demdummy, 2);
        m=rr-2;
        n=cc-1;
%           W-NW and NW-N = NW sector or rotated/mirrored other sectors to NW sector
        while ( n>0 ) % cols
            if ( m>0 ) % rows
                if ( rr-m<cc-n ) % "south of NW-diagonal"
                    Z=(rr-m)*(Rdummy(m,n+1)-Rdummy(m+1,n+1)) ...
                        +(cc-n)*((rr-m)*(Rdummy(m+1,n+1)-Rdummy(m,n+1))+Rdummy(m,n+1)) ...
                          /(cc-n-1);
                elseif ( rr-m>cc-n ) % "north of NW-diagonal"
                    Z=(cc-n)*(Rdummy(m+1,n)-Rdummy(m+1,n+1)) ...
                        +(rr-m)*((cc-n)*(Rdummy(m+1,n+1)-Rdummy(m+1,n))+Rdummy(m+1,n)) ...
                          /(rr-m-1);
                elseif ( rr-m==cc-n ) % "NW-diagonal"
                    Z=Rdummy(m+1,n+1)*(rr-m)/(rr-m-1);
                end
                if isfinite(demdummy(m,n))
                    if demdummy(m,n)>Z
                        azidummy=abs(atan((rr-m)/(cc-n))+secdummy);
                        vertdummy=atan(demdummy(m,n)/sqrt((rr-m)^2+(cc-n)^2));
                        if vertdummy>Vert(1,1) || vertdummy<Vert(2,1)
                            viewdummy(m,n)=-1; % invisible
                        else
                            viewdummy(m,n)=azidummy; % vertically visible
                        end
                        Rdummy(m,n)=demdummy(m,n);
                    else
                        viewdummy(m,n)=-1; % invisible
                        Rdummy(m,n)=Z;
                    end
                else
                    viewdummy(m,n)=NaN; % like invisible, but NaN
                    Rdummy(m,n)=Z;
                end
                m=m-1;
            else
                m=rr-1;
                n=n-1;
            end
        end
        clear m n
        viewW(indexdummy)=viewdummy;
        R(indexdummy)=Rdummy;
        clear demdummy viewdummy Rdummy indexdummy N0_rc N0_rc_hor ...
         N0_rc_vert
    end
    clear rr cc secdummy azidummy vertdummy Z dir
    clear R
%       visibility decision using azimuth of vertically visible pixels
    if Azi(1,1)>Azi(2,1)
        viewW(viewW<Azi(1,1) & viewW>Azi(2,1))=-1;
    else
        viewW(viewW<Azi(1,1) | viewW>Azi(2,1))=-1;
    end
    viewW(buffer_cam==1)=-1;
    viewW(viewW>=0)=1;
    viewW(viewW<0)=0;
    viewAzi=Azi;
    viewAziSec=azi_sec;
    clear demWvs R index buffer_cam Azi azi_sec N
end % if vs==1
%%%%%%%%%%%%%%%%%%%%%%% projection + classification %%%%%%%%%%%%%%%%%%%%%%%
toc
disp('Starting projection')
% prepare DEM dependent on the viewshed
%   initialization
demW=NaN(6, headerv_W(2,1)*headerv_W(1,1)); % if row, col included: 6
k=0;
%   DEM grid point coordinates (middle of pixel)
for i=1:headerv_W(1,1) % ncols
    for j=1:headerv_W(2,1) % nrows
        if viewW(j,i)==1 % visible
            k=k+1;
            demW(1,k)=headerv_W(3,1)-(0.5*headerv_W(5,1))+(i*headerv_W(5,1)); % longitude
            demW(2,k)=headerv_W(4,1)+(headerv_W(2,1)*headerv_W(5,1)) ...
                +(0.5*headerv_W(5,1))-(j*headerv_W(5,1)); % latitude
            demW(3,k)=demWrcz(j,i); %altitude
            demW(4,k)=NaN;%varWrcz(j,i); space for an additional variable
            demW(5,k)=j; %ascii row
            demW(6,k)=i; %ascii col
        end
    end
end
demW=demW(1:6,isfinite(demW(1,:))); % if row, col included: 1:6
clear i j k
% prepare projection parameters
%   in case that optimisation off
toc
if os<2
    Xopt=X0;
end
% project and scale the world coordinates (xyzW in m) to the photo plane
% (crP in pixels)
%   1) translate world reference system W to reference system T
%      -> set coordinate system origin = viewpoint
%   2) rotate reference system T to camera reference system C
%      -> set coordinate system axis in the directions right (U), up (V)
%      and depth (N) of image
%   3) project and scale camera reference system C to perspective
%   projection system P
%      -> project to metric cam image using central projection
%      (at depth of focal length)
%      -> scale metric to pixel based image size
%      -> translate to the origin of the 2D coordinate system of
%      the photograph
%      -> ceil to integer pixel values
%   4) extract RGB pixel values and classify
%      -> keep rows and columns of DEM (Arc/Info ASCII Grid)
%      -> keep pixel rows and columns of photograph
%      -> extract RGB values for pixel rows and columns
%      -> classify with one of the available methods
%
%   call steps 1) to 3) of DEM projection
%
[demP]=Proj_PRACTISE(demW, Xopt, pix_c, pix_r, demWrcz, headerv_W);
% loop for classification of multiple photographs
for photoloop=1:N_images
    if is==1 && N_images>1 && photoloop>1
%   define the respective photo file
        fin_image=char(fin_images(photoloop));
%   read photo
        photo = imread([fin_imagepath, fin_image]);
        fin_imageDot=regexp(fin_image, '\.');
        f_name=fin_image(1:fin_imageDot-1);
        clear fin_imageDot
%       check photograph size
        if size(photo, 1)~=pix_r || size(photo, 2)~=pix_c
            disp(['The current photograph has a different dimension ', ...
	        '(rows and/or columns) than all processed photographs before.'])
            error(['Please process the photograph ''', fin_image, ...
	        ''' separately to avoid projection errors.'])
        end
        if cs3==3
            cs=3;
        end
    end
%   classification
    disp(['Load and start classification of the photograph ''', ...
     fin_image, ''''])
    if photoloop==1
%       initialize for 4)
        if rs==0
%           DEMrow, ~col, PHOTOcol, ~row, RGB, PHOTOclass
            rgbimage=NaN(8, size(demW, 2));
        else
%           ..., additional variable, DEMlongitude, ~latitude
            rgbimage=NaN(11, size(demW, 2));
        end
%       4) extract RGB pixel values and classify
%           extract pixel positions inside photograph boundaries
        i=find(demP(1,:)<pix_c & demP(1,:)>0 & demP(2,:)<pix_r & demP(2,:)>0);
        rgbimage(1:2,i)=demW(5:6,i); % DEM: row col
        rgbimage(3:4,i)=demP(:,i); % image: col row
        if rs>0
            rgbimage(10:11,i)=demW(1:2,i); % DEM: longitude latitude
        end
%           delete pixels outside of the photograph boundaries
        demW=demW(1:3,i);
        demP=demP(1:2,i);
%           delete NaNs
        i=find(isnan(rgbimage(1,:)));
        rgbimage(:,i)=[];
%           extract pixel positions outside of the excluded photograph areas
%         i=find((demP(1,:)>800 | demP(1,:)<400) & demP(2,:)<pix_r & demP(2,:)>0);
%         demW(:,i)=[];
%         demP(:,i)=[];
%         rgbimage(:,i)=[];
    end
%           extract RGB values
    for i=1:length(demP)
        rgbimage(5:7,i)=photo(rgbimage(4,i),rgbimage(3,i),1:3); % RGB values
    end
    clear i
%       first assumption: interactive classification (if cs~=3 will be changed to non-interactive)
    cs_intact='y';
    cs3=cs;
%       interactive loop
    while cs_intact=='y'
        if cs==3
            cs=1;
        elseif cs3~=3
            cs_intact='n';
        end
%           classification switch
        if cs==0
            disp('Manual classification of snow: All visible pixels in investigation area')
%               initialisation
            rgbimage(8,:)=0; % no snow
%               classification rule 1
            disp('with RGB values greater or equal than the specified thresholds in the bands')
            disp(['Red ', num2str(thres_rgb(1,1)), ','])
            disp(['Green ', num2str(thres_rgb(2,1)), ','])
            disp(['Blue ', num2str(thres_rgb(3,1)), ','])
            i=find(rgbimage(5,:)>=thres_rgb(1,1) & rgbimage(6,:)>=thres_rgb(2,1) & ...
                rgbimage(7,:)>=thres_rgb(3,1));
%               classification rule 2
            disp('and a maximum delta value (maximum to minimum RGB value) of a pixel')
            disp(['less than ', num2str(delta_rgb), '.'])
            j=find((max(rgbimage(5:7,i))-min(rgbimage(5:7,i)))<=delta_rgb);
%               value assignment for complete investigation area
            rgbimage(8,i(j))=1; % snow
            clear i j
            close(figure(2));
        elseif cs>0
            if cs==1
                disp('Automatic classification of snow: All visible pixels ')
            elseif cs==2
                disp('Automatic PCA-based classification: All visible pixels ')
            end
            disp('in complete investigation area')
%               initialisation
            if cs==1
                rgbimage(8,:)=0; % no snow
            end
            % moving average window size is not an odd number? (maw+1)
            if mod(movavgwindow, 2)==0
                movavgwindow=movavgwindow+1;
                maw_dummy=1;
            else
                maw_dummy=0;
            end
            maw=movavgwindow;
            tbo=thres_b_orig;
            if cs==2
                tbl=thres_b_low;
            end
%               blue band histogram
            factor=1;
            rgbhist1=hist(rgbimage(7,:), 0:255);
            rgbhistmean1=filter(ones(1,maw)/maw,1,rgbhist1')';
            rgbhistmean1(1:maw-1)=NaN;
            rgbhistmean1=rgbhistmean1([ceil(maw/2):end, 1:ceil(maw/2)-1]);
%                   if Financial Timeseries Toolbox is available, e.g.
%            rgbhistmean1=tsmovavg(rgbhist1, 's', maw);
%            rgbhistmean1=rgbhistmean1([(ceil(maw/2)):end,1:(ceil(maw/2))-1]);
%               slope of histogram
            rgbhistmeandiff=[diff(rgbhistmean1),NaN];
%               first negative slope greater than or equal to threshold
            i=find(rgbhistmeandiff(round(tbo/factor):round(max(rgbimage(7,:))/factor))<0);
            i=i(1)+round(tbo/factor)-1;
%               local minimum (first positive slope after negative slope)
            j=find(rgbhistmeandiff(i:round(max(rgbimage(7,:))/factor))>0);
            if j<max(rgbimage(7,:))
                thres_b_image1=(j(1)+i-1-1)*factor;
            else
                thres_b_image1=tbo;
            end
            clear i j rgbhistmeandiff
            disp(['with blue band values greater than or equal to ', num2str(thres_b_image1)])
            if cs==2
                disp('are classified as snow')
            end
            if thres_b_image1==tbo
                disp(['as no local minima in the blue band histogram (moving average window size ', ...
		         num2str(maw), ') was found above'])
            end
            % moving average window size is an odd number? (display)
            if maw_dummy==1
                disp(['Warning: The moving average window size was', ...
                 'not odd, hence the window size was increased by one'])
                clear maw_dummy
            end
            if cs==2
                disp(['and all visible pixels with blue band values smaller than ', ...
                 num2str(tbl), ' are assumed as rock'])
            end
            if cs==1
%               classification rule: snow
                isnow1=find(rgbimage(7,:)>=thres_b_image1);
            elseif cs==2
                rc_rgbpca=rgbimage';
                rgb=rc_rgbpca(:,5:7); %rgb
                clear rc_rgbpca
                rgb_centnorm=(rgb-repmat(mean(rgb), size(rgb, 1), ...
                 1))./repmat(std(rgb), size(rgb, 1), 1);
                [dummyU dummyS rgb_pc]=svd(rgb_centnorm, 0);
                rgb_sc=rgb_centnorm*rgb_pc;
                clear rgb_centnorm dummyU dummyS
%                   if Statistics Toolbox is available, e.g.
%                 [pc, sc]=princomp(zscore(rgb));
%               normalise values
                pca=(rgb_sc-min(min(rgb_sc)))/(max(max(rgb_sc))-min(min(rgb_sc)));
%                 pca=(sc-min(min(sc(:,1))))/(max(max(sc(:,1)))-min(min(sc(:,1))));
                clear rgb_sc % pc sc
%               classification rules
%                   snow (sun)
                isnow1=find(rgb(:,3)>=thres_b_image1);
%                   snow (shadow)
                isnow2=find((pca(:,3)<pca(:,2) & rgb(:,3)>=tbl) & rgb(:,3)<thres_b_image1);
%                   rock (sun)
                irock=find(~((pca(:,3)<pca(:,2) & rgb(:,3)>=tbl) | rgb(:,3)>=thres_b_image1) & rgb(:,3)<rgb(:,1));
%                   unsure
                i5050=find(~((pca(:,3)<pca(:,2) & rgb(:,3)>=tbl) | rgb(:,3)>=thres_b_image1) & rgb(:,3)>=rgb(:,1));
% %               classification rules
% %                   snow (sun & shadow)
%                 isnow=find(pca(:,3)<pca(:,2) | rgb(:,3)>=thres_b_image1);
% %                   rock (sun)
%                 irock=find(~((pca(:,3)<pca(:,2) | rgb(:,3)>=thres_b_image1)) & rgb(:,3)<rgb(:,1));
% %                   unsure
%                 i5050=find(~((pca(:,3)<pca(:,2) | rgb(:,3)>=thres_b_image1)) & rgb(:,3)>=rgb(:,1));
                clear rgb
            end
%               value assignment for complete investigation area
            rgbimage(8,isnow1)=1; % snow
            if cs==2
                rgbimage(8,isnow2)=2; % snow
                rgbimage(8,irock)=-1; % rock
                maxi=thres_b_image1;
                mini=max([min(rgbimage(7,i5050)-1), tbl-1]);
                rgbimage(8,i5050)=(1-0)/(maxi-mini)*(rgbimage(7,i5050) ...
                 -mini); % unsure
                rgbimage(8,rgbimage(8,:)>1)=1;
                rgbimage(8,rgbimage(8,:)<0)=0;
                clear maxi mini
            end
%           plot blue band (DN frequency) histogram
%            close(figure(2));
%            if IsOctave_PRACTISE()
%                figure2=figure(2,'position',get(0,'screensize')+[100 50 -150 -150]);
%            else
%                figure2=figure(2);
%                set(figure2, 'units', 'normalized', 'outerposition', [0 0 1 1], ...
%                  'PaperUnits', 'centimeters', 'PaperPosition', [0 0 50 35], ...
%                  'PaperSize', [50 35], 'PaperType', '<custom>');
%            end
%            set(figure2, 'NumberTitle', 'off', 'Name', 'Blue band histogram');
%            set(figure2, 'Color', [1 1 1], 'InvertHardcopy', 'off');
%			axes1 = axes('FontSize', 15.91, 'TickDir', 'out', ...
%			  'Parent', figure2);
%			if ~IsOctave_PRACTISE()
%				set(axes1, 'FontWeight', 'demi')
%			end
%            xlabel(axes1, 'Digital Number', 'FontSize', 18.185);
%            ylabel(axes1, 'Number of pixels', 'FontSize', 18.185);
%            box(axes1, 'on');
%            hold(axes1, 'all');
%            %
%            hist(rgbimage(7,:), 0:255)
%            figure2hist = findobj(gca, 'Type', 'patch');
%            set(figure2hist, 'EdgeColor', 'none', 'FaceColor', [0 0.749 0.749])
%            axis1=axis;
%            axis([0 255 0 axis1(4)])
%            plot(0:factor:255, rgbhistmean1, 'LineWidth', 3, 'Color', [0.0784 0.1686 0.5490])
%            plot([tbo tbo], [0 axis1(4)], 'r', 'LineWidth', 2)
%            plot([thres_b_image1 thres_b_image1], [0 axis1(4)], 'g', 'LineWidth', 2)
%            legend('Blue band DN frequency', 'Smoothed blue band DN frequency', ...
%                'Snow threshold (seed)', 'Snow threshold (optimum)', 'Location', 'NorthWest');
%            %
%            if IsOctave_PRACTISE()
%                fig2_leg=get(figure2,'children');
%                fig2_leg_pos=get(fig2_leg(1), 'position');
%                set(fig2_leg(1), 'position', fig2_leg_pos+[0, 0, 0.03, 0]);
%            end
%            %
            clear figure2 axes1 figure2hist axis1 fig2_leg fig2_leg_pos
            clear tbo maw
        end
%       plot photograph, classification (projection) and GCP
        close(figure(4));
        if IsOctave_PRACTISE()
            figure4=figure(4,'position',get(0,'screensize')+[100 50 -150 -150]);
        else
            figure4=figure(4);
            set(figure4, 'units', 'normalized', 'outerposition', [0 0 1 1], ...
            'PaperUnits', 'centimeters', 'PaperPosition', [0 0 50 35], ...
            'PaperSize', [50 35], 'PaperType', '<custom>');
        end
        if cs==0
            set(figure4, 'NumberTitle', 'off', 'Name', 'Manual snow cover classification');
        elseif cs==1
            set(figure4, 'NumberTitle', 'off', 'Name', 'Projection Result');
        elseif cs==2
            set(figure4, 'NumberTitle', 'off', 'Name', 'Automatic blue band + pca snow cover classification');
        end
        set(figure4, 'Color', [1 1 1], 'InvertHardcopy', 'off');
        axes1 = axes('FontSize', 15.91, 'TickDir', 'out', ...
          'XAxisLocation', 'top', 'YDir', 'reverse', ...
          'DataAspectRatio', [1 1 1], 'Parent', figure4);
		if ~IsOctave_PRACTISE()
			set(axes1, 'FontWeight', 'demi')
		end
        axis(axes1, [0.5 pix_c+0.5 0.5 pix_r+0.5]);
        xlabel(axes1, 'Pixel column', 'FontSize', 18.185);
        ylabel(axes1, 'Pixel row', 'FontSize', 18.185);
        box(axes1, 'on');
        hold(axes1, 'all');
        %
        image(photo, 'Parent', axes1);
        if cs==0 || cs==1
            plot(rgbimage(3,(rgbimage(8,:)>=0)), rgbimage(4,(rgbimage(8,:)>=0)), '.', 'Color', [1 0 0], 'MarkerSize', cs_MarkSiz)
%            plot(rgbimage(3,(rgbimage(8,:)==1)), rgbimage(4,(rgbimage(8,:)==1)), '.', 'Color', [1 0 0], 'MarkerSize', cs_MarkSiz)
            Leg=legend('DEM point');
        elseif cs==2
            if IsOctave_PRACTISE()
                cs_MarkSiz2=cs_MarkSiz;
            else
                cs_MarkSiz2=cs_MarkSiz+2;
            end
            plot(rgbimage(3,(rgbimage(8,:)==0)), rgbimage(4,(rgbimage(8,:)==0)), '.', 'Color', [0 0 1], 'MarkerSize', cs_MarkSiz)
            plot(rgbimage(3,(rgbimage(8,:)<=0.25 & rgbimage(8,:)>0)), rgbimage(4,(rgbimage(8,:)<=0.25 & rgbimage(8,:)>0)), ...
	        '.', 'Color', [0 0.5 1], 'MarkerSize', cs_MarkSiz2)
            plot(rgbimage(3,(rgbimage(8,:)<0.75 & rgbimage(8,:)>0.25)), rgbimage(4,(rgbimage(8,:)<0.75 & rgbimage(8,:)>0.25)), ...
	        '.', 'Color', [1 1 0], 'MarkerSize', cs_MarkSiz2)
            plot(rgbimage(3,(rgbimage(8,:)<1 & rgbimage(8,:)>=0.75)), rgbimage(4,(rgbimage(8,:)<1 & rgbimage(8,:)>=0.75)), ...
	        '.', 'Color', [1 0.5 0], 'MarkerSize', cs_MarkSiz2)
            plot(rgbimage(3,(rgbimage(8,:)==1)), rgbimage(4,(rgbimage(8,:)==1)), '.', 'Color', [1 0 0], 'MarkerSize', cs_MarkSiz)
            Leg=legend('No snow', 'No snow (>=75%)', 'No snow/Snow (<75%)', 'Snow (>=75%)', 'Snow');
        end
        if IsOctave_PRACTISE()
            set(Leg, 'location', 'southoutside', 'color', [0.75 0.75 0.75])
            fig4_leg=get(figure4,'children');
            fig4_leg_pos=get(fig4_leg(1), 'position');
            set(fig4_leg(1), 'position', fig4_leg_pos+[0, 0, 0.03, 0]);
        else
            LegProp=get(Leg, 'children');
            set(LegProp([1:3:end]), 'markersize', 10);
        end
        %
        clear figure4 axes1 Leg LegProp fig4_leg fig4_leg_pos
        % interactive questions
        if cs_intact=='y'
            if cs==1
                cs_quest=['Do you want to try a classification method ', ...
                 'other than the blue band classification?'];
                cs_answ1='Yes using the manual classification';
                cs_answ2='Yes using the combined blue band and pca classification';
            elseif cs==0
                button=['Do you want to try a classification method ', ...
                 'other than the manual classification?'];
                cs_answ1='Yes using the blue band classification';
                cs_answ2='Yes using the combined blue band and pca classification';
            elseif cs==2
                cs_quest=['Do you want to try a classification method ', ...
                 'other than the combined blue band and pca classification?'];
                cs_answ1='Yes using the manual classification';
                cs_answ2='Yes using the blue band classification';
            end
            cs_answ3='No';
            cs_title='Snow classification';
            if photoloop==1 && os~=3 && ~exist('oscs_ok1_dummy', 'var')
%                 uisetpref('clearall'); % if you want to reset the dialog box properties
                if ~IsOctave_PRACTISE()
                    oscs_ok1=uigetpref('userhelp', 'ctrlc', 'Hint', ...
                    ['If switching between different windows ', ...
                    'while a message box is open is not possible. ', ...
                    'Use the keyboard combination ''Ctrl + C'' when the ', ...
                    'message box is active to allow this option.'], ...
                    'Ok');
                end
                oscs_ok1_dummy=1;
            end
            button=questdlg(cs_quest, cs_title, cs_answ1, cs_answ2, cs_answ3, cs_answ3);
            switch button
                case ''
                    disp('PRACTISE has been terminated by user.')
                    diary off
                    return
                case 'Yes using the manual classification'
                    cs=0;
                case 'Yes using the blue band classification'
                    cs=1;
                case 'Yes using the combined blue band and pca classification'
                    cs=2;
                case 'No'
                    cs_intact='n';
            end
            if cs==1
                cs_enter={'Threshold blue band'; 'Moving average window'};
                cs_def={'127'; '5'};
                cs_title='Blue band classification';
            elseif cs==2
                cs_enter={'Threshold blue band (snow)'; 'Moving average window (snow)'; ...
                    'Threshold blue band (rock)'};
                cs_def={'127'; '5'; '63'};
                cs_title='Combined blue band and PCA classification';
            elseif cs==0
                cs_enter={'Threshold red band'; 'Threshold green band'; ...
                    'Threshold blue band'; ...
                    'Maximum delta between maximum and minimum RGB value of a pixel'};
                cs_def={'127'; '127'; '127'; '10'};
                cs_title='Manual classification';
            end
            if cs_intact=='y'
                if IsOctave_PRACTISE()
                    cs_answer=inputdlg(cs_enter, cs_title, 1, cs_def);
                else
                    cs_answer=inputdlg(cs_enter, cs_title, 1, cs_def, 'on');
                end
            end
            if ~exist('cs_answer', 'var')
            elseif isempty(cs_answer)
                disp('PRACTISE has been terminated by user.')
                diary off
                return
            else
                if cs>0
                    thres_b_orig=str2num(cs_answer{1});
                    movavgwindow=str2num(cs_answer{2});
					if cs==2
                        thres_b_low=str2num(cs_answer{3});
                    end
                elseif cs==0
                    thres_rgb=[str2num(cs_answer{1});str2num(cs_answer{2}); ...
                    str2num(cs_answer{3})];
                    delta_rgb=str2num(cs_answer{4});
                end
            end
            clear cs_answer
        end
    end
    clear cs_intact button cs_quest cs_title cs_answ1 cs_answ2 cs_answ3 ...
        cs_enter cs_def maw_dummy oscs_ok1 oscs_ok1_dummy cs_MarkSiz2
%   remote sensing classification
    if rs>0
        toc
        disp('Start satellite image classification')
        if rs==1
            disp('of a raw Landsat scene using a photo-calibrated NDSI threshold')
        elseif rs==2
            disp('of a preprocessed NDSI map using a photo-calibrated threshold')
        end
%       A) control that rgbimage/DEM pixel size is smaller or equal to the Landsat pixel
%          size, in a optimum case at least 10 times smaller
        if satMETA.cellsize<headerv_W(5)
            error(['Nearest neighbour/majority approach using pixel centres seems to be not valid: ', ...
                'DEM pixel size is larger than the satellite pixel size of ', num2str(satMETA.cellsize), 'm.'])
        elseif satMETA.cellsize<10*headerv_W(5)
            disp(['Warning: Check if the nearest neighbour/majority approach using pixel centres is valid: ', ...
                'satellite pixel size (', num2str(satMETA.cellsize), 'm) and DEM pixel size (', ...
                num2str(headerv_W(5)), 'm).'])
        end
%           control the timeshift between photo and satellite image
        dummy=satJDAY-datenum(time.year, time.month, time.day, ...
            time.hour-time.UTC, time.min, time.sec);
        disp('The timeshift between photo and satellite image is ')
        if abs(dummy)<1
            disp([num2str(abs(dummy)*24), ' hours.'])
        elseif abs(dummy)>=1
            disp([num2str(abs(dummy)), ' day(s).'])
        end
        if abs(dummy)>=1
            disp(['Warning: Check if the timeshift between satellite ', ...
             'image and photograph is valid for a calibration.'])
        end
        clear dummy
%       B) calculate NDSI of Landsat image (if necessary and maybe use given bounding box)
%           Bounding box definitions
        if strcmp(satBB, 'full')
            disp('Satellite image classification will be calculated for the full extent of the image:')
%               define bounding box (N, E, S, W)
            satBBin_ext=satBB;
            clear satBB
            satBB=[satMETA.y_ul+0.5*satMETA.cellsize, ...
                satMETA.x_ul+satMETA.ncols*satMETA.cellsize-0.5*satMETA.cellsize, ...
                satMETA.y_ul-satMETA.nrows*satMETA.cellsize+0.5*satMETA.cellsize, ...
                satMETA.x_ul-0.5*satMETA.cellsize];
        elseif all(isnumeric(satBB))
            if length(satBB)==1
                satBBin_ext=satBB;
                satBB=[max(rgbimage(11,:))+satBBin_ext, max(rgbimage(10,:))+satBBin_ext, ...
                    min(rgbimage(11,:))-satBBin_ext, min(rgbimage(10,:))-satBBin_ext];
            elseif length(satBB)>1 && length(satBB)~=4
                error(['Check bounding box definition, either use four coordinates (N, E, S, W), ', ...
                'create a rectangle buffer area (1 value) around the maximum photo snow map extent', ...
                'or define the string ''full'' for full extent.'])
            end
            disp('Satellite image classification will be calculated for a bounding box:')
        else
            error(['Check bounding box definition, either use four coordinates (N, E, S, W), ', ...
                'create a rectangle buffer area (1 value) around the maximum photo snow map extent', ...
                'or define the string ''full'' for full extent.'])
        end
%               test bounding box coordinates
        if (satBB(2)<satBB(4) || satBB(1)<satBB(3))
            error(['Check if the bounding box coordinates and the order ', ...
                '(N, E, S, W) are correct .'])
        end
%               extract satellite image pixel centres inside bounding box
%                   latitude and longitude
        satLAT=(satMETA.y_ul:-satMETA.cellsize:satMETA.y_ul-(satMETA.nrows-1)*satMETA.cellsize)'; % N to S
        i=find(satLAT<=satBB(1) & satLAT>=satBB(3));
        satBBlat=satLAT(i);
        satLONG=satMETA.x_ul:satMETA.cellsize:satMETA.x_ul+(satMETA.ncols-1)*satMETA.cellsize; % W to E
        j=find(satLONG>=satBB(4) & satLONG<=satBB(2));
        satBBlong=satLONG(j);
        clear satLAT satLONG
        satBBext=[max(satBBlat)+0.5*satMETA.cellsize, ...
            max(satBBlong)+0.5*satMETA.cellsize, ...
            min(satBBlat)-0.5*satMETA.cellsize, ...
            min(satBBlong)-0.5*satMETA.cellsize];
        if satBB(1)>satBBext(1)+0.5*satMETA.cellsize || ...
         satBB(2)>satBBext(2)+0.5*satMETA.cellsize || ...
         satBB(3)<satBBext(3)-0.5*satMETA.cellsize || ...
         satBB(4)<satBBext(4)-0.5*satMETA.cellsize
            disp(['Warning: The extent of the desired bounding box ', ...
             'is larger than the used satellite image.'])
        end
        disp(['Used borders (in m) are ', ...
            'in the North: ', num2str(satBBext(1)), ...
            ', in the East: ', num2str(satBBext(2)), ...
            ', in the South: ', num2str(satBBext(3)), ...
            ', and in the West: ', num2str(satBBext(4))])
        satBBlat=repmat(satBBlat, 1, length(j));
        satBBlong=repmat(satBBlong, length(i), 1);
%           calculate NDSI
        if rs==1
%               green and mid-infrared band
            satBBg_dn=satG(i,j);
            satBBmir_dn=satMIR(i,j);
            satBBnir_dn=satNIR(i,j);
            if exist('satLook', 'var')
                satBBlook=satLook(i,j,:);
            end
            clear i j satG satMIR satNIR satLook
%               set no data pixel to NaN
            i=find(satBBg_dn==satMETA.nodata | satBBmir_dn==satMETA.nodata | satBBnir_dn==satMETA.nodata);
            satBBg_dn=double(satBBg_dn);
            satBBmir_dn=double(satBBmir_dn);
            satBBnir_dn=double(satBBnir_dn);
            satBBg_dn(i)=NaN;
            satBBmir_dn(i)=NaN;
            satBBnir_dn(i)=NaN;
            clear i
%               convert Digital Numbers to Reflectances
            if strcmp(fin_satname_meta(1:3), 'LT5') || strcmp(fin_satname_meta(1:3), 'LE7')
%                   for Landsat 7
%                       convert Digital Numbers to Radiances
                satBBg_rad=((satMETA.bGlmax-satMETA.bGlmin)/(satMETA.bGqcalmax-satMETA.bGqcalmin))* ...
                    (satBBg_dn-satMETA.bGqcalmin)+satMETA.bGlmin;
                satBBmir_rad=((satMETA.bMIRlmax-satMETA.bMIRlmin)/(satMETA.bMIRqcalmax-satMETA.bMIRqcalmin))* ...
                    (satBBmir_dn-satMETA.bMIRqcalmin)+satMETA.bMIRlmin;
                satBBnir_rad=((satMETA.bNIRlmax-satMETA.bNIRlmin)/(satMETA.bNIRqcalmax-satMETA.bNIRqcalmin))* ...
                    (satBBnir_dn-satMETA.bNIRqcalmin)+satMETA.bNIRlmin;
%                       set negative Radiances (due to Landsat 7 metadata parameters) to Zero
                satBBg_rad(satBBg_rad<0)=0;
                satBBmir_rad(satBBmir_rad<0)=0;
                satBBnir_rad(satBBnir_rad<0)=0;
%                       set DN NaN to no data pixel and convert to 8bit uint again
                i=find(isnan(satBBg_dn) | isnan(satBBmir_dn) | isnan(satBBnir_dn));
                satBBg_dn(i)=satMETA.nodata;
                satBBmir_dn(i)=satMETA.nodata;
                satBBnir_dn(i)=satMETA.nodata;
                clear i
                satBBg_dn=uint8(satBBg_dn);
                satBBmir_dn=uint8(satBBmir_dn);
                satBBnir_dn=uint8(satBBnir_dn);
%                       calculate the Day Of Year (DOY) of the Landsat 7 image
                dummy=datevec(satJDAY);
                satDOY=floor(satJDAY-datenum(dummy(1)-1, 12, 31));
                clear dummy
%                       convert Radiances to Reflectances (according to Landsat 7 handbook)
                satBBg_refl=pi*satBBg_rad*(satDistEtoS(satDOY,2)^2)/(satRadSol(satRadSol(:,1)==2,2) ...
                    *cosd(90-satMETA.sunelev));
                satBBmir_refl=pi*satBBmir_rad*(satDistEtoS(satDOY,2)^2)/(satRadSol(satRadSol(:,1)==5,2) ...
                    *cosd(90-satMETA.sunelev));
                satBBnir_refl=pi*satBBnir_rad*(satDistEtoS(satDOY,2)^2)/(satRadSol(satRadSol(:,1)==4,2) ...
                    *cosd(90-satMETA.sunelev));
            elseif strcmp(fin_satname_meta(1:3), 'LC8')
%                   for Landsat 8
                satBBg_refl=(satMETA.bGreflmult*satBBg_dn+satMETA.bGrefladd)/cosd(90-satMETA.sunelev);
                satBBmir_refl=(satMETA.bMIRreflmult*satBBmir_dn+satMETA.bMIRrefladd)/cosd(90-satMETA.sunelev);
                satBBnir_refl=(satMETA.bNIRreflmult*satBBnir_dn+satMETA.bNIRrefladd)/cosd(90-satMETA.sunelev);
%                       set negative Reflectances (due to Landsat 8
%                       metadata parameters) to Zero
                satBBg_refl(satBBg_refl<0)=0;
                satBBmir_refl(satBBmir_refl<0)=0;
                satBBnir_refl(satBBnir_refl<0)=0;
%                       set DN NaN to no data pixel and convert to 8bit uint again
                i=find(isnan(satBBg_dn) | isnan(satBBmir_dn) | isnan(satBBnir_dn));
                satBBg_dn(i)=satMETA.nodata;
                satBBmir_dn(i)=satMETA.nodata;
                satBBnir_dn(i)=satMETA.nodata;
                clear i
                satBBg_dn=uint16(satBBg_dn);
                satBBmir_dn=uint16(satBBmir_dn);
                satBBnir_dn=uint16(satBBnir_dn);
            end
%               calculate NDSI from reflectances (band 2 and 5)
            satBBndsi=(satBBg_refl-satBBmir_refl)./(satBBg_refl+satBBmir_refl);
%                   no classification decision if NIR reflectances <= 0.11 (presumed to be water)
            i=find(satBBnir_refl<=0.11);
            satBBndsi(i)=2;
            clear i
        elseif rs==2
            satBBndsi=satNDSI(i,j);
            clear i j
        end
%           use mask to remove not clear pixels from NDSI map
        if rsms>0
%               extent of mask (N, E, S and W)
            satmaskEXT(1)=headerv_satmask(4,1)+ ...
             (headerv_satmask(2,1)*headerv_satmask(5,1)); %N
            satmaskEXT(2)=headerv_satmask(3,1)+ ...
             (headerv_satmask(1,1)*headerv_satmask(5,1)); % E
            satmaskEXT(3)=headerv_satmask(4,1); % S
            satmaskEXT(4)=headerv_satmask(3,1); % W
            dummy_ij(1:2)=(satmaskEXT(1:2)-satBBext(1:2))./headerv_satmask(5,1);
            dummy_ij(3:4)=(satBBext(3:4)-satmaskEXT(3:4))./headerv_satmask(5,1);
            ij=dummy_ij;
            if all(dummy_ij<0)
                disp(['Warning: Satellite image (ROI)) and mask data do ', ...
                 'not overlap, assuming clear pixel in the ROI'])
            elseif any(dummy_ij>=0) && any(dummy_ij<0)
                dummynr=find(dummy_ij<0)
                dummydir={'North', 'East', 'South', 'West'}
                disp('Warning: The satellite image ROI extent to the ')
                for i=1:length(dummynr)
                    disp(dummydir{dummynr(i)})
                end
                disp(['is larger than the available mask data, ', ...
                 'assuming clear pixel in the area not covered'])
                ij(dummynr)=0;
            end
            satmaskBB=satmask(ij(1)+1:headerv_satmask(2,1)-ij(3), ...
             ij(4)+1:headerv_satmask(1,1)-ij(2));
            if exist('dummynr', 'var')
                for i=1:length(dummynr)
                    if i==1
                        satmaskBB=[zeros(abs(dummy_ij(i)), length(satmaskBB(1,:))); ...
                         satmaskBB];
                    elseif i==2
                        satmaskBB=[satmaskBB, ...
                         zeros(length(satmaskBB(:,1)), abs(dummy_ij(i)))];
                    elseif i==3
                        satmaskBB=[satmaskBB; ...
                         zeros(abs(dummy_ij(i)), length(satmaskBB(1,:)))];
                    elseif i==4
                        satmaskBB=[zeros(length(satmaskBB(:,1)), abs(dummy_ij(i))), ...
                         satmaskBB];
                    end
                end
            end
            clear ij dummy_ij
%               apply mask on NDSI map
            satBBndsi(satmaskBB)=2;
        end
%       C) find out which rgbimage pixel is belonging to which Landsat pixel
        if (max(rgbimage(11,:))<satBB(3)) || (max(rgbimage(10,:))<satBB(4)) || ...
         (min(rgbimage(11,:))>satBB(1)) || (min(rgbimage(10,:))>satBB(2))
            error(['The photographed area is not covered by the extent of the ', ...
             'bounding box of the Landsat 7 image.'])
        elseif (max(rgbimage(11,:))>satBB(1)) || (max(rgbimage(10,:))>satBB(2)) || ...
         (min(rgbimage(11,:))<satBB(3)) || (min(rgbimage(10,:))<satBB(4))
            disp(['The photographed area is not completely covered by the extent ', ...
             'of the bounding box of the Landsat 7 image, only the overlapping area ', ...
             'will be used in the optimisation.'])
        end
%           initialise indices matrix
        index=reshape(1:size(satBBlat, 1)*size(satBBlat, 2), size(satBBlat));
%           find overlap of photographed DEM pixel and Landsat bounding box
        indBB=find(rgbimage(10,:)>=min(min(satBBlong))-0.5*satMETA.cellsize & ...
            rgbimage(10,:)<=max(max(satBBlong))+0.5*satMETA.cellsize & ...
            rgbimage(11,:)>=min(min(satBBlat))-0.5*satMETA.cellsize & ...
            rgbimage(11,:)<=max(max(satBBlat))+0.5*satMETA.cellsize);
        indBBnear=griddata(satBBlong/1e5, satBBlat/1e6, index, rgbimage(10,:)/1e5, rgbimage(11,:)/1e6, 'nearest');
%           merge overlapping photo and satellite pixel
        PhRs=NaN(2, length(rgbimage(8,:)));
%               fill with photo classification
        PhRs(1,indBB)=rgbimage(8,indBB);
%               remove optionally photo pixel with classification 'unsure' from optimisation
        if rsps==0
            i=find(rgbimage(8,:)<1 & rgbimage(8,:)>0);
            PhRs(1,i)=NaN;
            clear i
        end
%               fill with NDSI values
        PhRs(2,:)=satBBndsi(indBBnear);
        clear indBB indBBnear index
        PhRs(2,isnan(PhRs(1,:)))=NaN;
        PhRs(:,PhRs(2,:)==2)=NaN;
%           check NDSI threshold seed
        UB_NDSI=max(PhRs(2,PhRs(2,:)<=1));
        LB_NDSI=min(PhRs(2,:));
        if NDSIthres0<LB_NDSI || NDSIthres0>UB_NDSI
            disp(['The initial estimated NDSI threshold (', num2str(NDSIthres0), ...
                ') is out of bounds of the existing NDSI values in the ', ...
                'photographed area.'])
            NDSIthres0=(UB_NDSI+LB_NDSI)/2;
            disp(['The initial estimated NDSI threshold has been changed to ', ...
                num2str(NDSIthres0)])
        end
%           optimise NDSI threshold
        [NDSIbct10]=DDS_MBCT12prob_PRACTISE(NDSIthres0, PhRs(1,:), PhRs(2,:), 'BCT1', rsps);
%         [NDSIbct20]=DDS_MBCT12prob_PRACTISE(NDSIthres0, PhRs(1,:), PhRs(2,:), 'BCT2', rsps);
%         [NDSImbct0]=DDS_MBCT12prob_PRACTISE(NDSIthres0, PhRs(1,:), PhRs(2,:), 'MBCT', rsps);
        [NDSIthresOpt_detail, NDSIOpt]= ...
         DDS_PRACTISE('DDS_MBCT12prob_PRACTISE', NDSIthres0, LB_NDSI, ...
         UB_NDSI, DDS_R_rs, DDS_MaxEval_rs, PhRs(1,:), PhRs(2,:), ...
         'BCT1', rsps); % 'BCT1', rsps, 2);
        [NDSIbct1Opt]=DDS_MBCT12prob_PRACTISE(NDSIthresOpt_detail, PhRs(1,:), PhRs(2,:), 'BCT1', rsps);
%         [NDSIbct2Opt]=DDS_MBCT12prob_PRACTISE(NDSIthresOpt_detail, PhRs(1,:), PhRs(2,:), 'BCT2', rsps);
%         [NDSImbctOpt]=DDS_MBCT12prob_PRACTISE(NDSIthresOpt_detail, PhRs(1,:), PhRs(2,:), 'MBCT', rsps);
        NDSIthresOpt=round(10^2*NDSIthresOpt_detail)/10^2;;
        disp(['The optimised NDSI threshold is ', num2str(NDSIthresOpt)])
%       F) classify satellite image pixels (in bounding box) with the optimised NDSI threshold
        satBBmap=satBBndsi;
        satBBmap(satBBndsi>NDSIthresOpt_detail & satBBndsi<=1)=1;
        satBBmap(satBBndsi<=NDSIthresOpt_detail & satBBndsi>-1)=0;
%           plot LandsatLook image (or NDSI map) and the classified maps of photograph and satellite image
        close(figure(5));
        if IsOctave_PRACTISE()
            figure5=figure(5,'position',get(0,'screensize')+[100 50 -150 -150]);
        else
            figure5=figure(5);
            set(figure5, 'units', 'normalized', 'outerposition', [0 0 1 1], ...
            'PaperUnits', 'centimeters', 'PaperPosition', [0 0 50 35], ...
            'PaperSize', [50 35], 'PaperType', '<custom>');
        end
        if rs==1
            set(figure5, 'NumberTitle', 'off', 'Name', 'Optimised Landsat snow cover map');
        elseif rs==2
            set(figure5, 'NumberTitle', 'off', 'Name', 'Optimised snow cover map of preprocessed NDSI data');
        end
        set(figure5, 'Color', [1 1 1], 'InvertHardcopy', 'off');
        axes1 = axes('FontSize', 15.91, 'TickDir', 'out', ...
          'XAxisLocation', 'top', 'YDir', 'reverse', ...
          'DataAspectRatio', [1 1 1], 'Parent', figure5);
		if ~IsOctave_PRACTISE()
			set(axes1, 'FontWeight', 'demi')
		end
        axis(axes1, [min(satBBlong(1,:))-0.5 max(satBBlong(1,:))+0.5 ...
            min(satBBlat(:,1))-0.5 max(satBBlat(:,1))+0.5]);
        xlabel(axes1, 'Longitude', 'FontSize', 18.185);
        ylabel(axes1, 'Latitude', 'FontSize', 18.185);
        box(axes1, 'on');
        hold(axes1, 'all');
        %
        if exist('satBBlook', 'var')
            if IsOctave_PRACTISE()
                image(satBBlong(1,:), flipud(satBBlat(:,1)), flipud(satBBlook))
            else
                image(satBBlong(1,:), satBBlat(:,1), satBBlook)
            end
        else
            satBBndsi_dummy=satBBndsi;
            satBBndsi_dummy(satBBndsi==2)=NaN;
            if IsOctave_PRACTISE()
                imagesc(satBBlong(1,:), flipud(satBBlat(:,1)), flipud(satBBndsi_dummy))
            else
                imagesc(satBBlong(1,:), satBBlat(:,1), satBBndsi_dummy)
            end
            clear satBBndsi_dummy
        end
        set(gca, 'YDir', 'normal');
        hold on
        if cs==2 && rsps==1
            plot(rgbimage(10,(rgbimage(8,:)<1 & rgbimage(8,:)>0)), ...
             rgbimage(11,(rgbimage(8,:)<1 & rgbimage(8,:)>0)), 'y.', ...
             'MarkerSize', rs_MarkSiz)
        end
        plot(rgbimage(10,(rgbimage(8,:)==1)), ...
         rgbimage(11,(rgbimage(8,:)==1)), 'r.', 'MarkerSize', rs_MarkSiz)
        plot(rgbimage(10,(rgbimage(8,:)==0)), ...
         rgbimage(11,(rgbimage(8,:)==0)), 'b.', 'MarkerSize', rs_MarkSiz)
%         plot(satBBlong(satBBmap==0), satBBlat(satBBmap==0), '+', ...
%             'Color', [0 0.5 1], 'MarkerSize', rs_MarkSiz)
        plot(satBBlong(satBBmap==2), satBBlat(satBBmap==2), 'k+', 'MarkerSize', rs_MarkSiz+2)
        plot(satBBlong(satBBmap==1), satBBlat(satBBmap==1), 'w+', 'MarkerSize', rs_MarkSiz+1)
        if cs==2 && rsps==1
            Leg=legend('No snow/Snow (<100%, photo)', 'Snow (photo)', ...
            'No snow (photo)', 'Satellite image mask', ...
            'Snow (satellite image)');
        else
            Leg=legend('Snow (photo)', 'No snow (photo)', ...
            'Satellite image mask', 'Snow (satellite image)');
        end
        if IsOctave_PRACTISE()
            set(Leg, 'Color', [0.75 0.75 0.75], 'Location', 'SouthOutside')
            fig5_leg=get(figure5,'children');
            fig5_leg_pos=get(fig5_leg(1), 'position');
            set(fig5_leg(1), 'position', fig5_leg_pos+[0, 0, 0.03, 0]);
        else
            set(Leg, 'Color', [0.75 0.75 0.75], 'Location', 'SouthWest')
            LegProp=get(Leg, 'children');
            set(LegProp([1:3:end]), 'markersize', 6);
        end
        %
        clear figure5 axes1 Leg LegProp rs_MarkSiz fig5_leg fig5_leg_pos
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save & write %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc
    disp('start saving and writing files')
%   define output file names
%       variables (MAT-file)
    fout_mat=[f_name, '_proj.mat'];
%       photo
%           Viewshed (Arc/Info ASCII Grid)
    if vs==1 && photoloop==1
        fout_viewW_asc=[f_name, '.view.asc'];
    end
%           Classification (Arc/Info ASCII Grid)
    fout_classW_asc=[f_name, '_class.asc'];
%       satellite image
    if rs==1
%           NDSI (Arc/Info ASCII Grid)
        fout_ndsi_asc=[lower(fin_satname_meta(1:3)), '_',  f_name, ...
         '_ndsi.asc'];
%           Classification (Arc/Info ASCII Grid)
        fout_satmap_asc=[lower(fin_satname_meta(1:3)), '_', f_name, ...
         '_class.asc'];
    elseif rs==2
        fout_satmap_asc=['sat_', f_name, '_class.asc'];
    end
%       figures
%           GCPs (and DDS optimisation)
    if os==1 && photoloop==1
        fout_figname1=[f_name, '_GCPs'];
    elseif os==2 || os==3
        fout_figname1=[f_name, '_GCPs_DDS_opti'];
    end
%           photograph
    if cs==0
%               manual
        fout_figname4=[f_name, '_', num2str(thres_rgb(1,1)), ...
         '_', num2str(delta_rgb), '_manu'];
    elseif cs>0
%               automatic
%                   blue band DN frequency histogram
        fout_figname2=[f_name, '_DN_freq_hist'];
        if cs==1
%                   blue band classification
            fout_figname4=[f_name, '_', num2str(thres_b_image1), ...
             '_auto'];
        elseif cs==2
%                   blue band + pca classification
            fout_figname4=[f_name, '_', num2str(thres_b_image1), ...
             '_auto_pca'];
        end
    end
%           satellite image classification
    if rs>0
        if rs==1
            fout_figname5=[lower(fin_satname_meta(1:3)), '_', f_name, ...
             '_', regexprep(num2str(NDSIthresOpt), ...
             '\.', '_')];
        elseif rs==2
            fout_figname5=['sat_', f_name, '_', ...
             regexprep(num2str(NDSIthresOpt), '\.', ...
             '_')];
        end
    end
%   write
%       Arc/Info ASCII Grids
    if headerv_W(6,1)<=2 && headerv_W(6,1)>=-1
        header_W(6,:)=cellstr('NODATA_value  -9999');
        headerv_W(6,1)=-9999;
    end
%           photograph
%               map
    classW=NaN(headerv_W(2,1), headerv_W(1,1));
    classW(:,:)=headerv_W(6,1);
    for i=1:size(rgbimage, 2)
        classW(rgbimage(1,i),rgbimage(2,i))=rgbimage(8,i);
    end
    fid=fopen([fout_folder, fout_classW_asc], 'w');
    if fid==-1
        print=('Classification output file could not be created.');
    end
    for i=1:6
        fprintf(fid, '%s\n', char(header_W(i,:)));
    end
    fclose(fid);
    clear i fid
    dlmwrite([fout_folder, fout_classW_asc], classW, 'delimiter', ' ', '-append')
%               viewshed
    if vs==1 && photoloop==1
        viewW(classW(:,:)==headerv_W(6,1))=headerv_W(6,1);
        fid=fopen([fout_folder, fout_viewW_asc], 'w');
        if fid==-1
            print=('Viewshed output file could not be created.');
        end
        for i=1:6
            fprintf(fid, '%s\n', char(header_W(i,:)));
        end
        fclose(fid);
        clear i fid
        dlmwrite([fout_folder, fout_viewW_asc], viewW, 'delimiter', ' ', '-append')
    end
    if rs>0
        dummystr=[size(satBBlong,2), size(satBBlat,1), min(min(satBBlong))-0.5*satMETA.cellsize, ...
            min(min(satBBlat))-0.5*satMETA.cellsize, satMETA.cellsize, headerv_W(6,1)]';
        for i=1:6
            satHeader{i,1}=[headern_W{i,1}, ' ', num2str(dummystr(i,:))];
        end
        clear i dummystr
%           satellite image (bounding box)
%               NDSI
        if rs==1
            satBBndsi(isnan(satBBndsi))=headerv_W(6,1);
            fid=fopen([fout_folder, fout_ndsi_asc], 'w');
            if fid==-1
                print=('NDSI (bounding box) output file could not be created.');
            end
            for i=1:6
                fprintf(fid, '%s\n', char(satHeader(i,:)));
            end
            fclose(fid);
            clear i fid
            dlmwrite([fout_folder, fout_ndsi_asc], satBBndsi, 'delimiter', ' ', '-append')
        end
%               map
        satBBmap(isnan(satBBmap))=headerv_W(6,1);
        fid=fopen([fout_folder, fout_satmap_asc], 'w');
        if fid==-1
            print=('Classification (bounding box) output file could not be created.');
        end
        for i=1:6
            fprintf(fid, '%s\n', char(satHeader(i,:)));
        end
        fclose(fid);
        clear i fid
        dlmwrite([fout_folder, fout_satmap_asc], satBBmap, 'delimiter', ...
         ' ', '-append')
    end
    disp('Load workspace using the command:')
    disp(['load(''', fout_folder, fout_mat, ''')'])
%       Matlab/Octave figures
%           photograph
    if os>1 || (os==1 && photoloop==1)
%               GCPs (and DDS optimisation)
        if IsOctave_PRACTISE()
            hgsave(figure(1), [fout_folder, fout_figname1, '.ofig']);
              %  Load in Octave and
              %  get graphics, axes, object and legend handles for 'set'-options
              %   hgload([fout_folder, fout_figname1,'.ofig']);
              %   fig1_f=get(gcf);
              %   fig1_a=get(gca);
              %   fig1_o=get(gco);
              %   fig1_l=get(gcf,'children')
              %      to test if fig1_l(1),fig1_l(2),...:
              %   set(fig1_l(1),'color',[0 0 0]); set(fig1_l(1),'color',[1 1 1])
            %dummy_pos=get(figure(1),'position');
            %print (figure(1), [fout_folder, fout_figname1, '.png'], ...
            %  ['-S', num2str(dummy_pos(3)), ',', num2str(dummy_pos(4))])
            %clear dummy_pos
        else
            saveas(figure(1), [fout_folder, fout_figname1], 'fig');
            %print (figure(1), '-dtiff', [fout_folder, fout_figname1, '.tif'])
        end
    end
    if cs>0
%               DN frequency histogram
%        if IsOctave_PRACTISE()
%            hgsave(figure(2), [fout_folder, fout_figname2, '.ofig']);
%        else
%            saveas(figure(2), [fout_folder, fout_figname2], 'fig');
%        end
    end
%               classification
    if IsOctave_PRACTISE()
        hgsave(figure(4), [fout_folder, fout_figname4, '.ofig']);
    else
        saveas(figure(4), [fout_folder, fout_figname4], 'fig');
    end
%          satellite image and photo map
    if rs>0
        if IsOctave_PRACTISE()
            hgsave(figure(5), [fout_folder, fout_figname5, '.ofig']);
        else
            saveas(figure(5), [fout_folder, fout_figname5], 'fig');
        end
    end
%   save variables
%       tidy up
    clear ans c classW demWrcz dummy f_name fout_figname1 fout_figname2 ...
     fout_figname4 fout_figname5 headern_W l r
    if ~exist('fout_mat', 'var')
        error('Matlab output file name is not defined.')
    end
    save([fout_folder, fout_mat], 'cam', 'cam_rc', 'cam_off', ...
        'cam_rol', 'cam_foc', 'header_W', 'headerv_W', 'cam_hei', ...
        'cam_wid', 'pix_r', 'pix_c', 'X0', 'fin_demW', 'demW', ...
        'fin_folder', 'fin_imagepath', 'fin_image', ...
        'rgbimage', 'fout_folder', 'fout_classW_asc', 'fout_mat', ...
        'vs', 'os', 'cs', 'rs', 'rsms', 'rsps', 'is');
    if vs==0
        save([fout_folder, fout_mat], 'fin_viewWpath', 'fin_viewW', '-append');
    elseif vs==1
        save([fout_folder, fout_mat], 'fout_viewW_asc', 'buffer_radius', ...
            'viewAzi', 'viewAziSec', 'azi_sec_name', 'Vert', '-append');
    end
    if os>0
        save([fout_folder, fout_mat], 'gcpW_name', 'gcpW', 'gcpP', ...
            'gcpRMSE_orig', 'fin_gcpWpath', 'fin_gcpW', ...
            '-append');
        if os>1
            save([fout_folder, fout_mat], 'UB_gcp', 'LB_gcp', ...
                'DDS_MaxEval_os', 'DDS_R_os', 'Xopt', 'gcpRMSE_opt', ...
                'cam_orig', 'cam_rc_orig', 'gcpP_orig', '-append');
            if any(UB_gcp_orig~=UB_gcp) || any(LB_gcp_orig~=LB_gcp)
                save([fout_folder, fout_mat], 'UB_gcp_orig', ...
                 'LB_gcp_orig', '-append');
            end
            if gcpRMSE_orig~=gcpRMSE_0
                save([fout_folder, fout_mat], 'gcpRMSE_0', 'gcpP_0', ...
                 'X_0', 'cam_0', 'cam_rc_0', '-append');
            end
            if os==3
               save([fout_folder, fout_mat], 'gcpRMSE_optthres', ...
                'gcpRMSE_countthres', '-append');
            end
        end
    end
    if cs==0
        save([fout_folder, fout_mat], 'thres_rgb', 'delta_rgb', '-append');
    elseif cs>0
        save([fout_folder, fout_mat], 'thres_b_orig', 'movavgwindow', ...
         'thres_b_image1', 'rgbhist1', 'rgbhistmean1', ...
         'isnow1', '-append');
        if cs==2
            save([fout_folder, fout_mat], 'thres_b_low', 'isnow2', ...
             'irock', 'i5050', 'rgb_pc', 'pca', '-append');
        end
    end
    if cs3==3
        save([fout_folder, fout_mat], 'cs3', '-append');
    end
    if rs>0
        save([fout_folder, fout_mat], 'NDSIthres0', ...
         'UB_NDSI', 'LB_NDSI', 'DDS_MaxEval_rs', 'DDS_R_rs', ...
         'NDSIthresOpt', 'NDSIthresOpt_detail', 'PhRs', ...
         'NDSIbct10', 'NDSIbct1Opt', 'NDSIOpt', ...
         'fin_satpath', 'fin_satfolder_image', ...
         'fout_satmap_asc', 'satBB', 'satBBext', ...
         'satBBlat', 'satBBlong', ...
         'satBBndsi', 'satMETA', 'satHeader', ...
         'satJDAY', 'time', '-append');
        if exist('satBBin_ext', 'var')
            save([fout_folder, fout_mat], 'satBBin_ext', '-append');
        end
        if rs==1
            save([fout_folder, fout_mat], 'fin_satname_meta', ...
             'fout_ndsi_asc', 'satBBg_dn', 'satBBg_refl', ...
             'satBBnir_dn', 'satBBnir_refl', 'satBBmir_dn', ...
             'satBBmir_refl', 'satMTL', '-append');
            if exist('satBBlook', 'var')
                save([fout_folder, fout_mat], 'fin_satfolder_look', ...
                 'fin_satname_look', '-append');
            end
            if strcmp(fin_satname_meta(1:3), 'LT5') | strcmp(fin_satname_meta(1:3), 'LE7')
                save([fout_folder, fout_mat], 'fin_satname_DistEtoS', ...
                 'fin_satname_RadSol', 'satDistEtoS', 'satRadSol', ...
                 'satDOY', 'satBBg_rad', 'satBBmir_rad', 'satBBnir_rad', ...
                 '-append');
            end
        elseif rs==2
            save([fout_folder, fout_mat], 'fin_satname_ndsi', '-append');
        end
        if rsms>0
            save([fout_folder, fout_mat], 'satmask_code', 'satmaskEXT', ...
             'header_satmask', 'headerv_satmask', 'fin_satfolder_mask', ...
             'fin_satname_mask', 'satmaskBB', '-append');
        end
    end
    if N_images>1
        save([fout_folder, fout_mat], 'fin_images', '-append');
    end
    toc
end
disp('PRACTISE was executed successfully.')
diary off
clear all
                                                                            