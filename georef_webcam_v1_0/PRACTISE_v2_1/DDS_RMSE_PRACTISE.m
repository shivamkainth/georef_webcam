% Haerer, Bernhardt and Schulz (2016)
% "PRACTISE - Photo Rectification And ClassificaTIon SoftwarE (V.2.1)"
%
%   written by
%   Stefan Haerer (LMU Munich)
%   08/2012
%   contact: stefan.haerer@boku.ac.at
%       updated by Stefan Haerer (BOKU Vienna, 30/11/2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ObjVal,gcpP] = DDS_RMSE_PRACTISE(X0, dem, header, gcpW, ...
    pix_c, pix_r)
%   Name:       DDS_RMSE_PRACTISE
%   Purpose:    DDS calls RMSE as objective function to
%               1) call Proj_PRACTISE to calculate projected ground control 
%               points (GCP),
%               2) calculate the root mean square error (RMSE) for the
%               distance (euclidean length of pixel rows and columns to the 
%               center of the photograph) between the observed and
%               projected GCP
%   
%   Output:     ObjVal (calculated RMSE)
%   Input:      gcp_obs (=[rows; cols] of observed GCP)
%               gcp_proj (=[rows; cols]) of projected GCP)
%               pix_r (number of pixel rows of photograph)
%               pix_c (number of pixel columns of photograph)
%               dem (DEM raster from ASCII-file [cf. ESRI ArcGIS])
%               header (values od DEM ASCII-file header [cf. ESRI ArcGIS]:
%                 ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[gcpP]=Proj_PRACTISE(gcpW, X0, pix_c, pix_r, dem, header);
% observed and projected distance to the center of the photograph 
rdistW=gcpW(5,:)-pix_r/2;
cdistW=gcpW(4,:)-pix_c/2;
rdistP=gcpP(2,:)-pix_r/2;
cdistP=gcpP(1,:)-pix_c/2;
% RMSE
ObjVal=sqrt((sum((rdistW-rdistP).^2+(cdistW-cdistP).^2))/(length(rdistW)));