

% https://in.mathworks.com/videos/plot-geographic-data-on-a-map-in-matlab-1545831202291.html
clear all, close all, clc
Lat = xlsread('AE4423_Datasheets.xls','Group 31','C6:V6'); % Latitude of each airport [degrees]
Lon = xlsread('AE4423_Datasheets.xls','Group 31','C7:V7'); % Longitude of each airport [degrees]

% plot(Lon, Lat,'r', 'LineWidth',2)
title('Results')
geoplot(Lat,Lon,'.','MarkerSize',10)
% geobasemap('usgsimageryonly')
% geobasemap('landcover')
% geobasemap('darkwater')
% geobasemap('grayland')
% geobasemap('bluegreen')
% geobasemap('grayterrain')
% geobasemap('colorterrain')
% geobasemap('none')