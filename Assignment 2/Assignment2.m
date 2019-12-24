%% Assignment 2 - (AE4423-19 Airline Planning and Optimization)
% Authors (Group 31):  Adri�n Chaves Garc�a-Mascaraque (5077982)
%                      Guillermo Lara Juan (5169860)
%                      Miguel �ngel S�ez Ortuno (4541235)
clc, clear all, close all

%% DATA

Demand_2014 = xlsread('AE4423_Datasheets.xls','Group 31','C13:V32');
Hour_coeff = xlsread('AE4423_Ass2_APO.xlsx','Hour Coefficients','D3:AA22');
f = 1.42; % [USD/gallon] 2019 Assuming the same value for the 2014

%% Airport Data

% Airport data
[~,Cities,~] = xlsread('AE4423_Ass2_APO.xlsx','Airport','B1:U1'); % Name of cities
[~,Airports,~] = xlsread('AE4423_Ass2_APO.xlsx','Airport','B2:U2'); % Name of airports (ICAO Code)
Lat = xlsread('AE4423_Ass2_APO.xlsx','Airport','B3:U3'); % Latitude of each airport [degrees]
Lon = xlsread('AE4423_Ass2_APO.xlsx','Airport','B4:U4'); % Longitude of each airport [degrees]
Runway_length = xlsread('AE4423_Ass2_APO.xlsx','Airport','B5:U5'); % Runway length [m]

% Hub
[~,Hub,~] = xlsread('AE4423_Datasheets.xls','Group 31','C2');
pos_hub = find(strcmp(Cities,Hub)); % This variable indicates the position of the hub in future vectors

N = length(Airports); % Number of airports
%% Aircraft Data - Aircraft Characteristics

% Aircraft 1: Regional turboprop
% Aircraft 2: Regional jet
% Aircraft 3: Single aisle twin engine jet

Speed = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B2:D2'); % [km/h]
Seats = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B3:D3'); 
Ave_TAT = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B4:D4'); % Average TAT [min]
Max_Range = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B5:D5'); % Maximum Range [km]
Req_RWY = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B6:D6'); % Runway Required [m]
Lease_cost = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B7:D7'); % [�]
CX = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B8:D8'); % Fixed operating cost CX [�]
CT = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B9:D9'); % Time cost parameter CT [�/hr]
CF = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B10:D10'); % Fuel cost parameter CF
Fleet = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B11:D11'); % Units

%% Pre-processing

% Distance
RE = 6371; % Radius of the Earth [km]
Distance = zeros(N);
Dsigma = zeros(N);

for i = 1:N
    for j = 1:N
        Dsigma(i,j) = 2*asin(sqrt((sind((Lat(i)-Lat(j))/2))^2+cosd(Lat(i))*...
            cosd(Lat(j))*(sind((Lon(i)-Lon(j))/2))^2));
        Distance(i,j) = RE*Dsigma(i,j); % [km]
    end
end

% To compute the Demand per hour, I think that the best idea would be to
% compute it when required, not the whole matrix because it would be huge
% (20x20x24 = 9600)