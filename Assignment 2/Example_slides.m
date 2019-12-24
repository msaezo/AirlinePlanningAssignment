%% Assignment 2 - (AE4423-19 Airline Planning and Optimization)
% Authors (Group 31):  Adrián Chaves García-Mascaraque (5077982)
%                      Guillermo Lara Juan (5169860)
%                      Miguel Ángel Sáez Ortuno (4541235)
clc, clear all, close all
%% TRIAL VERSION %%
% This code is used to check the results obtained for the exercise in class

%% DATA 
Demand = xlsread('Example_slides.xlsx','Example','B12:H18');
Hour_coeff = xlsread('Example_slides.xlsx','Example','B3:Y9');
f = 1.42; % [USD/gallon] 2019 Assuming the same value for the 2014

%% Airport Data
Lat = xlsread('Example_slides.xlsx','Example','B24:H24'); % Latitude of each airport [degrees]
Lon = xlsread('Example_slides.xlsx','Example','B25:H25'); % Longitude of each airport [degrees]
Runway_length = xlsread('Example_slides.xlsx','Example','B26:H26'); % Runway length [m]

%% Aircraft Data - Aircraft Characteristics
Speed = [810 870]; % [km/h]
Seats = [50 120]; 
Ave_TAT = [30 35]/60; % Average TAT [min --> h]
Max_Range = [3800 4500]; % Maximum Range [km]
Req_RWY = [1500 1600]; % Runway Required [m]
Lease_cost = [4540 7340]; % [€]
CX = [620 920]; % Fixed operating cost CX [€]
CT = [710 1000]; % Time cost parameter CT [€/hr]
CF = [2.4 2.9]; % Fuel cost parameter CF
Fleet = [3 2]; % Units

% Hub
pos_hub = 1; % This variable indicates the position of the hub in future vectors (Amsterdam - EHAM)

N = length(Lat); % Number of airports
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

% Costs (Matrix 20x20 per aircraft) (Smaller in this case, less airports)
for k = 1:length(Fleet)
    for i = 1:N
        for j = 1:N
%             if i == j
%                 Fixed_costs = 0;
%             else
%                 Fixed_costs = CX(k);
%             end
                Time_costs(i,j,k) = CT(k)*Distance(i,j)/Speed(k);
                Fuel_costs(i,j,k) = CF(k)*f/1.5*Distance(i,j);
%                 Operating_cost = Fixed_costs+Time_costs+Fuel_costs;
            
%             if i == pos_hub || j == pos_hub
%                 Operating_cost = 0.7*Operating_cost;
%             end

        end
    end
end

% Demand per hour (3D Matrix --> 20x20x24) (Not in this example)
Dem_ph = zeros(N,N,24);
for i = 1:N % Origin
    for j = 1:N % Destination
        Dem_ph(i,j,:) = Demand(i,j)*Hour_coeff(i,:);
    end
end

% Revenue
Revenues = zeros(N,N,24);
for i = 1:N
    for j = 1:N
        if i == j 
            Yield = 0;
        else
            Yield = 5.9*Distance(i,j)^(-0.76)+0.043;
        end
        Revenues(i,j,:) = Yield*Distance(i,j)*Dem_ph(i,j,:);
    end
end

% Operational Constraints
% - Range

% - Runway



time_step = 6/60; % Time step [min --> h]
time_stages = 0:time_step:24-time_step;