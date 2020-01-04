%% Assignment 2 - (AE4423-19 Airline Planning and Optimization)
% Authors (Group 31):  Adrián Chaves García-Mascaraque (5077982)
%                      Guillermo Lara Juan (5169860)
%                      Miguel Ángel Sáez Ortuno (4541235)
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
h_div = 240;
%% Aircraft Data - Aircraft Characteristics

% Aircraft 1: Regional turboprop
% Aircraft 2: Regional jet
% Aircraft 3: Single aisle twin engine jet

Speed = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B2:D2'); % [km/h]
Seats = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B3:D3'); 
Ave_TAT = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B4:D4'); % Average TAT [min]
Max_Range = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B5:D5'); % Maximum Range [km]
Req_RWY = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B6:D6'); % Runway Required [m]
Lease_cost = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B7:D7'); % [€]
CX = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B8:D8'); % Fixed operating cost CX [€]
CT = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B9:D9'); % Time cost parameter CT [€/hr]
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

% Demand - The matrix obtained "Demand_h" gives the demand in each flight
% depending on the hour. The final matrix has (20x20x240)
% Demand_h = zeros(N,N,h_div);
indx = 1;
% for i = 1:N
%     for j = 1:N
%         indx = 1;
%         for h = 1:h_div/10
%             for p = 1:10
%                 Demand_h(i,j,indx) = Demand_2014(i,j)*Hour_coeff(i,h);
%                 indx = indx + 1;
%             end
%         end
%     end
% end

for i = 1:N
    for j = 1:length(Hour_coeff)
        Demand_h(i,:,j) = Demand_2014(i,:)*Hour_coeff(i,j);
    end
end


%% Dynamic programming
time_stages = linspace(24-1/10,0,240);
% We have to compute the following for each aircraft, maybe we can put a "while" to do it until no aircraft in our fleet is left.
for aircraft_type = 1:length(Fleet) % To perform the dynamic programming for each aircraft type.
    possible_routes = 1; % This variable will be used to count the number of possible routes.
    for t = 1:240/10
        t_stage = time_stages(t);
        while t_stage ~= 0 % As we don't know when the aircraft route should finish, we will study all the possible routes finishing in each of the different time stages.
        % The assumption of the hub-and-spoke airline makes that the first
        % studied airport, in this case we are starting at the end of the
        % day so it should be the last one, has to be the hub.
        
        
        % Hacer un loop que siga haciendo lo de dynamic programming para
        % todos los tiempos y en cada tiempo para todos los aeropuertos
        
        % Hacer una matriz que vaya guardando el coste, la secuencia de los
        % aeropuertos, la secuencia del tiempo, la demanda.
        
        % Tal vez sea mejor hacer 
        end
    end
end
