%% Assignment 1 - (AE4423-19 Airline Planning and Optimization)
clc, clear all, close all
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Community129\cplex\matlab\x64_win64')

%% Data from excel

% Population
pop_2010 = xlsread('AE4423_Datasheets.xls','General','B4:B23');
pop_2017 = xlsread('AE4423_Datasheets.xls','General','C4:C23');

% Gross Domestic Product
GDP_2010 = xlsread('AE4423_Datasheets.xls','General','F4:F23');
GDP_2017 = xlsread('AE4423_Datasheets.xls','General','G4:G23');

[~,Airports,~] = xlsread('AE4423_Datasheets.xls','Group 31','C5:V5');
Lat = xlsread('AE4423_Datasheets.xls','Group 31','C6:V6');
Lon = xlsread('AE4423_Datasheets.xls','Group 31','C7:V7');
Runway_length = xlsread('AE4423_Datasheets.xls','Group 31','C8:V8');
Slots = xlsread('AE4423_Datasheets.xls','Group 31','C8:V8');

Demand_2014 = xlsread('AE4423_Datasheets.xls','Group 31','C13:V32');

N = length(Airports); % Number of airports

%% Pre-processing

% Population --> Linear regression y=m*Year+n
m_pop = zeros(1,N);
n_pop = zeros(1,N);

for i = 1:N
    m_pop(i) = (pop_2017(i)-pop_2010(i))/(2017-2010);
    n_pop(i) = pop_2017(i)-m_pop(i)*2017;
end

pop_2014 = m_pop*2014+n_pop;
pop_2019 = m_pop*2019+n_pop;

% Gross Domestic Product --> Linear regression y=m*Year+n
m_GDP = zeros(1,N);
n_GDP = zeros(1,N);

for i = 1:N
    m_GDP(i) = (GDP_2017(i)-GDP_2010(i))/(2017-2010);
    n_GDP(i) = GDP_2017(i)-m_GDP(i)*2017;
end

GDP_2014 = m_GDP*2014+n_GDP;
GDP_2019 = m_GDP*2019+n_GDP;

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

% Ordinary Least Squares
logDemand_2014 = log(Demand_2014);
logpop_2014 = log(pop_2014);
logGDP_2014 = log(GDP_2014);
logDistance = log(Distance);



%% Aircraft Data

% Aircraft Characteristics

% Aircraft 1: Regional turboprop
% Aircraft 2: Regional jet
% Aircraft 3: Single aisle twin engine jet
% Aircraft 4: Twin aisle, twin engine jet

Speed = [550 820 850 870]; % [km/h]
Seats = [45 70 150 320]; 
Average_TAT = [25 35 45 60]; % [mins]
Max_range = [1500 3300 6300 12000]; % Maximum Range [km]
Rwy_req = [1400 1600 1800 2600]; % Runway required [m]

% Cost
Lease_cost = [15000 34000 80000 190000]; % Weekly lease cost [€]
CX =[300 600 1250 2000]; % Fixed operating cost CX [€]
CT = [750 775 1400 2800]; % Time cost parameter CT [€/hr]
CF = [1.0 2.0 3.75 9.0]; % Fuel cost parameter CF

%% Initialize model

model = 'Assignment1_group31'; % Name of the model
cplex = Cplex(model); % Initialize Cplex
cplex.Model.sense = 'maximize'; % Maximize profit


