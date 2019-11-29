%% Assignment 1 - (AE4423-19 Airline Planning and Optimization)
clc, clear all, close all
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Community129\cplex\matlab\x64_win64')

%% Data

f = 1.42; % [USD/gallon] 2019
LF = 0.75; % Load Factor

%% Data from excel

% Population: It has to be multiplied by 1000, because those values are 
% given per 1000 inhabitants
pop_2014 = xlsread('AE4423_Datasheets.xls','General','B4:B23')*1000; % Population 2014
pop_2019 = xlsread('AE4423_Datasheets.xls','General','C4:C23')*1000; % Population 2019

% Gross Domestic Product
GDP_2014 = xlsread('AE4423_Datasheets.xls','General','F4:F23'); % GDP 2014 [USD]
GDP_2019 = xlsread('AE4423_Datasheets.xls','General','G4:G23'); % GDP 2019 [USD]

% Airport data
[~,Cities,~] = xlsread('AE4423_Datasheets.xls','Group 31','C4:V4'); % Name of cities
[~,Airports,~] = xlsread('AE4423_Datasheets.xls','Group 31','C5:V5'); % Name of airports
Lat = xlsread('AE4423_Datasheets.xls','Group 31','C6:V6'); % Latitude of each airport [degrees]
Lon = xlsread('AE4423_Datasheets.xls','Group 31','C7:V7'); % Longitude of each airport [degrees]
Runway_length = xlsread('AE4423_Datasheets.xls','Group 31','C8:V8'); % Runway length []
Slots = xlsread('AE4423_Datasheets.xls','Group 31','C8:V8');

% Hub
[~,Hub,~] = xlsread('AE4423_Datasheets.xls','Group 31','C2');
pos_hub = find(strcmp(Cities,Hub)); % This variable indicates the position of the hub in future vectors

Demand_2014 = xlsread('AE4423_Datasheets.xls','Group 31','C13:V32');

N = length(Airports); % Number of airports

%% Aircraft Data

% Aircraft Characteristics

% Aircraft 1: Regional turboprop
% Aircraft 2: Regional jet
% Aircraft 3: Single aisle twin engine jet
% Aircraft 4: Twin aisle, twin engine jet

Speed = [550 820 850 870]; % [km/h]
Seats = [45 70 150 320]; 
Average_TAT = [25 35 45 60]/60; % Turn Around Time [mins --> h]
Max_range = [1500 3300 6300 12000]; % Maximum Range [km]
Rwy_req = [1400 1600 1800 2600]; % Runway required [m]

% Cost
Lease_cost = [15000 34000 80000 190000]; % Weekly lease cost [€]
CX =[300 600 1250 2000]; % Fixed operating cost CX [€]
CT = [750 775 1400 2800]; % Time cost parameter CT [€/hr]
CF = [1.0 2.0 3.75 9.0]; % Fuel cost parameter CF

AC_number = length(Speed); % Number of aircraft

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

% Demand 2019
% Beta = inv(X^T X) X^T y; y is the demand and X is the matrix defined by
% the coefficients of the gravity model

X = zeros(N^2,4);
y = zeros(N^2,1);
B = zeros(4,1);

% Be aware that ln(0) is -inf.
for i = 1:N
    for j = 1:N
        k = i+N*(j-1);
        X(k,1) = 1;
        X(k,2) = log(pop_2014(i)*pop_2014(j));
        X(k,3) = log(GDP_2014(i)*GDP_2014(j));
        if i ~= j
%             y(k) = 0;
%             X(k,4) = 0;
%         else
            y(k) = log(Demand_2014(i,j)); 
            X(k,4) = -log(f*Distance(i,j));
        end
    end
end

B = inv(X'*X)*X'*y; % Vector that contains the coefficients of the linear 
% function obtained when taking logarithms of the gravity model that 
% defines the demand. (k, b1, b2, b3)

Demand_2019 = zeros(N); % Demand 2019 matrix definition

for i = 1:N
    for j = 1:N
        if i ~= j
%             Demand_2019(i,j) = 0; % We do not need to put this because it is already zero
%         else
            logDemand_2019 = B(1)+B(2)*log(pop_2019(i)*pop_2019(j))+...
                B(3)*log(GDP_2019(i)*GDP_2019(j))-B(4)*log(f*Distance(i,j));
            Demand_2019(i,j) = exp(logDemand_2019);
        end
    end
end

% We need to have integer numbers in the demand so we can choose the way of
% rounding it: ceil, fix, floor, round
Demand_2019 = round(Demand_2019);

% Vector g
g = ones(1,N);
g(pos_hub) = 0;

%% Initialize model

model = 'Assignment1_group31'; % Name of the model
cplex = Cplex(model); % Initialize Cplex
cplex.Model.sense = 'maximize'; % Maximize profit

DV = N^2*2+N^2*AC_number+AC_number; % Decision variables 
% There are 2404 decision variables splitted in:
% - wij that indicate the flow between i and j passing through the hub. (400)
% - xij that indicate the direct flow between i and j (not going through the hub) (400)
% - zij^k that indicate the number flight performed by AC^k from i to k. (1600)
% - AC^k that indicate the number of aircraft of each type.

obj = zeros(1,DV); % Objective function definition
% The Decision Variables are organized:
% First -> wij 
% Second -> xij
% Third -> zij^k
% Fourth -> AC^k

%% Objective function

for i = 1:N
    for j = 1:N
        % Revenues
        if i == j 
            Yield = 0;
        else
            Yield = 5.9*Distance(i,j)^(-0.76)+0.43;
        end
        
        obj(varindex_3(1,i,j,k,N,AC_number)) = Yield*Distance(i,j);
        obj(varindex_3(2,i,j,k,N,AC_number)) = Yield*Distance(i,j);
        
        for k = 1:AC_number
        % Costs
            if i == j
                Fixed_costs = 0;
            else
                Fixed_costs = CX(k);
            end
            Time_costs = CT(k)*Distance(i,j)/Speed(k);
            Fuel_costs = CF(k)*f/1.5*Distance(i,j);
            Operating_cost = Fixed_costs+Time_costs+Fuel_costs;
            
            if i == pos_hub || j == pos_hub
                Operating_cost = 0.7*Operating_cost;
            end
            
            obj(varindex_3(3,i,j,k,N,AC_number)) = -Operating_cost;
            obj(varindex_3(4,i,j,k,N,AC_number)) = -Lease_cost(k);
        end
    end
end
% NaNvector = find(isnan(obj))--> Used to check if there were any cell with
% NaN value. Found when 0^(-0.76) when computing the Yield. I think that
% that problem is solved.

lb = zeros(DV,1);
ub = inf(DV,1);
obj = obj';
ctype1 = char(ones(1, (DV)) * ('I'));  %I=integer. Other options, C=continous, B=binary.
ctype = strcat(ctype1);

cplex.addCols(obj, [], lb, ub, ctype,NameDV); %This can also be done without NameDV

%% Constraint 1: Demand verification


for i=1:N
    for j=1:N
        C1 = zeros(1,DV);
        C1(varindex_3(1,i,j,k,N,AC_number)) = 1;
        C1(varindex_3(2,i,j,k,N,AC_number)) = 1;
        RHS = Demand(i,j);
        cplex.addRows(-inf, C1, Demand(i,j),sprintf('Constraint1%d_%d',i,j));
    end
end


%% Solve the model

cplex.writeModel([model '.lp']) %Store the model to an .lp file for debugging
cplex.Param.timelimit.Cur = 10; %Timelimit of the solver in seconds, more useful for larger models
cplex.solve();

% for i=1:Nodes
%     for j=1:Nodes
%         flow_this_leg = cplex.Solution.x(varindex_3(1,i,j));
%         if flow_this_leg>0
%             sprintf('%d passengers from %d to %d', flow_this_leg, i, j)
%         end
%     end
% end
            
