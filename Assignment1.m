
%% Assignment 1 - (AE4423-19 Airline Planning and Optimization)
clc, clear all, close all
% addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Community129\cplex\matlab\x64_win64')
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab\x64_win64') % Guille
% addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Community129\cplex\matlab\x64_win64')
set(0,'defaulttextinterpreter','latex')

%% Data

f = 1.42; % [USD/gallon] 2019 Assuming the same value for the 2014
LF = 0.75; % Load Factor
BT = 10*7; % [h/week]

%% Data from excel 
%The data is extracted from Excel and it is written into variables. 

% Population: It has to be multiplied by 1000, because those values are 
% given per 1000 inhabitants
pop_2014 = xlsread('AE4423_Datasheets.xls','General','B4:B23');%/1000; % Population 2014
pop_2017 = xlsread('AE4423_Datasheets.xls','General','C4:C23');%/1000; % Population 2019

% Gross Domestic Product
GDP_2014 = xlsread('AE4423_Datasheets.xls','General','F4:F23'); % GDP 2014 [USD]
GDP_2017 = xlsread('AE4423_Datasheets.xls','General','G4:G23'); % GDP 2019 [USD]

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
TAT = [25 35 45 60]/60; % Turn Around Time [mins --> h]
Max_range = [1500 3300 6300 12000]; % Maximum Range [km]
Rwy_req = [1400 1600 1800 2600]; % Runway required [m]

% Cost
Lease_cost = [15000 34000 80000 190000]; % Weekly lease cost [€]
CX =[300 600 1250 2000]; % Fixed operating cost CX [€]
CT = [750 775 1400 2800]; % Time cost parameter CT [€/hr]
CF = [1.0 2.0 3.75 9.0]; % Fuel cost parameter CF

AC_number = length(Speed); % Number of aircraft

%% Pre-processing

% Population --> Linear regression y=m*Year+n
m_pop = zeros(1,N);
n_pop = zeros(1,N);

for i = 1:N
    m_pop(i) = (pop_2017(i)-pop_2014(i))/(2017-2014);
    n_pop(i) = pop_2017(i)-m_pop(i)*2017;
end

pop_2019 = m_pop*2019+n_pop;

% Gross Domestic Product --> Linear regression y=m*Year+n
m_GDP = zeros(1,N);
n_GDP = zeros(1,N);

for i = 1:N
    m_GDP(i) = (GDP_2017(i)-GDP_2014(i))/(2017-2014);
    n_GDP(i) = GDP_2017(i)-m_GDP(i)*2017;
end

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

% Demand 2019
% Beta = inv(X^T X) X^T y; y is the demand and X is the matrix defined by
% the coefficients of the gravity model

X = zeros(N^2-N,4);
y = zeros(N^2-N,1);
B = zeros(4,1);
c=1;
% Be aware that ln(0) is -inf.
for i = 1:N
    for j = 1:N
        if i ~= j
%             k = i+N*(j-1)
            X(c,1) = 1;
            X(c,2) = log(pop_2014(i))+log(pop_2014(j));
            X(c,3) = log(GDP_2014(i))+log(GDP_2014(j));
            X(c,4) = -log(f*Distance(i,j));
            y(c) = log(Demand_2014(i,j)); 
            c=c+1;
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



Demand_2014_trial = zeros(N); % Demand 2014 matrix definition
% Demand_2014_trial_2 = zeros(N); % Demand 2014 matrix definition
for i = 1:N
    for j = 1:N
        if i ~= j
%             Demand_2019(i,j) = 0; % We do not need to put this because it is already zero
%         else
            logDemand_2014_trial = B(1)+B(2)*log(pop_2014(i)*pop_2014(j))+...
                B(3)*log(GDP_2014(i)*GDP_2014(j))-B(4)*log(f*Distance(i,j));
            Demand_2014_trial(i,j) = exp(logDemand_2014_trial);
        end
    end
end

Demand_2014_trial = round(Demand_2014_trial);

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
% - zij^k that indicate the number fl ight performed by AC^k from i to k. (1600)
% - AC^k that indicate the number of aircraft of each type.

obj = zeros(1,DV); % Objective function definition
% The Decision Variables are organized:
% First -> wij 
% Second -> xij
% Third -> zij^k
% Fourth -> AC^k

%% Objective function
k=1;
for i = 1:N
    for j = 1:N
        % Revenues
        if i == j 
            Yield = 0;
        else
            Yield = 5.9*Distance(i,j)^(-0.76)+0.043;
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

lb = zeros(DV,1);
ub = inf(DV,1);
obj = obj';
ctype1 = char(ones(1, (DV)) * ('I'));  %I=integer. Other options, C=continous, B=binary.
ctype = strcat(ctype1);

cplex.addCols(obj, [], lb, ub, ctype); %This can also be done without NameDV

%% Constraints
Constraints(cplex,N,DV,AC_number,Demand_2019,g,Seats,LF,Speed,Distance,TAT,BT,Max_range,Runway_length,Rwy_req)

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

obj_matrix=cplex.Model.obj;
A_total=cplex.Model.A;
RHS_total=cplex.Model.rhs;
solution=cplex.Solution.x;
%% Post-Processing
%In this section the KPI performance parameters will be calculated. 
ASK = 0;

Operating_cost_total = 0;
Revenue_tot = 0;
for i = 1:N
    for j = 1:N
        Lease_cost_tot = 0;
        for k = 1:AC_number
            %CASK
            z_sol = (cplex.Solution.x(varindex_3(3,i,j,k,N,AC_number)));
            ASK = ASK + z_sol*Distance(i,j)*Seats(k);
            Lease_cost_tot = Lease_cost_tot + Lease_cost(k)*cplex.Solution.x(varindex_3(4,i,j,k,N,AC_number));
            
            if i == j
                Fixed_costs = 0;
            else
                Fixed_costs = CX(k);
            end
            Time_costs = CT(k)*Distance(i,j)/Speed(k);
            Fuel_costs = CF(k)*f/1.5*Distance(i,j);
            Operating_cost = (Fixed_costs+Time_costs+Fuel_costs)*cplex.Solution.x(varindex_3(3,i,j,k,N,AC_number));
            
            if i == pos_hub || j == pos_hub
                Operating_cost = 0.7*Operating_cost;
            end

            Operating_cost_total = Operating_cost_total + Operating_cost;
            
          
        end
        %Revenue Calculation
        if i == j 
            Yield = 0;
        else
            Yield = 5.9*Distance(i,j)^(-0.76)+0.043;
        end
        Revenue_tot = Revenue_tot + Yield*Distance(i,j)*(cplex.Solution.x(varindex_3(1,i,j,k,N,AC_number))+cplex.Solution.x(varindex_3(2,i,j,k,N,AC_number)));
        
    end
end

%RPK
RPK_x = 0;
RPK_w = 0;
for i = 1:N
    for j = 1:N
        RPK_x = RPK_x + Distance(i,j)*cplex.Solution.x(varindex_3(2,i,j,k,N,AC_number));
        RPK_w = RPK_w + (Distance(i,pos_hub)+Distance(pos_hub,j))*cplex.Solution.x(varindex_3(1,i,j,k,N,AC_number));
    end
end

Cost_tot = Lease_cost_tot + Operating_cost_total;
CASK = Cost_tot/ASK;
RASK = Revenue_tot/ASK;
RPK = RPK_x + RPK_w;
LoadFactor = RPK/ASK;
Yield_tot = (Revenue_tot/RPK);
BELF = CASK/Yield_tot;

%% Graphs

% https://in.mathworks.com/videos/plot-geographic-data-on-a-map-in-matlab-1545831202291.html
% COLOUR SPECIFICATION

blue = [0, 0.4470, 0.7410];
orange = [0.8500, 0.3250, 0.0980];
yellow = [0.9290, 0.6940, 0.1250];
purple = [0.4940, 0.1840, 0.5560];
green = [0.4660, 0.6740, 0.1880];
cyan = [0.3010, 0.7450, 0.9330];
red = [0.6350, 0.0780, 0.1840];

Pax = zeros(N);
Direct_pax = zeros(N);
Nodirect_pax = zeros(N);
Tot_ndp = 0;
Tot_dp = 0;
flighs_ac1 = zeros(N);
flighs_ac2 = zeros(N);
flighs_ac3 = zeros(N);
flighs_ac4 = zeros(N);

for i = 1:N
    for j = 1:N
        Pax(i,j) = Pax(i,j)+solution(varindex_3(2,i,j,k,N,AC_number));
        Direct_pax(i,j) = solution(varindex_3(2,i,j,k,N,AC_number));
        Tot_dp = Tot_dp + Direct_pax(i,j);
        
        Pax_through_hub = solution(varindex_3(1,i,j,k,N,AC_number));
        Tot_ndp = Tot_ndp + Pax_through_hub;
        
        Pax(i,pos_hub) = Pax(i,pos_hub)+Pax_through_hub;
        Pax(pos_hub,j) = Pax(pos_hub,j)+Pax_through_hub;
        
        Nodirect_pax(i,pos_hub) = Nodirect_pax(i,pos_hub)+Pax_through_hub;
        Nodirect_pax(pos_hub,j) = Nodirect_pax(pos_hub,j)+Pax_through_hub;
        
        flighs_ac1(i,j) = flighs_ac1(i,j)+solution(varindex_3(3,i,j,1,N,AC_number));
        flighs_ac2(i,j) = flighs_ac2(i,j)+solution(varindex_3(3,i,j,2,N,AC_number));
        flighs_ac3(i,j) = flighs_ac3(i,j)+solution(varindex_3(3,i,j,3,N,AC_number));
        flighs_ac4(i,j) = flighs_ac4(i,j)+solution(varindex_3(3,i,j,4,N,AC_number));
        
    end
end

Pax_norm = Pax/max(max(Pax));
Nodirect_pax_norm = Nodirect_pax/max(max(Pax));
Direct_pax_norm = Direct_pax/max(max(Pax));

Tot_dp = sum(sum(Direct_pax));
Tot_ndp = sum(sum(Nodirect_pax));
Tot_pax = sum(sum(Pax));

Tot_fac1 = sum(sum(flighs_ac1));
Tot_fac2 = sum(sum(flighs_ac2));
Tot_fac3 = sum(sum(flighs_ac3));
Tot_fac4 = sum(sum(flighs_ac4));
Tot_f = Tot_fac1 + Tot_fac2 + Tot_fac3 + Tot_fac4;

% PAX = [Direct_pax ]

figure()
pie([Tot_fac1/Tot_f Tot_fac2/Tot_f Tot_fac3/Tot_f Tot_fac4/Tot_f]);
legend('Regional turboprop','Regional jet','Single aisle twin engine jet',...
'Twin aisle, twin engine jet','Location','bestoutside','Orientation',...
'vertical')
title('Flights flown by each type of AC')

figure()
pie([Tot_dp/Tot_pax Tot_ndp/Tot_pax]);
legend('Direct','No direct','Location','bestoutside','Orientation',...
'horizontal')
title('Type of passengers')

figure()
% geoplot([Lat(1) Lat(2)],[Lon(1) Lon(2)],'Linewidth',2)
geoplot(Lat,Lon,'.','MarkerSize',20)
hold on
for i = 1:N
    for j = 1:N
        if Pax_norm(i,j) > 0
            geoplot([Lat(i) Lat(j)],[Lon(i) Lon(j)],'k','Linewidth',5*Pax_norm(i,j))
            
        end
    end
%     text(Lat(i),Lon(i),Airports(i))
end
geobasemap('bluegreen')
% geobasemap('usgsimageryonly')
% geobasemap('landcover')
% geobasemap('darkwater')
% geobasemap('grayland')
% geobasemap('grayterrain')
% geobasemap('colorterrain')
% geobasemap('none')
title('Flow of passengers')


figure()
geoplot([Lat(1) Lat(2)],[Lon(1) Lon(2)],'color',red,'Linewidth',15*flighs_ac1(1,2)/Tot_f)            
hold on
geoplot([Lat(1) Lat(2)],[Lon(1) Lon(2)],'--','color',blue,'Linewidth',15*flighs_ac1(1,2)/Tot_f)
geoplot([Lat(1) Lat(2)],[Lon(1) Lon(2)],':','color',green,'Linewidth',15*flighs_ac1(1,2)/Tot_f)
geoplot([Lat(1) Lat(2)],[Lon(1) Lon(2)],'-.','color',yellow,'Linewidth',15*flighs_ac1(1,2)/Tot_f)
geoplot(Lat,Lon,'b.','MarkerSize',20)
for i = 1:N
    for j = 1:N
        if flighs_ac1(i,j) > 0
            geoplot([Lat(i) Lat(j)],[Lon(i) Lon(j)],'color',red,'Linewidth',4)            
        end
        if flighs_ac2(i,j) > 0
            geoplot([Lat(i) Lat(j)],[Lon(i) Lon(j)],'--','color',blue,'Linewidth',3)
        end
        if flighs_ac3(i,j) > 0
            geoplot([Lat(i) Lat(j)],[Lon(i) Lon(j)],':','color',green,'Linewidth',2)            
        end
        if flighs_ac4(i,j) > 0
            geoplot([Lat(i) Lat(j)],[Lon(i) Lon(j)],'-.','color',yellow,'Linewidth',2)
        end
        
    end
%     text(Lat(i),Lon(i),Airports(i))
end
geobasemap('bluegreen')
title('Aircraft Routes')
legend('Regional turboprop','Regional jet','Single aisle twin engine jet',...
'Twin aisle, twin engine jet','Location','northeast','NumColumns',2)



figure()
% geoplot([Lat(1) Lat(2)],[Lon(1) Lon(2)],'Linewidth',2)
geoplot(Lat,Lon,'.','MarkerSize',20)
hold on
for i = 1:N
    for j = 1:N
        if Direct_pax_norm(i,j) > 0
            geoplot([Lat(i) Lat(j)],[Lon(i) Lon(j)],'r','Linewidth',5*Direct_pax_norm(i,j))            
        end
        if Nodirect_pax_norm(i,j) > 0
            geoplot([Lat(i) Lat(j)],[Lon(i) Lon(j)],'b--','Linewidth',5*Nodirect_pax_norm(i,j))
        end
    end
%     text(Lat(i),Lon(i),Airports(i))
end
geobasemap('bluegreen')
% geobasemap('usgsimageryonly')
% geobasemap('landcover')
% geobasemap('darkwater')
% geobasemap('grayland')
% geobasemap('grayterrain')
% geobasemap('colorterrain')
% geobasemap('none')
title('Flow of passengers')
legend('Airports','Direct pax','No direct pax')