%% Assignment 1 - (AE4423-19 Airline Planning and Optimization)
clc, clear all, close all
% addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Community129\cplex\matlab\x64_win64')
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab\x64_win64')

%% Data

f = 1.42; % [USD/gallon] 2019 Assuming the same value for the 2014
LF = 0.75; % Load Factor
BT = 10*7; % [h/week]

%% Data from excel

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
            
%             Demand_2014_trial_2(i,j) = exp(B(1))*(((pop_2014(i)*pop_2014(j))^B(2)*(GDP_2014(i)*GDP_2014(j))^B(3))/((f*Distance(i,j))^B(4)));
        end
    end
end

Demand_2014_trial = round(Demand_2014_trial);
% Demand_2014_trial_2 = round(Demand_2014_trial_2);
% Vector g
g = ones(1,N);
g(pos_hub) = 0;


% Demand_2019=Demand_2014;
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

%% Constraint 1: Demand verification
C1_matrix = zeros(N^2,DV);
row_c1 = 1;
for i = 1:N
    for j = 1:N
        C1 = zeros(1,DV);
        C1(varindex_3(1,i,j,k,N,AC_number)) = 1;
        C1(varindex_3(2,i,j,k,N,AC_number)) = 1;
        RHS = Demand_2019(i,j);
        cplex.addRows(-inf,C1,RHS,sprintf('Constraint1%d_%d',i,j));
        C1_matrix(row_c1,:) = C1; 
        row_c1 = row_c1 + 1;
    end
end

%% Constraint 2: Making sure that w is 0 if the origin or the destination 
% is the hub
C2_matrix = zeros(N^2,DV);
row_c2 = 1;
for i = 1:N
    for j = 1:N
        C2 = zeros(1,DV);
        C2(varindex_3(1,i,j,k,N,AC_number)) = 1;
        RHS = Demand_2019(i,j)*g(i)*g(j);
        cplex.addRows(-inf,C2, RHS,sprintf('Constraint2%d_%d',i,j));
        C2_matrix(row_c2,:) = C2; 
        row_c2 = row_c2 + 1;
    end
end

%% Constraint 3: Capacity

C3_matrix = zeros(N^2,DV);
row_c3 = 1;
for i = 1:N
    for j = 1:N
        C3 = zeros(1,DV);
        C3(varindex_3(2,i,j,k,N,AC_number)) = 1;
        
        for m = 1:N
            C3(varindex_3(1,i,m,k,N,AC_number)) = C3(varindex_3(1,i,m,k,N,AC_number))+(1-g(j));
            C3(varindex_3(1,m,j,k,N,AC_number)) = C3(varindex_3(1,m,j,k,N,AC_number))+(1-g(i));
        end
        
        
        for k = 1:AC_number
            C3(varindex_3(3,i,j,k,N,AC_number)) = -Seats(k)*LF;
        end
        C3_matrix(row_c3,:) = C3; 
        RHS = 0;
        cplex.addRows(-inf, C3, RHS,sprintf('Constraint3%d_%d',i,j));
        row_c3 = row_c3 + 1;
    end
end

%% Constraint 4: Continuity

C4_matrix = zeros(N*AC_number,DV);
row_c4 = 1;
for k = 1:AC_number
    for i = 1:N
        C4 = zeros(1,DV);
        for j = 1:N
            if i ~= j
                C4(varindex_3(3,i,j,k,N,AC_number)) = 1;
                C4(varindex_3(3,j,i,k,N,AC_number)) = -1;
            end
        end
        C4_matrix(row_c4,:) = C4; 
        cplex.addRows(0,C4,0,sprintf('Constraint4%d_%d',k,i));
        row_c4 = row_c4 + 1;
    end
end

%% Constraint 5: Aircraft productivity

C5_matrix = zeros(AC_number,DV);
row_c5 = 1;
for k = 1:AC_number
    C5 = zeros(1,DV);
    for i = 1:N
        for j = 1:N % if i=j the distance is 0 so the first variable should be zero? 
            if i ~= j
            C5(varindex_3(3,i,j,k,N,AC_number)) = Distance(i,j)/Speed(k) + TAT(k)*(1+(1-g(j)));
            end
            C5(varindex_3(4,i,j,k,N,AC_number)) = -BT;
        end
    end
    C5_matrix(row_c5,:) = C5;
    cplex.addRows(-inf,C5,0,sprintf('Constraint5%d',k));
    row_c5 = row_c5 + 1;
end

%% Constraint 6: Range verification
C6_matrix = zeros(N^2*AC_number,DV);
row_c6 = 1;
for k = 1:AC_number 
    for i = 1:N
        for j = 1:N
            C6 = zeros(1,DV);
            C6(varindex_3(3,i,j,k,N,AC_number)) = 1;
            RHS = 0;
            if Distance(i,j) <= Max_range(k)
                RHS = 1000000;                
            end
            cplex.addRows(-inf, C6, RHS,sprintf('Constraint6%d_%d_%d',k,i,j));
            C6_matrix(row_c6,:) = C6;
            row_c6 = row_c6 + 1;
        end
    end
end

%% Constraint 7: Runway verification
C7_matrix = zeros(N^2*AC_number,DV);
row_c7 = 1;
for k = 1:AC_number 
    for i = 1:N
        for j = 1:N
            C7 = zeros(1,DV);
            C7(varindex_3(3,i,j,k,N,AC_number)) = 1;
            RHS = 0;
            if Runway_length(i) >= Rwy_req(k) && Runway_length(j) >= Rwy_req(k)
                RHS = 1000000;                
            end
            cplex.addRows(-inf, C7, RHS,sprintf('Constraint7%d_%d_%d',k,i,j));
            C7_matrix(row_c7,:) = C7;
            row_c7 = row_c7 + 1;
        end
    end
end

% %% Constraint 8
% C8_matrix = zeros(1,DV);
% row_c8 = 1;
% C8 = zeros(1,DV);
% for k = 1:AC_number 
%     C8(varindex_3(4,i,j,k,N,AC_number)) = 1;
% end
% cplex.addRows(-inf, C8, inf,sprintf('Constraint8'));

%% Solve the model

cplex.writeModel([model '.lp']) %Store the model to an .lp file for debugging
cplex.Param.timelimit.Cur = 300; %Timelimit of the solver in seconds, more useful for larger models
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