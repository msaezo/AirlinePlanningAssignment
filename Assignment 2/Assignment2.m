%% Assignment 2 - (AE4423-19 Airline Planning and Optimization)
% Authors (Group 31):  Adrián Chaves García-Mascaraque (5077982)
%                      Guillermo Lara Juan (5169860)
%                      Miguel Ángel Sáez Ortuno (4541235)
clc, clear all, close all

%% DATA

Demand_2014 = xlsread('AE4423_Datasheets.xls','Group 31','C13:V32');
Hour_coeff = xlsread('AE4423_Ass2_APO.xlsx','Hour Coefficients','D3:AA22');
f = 1.42; % [USD/gallon] 2019 Assuming the same value for the 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Aircraft Data - Aircraft Characteristics

% Aircraft 1: Regional turboprop
% Aircraft 2: Regional jet
% Aircraft 3: Single aisle twin engine jet

Speed = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B2:D2'); % [km/h]
Seats = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B3:D3'); % Capacity
Ave_TAT = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B4:D4'); % Average TAT [min]
Max_Range = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B5:D5'); % Maximum Range [km]
Req_RWY = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B6:D6'); % Runway Required [m]
Lease_cost = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B7:D7'); % [€] This is already given per day
CX = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B8:D8'); % Fixed operating cost CX [€] 
CT = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B9:D9'); % Time cost parameter CT [€/hr]
CF = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B10:D10'); % Fuel cost parameter CF
Fleet = xlsread('AE4423_Ass2_APO.xlsx','Fleet Type','B11:D11'); % Units

AC_types = length(Fleet); % This variable indicates the number of aircraft types

%% PRE-PROCESSING

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Passengers per hour 

% Demand - The matrix obtained "Demand_h" gives the demand in each flight
% depending on the hour. The final matrix has (20x20x240)
% Demand_h = zeros(N,N,h_div);
% indx = 1;
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

Demand_h = zeros(N,N,24);
for i = 1:N % Departure airport
    for t = 1:length(Hour_coeff)
%         Demand_h(Departure airport, Arrival airport, Time)
        Demand_h(i,:,t) = Demand_2014(i,:)*Hour_coeff(i,t);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Costs (The leasing cost will be considered when the profit is computed)

Cost = zeros(N,N,AC_types);
for i = 1:N
    for j = 1:N
        for k = 1:AC_types
            
            if i == j
            	Fixed_costs = 0;
            else
            	Fixed_costs = CX(k);
            end
                
            Time_costs = CT(k)*Distance(i,j)/Speed(k);
            Fuel_costs = CF(k)*f/1.5*Distance(i,j);
%           Operating cost is the sum of the fixed, time and fuel costs
            Operating_cost = Fixed_costs + Time_costs + Fuel_costs;
            
%           It is said in Assignment 1 that: 
%           Flights departing or arriving at your hub airport all operating
%           costs can be assumed to be 30% lower due to economies of scale.
            if i == pos_hub || j == pos_hub
            	Operating_cost = 0.7*Operating_cost;
            end
                Cost(i,j,k) = Operating_cost;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Available aiports for each aircraft depending on the range and the runway
% length.

% As the only routes that are going to be performed are the ones that go
% and come from the hub, the distances from the hub to other airports are
% the ones evaluated.

av_airports_ac1 = av_airport(Runway_length, Req_RWY(1), Max_Range(1), Distance, N, pos_hub);
av_airports_ac2 = av_airport(Runway_length, Req_RWY(2), Max_Range(2), Distance, N, pos_hub);
av_airports_ac3 = av_airport(Runway_length, Req_RWY(3), Max_Range(3), Distance, N, pos_hub);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_division = 6; % [min]
min_per_day = 24*60; % [min] Minutes in a day
n_stages = min_per_day/time_division; % Number of stages in which the state space is divided

time_stages = linspace(24-1/10,0,n_stages); % It starts in 23:54 and ends in 0
% time_stages = linspace(24,1/10,n_stages); % It starts in 24:00 and ends in 0:06  --> we dont need this vector

% Minimum and maximum posible time stages for each airport (It depends on
% the aircraft since the velocity and the TAT are different)
airport_lim_tstage = zeros(N,2,AC_types);
flight_time = zeros(N,AC_types); % Vector in which is saved the time that each aircraft spends flying from hub to any other aiport

for k = 1:AC_types
    flight_time(:,k) = ceil((15 + 15 + Distance(:,2)/Speed(k)*60 + Ave_TAT(k))/6);
    airport_lim_tstage(:,:,k) = [1+flight_time(:,k) 240-flight_time(:,k)];
end

flight_time(pos_hub,:)=0;

%% 
% The following two variables are the controllers that indicate if there 
% are aircraft available or if new routes are still profitable.

profitable = 1;
AC_available = 1;
num_route = 1;
ROUTE_SCHEDULE = zeros(15,10,sum(Fleet));
while profitable == 1 && AC_available == 1
    
%     Now, passengers per route and per hour are computed assumming that 
%     we can capture the demand estimated for hours t, t-1, t+1 and t-2.
%     Then, the revenues from those passengers can be computed:

    [Act_pax_ac1, Rev_ac1, Act_pos_pax_ac1] = act_demand_and_rev(Demand_h, Seats(1), Distance, N);
    [Act_pax_ac2, Rev_ac2, Act_pos_pax_ac2] = act_demand_and_rev(Demand_h, Seats(2), Distance, N);
    [Act_pax_ac3, Rev_ac3, Act_pos_pax_ac3] = act_demand_and_rev(Demand_h, Seats(3), Distance, N);

%     Dynamic programming

% Aircraft 1
if Fleet(1)>0

[Profit_ac1, Nodes_ac1, profit_check1, node_schedule1, t_stages_schedule1] = Dynamic_Programming(Rev_ac1, Cost(:,:,1), airport_lim_tstage(:,:,1), av_airports_ac1, flight_time(:,1), pos_hub, n_stages);
% - Compute actual profit
number_pax_taken_ac1 = zeros(length(node_schedule1),1);
for p = 1:length(node_schedule1)
    Dep_airport = node_schedule1(p,1);
    Arr_airport = node_schedule1(p,2);
    Dep_time = floor(t_stages_schedule1(p,1)/10)+1;
    number_pax_taken_ac1(p) = Act_pax_ac1(Dep_airport, Arr_airport, Dep_time);
end

pax_taken_ac1 = [node_schedule1 t_stages_schedule1(:,1) number_pax_taken_ac1]; % This is part of the SCHEDULE
org_pax_ac1 = org_pax_fun(Demand_h, pax_taken_ac1);

check_pax_ac1 = sum(sum(sum(org_pax_ac1>Demand_h)));
[Final_Demand_h_ac1, pax_taken_ac1] = new_demand(Demand_h, pax_taken_ac1, check_pax_ac1); % This will give the demand after substracting the transported passengers and the actrual route.

% Now, the profit and the schedule can be obtained
% SCHEDULE (Dep(1), Arr(2), Dep_time(3), Arr_time(4), Pax(5), LF(6), Distance(7), Rev(8), Cost(9), Profit(10))

SCHEDULE_ac1 = [pax_taken_ac1(:,1:3) zeros(length(pax_taken_ac1),1) pax_taken_ac1(:,end) pax_taken_ac1(:,end)/Seats(1)];
for r = 1:length(pax_taken_ac1)
    SCHEDULE_ac1(r,4) = t_stages_schedule1(r,2); % Arrival time (Given in time stages)
    SCHEDULE_ac1(r,7) = Distance(SCHEDULE_ac1(r,1),SCHEDULE_ac1(r,2)); % Distance
    SCHEDULE_ac1(r,8) = (5.9*SCHEDULE_ac1(r,7)^(-0.76)+0.043)*SCHEDULE_ac1(r,7)*SCHEDULE_ac1(r,5); % Rev
    SCHEDULE_ac1(r,9) = Cost(SCHEDULE_ac1(r,1),SCHEDULE_ac1(r,2),1); % Cost
    SCHEDULE_ac1(r,10) = SCHEDULE_ac1(r,8) - SCHEDULE_ac1(r,9); % Profit
end 
FINAL_PROFIT_ac1 = sum(SCHEDULE_ac1(:,10))- Lease_cost(1);
else
    FINAL_PROFIT_ac1 = -100000;
end


% Aircraft 2
if Fleet(2)>0

[Profit_ac2, Nodes_ac2, profit_check2, node_schedule2, t_stages_schedule2] = Dynamic_Programming(Rev_ac2, Cost(:,:,2), airport_lim_tstage(:,:,2), av_airports_ac2, flight_time(:,2), pos_hub, n_stages);
% - Compute actual profit
number_pax_taken_ac2 = zeros(length(node_schedule2),1);
for p = 1:length(node_schedule2)
    Dep_airport = node_schedule2(p,1);
    Arr_airport = node_schedule2(p,2);
    Dep_time = floor(t_stages_schedule2(p,1)/10)+1;
    number_pax_taken_ac2(p) = Act_pax_ac2(Dep_airport, Arr_airport, Dep_time);
end

pax_taken_ac2 = [node_schedule2 t_stages_schedule2(:,1) number_pax_taken_ac2]; % This is part of the SCHEDULE
org_pax_ac2 = org_pax_fun(Demand_h, pax_taken_ac2);

check_pax_ac2 = sum(sum(sum(org_pax_ac2>Demand_h)));
[Final_Demand_h_ac2, pax_taken_ac2] = new_demand(Demand_h, pax_taken_ac2, check_pax_ac2); % This will give the demand after substracting the transported passengers and the actrual route.

% Now, the profit and the schedule can be obtained
% SCHEDULE (Dep(1), Arr(2), Dep_time(3), Arr_time(4), Pax(5), LF(6), Distance(7), Rev(8), Cost(9), Profit(10))

SCHEDULE_ac2 = [pax_taken_ac2(:,1:3) zeros(length(pax_taken_ac2),1) pax_taken_ac2(:,end) pax_taken_ac2(:,end)/Seats(2)];
for r = 1:length(pax_taken_ac2)
    SCHEDULE_ac2(r,4) = t_stages_schedule2(r,2); % Arrival time (Given in time stages)
    SCHEDULE_ac2(r,7) = Distance(SCHEDULE_ac2(r,1),SCHEDULE_ac2(r,2)); % Distance
    SCHEDULE_ac2(r,8) = (5.9*SCHEDULE_ac2(r,7)^(-0.76)+0.043)*SCHEDULE_ac2(r,7)*SCHEDULE_ac2(r,5); % Rev
    SCHEDULE_ac2(r,9) = Cost(SCHEDULE_ac2(r,1),SCHEDULE_ac2(r,2),2); % Cost
    SCHEDULE_ac2(r,10) = SCHEDULE_ac2(r,8) - SCHEDULE_ac2(r,9); % Profit
end 
FINAL_PROFIT_ac2 = sum(SCHEDULE_ac2(:,10))- Lease_cost(2);
else
    FINAL_PROFIT_ac2 = -100000;
end


% Aircraft 3
if Fleet(3)>0

[Profit_ac3, Nodes_ac3, profit_check3, node_schedule3, t_stages_schedule3] = Dynamic_Programming(Rev_ac3, Cost(:,:,3), airport_lim_tstage(:,:,3), av_airports_ac3, flight_time(:,3), pos_hub, n_stages);
% - Compute actual profit
number_pax_taken_ac3 = zeros(length(node_schedule3),1);
for p = 1:length(node_schedule3)
    Dep_airport = node_schedule3(p,1);
    Arr_airport = node_schedule3(p,2);
    Dep_time = floor(t_stages_schedule3(p,1)/10)+1;
    number_pax_taken_ac3(p) = Act_pax_ac3(Dep_airport, Arr_airport, Dep_time);
end

pax_taken_ac3 = [node_schedule3 t_stages_schedule3(:,1) number_pax_taken_ac3]; % This is part of the SCHEDULE
org_pax_ac3 = org_pax_fun(Demand_h, pax_taken_ac3);

check_pax_ac3 = sum(sum(sum(org_pax_ac3>Demand_h)));
[Final_Demand_h_ac3, pax_taken_ac3] = new_demand(Demand_h, pax_taken_ac3, check_pax_ac3); % This will give the demand after substracting the transported passengers and the actrual route.

% Now, the profit and the schedule can be obtained
% SCHEDULE (Dep(1), Arr(2), Dep_time(3), Arr_time(4), Pax(5), LF(6), Distance(7), Rev(8), Cost(9), Profit(10))

SCHEDULE_ac3 = [pax_taken_ac3(:,1:3) zeros(length(pax_taken_ac3),1) pax_taken_ac3(:,end) pax_taken_ac3(:,end)/Seats(3)];
for r = 1:length(pax_taken_ac3)
    SCHEDULE_ac3(r,4) = t_stages_schedule3(r,2); % Arrival time (Given in time stages)
    SCHEDULE_ac3(r,7) = Distance(SCHEDULE_ac3(r,1),SCHEDULE_ac3(r,2)); % Distance
    SCHEDULE_ac3(r,8) = (5.9*SCHEDULE_ac3(r,7)^(-0.76)+0.043)*SCHEDULE_ac3(r,7)*SCHEDULE_ac3(r,5); % Rev
    SCHEDULE_ac3(r,9) = Cost(SCHEDULE_ac3(r,1),SCHEDULE_ac3(r,2),3); % Cost
    SCHEDULE_ac3(r,10) = SCHEDULE_ac3(r,8) - SCHEDULE_ac3(r,9); % Profit
end 
FINAL_PROFIT_ac3 = sum(SCHEDULE_ac3(:,10))- Lease_cost(3);
else
    FINAL_PROFIT_ac3 = -100000;
end

%   Choose the most profitable

PROFIT_per_aircraft = [FINAL_PROFIT_ac1 FINAL_PROFIT_ac2 FINAL_PROFIT_ac3];
max_PROFIT = max(PROFIT_per_aircraft);

if max_PROFIT > 0
ac_choosen = find(max_PROFIT==PROFIT_per_aircraft) ;
    if ac_choosen == 1
        Fleet(1) = Fleet(1)-1; % Remove one aircraft
        ROUTE_SCHEDULE(1:length(SCHEDULE_ac1(:,1)),:,num_route) = SCHEDULE_ac1; % The schedule is saved in the result matrix
        TOT_PROFIT(num_route) = FINAL_PROFIT_ac1;
        num_route = num_route + 1;
        Demand_h = Final_Demand_h_ac1; % The final demand is changed
        
    elseif ac_choosen == 2
        Fleet(2) = Fleet(2)-1; % Remove one aircraft
        ROUTE_SCHEDULE(1:length(SCHEDULE_ac2(:,1)),:,num_route) = SCHEDULE_ac2; % The schedule is saved in the result matrix
        TOT_PROFIT(num_route) = FINAL_PROFIT_ac2;
        num_route = num_route + 1;
        Demand_h = Final_Demand_h_ac2; % The final demand is changed

    elseif ac_choosen == 3
        Fleet(3) = Fleet(3)-1; % Remove one aircraft
        ROUTE_SCHEDULE(1:length(SCHEDULE_ac3(:,1)),:,num_route) = SCHEDULE_ac3; % The schedule is saved in the result matrix
        TOT_PROFIT(num_route) = FINAL_PROFIT_ac3;
        num_route = num_route + 1;
        Demand_h = Final_Demand_h_ac3; % The final demand is changed
    end
else
    profitable = 0;
end

%     if another something --> This constraint is also considered with the
%     profit.
%         AC_available = 0;
%     end

end


%% FUNCTIONS

% Function to compute the demand and the revenue.
function [X, Rev, Act_pos_pax] = act_demand_and_rev(Demand_h, Seats, Distance, N)
X = zeros(N,N,24); % X represents actual flow of passengers at those hours
Rev = zeros(N,N,24); % Rev represents actual revenues at those hours

	for t = 1:24
        if t == 1 
        	X(:,:,t) = Demand_h(:,:,t) + Demand_h(:,:,t+1);
        elseif t == 2
        	X(:,:,t) = Demand_h(:,:,t)+Demand_h(:,:,t+1)+Demand_h(:,:,t-1);
        elseif t == 24
        	X(:,:,t) = Demand_h(:,:,t)+Demand_h(:,:,t-1)+Demand_h(:,:,t-2);
        else
            X(:,:,t) = Demand_h(:,:,t)+Demand_h(:,:,t-1)+Demand_h(:,:,t-2)+Demand_h(:,:,t+1);
        end
        
        % Let's check if the actual demand is higher that the actual
        % capacity. If that is the case, the actual demand is changed by
        % the capacity of the aircraft.
        
        Act_pos_pax = X; % Actual possible passenges
        X(X>Seats) = Seats;
        
        % Having the actual capacity, the revenues can be computed.
        Rev(:,:,t) = (5.9*Distance.^(-0.76)+0.043).*Distance.*X(:,:,t);
    end
    
Rev(isnan(Rev)) = 0;
end

% Function used to eliminate does airports that cannot be operated due to
% the range and runway length constraints.
function available_airports = av_airport(RWY_length, Req_RWY, Max_Range, Distance, N, pos_hub)
p = 1;
for i = 1:N
    if Distance(i,pos_hub)<Max_Range && RWY_length(i) >= Req_RWY
        available_airports(p) = i;
        p = p + 1;
    end
end
end

% Function that is in charge of the dynamic programming
function [Profit, Nodes, profit_check, node_schedule_matrix, t_stages_schedule] = Dynamic_Programming(Rev, Cost, airport_lim_tstage, av_airport, flight_time, pos_hub, n_stages)

Profit = zeros(length(av_airport), n_stages);
Nodes = zeros(length(av_airport), n_stages); % Will contain the nodes to which is the most profitable to fly to
new_pos_hub = av_airport == pos_hub; % In our case, we would not need this variable since in all cases airport 1 is considered and the hub is in pos 2.

for t = 1:n_stages % Study each time stage

    pos = n_stages - t + 1;
    hour = floor(pos/10)+1; % 1 is added because hour 0 is saved in position 1.

    if t == 1 % Initial condition
        Profit(:,pos) = -100000;
        Nodes(:,pos) = av_airport;
        Profit(new_pos_hub,pos) = 0; % Since it has to end in the hub
        
    else % For the rest of time stages
        % We need to take into account the time the aircraft takes!!!!!!!
        for i = 1:length(av_airport) % Departure airport
            if av_airport(i) == pos_hub % If the airport studied is the hub, we have to evaluate all available airports

                for j = 1:length(av_airport) % Arrival airport
                    if pos >= airport_lim_tstage(av_airport(j),1) && pos <= airport_lim_tstage(av_airport(j),2) && j~=pos_hub % We need to specify because that vector contains all the airports
                    % We are getting the profit that we obtain when the aircraft departs from
                    % the hub to each airport at each time stage. Rev(Dep,Arr,Time) Cost(Sym)
                        ind_profit = Rev(av_airport(i),av_airport(j),hour) - Cost(av_airport(i),av_airport(j));
                        profit_from_each_airport(j,pos) = Profit(j,pos+flight_time(av_airport(j))) + ind_profit;
                        
                    elseif pos < airport_lim_tstage(av_airport(j),1) && j~=pos_hub
                        ind_profit = Rev(av_airport(i),av_airport(j),hour) - Cost(av_airport(i),av_airport(j));
                        profit_from_each_airport(j,pos) = Profit(j,pos+flight_time(av_airport(j))) + ind_profit;
                    else
                        profit_from_each_airport(j,pos)= Profit(j,pos+1);
                    end
                end

                    Profit(i,pos) = max(profit_from_each_airport(:,pos));
                    
                    posible_nodes = av_airport(Profit(i,pos)==profit_from_each_airport(:,pos));
                    Nodes(i,pos) = posible_nodes(end);
%                     Nodes(i,pos) = find(Profit(i,pos)==profit_from_each_airport);


            else % In case, it is not the hub, only the hub and the same airport are studied


                if pos >= airport_lim_tstage(av_airport(i),1) && pos <= airport_lim_tstage(av_airport(i),2)
                    profit_hub = Profit(new_pos_hub,pos+flight_time(av_airport(i))) + Rev(av_airport(i),pos_hub,hour) - Cost(av_airport(i),pos_hub);
                    profit_same_airport = Profit(i,pos+1); % This means that the aircraft does not fly
                else
                    profit_hub = -100000;
                    profit_same_airport = -100000;
                end

                Profit(i,pos) = max([profit_same_airport profit_hub]);
                if profit_hub == profit_same_airport
                    Nodes(i,pos) = av_airport(i);
                elseif Profit(i,pos) == profit_same_airport
                    Nodes(i,pos) = av_airport(i);
                elseif Profit(i,pos) == profit_hub
                    Nodes(i,pos) = pos_hub;
                end
            end
        end
    end
end

% Schedule  
% Now the schedule can be computed from the previous results.
% When the profit keeps the same, it means that the aircraft does not fly.
% When it changes, we observe the Nodes matrix and save where the aircraft
% has to fly. Then, we need to compute the time it takes to reach that
% second airport and the profit that we would get.
next_node = 1;
node_schedule(next_node) = pos_hub; % Initial node
n_flights = 0;
t = 1;
s = pos_hub;
while t < n_stages
%     max_profit1 = max(Profit(:,t));
%     max_profit2 = max(Profit(:,t+1));
    
    max_profit1 = Profit(s,t);
    max_profit2 = Profit(s,t+1);
    if max_profit1 ~= max_profit2
        
        next_node = next_node + 1;
        n_flights = n_flights + 1;
        node_schedule(next_node) = Nodes(max_profit1==Profit(:,t),t);
        
        profit_check(n_flights) = Rev(node_schedule(next_node-1),node_schedule(next_node),floor(t/10)+1)-Cost(node_schedule(next_node-1),node_schedule(next_node));
        t_stages_schedule(n_flights,1) = t; % Departure time stage
        
        var = sum(node_schedule(next_node-1)+node_schedule(next_node))-pos_hub; % This variable is used to take the time of that flight
        t = t + flight_time(var);
        t_stages_schedule(n_flights,2) = t; % Arrival time stage
    else
        t = t + 1;
    end
    s = find(node_schedule(next_node)==av_airport);
end

node_schedule_matrix = [node_schedule(1:end-1)' node_schedule(2:end)'];
end

% Function that organized the passengers in their corresponding hour.
function org_pax = org_pax_fun(Demand, pax_taken)
org_pax = zeros(20,20,24);
for t = 1:length(pax_taken)
    h = 1;
    rem_pax = pax_taken(t,end);
    dep_airport = pax_taken(t,1);
    arr_airport = pax_taken(t,2);
    dep_time = floor(pax_taken(t,3)/10)+1;
    while rem_pax > 0
        if h == 1
            
            if Demand(dep_airport,arr_airport,dep_time) >= rem_pax
                org_pax(dep_airport,arr_airport,dep_time) = org_pax(dep_airport,arr_airport,dep_time) + rem_pax;
                rem_pax = 0;
            else
                org_pax(dep_airport,arr_airport,dep_time) = org_pax(dep_airport,arr_airport,dep_time) + (rem_pax - Demand(dep_airport,arr_airport,dep_time));
                rem_pax = (rem_pax - Demand(dep_airport,arr_airport,dep_time));
            end
            h = h + 1;
            
        elseif h == 2
            
            if Demand(dep_airport,arr_airport,dep_time-1) >= rem_pax
                org_pax(dep_airport,arr_airport,dep_time-1) = org_pax(dep_airport,arr_airport,dep_time-1) + rem_pax;
                rem_pax = 0;
            else
                org_pax(dep_airport,arr_airport,dep_time-1) = org_pax(dep_airport,arr_airport,dep_time-1) + (rem_pax - Demand(dep_airport,arr_airport,dep_time-1));
                rem_pax = (rem_pax - Demand(dep_airport,arr_airport,dep_time-1));
            end
            
            h = h + 1;
            
        elseif h == 3
            
            if Demand(dep_airport,arr_airport,dep_time+1) >= rem_pax
                org_pax(dep_airport,arr_airport,dep_time+1) = org_pax(dep_airport,arr_airport,dep_time+1) + rem_pax;
                rem_pax = 0;
            else
                org_pax(dep_airport,arr_airport,dep_time+1) = org_pax(dep_airport,arr_airport,dep_time+1) + (rem_pax - Demand(dep_airport,arr_airport,dep_time+1));
                rem_pax = (rem_pax - Demand(dep_airport,arr_airport,dep_time+1));
            end
            h = h + 1;
            
        else
            
            if Demand(dep_airport,arr_airport,dep_time-2) >= rem_pax
                org_pax(dep_airport,arr_airport,dep_time-2) = org_pax(dep_airport,arr_airport,dep_time-2) + rem_pax;
                rem_pax = 0;
            else
                org_pax(dep_airport,arr_airport,dep_time-2) = org_pax(dep_airport,arr_airport,dep_time-2) + (rem_pax - Demand(dep_airport,arr_airport,dep_time-2));
                rem_pax = (rem_pax - Demand(dep_airport,arr_airport,dep_time-2));
            end
            
        end
    end
end
end

% Function recompute demand
function [Demand_check_ac1, pax_taken_ac1] = new_demand(Demand_h, pax_taken_ac1, check_pax_ac1)

Demand_check_ac1 = Demand_h;

if check_pax_ac1 > 0
    new_org_pax = zeros(20,20,24);
    for t = 1:length(pax_taken_ac1)    
        
        dep = pax_taken_ac1(t,1);
        arr = pax_taken_ac1(t,2);
        dep_t = floor(pax_taken_ac1(t,3)/10)+1;
        pax_req = pax_taken_ac1(t,end);
        
        % Here it is studied each route individually
        if dep_t == 1 
        	act_pos_pax = Demand_check_ac1(dep,arr,dep_t) + Demand_check_ac1(dep,arr,dep_t+1);
        elseif dep_t == 2
        	act_pos_pax = Demand_check_ac1(dep,arr,dep_t)+Demand_check_ac1(dep,arr,dep_t+1)+Demand_check_ac1(dep,arr,dep_t-1);
        elseif dep_t == 24
        	act_pos_pax = Demand_check_ac1(dep,arr,dep_t)+Demand_check_ac1(dep,arr,dep_t-1)+Demand_check_ac1(dep,arr,dep_t-2);
        else
            act_pos_pax = Demand_check_ac1(dep,arr,dep_t)+Demand_check_ac1(dep,arr,dep_t-1)+Demand_check_ac1(dep,arr,dep_t-2)+Demand_h(dep,arr,dep_t+1);
        end
        
        if act_pos_pax >= pax_req 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        h = 1;
            while pax_req > 0
                if h == 1
            
                    if Demand_check_ac1(dep,arr,dep_t) >= pax_req
                        new_org_pax(dep,arr,dep_t) = new_org_pax(dep,arr,dep_t) + pax_req;
                        pax_req = 0;
                    else
                        new_org_pax(dep,arr,dep_t) = new_org_pax(dep,arr,dep_t) + (pax_req - Demand_check_ac1(dep,arr,dep_t));
                        pax_req = (pax_req - Demand_check_ac1(dep,arr,dep_t));
                    end
                    h = h + 1;
            
                elseif h == 2
            
                    if Demand_check_ac1(dep,arr,dep_t) >= pax_req
                        new_org_pax(dep,arr,dep_t-1) = new_org_pax(dep,arr,dep_t-1) + pax_req;
                        pax_req = 0;
                    else
                        new_org_pax(dep,arr,dep_t-1) = new_org_pax(dep,arr,dep_t-1) + (pax_req - Demand_check_ac1(dep,arr,dep_t-1));
                        pax_req = (pax_req - Demand_check_ac1(dep,arr,dep_t-1));
                    end
                    h = h + 1;
            
                elseif h == 3
            
                    if Demand_check_ac1(dep,arr,dep_t) >= pax_req
                        new_org_pax(dep,arr,dep_t+1) = new_org_pax(dep,arr,dep_t+1) + pax_req;
                        pax_req = 0;
                    else
                        new_org_pax(dep,arr,dep_t+1) = new_org_pax(dep,arr,dep_t+1) + (pax_req - Demand_check_ac1(dep,arr,dep_t+1));
                        pax_req = (pax_req - Demand_check_ac1(dep,arr,dep_t+1));
                    end
                    h = h + 1;
            
                else
            
                    if Demand_check_ac1(dep,arr,dep_t-2) >= pax_req
                        new_org_pax(dep,arr,dep_t-2) = new_org_pax(dep,arr,dep_t-2) + pax_req;
                        pax_req = 0;
                    else
                        new_org_pax(dep,arr,dep_t-2) = new_org_pax(dep,arr,dep_t-2) + (pax_req - Demand_check_ac1(dep,arr,dep_t-2));
                        pax_req = (pax_req - Demand_check_ac1(dep,arr,dep_t-2));
                    end
            
                end
            end  
            
            Demand_check_ac1(dep,arr,dep_t) = Demand_check_ac1(dep,arr,dep_t)-new_org_pax(dep,arr,dep_t);
            Demand_check_ac1(dep,arr,dep_t-1) = Demand_check_ac1(dep,arr,dep_t-1)-new_org_pax(dep,arr,dep_t-1);
            Demand_check_ac1(dep,arr,dep_t-2) = Demand_check_ac1(dep,arr,dep_t-2)-new_org_pax(dep,arr,dep_t-2);
            Demand_check_ac1(dep,arr,dep_t+1) = Demand_check_ac1(dep,arr,dep_t+1)-new_org_pax(dep,arr,dep_t+1);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %             
        else
            
            % The loop will enter here when the number of passengers that 
            % are required initially are higher than the actual possible 
            % that can be taken. So, the passengers transported are changed
            % by the possible ones and the Demand is changed to zero.
            
            pax_taken_ac1(t,end) = act_pos_pax;
            Demand_check_ac1(dep,arr,dep_t) = 0;
            Demand_check_ac1(dep,arr,dep_t-1) = 0;
            Demand_check_ac1(dep,arr,dep_t-2) = 0;
            Demand_check_ac1(dep,arr,dep_t+1) = 0;
        end
        
    end
else
    Demand_check_ac1 = Demand_check_ac1 - org_pax_ac1;
end
end
