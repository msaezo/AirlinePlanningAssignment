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

time_stages = linspace(24-1/10,0,n_stages); % It starts in 23.54 and ends in 0

% Minimum and maximum posible time stages for each airport (It depends on
% the aircraft since the velocity and the TAT are different)
airport_lim_tstage = zeros(N,2,AC_types);

for k = 1:AC_types
    time = ceil((15 + 15 + Distance(:,2)/Speed(k)*60 + Ave_TAT(k))/6);
    airport_lim_tstage(:,:,k) = [1+time 240-time];
end

%% 
% The following two variables are the controllers that indicate if there 
% are aircraft available or if new routes are still profitable.

profitable = 1;
AC_available = 1;

while profitable == 1 && AC_available == 1
    
%     Now, passengers per route and per hour are computed assumming that 
%     we can capture the demand estimated for hours t, t-1, t+1 and t-2.
%     Then, the revenues from those passengers can be computed:

    [Act_pax_ac1, Rev_ac1] = act_demand_and_rev(Demand_h, Seats(1), Distance, N);
    [Act_pax_ac2, Rev_ac2] = act_demand_and_rev(Demand_h, Seats(2), Distance, N);
    [Act_pax_ac3, Rev_ac3] = act_demand_and_rev(Demand_h, Seats(3), Distance, N);

%     Dynamic programming
% Aircraft 1
[Profit_ac1, Nodes_ac1, Schedule_ac1] = Dynamic_Programming(Rev_ac1, Cost(:,:,1), airport_lim_tstage(:,:,1), av_airports_ac1, pos_hub, n_stages);
% Aircraft 2
[Profit_ac2, Nodes_ac2, Schedule_ac2] = Dynamic_Programming(Rev_ac2, Cost(:,:,2), airport_lim_tstage(:,:,2), av_airports_ac2, pos_hub, n_stages);
% Aircraft 3
[Profit_ac3, Nodes_ac3, Schedule_ac3] = Dynamic_Programming(Rev_ac3, Cost(:,:,3), airport_lim_tstage(:,:,3), av_airports_ac3, pos_hub, n_stages);



























    if something
        profitable = 0;
    end
    
    if another something
        AC_available = 0;
    end
    
%     We need also the constraint of the block time.
end


%% FUNCTIONS

% Function to compute the demand and the revenue.
function [X, Rev] = act_demand_and_rev(Demand_h, Seats, Distance, N)
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
function [Profit, Nodes , Schedule] = Dynamic_Programming(Rev, Cost, airport_lim_tstage, av_airport, pos_hub, n_stages)

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
        
        for i = 1:length(av_airport) % Departure airport
            if av_airport(i) == pos_hub % If the airport studied is the hub, we have to evaluate all available airports

                for j = 1:length(av_airport) % Arrival airport
                    if pos >= airport_lim_tstage(av_airport(j),1) && pos <= airport_lim_tstage(av_airport(j),2) % We need to specify because that vector contains all the airports
                    % We are getting the profit that we obtain when the aircraft departs from
                    % the hub to each airport at each time stage. Rev(Dep,Arr,Time) Cost(Sym)
                        ind_profit = Rev(av_airport(i),av_airport(j),hour) - Cost(av_airport(i),av_airport(j));
                        profit_from_each_airport(j) = Profit(j,pos+1) + ind_profit;
                    else
                        profit_from_each_airport(j) = Profit(j,pos+1);
                    end
                end

                    Profit(i,pos) = max(profit_from_each_airport);
                    Nodes(i,pos) = av_airport(Profit(i,pos)==profit_from_each_airport);

            else % In case, it is not the hub, only the hub and the same airport are studied

                profit_same_airport = Profit(i,pos+1); % This means that the aircraft does not fly

                if pos >= airport_lim_tstage(av_airport(i),1) && pos <= airport_lim_tstage(av_airport(i),2)
                    profit_hub = Profit(new_pos_hub,pos+1) + Rev(av_airport(i),pos_hub,hour) - Cost(av_airport(i),pos_hub);
                else
                    profit_hub = -100000;
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

% Schedule = 
% Profit_updated =
Schedule = 1;











end