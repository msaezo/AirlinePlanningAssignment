%% Assignment 1 - (AE4423-19 Airline Planning and Optimization)
clc, clear all, close all
% addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Community129\cplex\matlab\x64_win64')
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab\x64_win64') % Guille
% addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Community129\cplex\matlab\x64_win64')
tic

%% DATA
Cap = xlsread('Group_31.xlsx','Flights','F2:F133');%Capacity of each Flight
[~,Flight,~]  = xlsread('Group_31.xlsx','Flights','A2:A133');
Demand = xlsread('Group_31.xlsx','Itineraries','G2:G461'); %Demand of each itinerary
Fare = xlsread('Group_31.xlsx','Itineraries','F2:F461'); %Fare of each itinerary
[~,Itin_f1,~] = xlsread('Group_31.xlsx','Itineraries','D2:D461');
[~,Itin_f2,~] = xlsread('Group_31.xlsx','Itineraries','E1:E461');
Itin_f1_mod = [' ' ; Itin_f1];
Itin_f2_mod = [' ' ; Itin_f2(2:end)];
Itin = [string(Itin_f1_mod) string(Itin_f2_mod)];
Flight = string(Flight);

%addzeros = length(Itin_f1) - length(Itin_f2);
%Itin_f2 = [zeros(addzeros,1) ; Itin_f2];
%Itin = [Itin_f1 Itin_f2(2:end)];

%%
% Considering the ficticious itinerary
Demand = [0;Demand];
Fare = [0;Fare];

P = length(Fare); %Number of itineraries 
L = length(Flight); %Number of flights 

Itin_b_f1 = xlsread('Group_31.xlsx','Recapture','A2:A441');
Itin_b_f2 = xlsread('Group_31.xlsx','Recapture','B2:B441');
Recap = xlsread('Group_31.xlsx','Recapture','C2:C441');
Itin_b_f1 = Itin_b_f1 + 2; %2 because we need +1 for the fictitious and +1 for the starting itin = 0
Itin_b_f2 = Itin_b_f2 + 2;

b = zeros(P); %Recapture ratio

% The following loop is used to organize the itineraries as well as to
% construct the recapture matrix
for i = 1:length(Itin_b_f1)
        Itin_b_f1_i = Itin_b_f1(i);
        Itin_b_f1_j = Itin_b_f2(i);
        b(Itin_b_f1_i,Itin_b_f1_j) = Recap(i);
end

% The following variables are used to define that the recapture rate to
% itinerary itself is equal to 1 for all p.
diagonal = ones(length(Fare),1);
b_diagonal = diag(diagonal);
b = b + b_diagonal;

% The recapture rate to the 'fictitious itinerary is equal to 1 fo all p.
b(:,1) = 1;
b(1,:) = 0; 

delta = zeros(L,P);
for i=1:L
    for j=1:P
        if Flight(i)==Itin(j,1) || Flight(i)==Itin(j,2)
            delta(i,j) = 1;
        end
        if j == 1
            delta(i,j) = 1;
        end
    end
end
Q = (delta*Demand);

DV = P*P;
% DV = (P-1)*P;
number_columns = P; % Number of itineraries considered. At the beginning 
% only the ones that wanted to go in the normal ones and they are moved to 
% the fictitious itinerary.
columns = [1:DV]; % This vector indicates which columns are considered. 
% Starting with the first 8 that are refered to the pax moved to the
% fictitious itinerary.

tic;

%% Initialize model

model = 'Assignment1_p2'; % Name of the model
cplex = Cplex(model); % Initialize Cplex
cplex.Model.sense = 'minimize'; % Minimize spillage

%% Objective function

obj = zeros(1,DV); % This vector indicates the position of each decision 
% variable. The first 8 elements are going to be related to the fictitious
% itinerary.

for p = 2:P
    for r = 1:P
        if sum(varindex_4(1, r, p, P) == columns) == 1
            obj(varindex_4(1, r, p, P)) = Fare(p)-b(p,r)*Fare(r);
        end
    end
end

lb = zeros(DV,1);
ub = inf(DV,1);
obj = obj';

cplex.addCols(obj, [], lb, ub); 

%% Constraint 1
for i = 1:L
    C1 = zeros(1,DV);
    for p = 1:P
        for r = 1:P
            if sum(varindex_4(1, r, p, P) == columns) == 1
                C1(varindex_4(1, r, p, P)) = delta(i,p);
            end
            if sum(varindex_4(1, p, r, P) == columns) == 1
                C1(varindex_4(1, p, r, P)) = delta(i,p)*b(r,p);
            end
        end
    end
    LHS = Q(i) - Cap(i);

    cplex.addRows(LHS,C1,inf,sprintf('Constraint1_%d',i));
end
%% Constraint 2
for p = 1:P
%     if p ~= 1
    C2 = zeros(1,DV);
    for r = 1:P
        if sum(varindex_4(1, r, p, P) == columns) == 1
            C2(varindex_4(1, r, p, P)) = 1;
        end
    end
    RHS = Demand(p);
    cplex.addRows(-inf,C2,RHS,sprintf('C nstraint2_%d',j));
%     end
end
%% Solution
cplex.writeModel([model '.lp']) %Store the model to an .lp file for debugging
%cplex.Param.timelimit.Cur = 10; %Timelimit of the solver in seconds, more useful for larger models
cplex.solve();
toc;