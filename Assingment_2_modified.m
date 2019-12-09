%% Assignment 1 - (AE4423-19 Airline Planning and Optimization)
clc, clear all, close all
% addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Community129\cplex\matlab\x64_win64')
%addpath('C:\Program
%Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab\x64_win64') Guille
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Community129\cplex\matlab\x64_win64')
Cap = xlsread('Group_31.xlsx','Flights','F2:F133');%Capacity of each Flight
[~,Flight,~]  = xlsread('Group_31.xlsx','Flights','A2:A133');
Demand = xlsread('Group_31.xlsx','Itineraries','G2:G461'); %Demand of each itinerary
Fare = xlsread('Group_31.xlsx','Itineraries','F2:F461'); %Fare of each itinerary
[~,Itin_f1,~] = xlsread('Group_31.xlsx','Itineraries','D2:D461');
[~,Itin_f2,~] = xlsread('Group_31.xlsx','Itineraries','E1:E461');

%addzeros = length(Itin_f1) - length(Itin_f2);
%Itin_f2 = [zeros(addzeros,1) ; Itin_f2];
Itin = [Itin_f1 Itin_f2(2:end)];

%Considering the ficticious itinerary

%%
%Cap = [0;Cap];
%Flight = ['CC0000';Flight];
Demand = [0;Demand];
Fare = [0;Fare];
Itin = ['0' '0' ; Itin];



P = length(Fare); %Number of itineraries 9 = 8 + 1
L = length(Flight); %Number of flights 6

Itin_b_f1 = xlsread('Group_31.xlsx','Recapture','A2:A441');
Itin_b_f2 = xlsread('Group_31.xlsx','Recapture','B2:B441');
Recap = xlsread('Group_31.xlsx','Recapture','C2:C441');
Itin_b_f1 = Itin_b_f1 + 2; %2 because we need +1 for the fictitious and +1 for the starting itin = 0
Itin_b_f2 = Itin_b_f2 + 2;

b = zeros(P); %Recapture ratio

for i = 1:length(Itin_b_f1)
        Itin_b_f1_i = Itin_b_f1(i);
        Itin_b_f1_j = Itin_b_f2(i);
        b(Itin_b_f1_i,Itin_b_f1_j) = Recap(i);
end

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
Q = (delta*Demand')';

DV = P*P;
% DV = (P-1)*P;
number_columns = P; % Number of itineraries considered. At the beginning 
% only the ones that wanted to go in the normal ones and they are moved to 
% the fictitious itinerary.
columns = [1:P]; % This vector indicates which columns are considered. 
% Starting with the first 8 that are refered to the pax moved to the
% fictitious itinerary.

Red_flag = "False";
while Red_flag == "False"
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
        cplex.addRows(-inf,C2,RHS,sprintf('Constraint2_%d',j));
    %     end
    end

    %% Solution
    cplex.writeModel([model '.lp']) %Store the model to an .lp file for debugging
    cplex.Param.timelimit.Cur = 10; %Timelimit of the solver in seconds, more useful for larger models
    cplex.solve();

    dual_variables = cplex.Solution.dual; % Dual variables

    t = cplex.Solution.x;
    pis = cplex.Solution.dual(1:L);
    sigmas = cplex.Solution.dual(L+1:end);

    piprice = zeros(P,P);
    for j=2:P
        for k=2:P
            for i=1:L
                piprice(j,k) = piprice(j,k)-pis(i)*delta(i,j)+b(j,k)*pis(i)*delta(i,k);
                %tnew(j,k) = (Fare(j)-pis(i)*delta(i,j)) - b(j,k)*(Fare(k)-pis(i)*delta(i,k))-sigmas(j);
            end
            tnew(j,k) = Fare(j)-b(j,k)*Fare(k)+piprice(j,k)-sigmas(j);
        end
    end
    
    if any(any(tnew<0))  
        Position = find(tnew==min(tnew(tnew<0)));
        columns = [columns Position(1)];
        
    else
        Red_flag = "True";
    end
    
end