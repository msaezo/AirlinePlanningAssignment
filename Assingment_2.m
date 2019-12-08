%% Assignment 1 - (AE4423-19 Airline Planning and Optimization)
clc, clear all, close all
% addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Community129\cplex\matlab\x64_win64')
%addpath('C:\Program
%Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab\x64_win64') Guille
addpath('C:\Program Files\IBM\ILOG\CPLEX_Studio_Community129\cplex\matlab\x64_win64')





Cap = [108 144 144 108 108 108];%Capacity of each Flight
Demand = [90 40 90 50 100 80 120 100]; %Demand of each itinerary
Fare = [150 130 150 130 150 150 100 100]; %Fare of each itinerary
Flight  = [301 102 101 302 201 202];
Itin = [301 0; 301 101; 102 0; 102 302; 101 0; 302 0; 201 0; 202 0];
P = length(Fare); %Number of itineraries
L = length(Flight); %Number of flights
b = zeros(P); %Recapture ratio
b(2,8) = 0.2;
delta = zeros(L,P);
for i=1:L
    for j=1:P
        if Flight(i)==Itin(j,1) || Flight(i)==Itin(j,2)
            delta(i,j) = 1;
            %Q(i,j) = delta(i,j)*Demand(j);
        end
    end
end
Q = (delta*Demand')';





%%

% delta = [ones(6,1) delta];
%% Initialize model

model = 'Assignment1_p2'; % Name of the model
cplex = Cplex(model); % Initialize Cplex
cplex.Model.sense = 'minimize'; % Minimize spillage
%DV = P*P + P;
DV = P;
obj = zeros(1,DV);


%% Objective function
    for j = 1:P
        obj(j)=Fare(j);
    end
lb = zeros(DV,1);
ub = inf(DV,1);
obj = obj';
%ctype1 = char(ones(1, (DV)) * ('C'));  %I=integer. Other options, C=continous, B=binary.
%ctype = strcat(ctype1);

cplex.addCols(obj, [], lb, ub);
%cplex.addCols(obj, [], lb, ub, ctype); %This can also be done without NameDV

%% Constraint 1
C1_matrix = zeros(L,DV);
row_c1 = 1;
for i = 1:L
    C1 = zeros(1,DV);
    for j = 1:P
        C1(j) = delta(i,j);
    end
    LHS = Q(i) - Cap(i);
    cplex.addRows(LHS,C1,inf,sprintf('Constraint1%d_%d',i,j));
    C1_matrix(row_c1,:) = C1; 
    row_c1 = row_c1+1;
end

%% Constraint 2

for j =1:P
    C2 = zeros(1,DV);
    C2(j) = 1;
    RHS = Demand(j);
    cplex.addRows(-inf,C2,RHS,sprintf('Constraint2%d_%d',i,j));
end

%% Constraint 3 
for j=1:P
    C3 = zeros(1,DV);
    C3(j) = 1;
    cplex.addRows(0,C3,inf,sprintf('Constraint3%d_%d',i,j));
end

%%
cplex.writeModel([model '.lp']) %Store the model to an .lp file for debugging
cplex.Param.timelimit.Cur = 10; %Timelimit of the solver in seconds, more useful for larger models
cplex.solve();

t = cplex.Solution.x;
pis = cplex.Solution.dual(1:L);
sigmas = cplex.Solution.dual(L+1:L+P);

piprice = zeros(8,8);
for j=1:P
    for k=1:P
        for i=1:6
            piprice(j,k) = piprice(j,k)-pis(i)*delta(i,j)+b(j,k)*pis(i)*delta(i,k);
            %tnew(j,k) = (Fare(j)-pis(i)*delta(i,j)) - b(j,k)*(Fare(k)-pis(i)*delta(i,k))-sigmas(j);
        end
        tnew(j,k) = Fare(j)-b(j,k)*Fare(k)+piprice(j,k)-sigmas(j);
    end
end

