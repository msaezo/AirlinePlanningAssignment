%Demand forecast
clc
clear all

%Initial Read
 pop_2014 = xlsread('AE4423_Datasheets.xls','General','B4:B23');
    GDP_2014 = xlsread('AE4423_Datasheets.xls','General','F4:F23');
    Demand_2014 = xlsread('AE4423_Datasheets.xls','Group 31','C13:V32');
    Lat = xlsread('AE4423_Datasheets.xls','Group 31','C6:V6');
    Lon = xlsread('AE4423_Datasheets.xls','Group 31','C7:V7');
    
    f_cost = 1.42;
    N = length(Lat);
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

    
% Fmincon
x0 = [1,1,1,1];
lb = [0,0,0,0];
ub = [1,1,1,1];
f=@(k,b1,b2,b3) LSE_demand(k,b1,b2,b3,pop_2014,GDP_2014,Demand_2014,Distance,N,f_cost,0);

options=optimset('Display','iter','Algorithm','sqp');
fun=@(x) f(x(1),x(2),x(3),x(4)); % [k,b1,b2,b3]

% Options for the optimization
options.Display         = 'iter-detailed';
options.Algorithm       = 'sqp';
options.FunValCheck     = 'off';
options.DiffMinChange   = 1e-5;         % Minimum change while gradient searching
options.DiffMaxChange   = 1e-4;         % Maximum change while gradient searching
options.TolCon          = 1e-6;         % Maximum difference between two subsequent constraint vectors [c and ceq]
options.TolFun          = 1e-6;         % Maximum difference between two subsequent objective value
options.TolX            = 1e-6;         % Maximum difference between two subsequent design vectors
options.PlotFcns        = {@optimplotx,@optimplotfunccount,@optimplotfval,@optimplotconstrviolation,@optimplotstepsize,@optimplotfirstorderopt};
options.MaxIter         = 500;           % Maximum iterations

tic;
[x,fval,exitflag,output]=fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
toc;

Demand_est = LSE_demand(x(1),x(2),x(3),x(4),pop_2014,GDP_2014,Demand_2014,Distance,N,f_cost,1);