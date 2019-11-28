% Exercise 1 - Group 16
% Startup
clc
clear
close all

% ---- Constants ----
fuelc17 = 1.42;
fuelc22 = 1.6;
nodes_eur = 20; %Number of european airports
nodes_intern = 24; %Number of international and european airports
SF = 10; %Scaling Factor between european and international flights

% ---- Read Data From Excel ----
fprintf('Reading data from Excel file\n')
demnd = xlsread('input16_Ex2.xlsx',2,'C15:V34'); % Demand in 2017
gdp17 = xlsread('input16_Ex2.xlsx',1,'G4:G27');  % GDP in 2017
pop17 = xlsread('input16_Ex2.xlsx',1,'C4:C27')/1000;  % Population in 2017
gdp10 = xlsread('input16_Ex2.xlsx',1,'F4:F27');  % GDP in 2010
pop10 = xlsread('input16_Ex2.xlsx',1,'B4:B27')/1000;  % Population in 2010
%(in Gravity Model Population is expressed per 1000 inhabitants)

% Airport Data
apdat.lat = xlsread('input16_Ex2.xlsx',2,'C6:Z6'); % Latitude
apdat.lon = xlsread('input16_Ex2.xlsx',2,'C7:Z7'); % Longitude
fprintf('Data Read from Excel\n')


% ---- Cosntruct the model matrices ----
% Construct matrix with distances from i to j
distm = zeros(nodes_intern,nodes_intern); % Initialisation
for j = 1:nodes_intern
    for i = 1:(j-1) % The matrix is symmetric
        % Distance from arc formula
            distm(i,j) = 6371*2*asin(sqrt(...
                sind(0.5*(apdat.lat(i)-apdat.lat(j)))^2 +...
                cosd(apdat.lat(i)) * cosd(apdat.lat(j)) * ...
                sind(0.5*(apdat.lon(i)-apdat.lon(j)))^2));
            distm(j,i) = distm(i,j);
    end
end

%% Calibration of Demand
% Construct the min squares matrix
Aminsq = ones(nodes_eur^2 - nodes_eur,4); % Min Squares matrix - nodes^2 data, 4 unknowns
Bminsq = ones(nodes_eur^2 - nodes_eur,1); % Min Squares independent vector
count = 1;
for i = 1:nodes_eur
    for j = 1:nodes_eur
        if j ~= i
            Aminsq(count,2) = log(pop17(i)*pop17(j));
            Aminsq(count,3) = log(gdp17(i)*gdp17(j));
            Aminsq(count,4) = -log(fuelc17*distm(i,j));
            Bminsq(count,1) = log(demnd(i,j));
            count = count + 1;
        end
    end
end
% ---- Min Squares solution ----
solvec = Aminsq'*Aminsq\(Aminsq'*Bminsq); % A^T*A*x = A^T*b
k  = exp(solvec(1));
b1 = solvec(2);
b2 = solvec(3);
b3 = solvec(4);
fprintf('Model Parameters:\n k  = %g \n b1 = %g \n b2 = %g \n b3 = %g \n',...
    k,b1,b2,b3);

% Reconstruct the demand and evaluate the error
demrc = zeros(nodes_eur);
for i = 1:nodes_eur
    for j=1:nodes_eur
        if i ~= j
            demrc(i,j) = round(k * (pop17(i)*pop17(j))^b1 * (gdp17(i)*gdp17(j))^b2 ...
                / (fuelc17*distm(i,j))^b3);
        end
    end
end
errm= demnd - demrc;
err = norm(reshape(errm,nodes_eur^2,1));
errp= norm(reshape(errm,nodes_eur^2,1))/norm(reshape(demnd,nodes_eur^2,1));
fprintf('Error Value:     %g\nPercentage Error: %g %%\n',err,errp*100);


%% Forecast
Years = [0; 7]; %2010 and 2017 years
For_time = Years(end) + 5; %We forecast 5 years more

%X matrix to solve phi = (X^T*X)*(X^T*y)
x_years = [ones(length(Years),1) Years]; 

% ---- Linear Forecast of Population
%Prepare Y matrix to solve phi = (X^T*X)*(X^T*y)
y_pop = zeros(2, nodes_intern);
y_pop(1,:) = pop10;
y_pop(2,:) = pop17;
phi_pop=zeros(2,nodes_intern);

%Get Coefficients of each equation for Population
for i = 1:nodes_intern 
phi_pop(:,i) = (x_years'*x_years)\(x_years'*y_pop(:,i));
end


%Linear Forecasting at year 2022: (2022-2017=5years)
pop22 = zeros(nodes_intern,1); %Population forecast, linear 
for i = 1:nodes_intern
pop22(i,1) = phi_pop(1,i) + phi_pop(2,i)*(For_time);
end


% ---- Linear Forecast of GDP
%Prepare Y matrix to solve phi = (X^T*X)*(X^T*y)
y_GDP = zeros(2, nodes_intern);
y_GDP(1,:) = gdp10;
y_GDP(2,:) = gdp17;
phi_GDP=zeros(2,nodes_intern);

%Get Coefficients of each equation for Population
for i = 1:nodes_intern 
phi_GDP(:,i) = (x_years'*x_years)\(x_years'*y_GDP(:,i));
end

%Linear Forecasting at year 2022: (2022-2017= 5years)
gdp22 = zeros(nodes_intern,1); %Population forecast, linear 
for i = 1:nodes_intern
gdp22(i,1) = phi_GDP(1,i) + phi_GDP(2,i)*(For_time);
end


%% Future Demand
% Calculate future demand with GDP and population forecasts
demand22 = zeros(nodes_intern);
for i = 1:nodes_intern
    for j=1:nodes_intern
        if i ~= j
            demand22(i,j) = round(k * (pop22(i)*pop22(j))^b1 * (gdp22(i)*gdp22(j))^b2 ...
                / (fuelc22*distm(i,j))^b3);
        end
    end
end

% Apply Scaling Factor for US Flights
demand22(:    ,21:24) = demand22(:    ,21:24)*SF;
demand22(21:24, 1:20) = demand22(21:24, 1:20)*SF;

% Write data in Excel sheet
fprintf('Writing Data in Excel\n')
xlswrite('input16_Ex2.xlsx',demand22,4,'C2:Z25');
fprintf('Data written\n')






