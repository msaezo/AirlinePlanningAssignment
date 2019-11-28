% Exercise 1 - Group 16
% Startup
clc
clear
close all
plots_on = 1; % 1 to plot, 0 to not plot. MAPPING TOOLBOX NEEDED. 

%% Preprocess
% ---- Constants ----
nodes = 20;  % Number of airports
nacft = 3; %length(acrft.num);  % Number of aircraft types
fuelc = 1.42; % Fuel costs
loadf = 0.75; % Load factor
avgut = 70; % Average Utilisation time of the aircraft

% ---- Read Data From Excel ----
fprintf('Reading data from Excel file\n')
% Airport Data
apdat.lat = xlsread('input16_Ex1.xlsx',2,'C6:V6'); % Latitude
apdat.lon = xlsread('input16_Ex1.xlsx',2,'C7:V7'); % Longitude
apdat.rwl = xlsread('input16_Ex1.xlsx',2,'C8:V8'); % Runway length

% Aircraft Data
acrft.num = xlsread('input16_Ex1.xlsx',2,'B12:F12'); % Available aircraft
acrft.vel = xlsread('input16_Ex1.xlsx',3,'B2:F2'); % Cruise speed
acrft.sts = xlsread('input16_Ex1.xlsx',3,'B3:F3'); % Seats available
acrft.tat = xlsread('input16_Ex1.xlsx',3,'B4:F4'); % Turn-around time
acrft.rng = xlsread('input16_Ex1.xlsx',3,'B5:F5'); % Range
acrft.rwl = xlsread('input16_Ex1.xlsx',3,'B6:F6'); % Runway length required
acrft.cl  = xlsread('input16_Ex1.xlsx',3,'B7:F7'); % Weekly lease costs
acrft.cx  = xlsread('input16_Ex1.xlsx',3,'B8:F8'); % Fixed operating costs
acrft.ct  = xlsread('input16_Ex1.xlsx',3,'B9:F9'); % Time costs parameter
acrft.cf  = xlsread('input16_Ex1.xlsx',3,'B10:F10'); % Fuel costs parameter

% Demand matrix
demnd = xlsread('input16_Ex1.xlsx',2,'C15:V34');

fprintf('Excel file read\n')
% ---- Preprocess the data ----
% Construct matrix with distances from i to j
distm = zeros(nodes,nodes); % Initialisation
for j = 1:nodes
    for i = 1:(j-1) % The matrix is symmetric
        % Distance from arc formula
            distm(i,j) = 6371*2*asin(sqrt(...
                sind(0.5*(apdat.lat(i)-apdat.lat(j)))^2 +...
                cosd(apdat.lat(i)) * cosd(apdat.lat(j)) * ...
                sind(0.5*(apdat.lon(i)-apdat.lon(j)))^2));
            distm(j,i) = distm(i,j);
    end
end

% Construct Revenue matrix 
yrpkm = 5.9.*distm.^-0.76 + 0.043; % RPK matrix
rvnum = distm.*yrpkm; % Revenue per passenger from i to j

for i = 1:nodes % Symmetric matrix
    rvnum(i,i) = 0; % Passengers from i to i generate no revenue
end

% Construct Cost Matrix
costm = zeros(nodes, nodes, nacft); % Initialisation
for k = 1:nacft
    for j = 1:nodes
        for i = 1:(j-1) % Symmetric matrix
            % Fixed Operating Costs
            cfixd = acrft.cx(k);
            % Time-based costs
            ctime = acrft.ct(k)/acrft.vel(k)*distm(i,j);
            % Fuel costs
            cfuel = 0.666666667*acrft.cf(k)*distm(i,j)*fuelc;
            % Sum of the costs
            costm(i,j,k) = cfixd + cfuel + ctime;
            if i == 1 || j == 1
                costm(i,j,k) = 0.7*costm(i,j,k);
            end
            costm(j,i,k) = costm(i,j,k);
        end
        costm(j,j,k) = 1e6; % Prevent the optimizer from scheduling flights from A to A
    end
end

% Turn-around time matrix
ltotm = zeros(nodes, nodes, nacft); % Initialisation
for k = 1:nacft
    for j = 1:nodes
        for i = 1:nodes % Non-symmetric matrix
            if (j == 1)
                ltotm(i,j,k) = max(2*acrft.tat(k)/60, 1);
            else
                ltotm(i,j,k) = acrft.tat(k)/60;
            end
        end
    end
end

%% C-Plex Setup
% ---- Objective Function ----
nvars = 2*nodes^2 + nodes^2 * nacft; % Number of decission variables
% Structure of the decission variable vector
% [xij (direct pax) wij (connecting pax) zijk (flights offered)]

% function to translate i-j-k index to dv index
clear i j k
varindex = @(i, j, k) (i + (j-1)*nodes)*(k == -1) + ... % locates X_ij if k == -1
    (nodes^2 + i + (j-1)*nodes)*(k == 0) + ...  % locates W_ij if k == 0
    (2*nodes^2 + i + (j - 1)*nodes + (k - 1)*nodes^2)*(k > 0); % Locates Z_ijk if k > 0

% Define the CPLEX model
model             =   'MCF_Model';
cplex             =   Cplex(model);
cplex.Model.sense =   'maximize';
ctype             =   char(ones(1,nvars) * ('I')); % Integer

% Objective function and bounds
% Maximize YIELD - COST
obj = [reshape(rvnum, nodes^2, 1);
     reshape(rvnum, nodes^2, 1);
    -reshape(costm, nacft*nodes^2,1)];

% Array with DV names
namedv = zeros(nvars,9);                                     
for k = -1:nacft
    for j = 1:nodes
        for i = 1:nodes
            idx = varindex(i, j, k);
            if k == -1
                namedv(idx,:) =  ['X_' num2str(i,'%02d') ',' num2str(j,'%02d') '__'];
            elseif k == 0
                namedv(idx,:) =  ['W_' num2str(i,'%02d') ',' num2str(j,'%02d') '__'];
            else
                namedv(idx,:) =  ['Z_' num2str(i,'%02d') ',' num2str(j,'%02d') '_' num2str(k,'%01d')];
            end
        end
    end
end
namedv = char(namedv);
% Bounds
lb = zeros(nvars,1);
ub = inf(nvars,1);
% Add objective to cplex
cplex.addCols(obj, [], lb, ub, ctype, namedv);
% ---- Constraints ----
% List of constraints of the problem
% - Demand
% - Transfer passengers - not to or from the hub
% - Flight Capacity + Range
% - Utilization time
% - Flights in and out each node
% TO IMPLEMENT
% - Slots

% Demand constraint
for i = 1:nodes
    for j = 1:nodes
        C1 = zeros(1, nvars);       %Setting coefficient matrix with zeros
        % x_ij + w_ij <= dem_ij
        C1(varindex(i, j,-1)) = 1; % x_ij
        C1(varindex(i, j, 0)) = 1; % w_ij
        cplex.addRows(0, C1, demnd(i,j),sprintf('Demand%d_%d',i,j));
    end
end

% Transfer Passengers constraint
for i = 1:nodes
    for j = 1:nodes
        C2 = zeros(1, nvars);       %Setting coefficient matrix with zeros
        % w_ij <= dem_ij*hub!
        C2(varindex(i, j, 0)) = 1;  % w_ij
        g = 1 - 1*((i == 1) || (j == 1));
        cplex.addRows(0, C2, g*demnd(i,j),sprintf('IsHub%d_%d',i,j));
    end
end

% Capacity constraint + Range and Runway Length
for i = 1:nodes
    for j = 1:nodes
        C3 = zeros(1, nvars);       %Setting coefficient matrix with zeros
        % x_ij + sum(w_im,m)*(1-gj) + sum(w_mj,m)*(1-gi) - sum(z_ijk * s_k * lf,k) <= 0
        % if d_ij >= range_k, cap_ijk = 0
        C3(varindex(i, j,-1)) = 1;  % x_ij
        if (j == 1) % Only add hub capacity terms if one of the points is a hub
            for m = 1:nodes
                C3(varindex(i, m, 0)) = 1;  % w_im* (1 - gj)
            end
        end
        if (i == 1) % Only add hub capacity terms if one of the points is a hub
            for m = 1:nodes
                C3(varindex(m, j, 0)) = 1;  % w_mj* (1 - gi)
            end
        end
        for k = 1:nacft
            rwy_c = (apdat.rwl(i)>= acrft.rwl(k))*(apdat.rwl(j)>= acrft.rwl(k)); % Runway length constraint
            if (distm(i,j) <= acrft.rng(k)) && (rwy_c)% If the distance is in range
                C3(varindex(i, j, k)) = -loadf*acrft.sts(k);
            end
        end
        cplex.addRows(-Inf, C3, 0, sprintf('Capacity%d_%d',i,j));
    end
end

% Balance constraint
for i = 1:nodes
    for k = 1:nacft
        C4 = zeros(1, nvars);       %Setting coefficient matrix with zeros
        % sum(z_ijk,j) == sum(z_jik,j)
        for j = 1:nodes
            C4(varindex(i, j, k)) =  1;  % z_ij
            C4(varindex(j, i, k)) = -1;  % z_ij
        end
        cplex.addRows(0, C4, 0, sprintf('Balance%d_%d',i,k));
    end
end

% Utilisation Constraint
for k = 1:nacft
    C5 = zeros(1, nvars);       %Setting coefficient matrix with zeros
    % sum(d_ij/v_k + LTO_ijk, i,j) <= avgut*n_k
    for i = 1:nodes
        for j = 1:nodes
            C5(varindex(i, j, k)) = distm(i,j)/acrft.vel(k) + ltotm(i,j,k);
        end
    end
    cplex.addRows(0, C5, acrft.num(k)*avgut, sprintf('Utilisation%d',k));
end

% ---- Execute model
cplex.Param.mip.limits.nodes.Cur    = 1e+8;   % max number of nodes
cplex.Param.timelimit.Cur           = 30;   % max time in

%   Run CPLEX
cplex.solve();
cplex.writeModel([model '.lp']);
%% Postprocess
sol  = cplex.Solution.x;
[~, tags] = xlsread('input16_Ex1.xlsx',2,'C5:V5');
tags = tags(1:nodes);
% Small violation of integer constraints ~1e-9
dpax = round(reshape(sol(1:nodes^2),nodes,nodes)); % Direct Passengers matrix
cpax = round(reshape(sol(1+nodes^2: 2*nodes^2),nodes,nodes));   % Connecting Passengers matrix
flts = round(reshape(sol(1+2*nodes^2: end),nodes,nodes,nacft)); % Flights matrix

% Capacity matrix
cap  = zeros(nodes);
for k = 1:nacft
    cap = cap + round(acrft.sts(k) * flts(:,:,k));
end

% Passenger flows (connecting and direct)
fpax = zeros(nodes);
for j = 1:nodes
    for i = 1:nodes
        fpax(i,j) = dpax(i,j);
        if i == 1
            fpax(i,j) = fpax(i,j) + sum(cpax(:,j));
        elseif j == 1
            fpax(i,j) = fpax(i,j) + sum(cpax(i,:));
        end
    end
end
% Operative Results and KPIs
tcost = sum(costm.*flts,'all') + sum(acrft.num.*acrft.cl); % Total cost
trvnu = sum(rvnum.*(dpax+cpax),'all'); % Total Revenue
task  = sum(cap.*distm,'all'); % Total ASK
trpk  = sum(fpax.*distm,'all');% Total RPK
fprintf('\n\nKey Performance Indicators:\n');
fprintf('   ASK:  %8.0f ASK/week\n', task);
fprintf('   RPK:  %8.0f RPK/week\n', trpk);
fprintf('  CASK:    %5.4f €/ASK\n', tcost/task);
fprintf('  RASK:    %5.4f €/ASK\n', trvnu/task);
fprintf(' Yield:    %5.4f €/RPK\n', trvnu/trpk);
fprintf('  ANLF:     %4.2f %%\n', trpk/task*100);
fprintf('  BELF:     %4.2f %%\n', (tcost/task)/(trvnu/trpk)*100);

fprintf('\n\nOperative Results:\n');
fprintf('   Costs:    %10.2f €\n', tcost);
fprintf('   Yield:    %10.2f €\n', trvnu);
fprintf(' Balance:    %10.2f €\n', trvnu - tcost);

namelist = {'Turboprop', 'Regional Jet','Single-Aisle Jet','Twin-Aisle Jet', 'Long Range Jet'};
for k = 1:nacft
    fprintf('\n\n%s flights:\n',namelist{k});
    ut = 0;
    ask= 0;
    for i  = 1:nodes
        for j = 1:nodes
            if flts(i,j,k) ~= 0
%                 fprintf('%s to %s: %i flights\n',tags{i}, tags{j}, flts(i,j,k));
                ut = ut + flts(i,j,k)*(distm(i,j)/acrft.vel(k) + ltotm(i,j,k));
                ask = ask + flts(i,j,k)*distm(i,j)*acrft.sts(k);
            end
        end
    end
    fprintf('Total %s Utilisation:  %5.2f h\n',namelist{k},ut);
    fprintf(' Unit %s Utilisation:  %5.2f h\n',namelist{k},ut/acrft.num(k));
    fprintf('Total %s Productivity: %7.0f ASK/week\n',namelist{k},ask);
    fprintf(' Unit %s Productivity: %7.0f ASK/week\n',namelist{k},ask/acrft.num(k));
end

% ---- Write matrices in Excel
fprintf('Writing Data in Excel\n')
xlswrite('input16_Ex1.xlsx',dpax,4,'B2:U21');
xlswrite('input16_Ex1.xlsx',cpax,4,'B24:U43');
xlswrite('input16_Ex1.xlsx',flts(:,:,1),4,'B46:U65');
xlswrite('input16_Ex1.xlsx',flts(:,:,2),4,'B68:U87');
xlswrite('input16_Ex1.xlsx',flts(:,:,3),4,'B90:U109');
fprintf('Data written\n')

%% ---- Plots
if plots_on == 1
    fprintf('Plotting maps...\n')
    % Graph objects for easier visualization 
    % Directed Graphs
%     Gdpax = digraph(dpax,tags);% Direct Passenger flows
%     Gcpax = digraph(cpax,tags);% Connecting Passenger flows
%     Gfpax = digraph(fpax,tags);% Direct Passenger flows
%     Gcap = digraph(cap,tags(1:nodes));% Capacity Passenger flows
%     Gfl1 = digraph(flts(:,:,1),tags);% Prop flights
%     Gfl2 = digraph(flts(:,:,2),tags);% Rjet flights
%     Gfl3 = digraph(flts(:,:,3),tags);% Sjet flights
    % Undirected Graphs
    Gdpax = graph(dpax,tags);% Direct Passenger flows
    Gcpax = graph(cpax,tags);% Connecting Passenger flows
    Gfpax = graph(fpax,tags);% Direct Passenger flows
    Gcap = graph(cap,tags(1:nodes));% Capacity Passenger flows
    Gfl1 = graph(flts(:,:,1),tags);% Prop flights
    Gfl2 = graph(flts(:,:,2),tags);% Rjet flights
    Gfl3 = graph(flts(:,:,3),tags);% Sjet flights
    
    maxval = max(max(fpax)); % Max value to scale graph edges
    lat = apdat.lat(1:nodes);
    lon = apdat.lon(1:nodes);
    ax = [-25 30 30 65]; % Axes - latitude and longitude
    
    % Direct, Connencting and total passenger flows - 3 subplots
    figure('units','normalized','outerposition',[0 0 0.5 0.7])
    plot_graph_map({Gcpax,Gdpax}, lat, lon, maxval, ax, {'b:','r'}) 
    xlabel('lat (deg)'); ylabel('lon (deg)'); title('Passenger Flows');
    legend('Connecting','Direct','Location','NorthWest')
    saveas(gca, 'Figures\passengers.eps','epsc');
    
    
    % Total Capacity Offered
    figure('units','normalized','outerposition',[0 0 0.33 0.5])
    maxval = max(max(cap)); % Max value to scale graph edges
    plot_graph_map(Gcap, lat, lon, maxval, ax, 'r')
    xlabel('lat (deg)'); ylabel('lon (deg)'); title('Capacity Offered');
    saveas(gca, 'Figures\capacity.eps','epsc');
    
    % Flights operated - 1 plot only
    maxval = max(max(max(flts))); % Max value to scale graph edges
    figure('units','normalized','outerposition',[0 0 0.5 0.7])
    plot_graph_map({Gfl1,Gfl2,Gfl3}, lat, lon, maxval, ax, {'-','k','r'})
    xlabel('lat (deg)'); ylabel('lon (deg)'); title('Operated Flights');
    legend('Turborop','Regional Jet','Single Aisle Jet','Twin Aisle Jet','Location','NorthWest')
    saveas(gca, 'Figures\flights.eps','epsc');
    
    fprintf('Maps Done\n')
end