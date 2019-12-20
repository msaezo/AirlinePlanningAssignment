function[] = Constraints(cplex,N,DV,AC_number,Demand_2019,g,Seats,LF,Speed,Distance,TAT,BT,Max_range,Runway_length,Rwy_req)
% Authors (Group 31):  Adrián Chaves García-Mascaraque (5077982)
%                      Guillermo Lara Juan (5169860)
%                      Miguel Ángel Sáez Ortuno (4541235)
%% Constraint 1: Demand verification
k = 1;
for i = 1:N
    for j = 1:N
        C1 = zeros(1,DV);
        C1(varindex_3(1,i,j,k,N,AC_number)) = 1;
        C1(varindex_3(2,i,j,k,N,AC_number)) = 1;
        RHS = Demand_2019(i,j);
        cplex.addRows(-inf,C1,RHS,sprintf('Constraint1%d_%d',i,j));
    end
end

%% Constraint 2: Making sure that w is 0 if the origin or the destination 
% is the hub

for i = 1:N
    for j = 1:N
        C2 = zeros(1,DV);
        C2(varindex_3(1,i,j,k,N,AC_number)) = 1;
        RHS = Demand_2019(i,j)*g(i)*g(j);
        cplex.addRows(-inf,C2, RHS,sprintf('Constraint2%d_%d',i,j));
    end
end

%% Constraint 3: Capacity
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
        RHS = 0;
        cplex.addRows(-inf, C3, RHS,sprintf('Constraint3%d_%d',i,j));
    end
end

%% Constraint 4: Continuity

for k = 1:AC_number
    for i = 1:N
        C4 = zeros(1,DV);
        for j = 1:N
            if i ~= j
                C4(varindex_3(3,i,j,k,N,AC_number)) = 1;
                C4(varindex_3(3,j,i,k,N,AC_number)) = -1;
            end
        end
        
        cplex.addRows(0,C4,0,sprintf('Constraint4%d_%d',k,i));
    end
end

%% Constraint 5: Aircraft productivity

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
    cplex.addRows(-inf,C5,0,sprintf('Constraint5%d',k));
end

%% Constraint 6: Range verification
for k = 1:AC_number 
    for i = 1:N
        for j = 1:N
            
            if i ~=j
            C6 = zeros(1,DV);
            C6(varindex_3(3,i,j,k,N,AC_number)) = 1;
            RHS = 0;
            if Distance(i,j) <= Max_range(k)
                RHS = 1000000;                
            end
            cplex.addRows(-inf, C6, RHS,sprintf('Constraint6%d_%d_%d',k,i,j));
            end
        end
    end
end

%% Constraint 7: Runway verification
for k = 1:AC_number 
    for i = 1:N
        for j = 1:N
            
            if i ~=j
            C7 = zeros(1,DV);
            C7(varindex_3(3,i,j,k,N,AC_number)) = 1;
            RHS = 0;
            if Runway_length(i) >= Rwy_req(k) && Runway_length(j) >= Rwy_req(k)
                RHS = 1000000;                
            end
            cplex.addRows(-inf, C7, RHS,sprintf('Constraint7%d_%d_%d',k,i,j));
            end
        end
    end
end
end