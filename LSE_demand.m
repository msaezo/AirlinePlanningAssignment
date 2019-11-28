function[LSE] = LSE_demand(k,b1,b2,b3,pop_2014,GDP_2014,Demand_2014,Distance,N,f_cost,Mode)
    if Mode == 0
        %LSE = 0;
        for i = 1:N
            for j = 1:N
                if i == j 
                    Demand_est(i,j) = 0;
                else
                    Demand_est(i,j) = k*((pop_2014(i)*pop_2014(j))^b1*(GDP_2014(i)*GDP_2014(j))^b2)...
                        /(f_cost*Distance(i,j))^b3;
                end
                %LSE = LSE + sqrt(abs(Demand_2014(i,j)^2-Demand_est(i,j)^2));
            end
        end
       LSE_mat = sqrt(abs(Demand_est.^2 - Demand_2014^2));
       LSE = sum( LSE_mat , 'all' ); 
       LSE = LSE/N;

    else
        LSE = 0;
        for i = 1:N
            for j = 1:N
                if i == j 
                    Demand_est(i,j) = 0;
                else
                    Demand_est(i,j) = k*((pop_2014(i)*pop_2014(j))^b1*(GDP_2014(i)*GDP_2014(j))^b2)...
                        /(f_cost*Distance(i,j))^b3;
                end
                %LSE = LSE + sqrt(abs(Demand_2014(i,j)^2-Demand_est(i,j)^2));
            end
        end
       LSE_mat = abs(Demand_est.^2 - Demand_2014^2);
       LSE = sum(LSE_mat,'all'); 
       LSE = LSE/N;
       LSE = Demand_est;
    end
end