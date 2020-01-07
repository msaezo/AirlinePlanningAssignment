
for i=1:20
    for j=1:20
        for t=1:24
            for k=1:3
                if dist(i,j) > Range(k) || runway(i) < runway_ac(k) ||  runway(j) < runway_ac(k)
                    if k==1
                        profit_turboprop(t,i,j) = -10000;
                    elseif k==2
                        profit_regjet(t,i,j) = -10000;
                    elseif k==3
                        profit_twinjet(t,i,j) = -10000;
                    end
%                 else  
%                     if i == j
%                         Fixed_costs = 0;
%                     else
%                         Fixed_costs = CX(k);
%                     end
%                     cost(i,j) = CT(k)*dist(i,j)/Speed(k)+CF(k)*f/1.5*dist(i,j) +  Fixed_costs;
%                     
%                     Dem_hour(t,i,j) = Dem(i,j)*Coef(t,i);
%                     if t==1 %X represents actual flow of passengers at those hours
%                         X(t,i,j) = Dem_hour(t,i,j)+Dem_hour(t+1,i,j);
%                     elseif t==2
%                         X(t,i,j) = Dem_hour(t,i,j)+Dem_hour(t+1,i,j)+Dem_hour(t-1,i,j);
%                     elseif t==24
%                         X(t,i,j) = Dem_hour(t,i,j)+Dem_hour(t-1,i,j)+Dem_hour(t-2,i,j);
%                     else
%                         X(t,i,j) = Dem_hour(t,i,j)+Dem_hour(t-1,i,j)+Dem_hour(t+1,i,j)+Dem_hour(t-2,i,j);
%                     end
%                     if X(t,i,j) > Cap(k) %capacity constraint
%                         X(t,i,j) = Cap(k);
%                     end
                    revenue(t,i,j) = [5.9*dist(i,j)^-0.76+0.043]*dist(i,j)*X(t,i,j);
                    if k==1
                        profit_turboprop(t,i,j) = revenue(t,i,j)-cost(i,j);
                    elseif k==2
                        profit_regjet(t,i,j) = revenue(t,i,j)-cost(i,j);
                    elseif k==3
                        profit_twinjet(t,i,j) = revenue(t,i,j)-cost(i,j);
                    end
                end
            end
        end
    end
end

                        
                        