%Demand forecast

x0 = [];
lb = [];
ub = [];
f=@(k,b1,b2,b3) LSE_demand(k,b1,b2,b3);

options=optimset('Display','iter','Algorithm','sqp');
fun=@(x) f(x(1),x(2),x(3),x(4)); % [deg,deg,m,deg,m,m]

% Options for the optimization
options.Display         = 'iter-detailed';
options.Algorithm       = 'sqp';
options.FunValCheck     = 'off';
options.DiffMinChange   = 1e-4;         % Minimum change while gradient searching
options.DiffMaxChange   = 1e-3;         % Maximum change while gradient searching
options.TolCon          = 1e-6;         % Maximum difference between two subsequent constraint vectors [c and ceq]
options.TolFun          = 1e-5;         % Maximum difference between two subsequent objective value
options.TolX            = 1e-6;         % Maximum difference between two subsequent design vectors
options.PlotFcns        = {@optimplotx,@optimplotfunccount,@optimplotfval,@optimplotconstrviolation,@optimplotstepsize,@optimplotfirstorderopt};
options.MaxIter         = 500;           % Maximum iterations

tic;
[x,fval,exitflag,output]=fmincon(fun,x0,[],[],[],[],lb,ub,[],options);
toc;