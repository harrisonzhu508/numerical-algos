function yvals = adaptive_algorithm(T,N,y0,TOL)
%Initialisation
%Calculate fitted values with T/N
h = T/N;
y1_fitted_RK2 = [0,0];
y1_fitted_Euler = [0,0];
NORM = 3;
hopt = 0;
c = 1; %keep count of optimal steps
count = zeros(N,1); %keep track of optimal steps
yvals = zeros(N,2);
yvals(1,:) = y0;
y = y0;

%while loop for adaptive algorithm
while count(c) <= T
    if NORM <= TOL %start new time step with y1 as starting point
        %updates
        c = c + 1; %increase count of steps
        yvals(c,:) = y1_fitted_Euler; %update new values
        count(c) = hopt + count(c-1);
        %new start
        y = y1_fitted_Euler; %new starting point
        y1_fitted_RK2 = RK2(h,y);
        y1_fitted_Euler = Euler(h,y);
        NORM = norm(y1_fitted_RK2 - y1_fitted_Euler);
        hopt = h*(TOL/NORM)^(1/2);
        h = hopt;
    else %if not, repeat time step
        y1_fitted_RK2 = RK2(h,y);
        y1_fitted_Euler = Euler(h,y);
        NORM = norm(y1_fitted_RK2 - y1_fitted_Euler);
        hopt = h*(TOL/NORM)^(1/2);
        h = hopt;
    end
end
plot(count, transpose(yvals(:,1)));
end
        
        


