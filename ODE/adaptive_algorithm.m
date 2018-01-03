function adaptive_algorithm(T,N,y0,TOL)
%Initialisation
%Calculate fitted values with T/N
h =  T/N;
y1_fitted_RK2 = RK2(h,y0);
y1_fitted_Euler = Euler(h,y0);
NORM = norm(y1_fitted_RK2 - y1_fitted_Euler);
endpt = 0;
c = 0; %keep count of optimal steps
count = zeros(N,1); %keep track of optimal steps
hopt = h;
yvals = zeros(N,1)

%while loop for adative algorithm
while endpt <= T
    h = hopt;
    if NORM <= TOL %start new time step with y1 as starting point
        y = y1_fitted_Euler; %new starting point
        c = c + 1; %increase count of steps
        y1_fitted_RK2 = RK2(h,y);
        y1_fitted_Euler = Euler(h,y);
        NORM = norm(y1_fitted_RK2 - y1_fitted_Euler);
        hopt = h*(TOL/NORM)^(1/2);
    else %if not, repeat time step
        y1_fitted_RK2 = RK2(h,y);
        y1_fitted_Euler = Euler(h,y);
        NORM = norm(y1_fitted_RK2 - y1_fitted_Euler);
        hopt = h*(TOL/NORM)^(1/2);
        h = hopt;
    end
    yvals(1,:) = y1_fitted_Euler;
    endpt = endpt + hopt;
    count(c) = hopt;
end
        
        
    
    
end



end

