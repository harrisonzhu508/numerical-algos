%Euler 1 step approximation function
%------------------------------------------------------------------------
function yvals = Euler(h,y)
%initialisation
y1 = [0,0]

%calculating the approximation point
y1 = y + h*[y(2), (1 - y(1)^2)*y(2) - y(1)];

end
