%Euler 1 step approximation function
%------------------------------------------------------------------------
function y1 = Euler(h,y0)
%calculating the approximation point
y1 = y0 + h*[y0(2), (1 - y0(1)^2)*y0(2) - y0(1)];

end
