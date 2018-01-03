%RK2 Implementation on an ODE problem with Van der Pol Oscillator
%input h, this only does it for 1 time step
%--------------------------------------------------------------------------------
function y1 = RK2(h,y)
%initialisation
func = @(y) [y(2), (1-y(1)^(2))*y(2) - y(1)];
y1 = [0,0];

%Calculate coefficients and then calculate y(t1) = y(t0 + h) ~ y1
k1 = func(y);
k2 = func(y + h*0.5*k1);
y1 = y + h*(k2);


end