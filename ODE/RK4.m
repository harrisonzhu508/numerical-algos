%RK4 Implementation on an ODE problem with Van der Pol Oscillator
%--------------------------------------------------------------------------------
function yvals = RK4(T,N,y0)
%initialisation
h = T/N;
func = @(y) [y(2), (1-y(1)^(2))*y(2) - y(1)];
t_range = linspace(0,T,N);
yvals = zeros(N,2);

%Calculate coefficients and then calculate y(t1) = y(t0 + h) ~ y1
yvals(1,:) = y0;
for i = 2:N
    y = yvals(i-1,:);
    k1 = func(y);
    k2 = func(y + h*0.5*k1);
    k3 = func(y + h*0.5*k2);
    k4 = func(y + h*k3);
    yvals(i,:) = y + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

%%plot
plot(t_range, transpose(yvals(:,1)));

%uncomment for interactive plot
%curve = animatedline('Color','k');
%set(gca, 'XLim', [0 T], 'YLim', [-10,10]);
%grid on;
%for i = 1:N+1
%   addpoints(curve, t_range(i),yvals(i,1));
%   drawnow
%end



end