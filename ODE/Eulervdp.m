function yvals = Eulervdp(T, N, y1_0, y2_0)
%Euler(100,20000,2,0)
%initialisation
%yvals is the matrix where we store our approximation poitns
%range of t; divided by steps of h = T/N
yvals = zeros(N+1,2);
y = [y1_0, y2_0];
yvals(1,:) = y;
h = T/N;
t_range = linspace(0,T,N+1);
k = 0;
%calculating the approximation points
for n = 1:N
    y = y + h*[y(2), (1 - y(1)^2)*y(2) - y(1)];
    yvals(n+1,:) = y;
    k = k + 1;
end
%%plot
plot(t_range, yvals(:,1))

%Uncomment for interactive plot
%curve = animatedline('Color','k');
%set(gca, 'XLim', [0 T], 'YLim', [-10,10]);
%grid on;
%for i = 1:N+1
%   addpoints(curve, t_range(i),yvals(i,1));
%   drawnow
%end

end