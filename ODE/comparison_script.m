%comparison

hold on
length(adaptive_algorithm(10,100,[2,0],0.01))
Eulervdp(10,100,2,0)
RK4(10,100,[2,0])
legend('adapt','Euler','RK4')
hold off
