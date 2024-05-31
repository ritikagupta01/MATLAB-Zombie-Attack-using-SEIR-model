N = 500; beta = 0.5; alpha = 0.2;
w = 0.3; K = 10; tr = 10; gamma = 0.05;

S0=199;E0=300;I0=100;R0=0;

tspan = 0:0.1:30;

dydt = @(t,y) [(beta*N - (alpha+w)*y(1) - ((beta-alpha)/K)*(y(1)*y(1)));
    (((beta-alpha)/K)*(y(1)*y(1)) - y(2)*(tr*y(3)+alpha));
    (tr*y(2)*y(3) - (gamma+alpha)*y(3));
    (gamma*y(3) + w*y(1) - alpha*y(4))];

ode45(dydt,tspan,[S0 E0 I0 R0]);

legend('S(t)','E(t)','I(t)','R(t)')
title('Zombie Outbreak Model')
xlabel('Time in days')
ylabel('Population')



