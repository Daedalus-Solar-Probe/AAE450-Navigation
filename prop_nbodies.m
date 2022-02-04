function [t,y,f] = prop_nbodies(y0,tspan,mu)

opts = odeset('RelTol',1e-6);
[t,y] = ode45(@(t,y)odefun(t,y,mu),tspan,y0,opts);

f = figure();

% plot central body
plot3(0,0,0,'o','MarkerSize',10)
hold on

for n = 0:6:length(y0)-1

    Rn = y(:,(1:3)+n);

    % object end point
    plot3(Rn(end,1),Rn(end,2),Rn(end,3),'.','MarkerSize',10)
    
    % object trajectory
    plot3(Rn(:,1),Rn(:,2),Rn(:,3),'-')

end

hold off

lim = max(max(y))*1.25;

grid on
xlim([-lim lim])
ylim([-lim lim])
zlim([-lim lim])

xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")


%% odefun
function dydt = odefun(~,y,mu)

dydt = zeros(size(y));

for i = 0:6:length(y)-1
    r = y((1:3)+i);
    v = y((4:6)+i);
    a = -mu*r/norm(r)^3;

    dydt((1:3)+i) = v;
    dydt((4:6)+i) = a;
end

end
end