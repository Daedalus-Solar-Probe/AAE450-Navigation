%% Clear
clear
clc

%%

start = juliandate(datetime);
stop = start + 365.25;

T_eph = linspace(start,stop,1000);

r_ecl = zeros(3,8,length(T_eph));
r_icrf = zeros(size(r_ecl));

% tic
% for i = 1:length(T_eph)
%     [r_ecl(:,:,i),r_icrf(:,:,i)] = ephemeris(T_eph(i));
% end
% toc

tic
[x_ecl,y_ecl,z_ecl,x_icrf,y_icrf,z_icrf] = ephemeris(T_eph);
toc

%% Plot J2000 Ecliptic

r_mercury = [x_icrf(1,:); y_icrf(1,:); z_icrf(1,:)];
r_venus = [x_icrf(2,:); y_icrf(2,:); z_icrf(2,:)];
r_earth = [x_icrf(3,:); y_icrf(3,:); z_icrf(3,:)];
r_mars = [x_icrf(4,:); y_icrf(4,:); z_icrf(4,:)];

figure(1)
plot3(0,0,0,'ko','MarkerSize',20,'MarkerFaceColor','y')
hold on

plot3(r_mercury(1,:),r_mercury(2,:),r_mercury(3,:),'k')
plot3(r_venus(1,:),r_venus(2,:),r_venus(3,:),'g')
plot3(r_earth(1,:),r_earth(2,:),r_earth(3,:),'b')
plot3(r_mars(1,:),r_mars(2,:),r_mars(3,:),'r')
hold off

grid on
xlim([-1.5 1.5])
ylim([-1.5 1.5])
zlim([-1.5 1.5])


% %% J2000 Ecliptic
% 
% r_mercury = reshape(r_ecl(:,1,:),3,[]);
% r_venus = reshape(r_ecl(:,2,:),3,[]);
% r_earth = reshape(r_ecl(:,3,:),3,[]);
% r_mars = reshape(r_ecl(:,4,:),3,[]);
% 
% figure(1)
% plot3(0,0,0,'ko','MarkerSize',20,'MarkerFaceColor','y')
% hold on
% 
% plot3(r_mercury(1,:),r_mercury(2,:),r_mercury(3,:),'r')
% plot3(r_venus(1,:),r_venus(2,:),r_venus(3,:),'r')
% plot3(r_earth(1,:),r_earth(2,:),r_earth(3,:),'r')
% plot3(r_mars(1,:),r_mars(2,:),r_mars(3,:),'r')
% hold off
% 
% grid on
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])
% zlim([-1.5 1.5])
% 
% %% J2000 ICRF
% 
% r_mercury = reshape(r_icrf(:,1,:),3,[]);
% r_venus = reshape(r_icrf(:,2,:),3,[]);
% r_earth = reshape(r_icrf(:,3,:),3,[]);
% r_mars = reshape(r_icrf(:,4,:),3,[]);
% 
% figure(2)
% plot3(0,0,0,'ko','MarkerSize',20,'MarkerFaceColor','y')
% hold on
% 
% plot3(r_mercury(1,:),r_mercury(2,:),r_mercury(3,:),'r')
% plot3(r_venus(1,:),r_venus(2,:),r_venus(3,:),'r')
% plot3(r_earth(1,:),r_earth(2,:),r_earth(3,:),'r')
% plot3(r_mars(1,:),r_mars(2,:),r_mars(3,:),'r')
% hold off
% 
% grid on
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])
% zlim([-1.5 1.5])
% 
% 
% %% MATLAB Ephemeris
% 
% figure(3)
% 
% r_mercury1 = planetEphemeris(T_eph','Sun','Mercury','405','AU')';
% r_venus1 = planetEphemeris(T_eph','Sun','Venus','405','AU')';
% r_earth1 = planetEphemeris(T_eph','Sun','Earth','405','AU')';
% r_mars1 = planetEphemeris(T_eph','Sun','mars','405','AU')';
% 
% figure(3)
% plot3(0,0,0,'ko','MarkerSize',20,'MarkerFaceColor','y')
% hold on
% 
% plot3(r_mercury1(1,:),r_mercury1(2,:),r_mercury1(3,:),'r')
% plot3(r_venus1(1,:),r_venus1(2,:),r_venus1(3,:),'r')
% plot3(r_earth1(1,:),r_earth1(2,:),r_earth1(3,:),'r')
% plot3(r_mars1(1,:),r_mars1(2,:),r_mars1(3,:),'r')
% hold off
% 
% grid on
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])
% zlim([-1.5 1.5])
% 
% 
% %% difference
% 
% r_mercury1 - r_mercury1
% 
% 
% % r = coes2svec(a(1),e(1),omega(1),OMEGA(1),I(1),M(1),0,3.986004418e5)
