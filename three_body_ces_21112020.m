clear all;
clc, 
close all
%input--------


m1=80.916;    % mass of the first ion
m2=80.916;    % mass of the second ion
m3=93.9240;   % mass of the third ion

q1=1;
q2=1;
q3=1;

x1=-2.42796600;   y1=-0.59922600;   z1= 0.03354500 ;                %coordinate of the first ion
x2=-0.29933500;   y2= 0.99950400;  z2=   -0.12871700;        %coordinate of the second ion
x3=2.47595900;    y3=-0.53573400;   z3=-0.06509600;    %coordinate of the third ion
		 


vx1=0;  vy1=0; vz1=0;  %initial velocities of the first ion  
vx2=0;  vy2=0; vz2=0;  % initial velocities of the second ion 
vx3=0;  vy3=0; vz3=0;  % initial velocities of the third ion 

%dt=0.2418884; %time step size fs
tspan= 0:300000; 

%---------------------------------------
angs2au=1.889726;      % angstrom to atomic unit (bohr radius)
au2fs=0.02418884;      % atomic unit to femtosecond
mpvsme=1822.888486424682;   % ratio of atomic mass and electron mass
hatoev=  27.211396; % taken from ksu webpage/ twice of ionization potential of H
m=[m1*mpvsme m2*mpvsme m3*mpvsme];

x=[x1,x2,x3];    % initial position in x
y=[y1,y2,y3];    % initial position in y
z=[z1,z2,z3];    % initial position in z


x_cm=(x(1)*m(1)+x(2)*m(2)+x(3)*m(3))/(m(1)+m(2)+m(3));  % center of mass position: x-values
y_cm=(y(1)*m(1)+y(2)*m(2)+y(3)*m(3))/(m(1)+m(2)+m(3));  % center of mass position: y-values
z_cm=(z(1)*m(1)+z(2)*m(2)+z(3)*m(3))/(m(1)+m(2)+m(3));  % center of mass position: z-values

r1=[(x1-x_cm)*angs2au (y1-y_cm)*angs2au (z1-z_cm)*angs2au]; % initial position in CM coordiates of the first ion
r2=[(x2-x_cm)*angs2au (y2-y_cm)*angs2au (z2-z_cm)*angs2au]; % initial position in CM coordiates of the second ion
r3=[(x3-x_cm)*angs2au (y3-y_cm)*angs2au (z3-z_cm)*angs2au]; % initial position in CM coordiates of the third ion

ini=[r1(1) r1(2) r1(3) r2(1) r2(2) r2(3) r3(1) r3(2) r3(3) vx1 vy1 vz1 vx2 vy2 vz2 vx3 vy3 vz3]; % initial conditions


options= odeset('AbsTol',1e-12,'RelTol',1e-12);


F = @(t,r) [r(10);r(11);r(12);r(13);r(14);r(15);r(16);r(17);r(18);...
    
((q1*q2/m(1))/(((r(1)-r(4)).^2 + (r(2)-r(5)).^2 + (r(3)-r(6)).^2).^(3/2)))*(r(1)-r(4)) + ((q1*q3/m(1))/(((r(1)-r(7)).^2 + (r(2)-r(8)).^2 + (r(3)-r(9)).^2).^(3/2)))*(r(1)-r(7));...
((q1*q2/m(1))/(((r(1)-r(4)).^2 + (r(2)-r(5)).^2 + (r(3)-r(6)).^2).^(3/2)))*(r(2)-r(5)) + ((q1*q3/m(1))/(((r(1)-r(7)).^2 + (r(2)-r(8)).^2 + (r(3)-r(9)).^2).^(3/2)))*(r(2)-r(8));...
((q1*q2/m(1))/(((r(1)-r(4)).^2 + (r(2)-r(5)).^2 + (r(3)-r(6)).^2).^(3/2)))*(r(3)-r(6)) + ((q1*q3/m(1))/(((r(1)-r(7)).^2 + (r(2)-r(8)).^2 + (r(3)-r(9)).^2).^(3/2)))*(r(3)-r(9));...

((q1*q2/m(2))/(((r(1)-r(4)).^2 + (r(2)-r(5)).^2 + (r(3)-r(6)).^2).^(3/2)))*(r(4)-r(1)) + ((q3*q2/m(2))/(((r(4)-r(7)).^2 + (r(5)-r(8)).^2 + (r(6)-r(9)).^2).^(3/2)))*(r(4)-r(7));...
((q1*q2/m(2))/(((r(1)-r(4)).^2 + (r(2)-r(5)).^2 + (r(3)-r(6)).^2).^(3/2)))*(r(5)-r(2)) + ((q3*q2/m(2))/(((r(4)-r(7)).^2 + (r(5)-r(8)).^2 + (r(6)-r(9)).^2).^(3/2)))*(r(5)-r(8));...
((q1*q2/m(2))/(((r(1)-r(4)).^2 + (r(2)-r(5)).^2 + (r(3)-r(6)).^2).^(3/2)))*(r(6)-r(3)) + ((q3*q2/m(2))/(((r(4)-r(7)).^2 + (r(5)-r(8)).^2 + (r(6)-r(9)).^2).^(3/2)))*(r(6)-r(9));...

((q1*q3/m(3))/(((r(1)-r(7)).^2 + (r(2)-r(8)).^2 + (r(3)-r(9)).^2).^(3/2)))*(r(7)-r(1)) + ((q3*q2/m(3))/(((r(4)-r(7)).^2 + (r(5)-r(8)).^2 + (r(6)-r(9)).^2).^(3/2)))*(r(7)-r(4));...
((q1*q3/m(3))/(((r(1)-r(7)).^2 + (r(2)-r(8)).^2 + (r(3)-r(9)).^2).^(3/2)))*(r(8)-r(2)) + ((q3*q2/m(3))/(((r(4)-r(7)).^2 + (r(5)-r(8)).^2 + (r(6)-r(9)).^2).^(3/2)))*(r(8)-r(5));...
((q1*q3/m(3))/(((r(1)-r(7)).^2 + (r(2)-r(8)).^2 + (r(3)-r(9)).^2).^(3/2)))*(r(9)-r(3)) + ((q3*q2/m(3))/(((r(4)-r(7)).^2 + (r(5)-r(8)).^2 + (r(6)-r(9)).^2).^(3/2)))*(r(9)-r(6))];

[t,r]=ode45(F, tspan, ini, options); % solving runge kutta 4th order

x1  = r(:,1);   y1  = r(:,2);   z1  = r(:,3);      % position of 1st ion
vx1 = r(:,10);   vy1 = r(:,11);   vz1 = r(:,12);      % velocity of 1st ion
px1 = m(1)*vx1; py1 = m(1)*vy1; pz1 = m(1)*vz1;    % momentum of 1st ion
p1=sqrt(px1.*px1+py1.*py1+pz1.*pz1);  % momentum magnitude of 1st ion

x2  = r(:,4);   y2  = r(:,5);    z2  = r(:,6);    % position of 2nd ion
vx2 = r(:,13);   vy2 = r(:,14);   vz2 = r(:,15);    % velocity of 2nd ion
px2 = m(2)*vx2; py2 = m(2)*vy2;  pz2 = m(2)*vz2;    % momentum of 2nd ion
p2 =sqrt(px2.*px2+py2.*py2+pz2.*pz2);  % momentum magnitude of 2nd ion

x3  = r(:,7);   y3  = r(:,8);    z3  = r(:,9);    % position of 3rd ion
vx3 = r(:,16);   vy3 = r(:,17);   vz3 = r(:,18);    % velocity of 3rd  ion
px3 = m(3)*vx3; py3 = m(3)*vy3;  pz3 = m(3)*vz3;    % momentum of 3rd  ion
p3 =sqrt(px3.*px3+py3.*py3+pz3.*pz3);  % momentum magnitude of 2nd ion


r12=((x1-x2).^2+(y1-y2).^2+(z1-z2).^2).^(1/2); % internuclear distance in au
r23=((x2-x3).^2+(y2-y3).^2+(z2-z3).^2).^(1/2); % internuclear distance in au
r13=((x1-x3).^2+(y1-y3).^2+(z1-z3).^2).^(1/2); % internuclear distance in au

angtheta1=(180/pi)*acos((px2.*px3+py2.*py3+pz2.*pz3)./(p2.*p3)); % angle between P2,P3
angtheta2=(180/pi)*acos((px1.*px2+py1.*py2+pz1.*pz2)./(p1.*p2)); % angle between P1,P2
angtheta3=(180/pi)*acos((px1.*px3+py1.*py3+pz1.*pz3)./(p1.*p3)); % angle between P2,P3

ang1 = cos(angtheta1(end)*pi/180); % angle between P2, P3 (cosine)
ang2 = cos(angtheta2(end)*pi/180); % angle between P2, P3 (cosine)
ang3 = cos(angtheta3(end)*pi/180); % angle between P2, P3 (cosine)

abs_ang1 = (angtheta1(end)) % angle between P2, P3 (cosine)
abs_ang2 = (angtheta2(end)) % angle between P1, P2 (cosine)
abs_ang3 = (angtheta3(end)) % angle between P1, P3 (cosine)


angsum = abs_ang1+abs_ang2+abs_ang3 % angle sum to see if it peaks at 360 degrees...



ke1= hatoev*(0.5*px1.*px1/m(1)+0.5*py1.*py1/m(1)+0.5*pz1.*pz1/m(1)); %kinetic energy of the 1st ion
ke2= hatoev*(0.5*px2.*px2/m(2)+0.5*py2.*py2/m(2)+0.5*pz2.*pz2/m(2)); %kinetic energy of the 2nd ion
ke3= hatoev*(0.5*px3.*px3/m(3)+0.5*py3.*py3/m(3)+0.5*pz3.*pz3/m(3)); %kinetic energy of the 3rd ion
ke=ke1+ke2+ke3;

d=ke-ke1-ke2-ke3;    %must be zero
pe=hatoev*(1./r12+1./r23+1./r13); %potential energy 
peke=pe+ke;

t=t*au2fs;
plot(t,ke,'k',t,ke1,'r',t,ke2,'b:',t,ke3,'y',t,pe,'m',t,peke,'c:',t,d,'g:','LineWidth',3)
title(ker_filename)
xlabel('time delay / fs');
ylabel('ker / ev');
legend('ke1+ke2+ke3','ke1','ke2','ke3','pe','pe+ke','Location','best');
legend('boxoff');

 
 ke1(end)
 ke2(end)
 ke3(end)
 ke(end)
 pe(1)
