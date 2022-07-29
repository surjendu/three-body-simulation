clear all;
clc, 
close all
tic
%input--------

%please provide input in the matlab function also
% ker_filename='ch(br81)3-ces-3body_non_com.csv';

m1=80.916;    % mass of the first ion
m2=80.916;    % mass of the second ion
m3=93.9230;   % mass of the third ion

q1=1;
q2=1;
q3=1;

x1=0.914949;   y1=1.601844;   z1=-0.04539;               %coordinate of the first ion
x2=-1.844898;  y2=-0.008681;   z2=-0.04567 ;           %coordinate of the second ion
x3=0.929983;   y3=-1.593171;   z3=-0.045499;    %coordinate of the third ion
		

% vx1=0;  vy1=0; vz1=0;  %initial velocities of the first ion  
% vx2=0;  vy2=0; vz2=0;  % initial velocities of the second ion 
% vx3=0;  vy3=0; vz3=0;  % initial velocities of the third ion 

tspan= 0:10000:300000; 

%---------------------------------------
angs2au=1.889726;      % angstrom to atomic unit (bohr radius)
au2fs=0.02418884;      % atomic unit to femtosecond
mpvsme=1822.888486424682;   % ratio of atomic mass and electron mass
hatoev=  27.211396; % taken from ksu webpage/ twice of ionization potential of H
m=mpvsme*[m1 m2 m3];

x=[x1,x2,x3];    % initial position in x
y=[y1,y2,y3];    % initial position in y
z=[z1,z2,z3];    % initial position in z


x_cm=(x(1)*m(1)+x(2)*m(2)+x(3)*m(3))/(m(1)+m(2)+m(3));  % center of mass position: x-values
y_cm=(y(1)*m(1)+y(2)*m(2)+y(3)*m(3))/(m(1)+m(2)+m(3));  % center of mass position: y-values
z_cm=(z(1)*m(1)+z(2)*m(2)+z(3)*m(3))/(m(1)+m(2)+m(3));  % center of mass position: z-values

% r1=[(x1-x_cm)*angs2au (y1-y_cm)*angs2au (z1-z_cm)*angs2au]; % initial position in CM coordiates of the first ion
% r2=[(x2-x_cm)*angs2au (y2-y_cm)*angs2au (z2-z_cm)*angs2au]; % initial position in CM coordiates of the second ion
% r3=[(x3-x_cm)*angs2au (y3-y_cm)*angs2au (z3-z_cm)*angs2au]; % initial position in CM coordiates of the third ion

% r1=[(x(1)-x_cm)*angs2au (y(1)-y_cm)*angs2au (z(1)-z_cm)*angs2au]; % initial position in CM coordiates of the first ion
% r2=[(x(2)-x_cm)*angs2au (y(2)-y_cm)*angs2au (z(2)-z_cm)*angs2au]; % initial position in CM coordiates of the second ion
% r3=[(x(3)-x_cm)*angs2au (y(3)-y_cm)*angs2au (z(3)-z_cm)*angs2au]; % initial position in CM coordiates of the third ion

r1=[(x(1)-x_cm) (y(1)-y_cm) (z(1)-z_cm)]; % initial position in CM coordiates of the first ion
r2=[(x(2)-x_cm) (y(2)-y_cm) (z(2)-z_cm)]; % initial position in CM coordiates of the second ion
r3=[(x(3)-x_cm) (y(3)-y_cm) (z(3)-z_cm)]; % initial position in CM coordiates of the third ion

%%
%random number within a sphere
%[https://www.mathworks.com/help/matlab/math/numbers-placed-randomly-within-volume-of-sphere.html]
%[https://www.mathworks.com/matlabcentral/fileexchange/67384-uniform-spherical-distribution-generator]
%[https://www.mathworks.com/matlabcentral/answers/93554-how-can-i-rotate-a-set-of-points-in-a-plane-by-a-certain-angle-about-an-arbitrary-point]
%[https://www.mathworks.com/matlabcentral/answers/437098-given-a-plane-defined-by-three-points-in-space-how-do-i-get-the-roll-and-pitch-angles-euler-angles]
close all;
n=1000; % number of random points
rvals = 2*rand(n,1)-1;
elevation = asin(rvals);
azimuth = 2*pi*rand(n,1);

radii_1 = (sqrt(r1(1)^2+r1(2)^2+r1(3)^2)); % distance in Ang of 1st hit from COM
radii_1=0.1*radii_1*(rand(n,1).^(1/3)); % 0.1 stands for 10% of this radii_1 we are varying for random distribution
[x1,y1,z1] = sph2cart(azimuth,elevation,radii_1);

x11=r1(1)+x1;
y11=r1(2)+y1;
z11=r1(3)+z1;

% figure
% plot3(x1,y1,z1,'.')
% axis equal

radii_2 = (sqrt(r2(1)^2+r2(2)^2+r2(3)^2)); % distance in Ang of 1st hit from COM
radii_2=0.1*radii_2*(rand(n,1).^(1/3)); % 0.1 stands for 10% of this radii_2 we are varying for random distribution
[x2,y2,z2] = sph2cart(azimuth,elevation,radii_2);

x22=r2(1)+x2;
y22=r2(2)+y2;
z22=r2(3)+z2;

% figure
% plot3(x2,y2,z2,'.')
% axis equal

x3=((m(1)+m(2)+m(3))*x_cm - x1.*m(1)-x2.*m(2))./ m(3);
y3=((m(1)+m(2)+m(3))*y_cm - y1.*m(1)-y2.*m(2))./ m(3);
z3=((m(1)+m(2)+m(3))*z_cm - z1.*m(1)-z2.*m(2))./ m(3);

x33=r3(1)+x3;
y33=r3(2)+y3;
z33=r3(3)+z3;
% 
% figure
% plot3(x3,y3,z3,'.')
% axis equal

%to check if we go the COM correct for each case
x_cm_f=(x11.*m(1)+x22.*m(2)+x33.*m(3))/(m(1)+m(2)+m(3));  % center of mass position: x-values
y_cm_f=(y11.*m(1)+y22.*m(2)+y33.*m(3))/(m(1)+m(2)+m(3));  % center of mass position: y-values
z_cm_f=(z11.*m(1)+z22.*m(2)+z33.*m(3))/(m(1)+m(2)+m(3));  % center of mass position: z-values

figure
plot3(x11,y11,z11,'.',x22,y22,z22,'x',x33,y33,z33,'o',x_cm_f,y_cm_f,z_cm_f,'*')
axis equal

%%

% r1=angs2au*[x11(1) y11(1) z11(1)];
% r2=angs2au*[x22(1) y22(1) z22(1)];
% r3=angs2au*[x33(1) y33(1) z33(1)]; 

% r1=angs2au*[x11(end) y11(end) z11(end)];
% r2=angs2au*[x22(end) y22(end) z22(end)];
% r3=angs2au*[x33(end) y33(end) z33(end)]; 

r1=angs2au*[x11 y11 z11];
r2=angs2au*[x22 y22 z22];
r3=angs2au*[x33 y33 z33]; 


%%
close all;
ke1_out=[];
ke2_out=[];
ke3_out=[];
ke_out=[];
pe_out=[];
parfor i=1:n
    loop_no=i
vx1=0;  vy1=0; vz1=0;  %initial velocities of the first ion  
vx2=0;  vy2=0; vz2=0;  % initial velocities of the second ion 
vx3=0;  vy3=0; vz3=0;  % initial velocities of the third ion 


ini=[r1(i,1) r1(i,2) r1(i,3) r2(i,1) r2(i,2) r2(i,3) r3(i,1) r3(i,2) r3(i,3) vx1 vy1 vz1 vx2 vy2 vz2 vx3 vy3 vz3]; % initial conditions


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

%actual angle to see in the end
abs_ang1 = (angtheta1(end)); % angle between P2, P3 (cosine)
abs_ang2 = (angtheta2(end)); % angle between P1, P2 (cosine)
abs_ang3 = (angtheta3(end)); % angle between P1, P3 (cosine)


angsum = abs_ang1+abs_ang2+abs_ang3; % angle sum to see if it peaks at 360 degrees...



ke1= hatoev*(0.5*px1.*px1/m(1)+0.5*py1.*py1/m(1)+0.5*pz1.*pz1/m(1)); %kinetic energy of the 1st ion
ke2= hatoev*(0.5*px2.*px2/m(2)+0.5*py2.*py2/m(2)+0.5*pz2.*pz2/m(2)); %kinetic energy of the 2nd ion
ke3= hatoev*(0.5*px3.*px3/m(3)+0.5*py3.*py3/m(3)+0.5*pz3.*pz3/m(3)); %kinetic energy of the 3rd ion
ke=ke1+ke2+ke3;

d=ke-ke1-ke2-ke3;    %must be zero
pe=hatoev*(1./r12+1./r23+1./r13); %potential energy 
peke=pe+ke;

% t=t*au2fs;
% plot(t,ke,'k',t,ke1,'r',t,ke2,'b:',t,ke3,'y',t,pe,'m',t,peke,'c:',t,d,'g:','LineWidth',3)
% title(ker_filename)
% xlabel('time delay / fs');
% ylabel('ker / ev');
% legend('ke1+ke2+ke3','ke1','ke2','ke3','pe','pe+ke','Location','best');
% legend('boxoff');

ker=[t ke ke1 ke2 ke3 pe peke (r12./angs2au) ];  % t in fs, ke,ke1,ke2,pe in eV, r in Angstorm
% dlmwrite(ker_filename, ker);

 ke1_out=[ke1_out;ke1(end)];
 ke2_out=[ke2_out;ke2(end)];
 ke3_out=[ke3_out;ke3(end)];
 ke_out=[ke_out;ke(end)];
 pe_out=[pe_out;pe(1)];
 
 %pause(5)
end

%%
% need to rewite th
    k1=[0;ke1_out;20]; 
    k2=[0;ke2_out;20];
    k3=[0;ke3_out;20];
    k=[0;ke_out;20];
    
   ke_raw=[k1 k2 k3 k];
   dlmwrite('ke_raw.csv',ke_raw);
   nbins = 100;
   [countsk1,centersk1]=hist(k1,nbins);
   [countsk2,centersk2]=hist(k2,nbins);
   [countsk3,centersk3]=hist(k3,nbins);
   [countsk,centersk]=hist(k,nbins);
      
   ke_hist=[centersk1' countsk1' centersk2' countsk2' centersk3' countsk3' centersk' countsk'];
   dlmwrite('ke_hist.csv',ke_hist);
   
   figure 
   plot(centersk1,countsk1)
   hold on
   plot(centersk2,countsk2)
   plot(centersk3,countsk3)
   plot(centersk,countsk)
   
toc