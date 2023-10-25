clear; close all;
stamp = 'plotted by plot\_situations\_for\_successive\_crossings.m';
%% construct shapes of HCS
x_HCS = linspace(-10,10,2001);
y_HCS = linspace(-10,10,2001);
[xx_HCS,yy_HCS] = meshgrid(x_HCS,y_HCS);
% shape1: plane
zz_HCS1 = zeros(size(xx_HCS));
% shape2: curved plane with 2.5 period
lambda2_HCS = 8;
zz_HCS2 = 3*sin(2*pi/lambda2_HCS*xx_HCS);
% shape3: curved plane with 2 period
lambda3_HCS = 10;
zz_HCS3 = 3*cos(2*pi/lambda3_HCS*xx_HCS);
%% construct orbits of PSP
x_PSP = linspace(-10,10,2001);
y_PSP = zeros(size(x_PSP));
% orbit1: curved line with 2.5 period
lambda1_PSP = 8;
z_PSP1 = -4*sin(2*pi/lambda1_PSP*x_PSP);
% orbit2: straight line
z_PSP2 = zeros(size(x_PSP));
% orbit3: curved line with 2.5 period
lambda3_PSP = 8;
z_PSP3 = -4*sin(2*pi/lambda3_PSP*x_PSP);
%% determine crossing locations
% situation1
x_cross1 = [-8,-4,0,4,8];
% situation2
x_cross2 = [-8,-4,0,4,8];
% situation3
x_cross3 = [-8.68,-5.00,-0.86,3.45,7.82];
z_cross3 = [2.02,-3.00,2.57,-1.69,0.60];
%% plot figure
figure();
LineWidth = 3;
FontSize = 20;
QuivLength = 2;
HeadSize = 5;
Azimuth = -1;
Elevation = 5;
%% situation 1: flapping only
% subplot(1,3,1)
subplot(3,1,1)
surf(xx_HCS,yy_HCS,zz_HCS1,xx_HCS); hold on;
shading interp; axis equal; axis off; colormap winter;
plot3(x_PSP,y_PSP,z_PSP1,'k','LineWidth',LineWidth); hold on;
for i_cross = 1 : 5
    if i_cross == 1 || i_cross == 3 || i_cross == 5
        q = quiver3(x_cross1(i_cross),0,0,0,0,QuivLength,'LineWidth',LineWidth);
    else
        q = quiver3(x_cross1(i_cross),0,0,0,0,-QuivLength,'LineWidth',LineWidth);
    end
    set(q,'MaxHeadSize',HeadSize,'Color','r');
end
view([-1,5]);
% text(-5,0,8,'(a) Flapping only','FontSize',FontSize);
text(-5,0,5,'(a) Flapping only','FontSize',FontSize);
%% situation 2: folds only
% subplot(1,3,2)
subplot(3,1,2)
surf(xx_HCS,yy_HCS,zz_HCS2); hold on;
shading interp; axis equal; axis off; colormap winter;
plot3(x_PSP,y_PSP,z_PSP2,'k','LineWidth',LineWidth); hold on;
for i_cross = 1 : 5
    x_sub = x_cross2(i_cross);
    norm_quiv = sqrt((2*pi/lambda2_HCS*cos(2*pi/lambda2_HCS*x_sub))^2+1);
    x_quiv = -2*pi/lambda2_HCS*cos(2*pi/lambda2_HCS*x_sub)/norm_quiv;
    z_quiv = 1/norm_quiv;
    if i_cross == 1 || i_cross == 3 || i_cross == 5
        q = quiver3(x_sub,0,0,QuivLength*x_quiv,0,QuivLength*z_quiv,'LineWidth',LineWidth);
    else
        q = quiver3(x_sub,0,0,-QuivLength*x_quiv,0,-QuivLength*z_quiv,'LineWidth',LineWidth);
    end
    set(q,'MaxHeadSize',HeadSize,'Color','r');
end
view([-10,5]);
% text(-4,0,8,'(b) Folds only','FontSize',FontSize);
text(-4,0,5,'(b) Folds only','FontSize',FontSize*1.2);
%% situation 3: flapping and folds
% subplot(1,3,3)
subplot(3,1,3)
surf(xx_HCS,yy_HCS,zz_HCS3); hold on;
shading interp; axis equal; axis off; colormap winter;
plot3(x_PSP,y_PSP,z_PSP3,'k','LineWidth',LineWidth); hold on;
for i_cross = 1 : 5
    x_sub = x_cross3(i_cross);
    z_sub = z_cross3(i_cross);
    norm_quiv = sqrt((2*pi/lambda3_HCS*sin(2*pi/lambda3_HCS*x_sub))^2+1);
    x_quiv = 2*pi/lambda3_HCS*sin(2*pi/lambda3_HCS*x_sub)/norm_quiv;
    z_quiv = 1/norm_quiv;
    if i_cross == 1 || i_cross == 3 || i_cross == 5
        q = quiver3(x_sub,0,z_sub,QuivLength*x_quiv,0,QuivLength*z_quiv,'LineWidth',LineWidth);
    else
        q = quiver3(x_sub,0,z_sub,-QuivLength*x_quiv,0,-QuivLength*z_quiv,'LineWidth',LineWidth);
    end
    set(q,'MaxHeadSize',HeadSize,'Color','r');
end
view([-10,5]);
% text(-8,0,8,'(c) Flapping and Folds','FontSize',FontSize);
text(-8,0,5.5,'(c) Flapping and Folds','FontSize',FontSize);
set(gca,'LineWidth',LineWidth/2,'FontSize',FontSize);
%% stamp
% legend('Current Sheet','Orbit of Satellite','CS Normal Direction','NumColumns',3);
legend('Current Sheet','Orbit of Satellite','CS Normal Direction','NumColumns',1);
text(0,0,-12,stamp);



