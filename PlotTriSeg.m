% subroutine plot TriSeg circles as a function of the time
function PlotTriSeg(xm_LV1,xm_SEP1,xm_RV1,ym1,t)
 for i = 1: length(t)
 xm_LV  = xm_LV1(i) ; % LV heart geometry variable, cm
 xm_SEP = xm_SEP1(i); % septum heart geometry variable, cm
 xm_RV  = xm_RV1(i); % RV heart geometry variable, cm
 ym    = ym1(i); % Heart geometry variable, cm

% xm_LV  = -0.7825
% xm_SEP = 0.3872
% xm_RV  = 0.9671
% ym  = 0.649

R_LV = 1/(2*xm_LV./(xm_LV.^2 + ym.^2));
R_SEP = 1/(2*xm_SEP./(xm_SEP.^2 + ym.^2));
R_RV = 1/(2*xm_RV./(xm_RV.^2 + ym.^2));

%I should do that for 3 different circle.
x_center_LV = xm_LV - R_LV;
y_center_LV = 0;

x_center_SEP =  xm_SEP - R_SEP;
y_center_SEP = 0;

x_center_RV = xm_RV - R_RV;
y_center_RV = 0;

theta_LV = linspace( pi/2 + (pi/2 + asin(ym/R_LV)), -pi/2 - (pi/2 + asin(ym/R_LV)), 1000);
xx_LV = R_LV*cos(theta_LV) + x_center_LV;
yy_LV = R_LV*sin(theta_LV) + y_center_LV;

theta_SEP = linspace( pi/2 - (pi/2 - asin(ym/R_SEP)), -pi/2 + (pi/2 - asin(ym/R_SEP)), 1000);
xx_SEP = R_SEP*cos(theta_SEP) + x_center_SEP;
yy_SEP = R_SEP*sin(theta_SEP) + y_center_SEP;

theta_RV = linspace( pi/2 + (pi/2 - asin(ym/R_RV)), -pi/2 - (pi/2 - asin(ym/R_RV)), 1000);
xx_RV = R_RV*cos(theta_RV) + x_center_RV;
yy_RV = R_RV*sin(theta_RV) + y_center_RV;
figure(23)
hold on
plot(xx_LV,yy_LV); axis equal;
plot(xx_SEP,yy_SEP); axis equal;
plot(xx_RV,yy_RV); axis equal;
plot([0 0],[-ym ym])
legend('LV','SEP','RV','Cap')
set(gca,'xdir','reverse')
 end
end
