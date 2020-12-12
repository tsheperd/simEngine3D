clear;close all;clc;
SAVE_PLOTS = 1;

UpperPOS = table2array(readtable('ValidationData\DoublePendulum\UpperPOS.csv'))';
UpperVEL = table2array(readtable('ValidationData\DoublePendulum\UpperVEL.csv'))';
UpperACC = table2array(readtable('ValidationData\DoublePendulum\UpperACC.csv'))';
UpperOmega = table2array(readtable('ValidationData\DoublePendulum\UpperOmega.csv'))';

LowerPOS = table2array(readtable('ValidationData\DoublePendulum\LowerPOS.csv'))';
LowerVEL = table2array(readtable('ValidationData\DoublePendulum\LowerVEL.csv'))';
LowerACC = table2array(readtable('ValidationData\DoublePendulum\LowerACC.csv'))';
LowerOmega = table2array(readtable('ValidationData\DoublePendulum\LowerOmega.csv'))';

time = UpperPOS(1,:);

%% Origin Plots BODY 1
% O Position plot
figure;
subplot(3,1,1);
hold on;
plot(time,UpperPOS(2,:)/1000);
plot(time,UpperPOS(3,:)/1000);
plot(time,UpperPOS(4,:)/1000);
title(['Solidworks Double Pendulum: O Global Position Body 1']);
xlabel("t (s)");
ylabel("position (m)");
legend('x','y','z');
hold off;

% O Velocity plot
%figure;
subplot(3,1,2);
hold on;
plot(time,UpperVEL(2,:)/1000);
plot(time,UpperVEL(3,:)/1000);
plot(time,UpperVEL(4,:)/1000);
title(['Solidworks Double Pendulum: O Global Velocity Body 1']);
xlabel("t (s)");
ylabel("velocity (m/s)");
legend('x','y','z');
hold off;

% O Acceleration plot
%figure;
subplot(3,1,3);
hold on;
plot(time,UpperACC(2,:)/1000);
plot(time,UpperACC(3,:)/1000);
plot(time,UpperACC(4,:)/1000);
title(['Solidworks Double Pendulum: O Global Acceleration Body 1']);
xlabel("t (s)");
ylabel("acceleration (m/s^2)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,['Solidworks_O_Body_1_Plot.png']);
end

% Omega Plot per body
figure;
hold on;
plot(time,UpperOmega(2,:)*pi/180);
title(['Solidworks Double Pendulum: \omega Body 1']);
xlabel("t (s)");
ylabel("\omega (rad/s)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,['Solidworks_doublePendulum_omega_Body_1_Plot.png']);
end



%%%%%%% BODY 2
%% Origin Plots
% O Position plot
figure;
subplot(3,1,1);
hold on;
plot(time,LowerPOS(2,:)/1000);
plot(time,LowerPOS(3,:)/1000);
plot(time,LowerPOS(4,:)/1000);
title(['Solidworks Double Pendulum: O Global Position Body 2']);
xlabel("t (s)");
ylabel("position (m)");
legend('x','y','z');
hold off;

% O Velocity plot
%figure;
subplot(3,1,2);
hold on;
plot(time,LowerVEL(2,:)/1000);
plot(time,LowerVEL(3,:)/1000);
plot(time,LowerVEL(4,:)/1000);
title(['Solidworks Double Pendulum: O Global Velocity Body 2']);
xlabel("t (s)");
ylabel("velocity (m/s)");
legend('x','y','z');
hold off;

% O Acceleration plot
%figure;
subplot(3,1,3);
hold on;
plot(time,LowerACC(2,:)/1000);
plot(time,LowerACC(3,:)/1000);
plot(time,LowerACC(4,:)/1000);
title(['Solidworks Double Pendulum: O Global Acceleration Body 2']);
xlabel("t (s)");
ylabel("acceleration (m/s^2)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,['Solidworks_O_Body_2_Plot.png']);
end

% Omega Plot per body
figure;
hold on;
plot(time,LowerOmega(2,:)*pi/180);
title(['Solidworks Double Pendulum: \omega Body 2']);
xlabel("t (s)");
ylabel("\omega (rad/s)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,['Solidworks_doublePendulum_omega_Body_2_Plot.png']);
end