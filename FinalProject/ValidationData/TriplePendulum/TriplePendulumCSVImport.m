clear;close all;clc;
SAVE_PLOTS = 1;

Body3POS = table2array(readtable('ValidationData\TriplePendulum\Body3POS.csv'))';
Body3VEL = table2array(readtable('ValidationData\TriplePendulum\Body3VEL.csv'))';
Body3ACC = table2array(readtable('ValidationData\TriplePendulum\Body3ACC.csv'))';
Body3Omega = table2array(readtable('ValidationData\TriplePendulum\Body3Omega.csv'))';

time = Body3POS(1,:);

%% Origin Plots BODY 1
% O Position plot
figure;
subplot(3,1,1);
hold on;
plot(time,Body3POS(2,:)/1000);
plot(time,Body3POS(3,:)/1000);
plot(time,Body3POS(4,:)/1000);
title(['Solidworks Triple Pendulum: O Global Position Body 3']);
xlabel("t (s)");
ylabel("position (m)");
ylim([-10,10]);
legend('x','y','z');
hold off;

% O Velocity plot
%figure;
subplot(3,1,2);
hold on;
plot(time,Body3VEL(2,:)/1000);
plot(time,Body3VEL(3,:)/1000);
plot(time,Body3VEL(4,:)/1000);
title(['Solidworks Triple Pendulum: O Global Velocity Body 3']);
xlabel("t (s)");
ylabel("velocity (m/s)");
legend('x','y','z');
ylim([-20,20]);
hold off;

% O Acceleration plot
%figure;
subplot(3,1,3);
hold on;
plot(time,Body3ACC(2,:)/1000);
plot(time,Body3ACC(3,:)/1000);
plot(time,Body3ACC(4,:)/1000);
title(['Solidworks Triple Pendulum: O Global Acceleration Body 3']);
xlabel("t (s)");
ylabel("acceleration (m/s^2)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,['Solidworks_O_Body_3_Plot.png']);
end

% Omega Plot per body
figure;
hold on;
plot(time,Body3Omega(2,:)*pi/180);
title(['Solidworks Triple Pendulum: \omega Body 3']);
xlabel("t (s)");
ylabel("\omega (rad/s)");
legend('x','y','z');
hold off;
if SAVE_PLOTS
	saveas(gcf,['Solidworks_triplePendulum_omega_Body_3_Plot.png']);
end

