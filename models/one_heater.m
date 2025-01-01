close all; clear all; clc

model = 'model_1heater.slx';
open_system(model);

C2K = 273.15;
T0 = 23; 
T0 = T0 + C2K;
T_inf = 23; 
T_inf = T_inf + C2K;
alfa = 0.01; 
alfa2 = 0.0075; 
cp = 500;
A = 0.0012;
m = 4e-3; 
U = 10;
eps = 0.5; 
sig = 5.67e-8;
K1 = (U*A)/(m*cp);
K2 = (eps*sig*A)/(m*cp);
Ku1 = alfa/(m*cp);

simOut = sim(model);
figure;
plot(simOut.simout);
title('Teplota topného tělesa');
xlabel('Čas (s)');
ylabel('Teplota (C)');
%%
tclab;
figure(1);
t1s = [];
h1s = [];
ht1 = 0;
h1(ht1);
for i = 1:2000
    tic;
    if i==1
        ht1 = 50;
        h1(ht1);
    end
    if i==1000
        ht1 = 0;
        h1(ht1);
    end
    t1 = T1C();
    h1s = [h1s,ht1];
    t1s = [t1s,t1];
    n = length(t1s);
    time = linspace(0,n+1,n);
    clf
    plot(time,t1s,'r.','MarkerSize',10);
    hold on;
    plot(simOut.simout);
    yyaxis left;
    ylabel('Temperature (degC)');
    yyaxis right;
    plot(time,h1s,'o-','LineWidth',1);
    ylabel('Heater (%)');
    xlabel('Time (sec)');
    legend('Temperature','Temperature model','Heater', 'Location','NorthEast');
    drawnow;
    t = toc;
    pause(max(0.01,1.0-t));
end

%%
h1(0);