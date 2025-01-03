close all; clear all; clc

model = 'model_2heaters.slx';
open_system(model);

C2K = 273.15;
T0 = 23; 
T0 = T0 + C2K;
T_inf = 23; 
T_inf = T_inf + C2K;
alfa = 0.01; 
alfa2 = 0.0075; 
cp = 500; 
A = 7e-4; 
As = 5e-4; 
m = 4e-3; 
U = 10;
eps = 0.5; 
sig = 5.67e-8;
K1 = (U*A)/(m*cp);
K2 = (eps*sig*A)/(m*cp);
K3 = (U*As)/(m*cp);
K4 = (eps*sig*As)/(m*cp);
Ku1 = alfa/(m*cp);
Ku2 = alfa2/(m*cp);
simOut = sim(model);
figure;
plot(simOut.dualHeater);
title('Teplota topného tělesa');
xlabel('Čas (s)');
ylabel('Teplota (C)');
%%
tclab;
figure(1)
t1s = [];
t2s = [];
h1s = [];
h2s = [];
ht1 = 0;
ht2 = 0;
h1(ht1);
h2(ht2);
for i = 1:2500
    tic;
    if i==100
        ht1 = 35;
        h1(ht1);
    end
    if i==500
        ht2 = 70;
        h2(ht2);
    end
    if i==1000
        ht1 = 25;
        h1(ht1);
    end
    if i==1500
        ht2 = 0;
        h2(ht2);
    end
    if i==1800
        ht1 = 0;
        h1(ht1);
        h2(ht2);
    end
    t1 = T1C();
    t2 = T2C();
    h1s = [h1s,ht1];
    h2s = [h2s,ht2];
    t1s = [t1s,t1];
    t2s = [t2s,t2];
    n = length(t1s);
    time = linspace(0,n+1,n);
    clf
    plot(time,t1s,'r.','MarkerSize',10);
    hold on;
    plot(time,t2s,'b.','MarkerSize',10);
    plot(simOut.simout);
    yyaxis right;
    plot(time,h1s,'o-','LineWidth',1);
    plot(time,h2s,'c-','LineWidth',1);
    ylabel('Heater (%)');
    yyaxis left;
    ylabel('Temperature (degC)');
    xlabel('Time (sec)');
    legend('Temperature 1','Temperature 2','Temp 1 model', 'Temp 2 model','Heater 1','Heater 2', 'Location','SouthEast');
    drawnow;
    t = toc;
    pause(max(0.01,1.0-t));
end

%%
h1(0);
h2(0);



