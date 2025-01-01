%% Měření změny v okolí pracovního bodu
tclab;
figure(1)
t1s = [];
t2s = [];
h1s = [];
h2s = [];
ht1 = 70;
ht2 = 40;
h1(ht1);
h2(ht2);
for i = 1:2400
    tic;
    if i==1200
        ht1 = 45;
        h1(ht1);
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
    yyaxis right;
    plot(time,h1s,'o-','LineWidth',1);
    plot(time,h2s,'c-','LineWidth',1);
    ylabel('Heater (%)');
    yyaxis left;
    ylabel('Temperature (degC)');
    xlabel('Time (sec)');
    legend('Temperature 1','Temperature 2','Heater 1','Heater 2', 'Location','SouthEast');
    drawnow;
    t = toc;
    pause(max(0.01,1.0-t));
end

%%
h1(0);
h2(0);

%%
h1s = h1s(200:end);
h2s = h2s(200:end);
t1s = t1s(200:end);
t2s = t2s(200:end);
h1s = h1s(:);
h2s = h2s(:);
t1s = t1s(:);
t2s = t2s(:);

%%
figure;
hold on;
plot(h1s)
plot(h2s)
plot(t1s)
plot(t2s)

%% Vytvoření přesnějšího modelu
Td = 15;
tau1 = 810/5;
tau2 = 130/5;
K = 21.5/25;

new_F1 = tf([K], [tau1, 1], "InputDelay", Td);
new_F1_2 = tf([K], [tau1*tau2, tau1+tau2, 1], "InputDelay", Td);
new_F1
new_F1_2

%% Porovnání přesnějšího modelu s naměřenými daty
t_end = 1201;   % Délka simulace v sekundách
dt = 1;        % Vzorkovací perioda (s)
time = 0:dt:t_end; % Od 0 do t_end s krokem dt

% Vstupní signál (konstantní výkon topení)
Q1 = -25 * ones(size(time));


% Simulace modelu
T_sim = lsim(new_F1, Q1, time);
T_sim_2 = lsim(new_F1_2, Q1, time);

% Porovnání s daty
figure;
plot(time, t1s-t1s(1), 'b', 'LineWidth', 1.5); hold on;
plot(time, T_sim, 'r--', 'LineWidth', 1.5);
plot(time, T_sim_2, 'g--', 'LineWidth', 1.5);
xlabel('Čas (s)');
ylabel('Teplota (°C)');
legend('Naměřená data', 'Model 1. řádu', 'Model 2. řádu');
grid on;

%% Návrh parametrů pomocí GMK
controlSystemDesigner(new_F1);

%%
disp(cd_C);
cd_C_Z = cd_C.Z{1};
cd_C_P = cd_C.P{1};
cd_C_K = cd_C.K;

Kp = cd_C_K
Ki = cd_C_K * -cd_C_Z

%% Počítání kritického zesílení a mezního dopravního zpoždění
s = tf('s');
R = Kp + Ki/s;

open_loop = new_F1 * R;

margin(open_loop);
[gm, pm, omega_gm, ~] = margin(open_loop);

mag_crit = bode(open_loop, omega_gm);
K_crit = 1 / mag_crit;
disp(K_crit);

for T_delay = 23.6:0.0001:23.7
    [num_delay, den_delay] = pade(T_delay, 4);
    open_delayed = series(open_loop, tf(num_delay, den_delay));
    closed_delayed = feedback(open_delayed, 1);
    poles = pole(pade(closed_delayed, 4));
    if any(real(poles) >= 0)
        fprintf('Mez stability při T_delay = %.2f s\n', T_delay);
        Td_crit = T_delay
        break;
    end
end


%% Testování kritického zesílení a mezního dopravního zpoždění
open_loop_K_crit = open_loop * (K_crit);
closed_loop_K_crit = feedback(open_loop_K_crit, 1);
figure;
step(closed_loop_K_crit);

open_loop_Td_crit = open_loop * tf(1, 1, "InputDelay", Td_crit);
closed_loop_Td_crit = feedback(open_loop_Td_crit, 1);
figure;
step(closed_loop_Td_crit);

