%%
close all; clear all; clc

model = 'model_2heaters_linear.slx';
open_system(model);

C2K = 273.15;
T0 = 23; 
T0 = T0 + C2K;
T_inf = 23; T_inf = T_inf + C2K;
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

%% Nalezení T10 a T20 pomocí numerického řešení soustavy nelineárních rovnic
Q10 = 70;
Q20 = 40;

options = optimoptions('fsolve', 'Display', 'iter');
x0 = [0; 0];

[x, fval, exitflag] = fsolve(@(x) TEquationSystem(K1, K2, K3, K4, Ku1, Ku2, Q10, Q20, T_inf, x), x0, options);

if exitflag > 0
    fprintf('Solution found:\n');
    fprintf('T10 = %.4f, T20 = %.4f\n', x(1)-C2K, x(2)-C2K);
    T10 = x(1);
    T20 = x(2);
else
    fprintf('Solution not found. Exit flag: %d\n', exitflag);
end

%% Definice parametrů pro model v přírůstkových souřadnicích
mat_A = [-K1-(4*K2*T10^3)-(4*K4*T10^3)-K3, K3+(4*K4*T20^3); K3+(4*K4*T10^3), -K1-(4*K2*T20^3)-(4*K4*T20^3)-K3];
mat_B = [Ku1, 0; 0, Ku2];
mat_C = [1, 0; 0, 1];

simOut = sim(model);
figure;
plot(simOut.dualHeaterLinear);
title('Teplota topného tělesa');
xlabel('Čas (s)');
ylabel('Teplota (C)');

%% Definice přenosových funkcí F1 a F2

numerator1 = [mat_B(1,1), -mat_B(1,1)*mat_A(2,2)];
numerator2 = [mat_B(2,2)*mat_A(1,2)];
denominator = [1, -(mat_A(1,1)+mat_A(2,2)), det(mat_A)];
F1 = tf(numerator1, denominator);
F2 = tf(numerator2, denominator);
F1
F2

p1 = pole(F1);
z1 = zero(F1);
G1 = dcgain(F1);
p2 = pole(F2);
z2 = zero(F2);
G2 = dcgain(F2);

%% Přechodová a impulsní charakteristika F1

% Časová osa
t = 0:0.01:2000;

% Vypočtená přechodová odezva
y_step = (((p1(1)-mat_A(2,2))*mat_B(1,1))/(p1(1)*(p1(1)-p1(2))))*exp(p1(1)*t) + (((p1(2)-mat_A(2,2))*mat_B(1,1))/(p1(2)*(p1(2)-p1(1))))*exp(p1(2)*t) + ((-mat_A(2,2)*mat_B(1,1))/(p1(1)*p1(2)));

% Vypočtená impulsní odezva
y_impulse = (((p1(1)-mat_A(2,2))*mat_B(1,1))/(p1(1)-p1(2)))*exp(p1(1)*t) + (((p1(2)-mat_A(2,2))*mat_B(1,1))/(p1(2)-p1(1)))*exp(p1(2)*t);

% Simulace přechodové odezvy
[y_step_num, t_step_num] = step(F1, t);

% Simulace impulsní odezvy
[y_impulse_num, t_impulse_num] = impulse(F1, t);

% Zobrazení výsledků
figure;
hold on;
plot(t, y_step, 'b', 'DisplayName', 'Přechodová charakteristika');
plot(t, y_impulse, 'r', 'DisplayName', 'Impulsní charakteristika');
xlabel('čas (s)');
ylabel('Odezva');
legend;
grid;
figure;
hold on;
plot(t, y_step_num, 'b', 'DisplayName', 'Simulace přechodové odezvy');
plot(t, y_impulse_num, 'r', 'DisplayName', 'Simulace impulsní odezvy'); 
xlabel('čas (s)');
ylabel('Odezva');
legend;
grid;

%% Frobeniova a Jordanova Stavová forma

A_F = [0,1;denominator(3),denominator(2)];
B_F = [0;1];
C_F = [numerator1(2),numerator1(1)];
D_F = 0;
[V,J]=eig(A_F)
A_J = inv(V)*A_F*V;
B_J = inv(V)*B_F;
C_J = C_F*V;
D_J = D_F;

F1_F = ss(A_F, B_F, C_F, D_F); % Frobeniova forma
F1_J = ss(A_J, B_J, C_J, D_J); % Jordanova forma

t = 0:0.01:10;
u = sin(t); % Sinusový signál
[y_tf, t_tf] = lsim(F1, u, t);
[y_F, t_F] = lsim(F1_F, u, t);
[y_J, t_J] = lsim(F1_J, u, t);

figure;
hold on;
%plot(t, u, 'y', 'DisplayName', 'Vstup', 'LineWidth', 1);
plot(t_tf, y_tf, 'r', 'DisplayName', 'Vnější popis', 'LineWidth', 1);
plot(t_F, y_F, 'g--', 'DisplayName', 'Frobeniova forma', 'LineWidth', 1);
plot(t_J, y_J, 'b-.', 'DisplayName', 'Jordanova forma', 'LineWidth', 1);
title('Porovnání vstupně/výstupního chování systému');
xlabel('Čas (s)');
ylabel('Výstup');
legend('show');
grid on;

%% Bodeho a Nyquistova frekvenční charakteristika

t = 0:0.0005:10000;
amplitudes = [1, 1, 1, 1];
frequencies = [2*pi*0.0005, 2*pi*0.001, 2*pi*0.005, 2*pi*0.01];
gains = [];
phases = [];
phases_deg = [];

for i=1:4
    u = [amplitudes(i)*sin(frequencies(i)*t)]; % Generování harmonického vstupu
    u = u(t > 2000); % Odříznutí přechodového děje
    t = t(t > 2000); % Odříznutí přechodového děje
    [y, ~] = lsim(F1, u, t);

    %figure; hold on;
    %plot(t, u) % Vykreslení vstupu
    %plot(t, y) % Vykreslení odezvy
    %legend("Vstup", "Výstup")

    % Amplituda odezvy
    amplitude_out = max(abs(y)) - mean(y); % Odečtení průměru (pokud není signál vycentrovaný kolem nuly)
    gain = amplitude_out / amplitudes(i);

    % Fáze odezvy
    fft_u = fft(u);
    fft_y = fft(y);
    [~, idx] = max(abs(fft_u));
    phase_u = angle(fft_u(idx));
    phase_y = angle(fft_y(idx));
    phi = phase_y - phase_u;
    phi_deg = rad2deg(phi);

    % Výsledky
    disp(['Zisk (dB): ', num2str(20*log10(gain))]);
    disp(['Fáze (°): ', num2str(phi_deg)]);

    gains = [gains, gain];
    phases = [phases, phi];
    phases_deg = [phases_deg, phi_deg];
end


bode_gains = 20*log10(gains);
bode_phases = phases_deg;

nyquist_real = gains.*cos(phases);
nyquist_imag = gains.*sin(phases);

[mag, phase, omega] = bode(F1);
mag = squeeze(mag);
phase = squeeze(phase);
omega = squeeze(omega);

figure; % Bode
subplot(2, 1, 1)
hold on;
plot(omega, 20*log10(mag))
scatter(frequencies, bode_gains, "xr")
set(gca,'xscale','log')
title("Bodeho charakteristika")
legend("", "Odsimulované vstupy")
xlabel("Frekvence (rad/s)")
ylabel("Amplituda (dB)")
subplot(2, 1, 2)
hold on;
plot(omega, phase)
scatter(frequencies, bode_phases, "xr")
set(gca,'xscale','log')
legend("", "Odsimulované vstupy")
xlabel("Frekvence (rad/s)")
ylabel("Fáze (°)")

figure; % Nyquist
hold on;
nyquist(F1)
scatter(nyquist_real, nyquist_imag, "xr")
title("Nyquistova charakteristika")
legend("", "Odsimulované vstupy")
xlabel("Reálná složka")
ylabel("Imaginární složka")
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
