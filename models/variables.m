%% Parametry pro simulaci

C2K = 273.15;

T0 = 23; %poc. teplota ve stupnich Celsia
T0 = T0 + C2K;

T_inf = 23; %teplota okolniho prostredu
T_inf = T_inf + C2K;


alfa = 0.01; %skalovani inzenyrskych jednotek z procent na Watty dodavaneho vykonu
alfa2 = 0.0075; %druhe topeni ma mensi vykon 0.75W !

cp = 500; %merna tepelna kapacita topneho telesa
A = 7e-4; %povrch topeni v m^2 - POZOR zmena oproti single-mass modelu, zmenseni ztrat vedenim do prostredi a naopak zvetseni predani tepla mezi telesy
As = 5e-4; %povrch kontaktu mezi telesy

m = 4e-3; %hmotnost topneho telesa
U = 10; %koeficient prestupu tepla
eps = 0.5; %emisivita - snizeni z 0.9 od vyrobce, zmenseni ztrat vyzarovanim
sig = 5.67e-8; %SB konstanta

%% Definice koeficientů K

K1 = (U*A)/(m*cp);
K2 = (eps*sig*A)/(m*cp);
K3 = (U*As)/(m*cp);
K4 = (eps*sig*As)/(m*cp);
Ku1 = alfa/(m*cp);
Ku2 = alfa2/(m*cp);



