Q1_range = [0:100];
Q2 = 40; % konstantní hodnota výkonu druhého topení, použijeme například zadanou hodnotu Q20

%% Nelineární model

scNonLinearT1 = [];
scNonLinearT2 = [];

for q1 = Q1_range
    [x, fval, exitflag] = fsolve(@(x) TEquationSystem(K1, K2, K3, K4, Ku1, Ku2, q1, Q2, T_inf, x), x0, options);
    
    if exitflag > 0
        fprintf('Solution found:\n');
        fprintf('T1 = %.4f, T2 = %.4f\n', x(1)-C2K, x(2)-C2K);
        scNonLinearT1 = [scNonLinearT1, x(1)];
        scNonLinearT2 = [scNonLinearT2, x(2)];
    else
        fprintf('Solution not found. Exit flag: %d\n', exitflag);
    end
end

%% Lineární model

C1 = -(mat_A(1,:)*[T10;T20])-mat_B(1,1)*Q10;
C2 = mat_B(2,2)*Q2-(mat_A(2,:)*[T10;T20])-mat_B(2,2)*Q20;

scLinearT2 = (-mat_B(1,1)*Q1_range - C1 + (mat_A(1,1)*C2)/mat_A(2,1))/(mat_A(1,2) - (mat_A(1,1)*mat_A(2,2))/(mat_A(2,1)));
scLinearT1 = (-C2 - mat_A(2,2)*scLinearT2)/(mat_A(2,1));

%% Grafy statické charakteristiky
figure;
hold on;
plot(Q1_range,scNonLinearT1-C2K,'r','MarkerSize',10);
plot(Q1_range,scLinearT1-C2K,'b','MarkerSize',10);
ylabel('Temperature T1 (degC)');
xlabel('Heater Q1 (%)');
legend('Non-Linear Model','Linear Model', 'Location','SouthEast')

figure;
hold on;
plot(Q1_range,scLinearT2-C2K,'r','MarkerSize',10);
plot(Q1_range,scNonLinearT2-C2K,'b','MarkerSize',10);
ylabel('Temperature T2 (degC)');
xlabel('Heater Q1 (%)');
legend('Non-Linear Model','Linear Model', 'Location','SouthEast')