clear all
clc
syms K
omega_t = linspace(0, 2*pi, 100); % calculate omega
K = -0.33333;
eta_0 = 2 * cos(omega_t);

eta_1 = -2 * K * cos(omega_t) + 2 * K + 2 * cos(omega_t) + 2;

eta_2 = 4 * K * cos(omega_t) - 4 * K - 2 * cos(omega_t) + 2;

eta_3 = -2 * K * cos(omega_t) + 2 * K - 2 * cos(omega_t) + 2;

H2 = eta_2 .* eta_1 - eta_3 .* eta_0;
figure(1)

subplot(2, 2, 1); plot(omega_t, eta_0); title('\eta_0'); xlabel('\omega t'); hold on; grid on
subplot(2, 2, 2); plot(omega_t, eta_1); title('\eta_1'); xlabel('\omega t'); hold on; grid on
subplot(2, 2, 3); plot(omega_t, eta_2); title('\eta_2'); xlabel('\omega t'); hold on; grid on
subplot(2, 2, 4); plot(omega_t, eta_3); title('\eta_3'); xlabel('\omega t'); hold on; grid on

figure(2)
plot(omega_t, H2)
hold on
grid on
x_0=-100
% Solve for eta_0
fun1 = @(K) 2 * cos(omega_t) - 0;
sol_eta_0 = fsolve(fun1, x_0);
disp(sol_eta_0)

% Solve for eta_1
fun2 = @(K) -2 * K * cos(omega_t) + 2 * K + 2 * cos(omega_t) + 2 - 0;
sol_eta_1 = fsolve(fun2, x_0);
disp(sol_eta_1)

% Solve for eta_2
fun3 = @(K) 4 * K * cos(omega_t) - 4 * K - 2 * cos(omega_t) + 2 - 0;
sol_eta_2 = fsolve(fun3, x_0);
disp(sol_eta_2)

% Solve for eta_3
fun4 = @(K) -2 * K * cos(omega_t) + 2 * K - 2 * cos(omega_t) + 2 - 0;
sol_eta_3 = fsolve(fun4, x_0);
disp(sol_eta_3)

% Solve for H2
fun5 = @(K) eta_2 .* eta_1 - eta_3 .* eta_0 - 0;
sol_H2 = fsolve(fun5, x_0);
disp(sol_H2)
