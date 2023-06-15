clear all
clc

syms C S

omega_t = linspace(0, 2*pi, 100); % calculate omega

eta_0 = C*cos(omega_t) - C + (-S-omega_t).*sin(omega_t);
eta_1 = 2*C*cos(omega_t) - 2*C - 2*omega_t.*sin(omega_t) + 4;
eta_2 = 2*S*sin(omega_t) + 6;
eta_3 = -2*C*cos(omega_t) + 2*C + 2*omega_t.*sin(omega_t) + 4;
eta_4 = -C*cos(omega_t) + C + (-S+omega_t).*sin(omega_t) + 1;
H3 = 64*S*sin(omega_t) + 64;

% figure(1)
% subplot(3, 2, 1); plot(omega_t, eta_0); title('\eta_0'); xlabel('\omega t'); hold on; grid on
% subplot(3, 2, 2); plot(omega_t, eta_1); title('\eta_1'); xlabel('\omega t'); hold on; grid on
% subplot(3, 2, 3); plot(omega_t, eta_2); title('\eta_2'); xlabel('\omega t'); hold on; grid on
% subplot(3, 2, 4); plot(omega_t, eta_3); title('\eta_3'); xlabel('\omega t'); hold on; grid on
% subplot(3, 2, 5); plot(omega_t, eta_4); title('\eta_4'); xlabel('\omega t'); hold on; grid on
% 
% figure(2)
% plot(omega_t, H3)
% hold on
% grid on


omega_t = linspace(0, 2*pi, 100); % calculate omega

x_0 = -1;

% Solve for C only
x_0 = -1;
fun1 = @(C) double(subs(C*cos(omega_t) - C + (-S-omega_t).*sin(omega_t), S, x_0));
sol_C = fsolve(fun1, x_0);
disp('Solutions for C:')
disp(sol_C)

% Solve for S only
fun2 = @(S) double(subs(2*S*sin(omega_t) + 6, C, x_0));
sol_S = fsolve(fun2, x_0);
disp('Solutions for S:')
disp(sol_S)

% Solve for both C and S
fun3 = @(x) [double(subs(C*cos(omega_t) - C + (-x(2)-omega_t).*sin(omega_t), [C, S], [x(1), x(2)]));
             double(subs(-C*cos(omega_t) + C + (-x(2)+omega_t).*sin(omega_t) + 1, [C, S], [x(1), x(2)]));
             double(subs(64*x(2)*sin(omega_t) + 64, [C, S], [x(1), x(2)]))];
sol_both = fsolve(fun3, [x_0, x_0]);
disp('Solutions for both C and S:')
disp(sol_both)