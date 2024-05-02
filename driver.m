clc;
clear;
close all;

f = @(u) [-(u(1)^2 + u(2)^2) * u(2); (u(1)^2 + u(2)^2) * u(1)];
x0 = [1; 0];
T = 2 * pi;
dt_values = [0.05, 0.025, 0.0125, 0.005, 0.0001]; % time steps

% errors
errors = zeros(1, length(dt_values)); 
Rerrors = zeros(1, length(dt_values)); 

for i = 1:length(dt_values)
    dt = dt_values(i);
    allErrors = errorAnalysis(f, x0, T, dt, 0);
    errors(i) = max(allErrors); % Euler max error
    allRErrors = errorAnalysis(f, x0, T, dt, 1);
    Rerrors(i) = max(allRErrors); % RK4 max error
end

% plotting error
figure;
loglog(dt_values, errors, 'o-', 'DisplayName', 'Euler Method Errors');
hold on;
loglog(dt_values, Rerrors, 's-', 'DisplayName', 'Runge-Kutta Method Errors');
xlabel('Time step');
ylabel('L2 Error');
title('Error Analysis for Convergence');
legend('Location', 'Best');
legend show;

% linear regression; Euler
coeffs = polyfit(log(dt_values), log(errors), 1);
slope = coeffs(1);
fprintf('Estimated order of convergence for Euler method: %f\n', abs(slope));

% Plotting the fitted line for Euler's method
fit_line = exp(polyval(coeffs, log(dt_values)));
loglog(dt_values, fit_line, '--', 'Color', [0.5, 0.5, 0.5], 'DisplayName', sprintf('Euler Fit (Order %.2f)', abs(slope)));

% linear regression; RK4
Rcoeffs = polyfit(log(dt_values), log(Rerrors), 1);
slopeRK = Rcoeffs(1);
fprintf('Estimated order of convergence for Runge-Kutta method: %f\n', abs(slopeRK));

% Plotting the fitted line for RK4 method
Rfit_line = exp(polyval(Rcoeffs, log(dt_values)));
loglog(dt_values, Rfit_line, '--', 'Color', [0.5, 0.5, 0.5], 'DisplayName', sprintf('RK4 Fit (Order %.2f)', abs(slopeRK)));

% Plotting the solution curves
T = 20*pi;
figure;
% dt = 0.1
subplot(1, 2, 1); % Euler
[T_vals, Y] = odeSolver(f, x0, T, 0.1, 0);
X_exact = cos(T_vals);
Y_exact = sin(T_vals);
plot(Y(1,:), Y(2,:), 'b-', 'DisplayName', 'Numerical Solution (Euler)');
hold on;
plot(X_exact, Y_exact, 'r--', 'DisplayName', 'Exact Solution');
legend show;
xlabel('x');
ylabel('y');
title('Numerical vs. Exact Solution (Euler, dt = 0.1)');
axis equal;
axis square;

subplot(1, 2, 2); % RK4
[T_vals, RY] = odeSolver(f, x0, T, 0.1, 1);
X_exact = cos(T_vals);
Y_exact = sin(T_vals);
plot(RY(1,:), RY(2,:), 'b-', 'DisplayName', 'Numerical Solution (RK4)');
hold on;
plot(X_exact, Y_exact, 'r--', 'DisplayName', 'Exact Solution');
legend show;
xlabel('x');
ylabel('y');
title('Numerical vs. Exact Solution (RK4, dt = 0.1)');
axis equal;
axis square;

% dt = 0.0001
figure;
subplot(1, 2, 1); % Euler
[T_vals, Y] = odeSolver(f, x0, T, 0.0001, 0);
X_exact = cos(T_vals);
Y_exact = sin(T_vals);
plot(Y(1,:), Y(2,:), 'b-', 'DisplayName', 'Numerical Solution (Euler)');
hold on;
plot(X_exact, Y_exact, 'r--', 'DisplayName', 'Exact Solution');
legend show;
xlabel('x');
ylabel('y');
title('Numerical vs. Exact Solution (Euler, dt = 0.0001)');
axis equal;
axis square;

subplot(1, 2, 2); % RK4
[T_vals, RY] = odeSolver(f, x0, T, 0.0001, 1);
X_exact = cos(T_vals);
Y_exact = sin(T_vals);
plot(RY(1,:), RY(2,:), 'b-', 'DisplayName', 'Numerical Solution (RK4)');
hold on;
plot(X_exact, Y_exact, 'r--', 'DisplayName', 'Exact Solution');
legend show;
xlabel('x');
ylabel('y');
title('Numerical vs. Exact Solution (RK4, dt = 0.0001)');
axis equal;
axis square;

% command for running simulation:
% runSimulation(3, 50, 0.05, 10);