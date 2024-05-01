function [T_vals, Y] = odeSolver(f, x0, endTime, dt, method)
    numSteps = round(endTime / dt);
    T_vals = linspace(0, endTime, numSteps+1);  % Changed T to T_vals and T to endTime for clarity
    Y = zeros(length(x0), numSteps+1);
    Y(:,1) = x0;

    for i = 1:numSteps
        if method == 0  % Euler method
            Y(:,i+1) = Y(:,i) + dt * f(Y(:,i));
        elseif method == 1  % Runge-Kutta method
            k1 = f(Y(:,i));
            k2 = f(Y(:,i) + 0.5 * dt * k1);
            k3 = f(Y(:,i) + 0.5 * dt * k2);
            k4 = f(Y(:,i) + dt * k3);
            Y(:,i+1) = Y(:,i) + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        end
    end
end
