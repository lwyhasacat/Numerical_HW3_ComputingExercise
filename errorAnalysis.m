function errors = errorAnalysis(f, x0, endTime, dt, method)
    [~, Y] = odeSolver(f, x0, endTime, dt, method);  % Solving the ODE
    numSteps = size(Y, 2);  % Number of steps based on the number of columns in Y
    exactY = [cos(linspace(0, endTime, numSteps)); sin(linspace(0, endTime, numSteps))];  % Generating the exact solution
    errors = zeros(1, numSteps);  % Initialize error storage

    for i = 1:numSteps
        errors(i) = norm(Y(:,i) - exactY(:,i));  % Calculate the Euclidean distance (L2 norm) at each step
    end
end
