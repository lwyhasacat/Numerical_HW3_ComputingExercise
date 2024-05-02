function [U_next, V_next] = latticeStepRK4(U, V, dt, beta)
    rows = size(U, 1);
    cols = size(U, 2);
    U_next = zeros(size(U));
    V_next = zeros(size(V));
    
    function [dUdt, dVdt] = derivatives(U, V)
        dUdt = V;
        dVdt = zeros(size(U));
        for j = 2:rows-1
            for k = 2:cols-1
                laplacian = U(j+1, k) + U(j-1, k) + U(j, k+1) + U(j, k-1) - 4*U(j, k);
                dVdt(j, k) = laplacian - beta * U(j, k)^3;
            end
        end
    end

    [k1U, k1V] = derivatives(U, V);
    [k2U, k2V] = derivatives(U + 0.5 * dt * k1U, V + 0.5 * dt * k1V);
    [k3U, k3V] = derivatives(U + 0.5 * dt * k2U, V + 0.5 * dt * k2V);
    [k4U, k4V] = derivatives(U + dt * k3U, V + dt * k3V);

    U_next = U + dt * (k1U + 2*k2U + 2*k3U + k4U) / 6;
    V_next = V + dt * (k1V + 2*k2V + 2*k3V + k4V) / 6;
end