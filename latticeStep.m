function [U_next, V_next] = latticeStep(U, V, dt, beta)
    rows = size(U, 1);
    cols = size(U, 2);
    U_next = zeros(size(U));
    V_next = zeros(size(V));
    
    for j = 2:rows-1
        for k = 2:cols-1
            laplacian = U(j+1, k) + U(j-1, k) + U(j, k+1) + U(j, k-1) - 4*U(j, k);
            V_next(j, k) = V(j, k) + dt * (laplacian - beta * U(j, k)^3);
            U_next(j, k) = U(j, k) + dt * V_next(j, k);
        end
    end
end
