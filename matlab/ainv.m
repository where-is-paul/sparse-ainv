function [Z, D] = ainv(A)
    n = size(A, 1);
    Z = eye(n);
    D = eye(n);
    for i = 1:n
        for j = i:n
            D(i,i) = A(i,:) * Z(:,j);
        end
        for j = i+1:n
            Z(:, j) = Z(:, j) - D(j,j) / D(i, i) * Z(:, i);
        end
    end
end