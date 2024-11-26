clear; close all;

% Inputs
% Parameters
n = 100; % Grid lines in x, y, and z directions
region = [0 2*pi 0 2*pi 0 1]; % Domain [x_min x_max y_min y_max z_min z_max]
%otol = 1e-5; % Error tolerance
const = 3; % Value of a

a = region(1); b = region(2); c = region(3); d = region(4);
e = region(5); f = region(6);

h1 = (b-a)/n;   
h2 = (d-c)/n;  
h3 = (f-e)/n;  

x = a:h1:b;    
y = c:h2:d;
z = e:h3:f;

% Initialize variables
u2 = zeros(n+1, n+1, n+1); % The solution component
u1 = 100 * ones(n+1, n+1, n+1);
f = zeros(n+1, n+1, n+1);

% Exact solution function
exactu = @(x, y, z) sin(x) .* cos(y) * (z)^2;

% Calculate the right-hand side function f
for i = 1:n+1
    for j = 1:n+1
        for k = 1:n+1
            f(i, j, k) = -(z(k)+2*const)*(sin(x(i)) * cos(x(i)) * cos(y(j))) / (const - z(k) * cos(x(i))) ...
                         -  sin(x(i)) .* cos(y(j)) ...
                         + (z(k)^2 * sin(x(i)) .* cos(y(j)))/(const - z(k) * cos(x(i)))^2 ...
                         + (sin(x(i)) .* cos(y(j)) .* (z(k))^2 );
        end
    end
end

% Apply boundary conditions
for i = 1:n+1
    for j = 1:n+1
        for k = 1:n+1
            if i == 1 || i == n+1 || j == 1 || j == n+1 || k == 1 || k == n+1
                u2(i, j, k) = exactu(x(i), y(j), z(k));
            end
        end
    end
end

% Jacobi iteration with error tracking
iter = 0;
max_iter = 1000; % Maximum number of iterations
L2_errors = zeros(1, max_iter);
Linf_errors = zeros(1, max_iter);

while iter < max_iter
    iter = iter + 1;
    u1 = u2;
    for i = 2:n
        for j = 2:n
            for k = 2:n
                u2(i, j, k) = (f(i, j, k) ...
                    + ((1 / ((z(k)) * h1)^2) + (sin(x(i))) / ((z(k)) * (const - (z(k)) * cos(x(i))) * 2 * h1)) * u1(i+1, j, k) ...
                    + ((1 / ((z(k)) * h1)^2) - (sin(x(i))) / ((z(k)) * (const - (z(k)) * cos(x(i))) * 2 * h1)) * u1(i-1, j, k) ...
                    + (1 / ((const - (z(k)) * cos(x(i)))^2 * h2^2)) * u1(i, j+1, k) ...
                    + (1 / ((const - (z(k)) * cos(x(i)))^2 * h2^2)) * u1(i, j-1, k) ...
                    + ((1 / (h3)^2) + (const * cos(x(i))) / ((z(k)) * (const - (z(k)) * cos(x(i))) * 2 * h3)) * u1(i, j, k+1) ...
                    + ((1 / (h3)^2) - (const * cos(x(i))) / ((z(k)) * (const - (z(k)) * cos(x(i))) * 2 * h3)) * u1(i, j, k-1)) ...
                    / ((2 / (z(k) * h1)^2) + (2 / ((const - z(k) * cos(x(i)))^2 * h2^2)) + (2 / (h3^2)) + 1);
            end
        end
    end

    % Exact solution for comparison
    usol = zeros(n+1, n+1, n+1);
    for i = 1:n+1
        for j = 1:n+1
            for k = 1:n+1
                usol(i, j, k) = exactu(x(i), y(j), z(k));
            end
        end
    end

    % Calculate the errors
    L2_errors(iter) = sqrt(sum(sum(sum((u2 - usol).^2))) * h1 * h2 * h3);
    Linf_errors(iter) = max(max(max(abs(u2 - usol))));
end

% Display final errors
disp('The final L2 error is:');
disp(L2_errors(iter));
disp('The final Linf error is:');
disp(Linf_errors(iter));

% Plot error vs iterations
figure;
plot(1:max_iter, L2_errors, '-o', 'DisplayName', 'L2 Error');
hold on;
plot(1:max_iter, Linf_errors, '-x', 'DisplayName', 'Linf Error');
xlabel('Number of Iterations');
ylabel('Error Magnitude');
title('Error vs. Number of Iterations');
legend('Location', 'northeast');
grid on;
