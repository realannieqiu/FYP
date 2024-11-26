clear; close all;

% Inputs
n = input('Enter the grid lines in x, y, and z directions: ');
region = input('Enter the domain [a b c d e f]: '); % [x_min x_max y_min y_max z_min z_max]
tol = input('Enter error tolerance, say 1e-5: ');
const = input('Enter value of a: ');

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
        for const = 1:n+1
            if i == 1 || i == n+1 || j == 1 || j == n+1 || k == 1 || k == n+1
                u2(i, j, k) = exactu(x(i), y(j), z(k));
            end
        end
    end
end

% Jacobi iteration
iter = 0;
while max(max(max(abs(u1 - u2)))) > tol
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
end

% Display number of iterations
disp('The number of iterations is:');
disp(iter);

% Exact solution for comparison
usol = zeros(n+1, n+1, n+1);
for i = 1:n+1
    for j = 1:n+1
        for k = 1:n+1
            usol(i, j, k) = exactu(x(i), y(j), z(k));
        end
    end
end

% Calculate the error
e = max(max(max(abs(u2 - usol))));

% Plot results
[X, Y, Z] = meshgrid(x, y, z);
figure(1); slice(X, Y, Z, u2, [], [], z); title('The estimated solution with FDM');
figure(2); slice(X, Y, Z, u2 - usol, [], [], z); title('The error plot');
figure(3); slice(X, Y, Z, usol, [], [], z); title('The Exact Solution');

% Parameters
torus_a = const; % Major radius
r_vals = linspace(region(5), region(6), n+1);
theta = linspace(0, 2*pi, n+1);
phi = linspace(0, 2*pi, n+1);

% Correct grid structure using ndgrid
[R, Theta, Phi] = ndgrid(r_vals, theta, phi);

% Convert to Cartesian coordinates
X1 = (torus_a - R .* cos(Theta)) .* cos(Phi);
Y1 = (torus_a - R .* cos(Theta)) .* sin(Phi);
Z1 = R .* sin(Theta);

disp(size(X1));
disp(size(Y1));
disp(size(Z1));


% Example: Reshape u2 to match the 3D grid dimensions
SolutionInVolume = reshape(u2, size(R));
disp(size(SolutionInVolume));


figure;
slice(X1, Y1, Z1, SolutionInVolume, [], [], linspace(r_vals(1), r_vals(end), 5)); % Example slicing along r-axis
% [X_test, Y_test, Z_test] = ndgrid(1:10, 1:10, 1:10);
% V_test = rand(size(X_test)); % Random 3D data
% figure;
% slice(X_test, Y_test, Z_test, V_test, [], [], [3, 6, 9]); % Test slicing planes
shading interp;
colormap jet;
colorbar;
title('Slice Plot of Solution in Toroidal Volume');
xlabel('\theta'); ylabel('\omega'); zlabel('r');

