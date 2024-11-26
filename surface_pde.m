clear; close all;

% Inputs
n = input('Enter the grid lines in x or y direction: ');
region = input('Enter the domain [a b c d]: ');
tol = input('Enter error tolerance, say 1e-5: ');
k = input('Enter value of a: ');
p = input('Enter value of R: ');

a = region(1); b = region(2); c = region(3); d = region(4);
h1 = (b-a)/n;   
h2 = (d-c)/n;  
x = a:h1:b;    
y = c:h2:d;

% Initialize variables
u2 = zeros(n+1, n+1); % The solution component
u1 = 100 * ones(n+1, n+1);
f = zeros(n+1, n+1);

% Exact solution function
exactu = @(x, y) sin(x) .* cos(y);

% Calculate the right-hand side function f
for i = 1:n+1
    for j = 1:n+1
        f(i, j) = -(sin(x(i)) .* cos(y(j)) .* cos(x(i))) / (k - p * cos(x(i))) ...
                  + 2 * sin(x(i)) .* cos(y(j)) ...
                  + (sin(x(i)) .* cos(y(j))) / (k - p * cos(x(i)))^2;
    end
end

% Apply boundary conditions
for i = 1:n+1
  u2(1,i) = exactu(a,y(i));
  u2(n+1,i) = exactu(b,y(i));
  u2(i,1) = exactu(x(i),c);
  u2(i,n+1) = exactu(x(i),d);
end

% Jacobi iteration
iter = 0; % Rename k to iter
while max(max(abs(u1 - u2))) > tol
    iter = iter + 1;
    u1 = u2;
    for i = 2:n
        for j = 2:n
            u2(i, j) = (f(i, j) ...
                + ((1 / (p * h1)^2) + (sin(x(i))) / (p * (k - p * cos(x(i))) * 2 * h1)) * u1(i+1, j) ...
                + ((1 / (p * h1)^2) - (sin(x(i))) / (p * (k - p * cos(x(i))) * 2 * h1)) * u1(i-1, j) ...
                + (1 / (k - p * cos(x(i)))^2 * h2^2) * u1(i, j+1) ...
                + (1 / ((k - p * cos(x(i))))^2 * h2^2) * u1(i, j-1)) ...
                / ((2 / (p * h1)^2) + (2 / ((k - p * cos(x(i))))^2 * h2^2) + 1);
        end
    end
end

% Display number of iterations
disp('The number of iterations is:');
disp(iter);

% Exact solution for comparison
usol = zeros(n+1, n+1);
for i = 1:n+1
    for j = 1:n+1
        usol(i, j) = exactu(x(i), y(j));
    end
end

% Calculate the error
e = max(max(abs(u2 - usol)));

% Plot results
figure(1); mesh(u2); title('The estimated solution with FDM');
figure(2); mesh(u2 - usol); title('The error plot');
figure(3); mesh(usol); title('The Exact Solution');



% Example: Torus surface
theta = linspace(0, 2*pi, n+1); % Angular parameter
phi = linspace(0, 2*pi, n+1);   % Angular parameter
[Theta, Phi] = meshgrid(theta, phi);

% Convert to Cartesian coordinates
X = (3 - 1 * cos(Theta)) .* cos(Phi);
Y = (3 - 1 * cos(Theta)) .* sin(Phi);
Z = 1 * sin(Theta);

% Assuming u2 is the solution on a (theta, phi) grid
SolutionOnSurface = reshape(u2, size(X));

figure;
surf(X, Y, Z, SolutionOnSurface);
shading interp; % Smooth shading
colormap jet; % Color map
colorbar; % Add color bar for solution magnitude
title('Solution on Curvilinear Surface');
xlabel('X'); ylabel('Y'); zlabel('Z');

