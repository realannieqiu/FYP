
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

% Parameters
n = 100; % Grid lines in x or y direction
region = [0 2*pi 0 2*pi]; % Domain [a b c d]
%tol = 1e-5; % Error tolerance
k = 3; % Value of k
p = 1; % Value of p

% Grid setup
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


% Jacobi iteration with error tracking
max_iter = 1000; % Fixed number of iterations
L2_errors = zeros(1, max_iter);
Linf_errors = zeros(1, max_iter);
iter = 0;


while iter < max_iter
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

       % Exact solution for comparison
    usol = zeros(n+1, n+1);

    for i = 1:n+1
        for j = 1:n+1
            usol(i, j) = exactu(x(i), y(j));
        end
    end

    % Calculate the errors
    L2_errors(iter) = sqrt(sum(sum((u2 - usol).^2)) * h1 * h2);
    Linf_errors(iter) = max(max(abs(u2 - usol)));
end

% Display errors and iterations
disp('The number of iterations is:');
disp(iter);

% Plot results
figure(1); mesh(x, y, u2); title('The Estimated Solution with FDM');
xlabel('x'); ylabel('y'); zlabel('u');

figure(2); mesh(x, y, u2 - usol); title('The Error Plot');
xlabel('x'); ylabel('y'); zlabel('Error');

figure(3); mesh(x, y, usol); title('The Exact Solution');
xlabel('x'); ylabel('y'); zlabel('u');

% Plot error vs iterations
figure(4);
plot(1:max_iter, L2_errors, '-o', 'DisplayName', 'L_{2} Error');
hold on;
plot(1:max_iter, Linf_errors, '-x', 'DisplayName', 'L_{\infty} Error');
xlabel('Number of Iterations');
ylabel('Error Magnitude');
title('Error Against the Number of Iterations');
legend('Location', 'northeast');
grid on;


% % Display errors and iterations
% disp('The number of iterations is:');
% disp(iter);
% disp('The L2 error is:');
% disp(L2_error);
% disp('The Linf error is:');
% disp(Linf_error);
% 
% % Plot results
% figure(1); mesh(x, y, u2); title('The Estimated Solution with FDM');
% xlabel('x'); ylabel('y'); zlabel('u');
% figure(2); mesh(x, y, u2 - usol); title('The Error Plot');
% xlabel('x'); ylabel('y'); zlabel('Error');
% figure(3); mesh(x, y, usol); title('The Exact Solution');
% xlabel('x'); ylabel('y'); zlabel('u');
% figure(4); bar([L2_error, Linf_error]); 
% set(gca, 'XTickLabel', {'L2 Error', 'Linf Error'});
% title('L2 and Linf Errors');
% ylabel('Error Magnitude');
