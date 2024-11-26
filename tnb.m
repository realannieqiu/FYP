% Parameters
t = linspace(0, 4 * pi, 100); % Parameter along the curve
a = 1; % Radius of the helix
b = 0.2; % Pitch of the helix

% Define the curve r_c(omega)
x = a * cos(t);
y = a * sin(t);
z = b * t;
curve = [x; y; z];

% Choose a point on the curve (index 50 here for demonstration)
index = 50;
r_c = curve(:, index); % Position vector at the chosen point

% Calculate Tangent (e1), Normal (e2), and Binormal (e3) at the chosen point
dr = gradient(curve, t(2) - t(1)); % Derivative along the curve
d2r = gradient(dr, t(2) - t(1)); % Second derivative along the curve

T = dr(:, index); % Tangent vector (not unit)
T = T / norm(T); % Normalize to get unit vector e1

N = d2r(:, index); % Normal vector (not unit)
N = N - dot(N, T) * T; % Ensure orthogonality
N = N / norm(N); % Normalize to get unit vector e2

B = cross(T, N); % Binormal vector (already unit if T and N are unit)
B = B / norm(B); % Normalize to get unit vector e3

% Plot the curve
figure;
plot3(x, y, z, 'k', 'LineWidth', 1.5); hold on;
plot3(r_c(1), r_c(2), r_c(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % Position vector point

% Plot TNB frame at the chosen point
quiver3(r_c(1), r_c(2), r_c(3), T(1), T(2), T(3), 0.5, 'r', 'LineWidth', 1.5, 'DisplayName', 'e_1 (Tangent)');
quiver3(r_c(1), r_c(2), r_c(3), N(1), N(2), N(3), 0.5, 'g', 'LineWidth', 1.5, 'DisplayName', 'e_2 (Normal)');
quiver3(r_c(1), r_c(2), r_c(3), B(1), B(2), B(3), 0.5, 'b', 'LineWidth', 1.5, 'DisplayName', 'e_3 (Binormal)');

% Labels and Legend
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('TNB Frame on r_{c}(\omega) Curve');
legend('Curve r_{c}(\omega)', 'Position Vector r_{c}', 'e_1 (Tangent)', 'e_2 (Normal)', 'e_3 (Binormal)');
grid on;
axis equal;
view(3); % 3D view

hold off;
