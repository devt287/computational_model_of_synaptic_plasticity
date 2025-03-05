a = 0.25;
b = 0.9;
pre=1.3;
% Parameters
K0 = 0.5; k1 = 2/10; k2 = 15/10; k3 = 1/10; k4 = a*120/10;
k11 = 2/10; k12 = 15/10; k13 = 1/10; k14 = a*80/10;
Km1 = 10; Km2 = 0.3; Km11 = 10; Km12 = b;
Km4 = 4; Ktot = 20; Ptot = 20; P0 = 0.5; Ca_basal = 0.1;
Km=3;
% mu and nu functions
mu_par = [1, 2, 40];
par2 = [0.9, 1/5];
v = 0.5; 

mu = @(Ca) (Ca^mu_par(1)) / ((Ca^mu_par(1)) + mu_par(2)^mu_par(1)) * mu_par(3);
nu = @(pK) v * par2(1) ./ (1 + par2(2) * pK) + 1 - par2(1);

ratio = 100;
pK_init = 0.018;
P_init = 0.085;
C_init = 5.708;
C_tot = (k11*(Ptot - P_init)*P_init) / (Km11 + Ptot - P_init);  
Km11 = 0.2;
%k11_new = C_tot * (Km11 + Ptot - P_init) / ((Ptot - P_init) * P_init);
k11_new=0.045;
init_val = [pK_init, P_init, C_init, Ca_basal];

% Define the nullcline equations
syms pK P

C_expr = nu(pK)*(Ktot - pK - K0)/(mu(Ca_basal));

dpK_dt = k1*((Ktot - pK - C_expr)/(Km1 + Ktot - pK - C_expr))*pK - k2*(pK/(Km2 + pK))*(P + P0) + k3*K0 + (k4*(Ca_basal^4)*(Ktot - pK - C_expr))/(Km4^4 + Ca_basal^4);
dP_dt = k11_new*((Ptot - P)/(Km11 + Ptot - P))*P - k12*(P/(Km12 + P))*(pre*pK + K0) + k13*P0 + (k14*(Ca_basal^3)/(Km^3 + Ca_basal^3))*(Ptot - P);

% Solve for nullclines
nullcline1 = solve(dpK_dt == 0, P);  % Solve dpK/dt = 0 for P
nullcline2 = solve(dP_dt == 0, pK);  % Solve dP/dt = 0 for pK

% Convert to MATLAB functions for plotting
nullcline1_func = matlabFunction(nullcline1, 'Vars', pK);
nullcline2_func = matlabFunction(nullcline2, 'Vars', P);

% Generate data for plotting
pK_vals = linspace(0, 20, 50000);
P_vals_nullcline1 = arrayfun(nullcline1_func, pK_vals); % P as a function of pK from dpK/dt = 0
P_vals = linspace(0, 20, 50000);
pK_vals_nullcline2 = arrayfun(nullcline2_func, P_vals); % pK as a function of P from dP/dt = 0

% (Optional) Vector field setup (can be commented out if not needed)
% [pK_grid, P_grid] = meshgrid(linspace(-1, 20, 20), linspace(-1, 20, 20));  % Create grid for vector field
% C_grid = arrayfun(@(x) nu(x), pK_grid) .* (Ktot - pK_grid - K0) ./ mu(Ca_basal);
% dpK_grid = k1 .* ((Ktot - pK_grid - C_grid) ./ (Km1 + Ktot - pK_grid - C_grid)) .* pK_grid - k2 .* (pK_grid ./ (Km2 + pK_grid)) .* (P_grid + P0) + k3 * K0 + (k4 .* (Ca_basal^4) .* (Ktot - pK_grid - C_grid)) ./ (Km^4 + Ca_basal^4);
% dP_grid = k11 .* ((Ptot - P_grid) ./ (Km11 + Ptot - P_grid)) .* P_grid - k12 .* (P_grid ./ (Km12 + P_grid)) .* (pK_grid + K0) + k13 * P0 + (k14 .* (Ca_basal^3) ./ (Km^3 + Ca_basal^3)) .* (Ptot - P_grid);
% vector_magnitude = sqrt(dpK_grid.^2 + dP_grid.^2);
% dpK_norm = dpK_grid ./ vector_magnitude;
% dP_norm = dP_grid ./ vector_magnitude;

% Plot 1: x from -1 to 20, y from -1 to 20 (without vector field)
figure;
% quiver(pK_grid, P_grid, dpK_norm, dP_norm, 'k'); % Vector field (removed)
hold on;
plot(pK_vals, P_vals_nullcline1, 'r', 'LineWidth', 2); % Red: dpK/dt = 0
plot(pK_vals_nullcline2, P_vals, 'b', 'LineWidth', 2); % Blue: dP/dt = 0
plot([-1, 20], [0, 0], 'k', 'LineWidth', 1); % x-axis
plot([0, 0], [-1, 20], 'k', 'LineWidth', 1); % y-axis
xlim([-1, 20]);
ylim([-1, 20]);
xlabel('pK');
ylabel('P');
title('Nullcline Plot for pK vs P');
legend('dpK/dt = 0', 'dP/dt = 0');
grid on;
hold off;
saveas(gcf, 'plot.png');

% Plot 2: x from -1 to 20, y from -1 to 1 (unchanged)
figure;
plot(pK_vals, P_vals_nullcline1, 'r', 'LineWidth', 2); % Red: dpK/dt = 0
hold on;
plot(pK_vals_nullcline2, P_vals, 'b', 'LineWidth', 2); % Blue: dP/dt = 0
plot([-1, 20], [0, 0], 'k', 'LineWidth', 1); % x-axis
plot([0, 0], [-1, 1], 'k', 'LineWidth', 1); % y-axis
xlim([-1, 20]);
ylim([-1, 1]);
xlabel('pK');
ylabel('P');
title('Nullcline Plot for pK vs P (Zoomed In: y from -1 to 1)');
legend('dpK/dt = 0', 'dP/dt = 0');
grid on;
hold off;
saveas(gcf, 'plot1.png');

% Plot 3: x from -1 to 1, y from -1 to 20 (unchanged)
figure;
plot(pK_vals, P_vals_nullcline1, 'r', 'LineWidth', 2); % Red: dpK/dt = 0
hold on;
plot(pK_vals_nullcline2, P_vals, 'b', 'LineWidth', 2); % Blue: dP/dt = 0
plot([-1, 1], [0, 0], 'k', 'LineWidth', 1); % x-axis
plot([0, 0], [-1, 20], 'k', 'LineWidth', 1); % y-axis
xlim([-1, 1]);
ylim([-1, 20]);
xlabel('pK');
ylabel('P');
title('Nullcline Plot for pK vs P (Zoomed In: x from -1 to 1)');
legend('dpK/dt = 0', 'dP/dt = 0');
grid on;
hold off;
saveas(gcf, 'plot2.png');








