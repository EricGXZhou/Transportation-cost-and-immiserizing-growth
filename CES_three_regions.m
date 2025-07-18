%% NOTE:
% This script analyzes a three-region general equilibrium model with CES preferences and transportation costs.
% It solves for equilibrium allocations and utilities as transportation productivity in Mexico varies.
% The model includes:
%   - Three regions (A, B, M) with goods and transportation sectors
%   - CES utility and transportation aggregation
%   - Endogenous labor allocation and transportation demand
%   - Parameter sweep over Z_TM (transportation productivity in Mexico)
%   - Plots utilities for US and Mexico as Z_TM changes
% The main function (CES) is defined in this file for system solving.

clear
clc


%% Parameters
% Productivity of good sector
Z_US = 5;
Z_A = Z_US; 
Z_B = Z_US; 
Z_M = 5;

% Productivity of transportation sector
Z_TUS = 5;
Z_TA = Z_TUS; 
Z_TB = Z_TUS;

% Distance matrix (symmetric)
d = struct( ...
    'AA', 1, 'BB', 1, 'MM', 1, ...
    'AB', 1, 'BA', 1, ...
    'AM', 1, 'MA', 1, ...
    'BM', 1, 'MB', 1 ...
);

rho1 = 1; % distance exponent
rho2 = 1; % commuting cost exponent

% Commuting costs function
commute = @(x, y) max(d.(x)^rho2, d.(y)^rho2);

% Commuting costs for i in route od
tau = struct( ...
    'A_AB', commute('AA', 'AB'), ...
    'B_AB', commute('BA', 'BB'), ...
    'M_AB', commute('MA', 'MB'), ...
    'A_BA', commute('BA', 'AA'), ...
    'B_BA', commute('BB', 'BA'), ...
    'M_BA', commute('MB', 'MA'), ...
    'A_AM', commute('AA', 'AM'), ...
    'M_AM', commute('MA', 'MM'), ...
    'B_AM', commute('BA', 'BM'), ...
    'A_MA', commute('AM', 'AA'), ...
    'M_MA', commute('MM', 'MA'), ...
    'B_MA', commute('BM', 'BA'), ...
    'B_BM', commute('BB', 'BM'), ...
    'M_BM', commute('MB', 'MM'), ...
    'A_BM', commute('AB', 'AM'), ...
    'B_MB', commute('BM', 'BB'), ...
    'M_MB', commute('MM', 'MB'), ...
    'A_MB', commute('AM', 'AB') ...
);

% policy parameters for i in route od
tilde_tau = struct( ...
    'A_AB', 1, ...
    'B_AB', 1, ...
    'M_AB', 1, ...
    'A_BA', 1, ...
    'B_BA', 1, ...
    'M_BA', 1, ...
    'A_AM', 1, ...
    'M_AM', 1, ...
    'B_AM', 1, ...
    'A_MA', 1, ...
    'M_MA', 1, ...
    'B_MA', 1, ...
    'B_BM', 1, ...
    'M_BM', 1, ...
    'A_BM', 1, ...
    'B_MB', 1, ...
    'M_MB', 1, ...
    'A_MB', 1 ...
);

sigma = 1; % elasticity of substitution for consumption
mu = 0.5; % home bias parameter for consumption
chi = 1; % elasticity of substitution for transportation services
lambda_o = 0.5; % origin bias for transportation services
lambda_d = 0.25; % destination bias for transportation services

t_A = 0.1;
t_B = 0.1;
t_M = 0.1;
L_US_total = 2;
L_M_total = 1;


% The range of Z_TM
Z_TM_range = linspace(1,10,100); % Values of Z_TM from 1 to 10

% Results
solutions = cell(length(Z_TM_range), 1);
utilities = zeros(length(Z_TM_range), 3); % Columns: U_A, U_B, U_M

% Diagnostic storage
wages = zeros(length(Z_TM_range), 3); % W_A, W_B, W_M
labors = zeros(length(Z_TM_range), 3); % L_A, L_B, L_M
trans_M_AB = zeros(length(Z_TM_range), 1);
trans_M_AM = zeros(length(Z_TM_range), 1);
trans_M_BM = zeros(length(Z_TM_range), 1);

%% Solutions
% The CES function is now in a separate file: CES.m

for i = 1:length(Z_TM_range)
    % Normalization
    W_M = 1;
    P_M = 1 / Z_M;
    mc_M = struct( ...
        'AB', d.AB ^ rho1 * tau.M_AB * tilde_tau.M_AB * W_M / Z_TM_range(i), ...
        'BA', d.BA ^ rho1 * tau.M_BA * tilde_tau.M_BA * W_M / Z_TM_range(i), ...
        'AM', d.AM ^ rho1 * tau.M_AM * tilde_tau.M_AM * W_M / Z_TM_range(i), ...
        'MA', d.MA ^ rho1 * tau.M_MA * tilde_tau.M_MA * W_M / Z_TM_range(i), ...
        'BM', d.BM ^ rho1 * tau.M_BM * tilde_tau.M_BM * W_M / Z_TM_range(i), ...
        'MB', d.MB ^ rho1 * tau.M_MB * tilde_tau.M_MB * W_M / Z_TM_range(i) ...
        );
    P_MM = P_M;

    % Set parameters
    params.Z_A = Z_A;
    params.Z_B = Z_B;
    params.Z_M = Z_M;
    params.Z_TA = Z_TA;
    params.Z_TB = Z_TB;
    params.Z_TM = Z_TM_range(i);
    params.d = d;
    params.rho1 = rho1;
    params.rho2 = rho2;
    params.tau = tau;
    params.tilde_tau = tilde_tau;
    params.sigma = sigma;
    params.mu = mu;
    params.chi = chi;
    params.lambda_o = lambda_o;
    params.lambda_d = lambda_d;
    params.t_A = t_A;
    params.t_B = t_B;
    params.t_M = t_M;
    params.L_US_total = L_US_total;
    params.L_M_total = L_M_total;

    params.W_M = W_M;
    params.P_M = P_M;
    params.mc_M = mc_M;
    params.P_MM = P_MM;

    % Initial guesses for the variables
    x0 = log([ones(9, 1); % Consumptions
          ones(3, 1); % Labor allocations in good sector
          ones(3, 1); % Labor allocations in transportation sector
          ones(18, 1); % Demand for transportation services
          ones(2, 1)]); % A and B's wage

    % Solve the system of equations
    options = optimoptions('fsolve', ...
    'Display', 'iter', ...
    'TolFun', 1e-16, ...
    'TolX', 1e-16, ...
    'MaxIterations', 10000, ...
    'MaxFunctionEvaluations', 50000);
    x_sol = fsolve(@(x) CES(x, params), x0, options);
    x = exp(x_sol);
    % Store the solution
    solutions{i} = [x(1:35); W_M];
    
    % Store key diagnostics for debugging
    wages(i, :) = [x(34), x(35), W_M]; % W_A, W_B, W_M
    labors(i, :) = [x(10), x(11), x(12)]; % L_A, L_B, L_M
    trans_M_AB(i) = x(18); % T_M_AB (Mexico transport on AB route)
    trans_M_AM(i) = x(21); % T_M_AM (Mexico transport on AM route)
    trans_M_BM(i) = x(24); % T_M_BM (Mexico transport on BM route)

    
    % Calculate utilities
    if sigma == 1
        U_A = x(1) ^ mu * x(2) ^ ((1 - mu) / 2) * x(3) ^ ((1 - mu) / 2);
        U_B = x(4) ^ ((1 - mu) / 2) * x(5) ^ mu * x(6) ^ ((1 - mu) / 2);
        U_M = x(7) ^ ((1 - mu) / 2) * x(8) ^ ((1 - mu) / 2) * x(9) ^ mu;
    else
        U_A = (mu ^ (1 / sigma) * x(1) ^ ((sigma - 1) / sigma) + ((1 - mu) / 2) ^ (1 / sigma) * x(2) ^ ((sigma - 1) / sigma) + ((1 - mu) / 2) ^ (1 / sigma) * x(3) ^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1));
        U_B = (((1 - mu) / 2) ^ (1 / sigma) * x(4) ^ ((sigma - 1) / sigma) + mu ^ (1 / sigma) * x(5) ^ ((sigma - 1) / sigma) + ((1 - mu) / 2) ^ (1 / sigma) * x(6) ^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1));
        U_M = (((1 - mu) / 2) ^ (1 / sigma) * x(7) ^ ((sigma - 1) / sigma) + ((1 - mu) / 2) ^ (1 / sigma) * x(8) ^ ((sigma - 1) / sigma) + mu ^ (1 / sigma) * x(9) ^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1));
    end

    utilities(i, :) = [U_A, U_B, U_M];
end

% Label the solutions
label = ["C_AA";"C_BA";"C_MA";"C_AB";"C_BB";"C_MB";"C_AM";"C_BM";"C_MM";"L_A";"L_B";"L_M";"L_TA";"L_TB";"L_TM";...
    "T_A_AB";"T_B_AB";"T_M_AB";"T_A_AM";"T_B_AM";"T_M_AM";"T_A_BM";"T_B_BM";"T_M_BM";"T_A_BA";"T_B_BA";"T_M_BA";"T_A_MA";"T_B_MA";"T_M_MA";"T_A_MB";"T_B_MB";"T_M_MB";...
    "W_A";"W_B";"W_M"];
for i = 1:length(Z_TM_range)
    solutions{i} = [label, solutions{i}];
end



% Extract utilities for each region
U_A = utilities(:, 1)'; % Utilities in A
U_B = utilities(:, 2)'; % Utilities in B
U_M = utilities(:, 3)'; % Utilities in Mexico

% Extract wages and labor allocations
W_A = wages(:, 1)';
W_B = wages(:, 2)';
W_M = wages(:, 3)';
L_A = labors(:, 1)';
L_B = labors(:, 2)';
L_M = labors(:, 3)';
T_M_AB = trans_M_AB;
T_M_AM = trans_M_AM;
T_M_BM = trans_M_BM;

%% Plotting utiltities and diagnostics
figure;
tiledlayout(2,2);

% 1. Utilities
nexttile;
hold on;
plot(Z_TM_range, U_A, 'b-', 'DisplayName', 'Utility in Texas', 'LineWidth', 1.5);
plot(Z_TM_range, U_B, 'g-', 'DisplayName', 'Utility in California', 'LineWidth', 1.5);
plot(Z_TM_range, U_M, 'r-', 'DisplayName', 'Utility in Mexico', 'LineWidth', 1.5);
xlabel('Z_{TM}'); ylabel('Utility'); title('Utilities'); legend; grid on; hold off;

% 2. Wages
nexttile;
hold on;
plot(Z_TM_range, W_A, 'b-', 'DisplayName', 'W_A (Texas)', 'LineWidth', 1.5);
plot(Z_TM_range, W_B, 'g-', 'DisplayName', 'W_B (California)', 'LineWidth', 1.5);
plot(Z_TM_range, W_M, 'r-', 'DisplayName', 'W_M (Mexico)', 'LineWidth', 1.5);
xlabel('Z_{TM}'); ylabel('Wage'); title('Wages'); legend; grid on; hold off;

% 3. Labor allocations
nexttile;
hold on;
plot(Z_TM_range, L_A, 'b-', 'DisplayName', 'L_A (Texas)', 'LineWidth', 1.5);
plot(Z_TM_range, L_B, 'g-', 'DisplayName', 'L_B (California)', 'LineWidth', 1.5);
plot(Z_TM_range, L_M, 'r-', 'DisplayName', 'L_M (Mexico)', 'LineWidth', 1.5);
xlabel('Z_{TM}'); ylabel('Labor in Goods'); title('Labor Allocations'); legend; grid on; hold off;

% 4. Mexico transport flows
nexttile;
hold on;
plot(Z_TM_range, T_M_AB, 'k-', 'DisplayName', 'T_{M,AB}', 'LineWidth', 1.5);
plot(Z_TM_range, T_M_AM, 'm-', 'DisplayName', 'T_{M,AM}', 'LineWidth', 1.5);
plot(Z_TM_range, T_M_BM, 'c-', 'DisplayName', 'T_{M,BM}', 'LineWidth', 1.5);
xlabel('Z_{TM}'); ylabel('Mexico Transport Flows'); title('Mexico Transport Flows'); legend; grid on; hold off;

% Save the figure
filename = sprintf('sigma_%.1f_mu_%.1f_chi_%.1f_lambda_o_%.1f_lambda_d_%.1f.png', sigma, mu, chi, lambda_o, lambda_d);

% Save the current figure
%saveas(gcf, filename);