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

%% Solutions
function F = CES(x_log, params)
    x = exp(x_log);  % Log-to-level transformation
    % Unpack variables
    C = x(1:9); % Consumptions
    % C(1): Consumption of goods from region A in region A
    % C(2): Consumption of goods from region B in region A
    % C(3): Consumption of goods from region M in region A
    % C(4): Consumption of goods from region A in region B
    % C(5): Consumption of goods from region B in region B
    % C(6): Consumption of goods from region M in region B
    % C(7): Consumption of goods from region A in region M
    % C(8): Consumption of goods from region B in region M
    % C(9): Consumption of goods from region M in region M

    L = x(10:12); % Labor allocations in good sector
    % L(1): Labor allocated to the good sector in region A
    % L(2): Labor allocated to the good sector in region B
    % L(3): Labor allocated to the good sector in region M

    LT = x(13:15); % Labor allocations in transportation sector
    % LT(1): Labor allocated to the transportation sector in region A
    % LT(2): Labor allocated to the transportation sector in region B
    % LT(3): Labor allocated to the transportation sector in region M

    T = x(16:33); % Demand for transportation services
    % T(1): Demand for A's services in route AB
    % T(2): Demand for B's services in route AB
    % T(3): Demand for M's services in route AB
    % T(4): Demand for A's services in route AM
    % T(5): Demand for B's services in route AM
    % T(6): Demand for M's services in route AM
    % T(7): Demand for A's services in route BM
    % T(8): Demand for B's services in route BM
    % T(9): Demand for M's services in route BM
    % T(10): Demand for A's services in route BA
    % T(11): Demand for B's services in route BA
    % T(12): Demand for M's services in route BA
    % T(13): Demand for A's services in route MA
    % T(14): Demand for B's services in route MA
    % T(15): Demand for M's services in route MA
    % T(16): Demand for A's services in route MB
    % T(17): Demand for B's services in route MB
    % T(18): Demand for M's services in route MB

    % Normalize W_M = 1
    W_A = x(34); % Wage in A
    W_B = x(35); % Wage in B

    % Unpack parameters
    Z_A = params.Z_A;
    Z_B = params.Z_B;
    Z_M = params.Z_M;
    Z_TA = params.Z_TA;
    Z_TB = params.Z_TB;
    Z_TM = params.Z_TM;
    d = params.d;
    rho1 = params.rho1;
    rho2 = params.rho2;
    tau = params.tau;
    tilde_tau = params.tilde_tau;
    sigma = params.sigma;
    mu = params.mu;
    chi = params.chi;
    lambda_o = params.lambda_o;
    lambda_d = params.lambda_d;
    t_A = params.t_A;
    t_B = params.t_B;
    t_M = params.t_M;
    L_US_total = params.L_US_total;
    L_M_total = params.L_M_total;

    W_M = params.W_M;
    P_M = params.P_M;
    mc_M = params.mc_M;
    P_MM = params.P_MM;

    % Delivered prices
    P_A = W_A / Z_A;
    P_AA = P_A;
    P_B = W_A / Z_B;
    P_BB = P_B;
    mc_A = struct( ...
        'AB', d.AB ^ rho1 * tau.A_AB * tilde_tau.A_AB * W_A / Z_TA, ...
        'BA', d.BA ^ rho1 * tau.A_BA * tilde_tau.A_BA * W_A / Z_TA, ...
        'AM', d.AM ^ rho1 * tau.A_AM * tilde_tau.A_AM * W_A / Z_TA, ...
        'MA', d.MA ^ rho1 * tau.A_MA * tilde_tau.A_MA * W_A / Z_TA, ...
        'BM', d.BM ^ rho1 * tau.A_BM * tilde_tau.A_BM * W_A / Z_TA, ...
        'MB', d.MB ^ rho1 * tau.A_MB * tilde_tau.A_MB * W_A / Z_TA ...
        );
    mc_B = struct( ...
        'AB', d.AB ^ rho1 * tau.B_AB * tilde_tau.B_AB * W_B / Z_TB, ...
        'BA', d.BA ^ rho1 * tau.B_BA * tilde_tau.B_BA * W_B / Z_TB, ...
        'AM', d.AM ^ rho1 * tau.B_AM * tilde_tau.B_AM * W_B / Z_TB, ...
        'MA', d.MA ^ rho1 * tau.B_MA * tilde_tau.B_MA * W_B / Z_TB, ...
        'BM', d.BM ^ rho1 * tau.B_BM * tilde_tau.B_BM * W_B / Z_TB, ...
        'MB', d.MB ^ rho1 * tau.B_MB * tilde_tau.B_MB * W_B / Z_TB ...
        );
    P_T = struct( ...
        'AB', (lambda_o * mc_A.AB ^ chi + lambda_d * mc_B.AB ^ chi + (1 - lambda_o - lambda_d) * mc_M.AB ^ chi) ^ (1 / chi), ...
        'BA', (lambda_o * mc_B.BA ^ chi + lambda_d * mc_A.BA ^ chi + (1 - lambda_o - lambda_d) * mc_M.BA ^ chi) ^ (1 / chi), ...
        'AM', (lambda_o * mc_A.AM ^ chi + lambda_d * mc_M.AM ^ chi + (1 - lambda_o - lambda_d) * mc_B.AM ^ chi) ^ (1 / chi), ...
        'MA', (lambda_o * mc_M.MA ^ chi + lambda_d * mc_A.MA ^ chi + (1 - lambda_o - lambda_d) * mc_B.MA ^ chi) ^ (1 / chi), ...
        'BM', (lambda_o * mc_B.BM ^ chi + lambda_d * mc_M.BM ^ chi + (1 - lambda_o - lambda_d) * mc_A.BM ^ chi) ^ (1 / chi), ...
        'MB', (lambda_o * mc_M.MB ^ chi + lambda_d * mc_B.MB ^ chi + (1 - lambda_o - lambda_d) * mc_A.MB ^ chi) ^ (1 / chi) ...
        );
    P_AB = P_A + t_A * P_T.AB;
    P_BA = P_B + t_B * P_T.BA;
    P_AM = P_A + t_A * P_T.AM;
    P_MA = P_M + t_M * P_T.MA;
    P_BM = P_B + t_B * P_T.BM;
    P_MB = P_M + t_M * P_T.MB;


    % Equations
    % Utility maximization
    F(1) = (P_AA + ((1 - mu) / 2 / mu) * P_AA ^ sigma * P_BA ^ (1 - sigma) + ((1 - mu) / 2 / mu) * P_AA ^ sigma * P_MA ^ (1 - sigma)) * C(1) - W_A * (L(1) + LT(1));
    F(2) = (P_AA + ((1 - mu) / 2 / mu) * P_AA ^ sigma * P_BA ^ (1 - sigma) + ((1 - mu) / 2 / mu) * P_AA ^ sigma * P_MA ^ (1 - sigma)) * C(2) - ((1 - mu) / 2 / mu) * (P_AA / P_BA) ^ sigma * W_A * (L(1) + LT(1));
    F(3) = (P_AA + ((1 - mu) / 2 / mu) * P_AA ^ sigma * P_BA ^ (1 - sigma) + ((1 - mu) / 2 / mu) * P_AA ^ sigma * P_MA ^ (1 - sigma)) * C(3) - ((1 - mu) / 2 / mu) * (P_AA / P_MA) ^ sigma * W_A * (L(1) + LT(1));
    F(4) = (P_BB + ((1 - mu) / 2 / mu) * P_BB ^ sigma * P_AB ^ (1 - sigma) + ((1 - mu) / 2 / mu) * P_BB ^ sigma * P_MB ^ (1 - sigma)) * C(4) - ((1 - mu) / 2 / mu) * (P_BB / P_AB) ^ sigma * W_B * (L(2) + LT(2));
    F(5) = (P_BB + ((1 - mu) / 2 / mu) * P_BB ^ sigma * P_AB ^ (1 - sigma) + ((1 - mu) / 2 / mu) * P_BB ^ sigma * P_MB ^ (1 - sigma)) * C(5) - W_B * (L(2) + LT(2));
    F(6) = (P_BB + ((1 - mu) / 2 / mu) * P_BB ^ sigma * P_AB ^ (1 - sigma) + ((1 - mu) / 2 / mu) * P_BB ^ sigma * P_MB ^ (1 - sigma)) * C(6) - ((1 - mu) / 2 / mu) * (P_BB / P_MB) ^ sigma * W_B * (L(2) + LT(2));
    F(7) = (P_MM + ((1 - mu) / 2 / mu) * P_MM ^ sigma * P_AM ^ (1 - sigma) + ((1 - mu) / 2 / mu) * P_MM ^ sigma * P_BM ^ (1 - sigma)) * C(7) - ((1 - mu) / 2 / mu) * (P_MM / P_AM) ^ sigma * W_M * L_M_total;
    F(8) = (P_MM + ((1 - mu) / 2 / mu) * P_MM ^ sigma * P_AM ^ (1 - sigma) + ((1 - mu) / 2 / mu) * P_MM ^ sigma * P_BM ^ (1 - sigma)) * C(8) - ((1 - mu) / 2 / mu) * (P_MM / P_BM) ^ sigma * W_M * L_M_total;
    F(9) = (P_MM + ((1 - mu) / 2 / mu) * P_MM ^ sigma * P_AM ^ (1 - sigma) + ((1 - mu) / 2 / mu) * P_MM ^ sigma * P_BM ^ (1 - sigma)) * C(9) - W_M * L_M_total;

    % Demand for transportation services
    F(10) = lambda_o * (mc_A.AB / P_T.AB) ^ (-chi) * t_A * C(4) - T(1);
    F(11) = lambda_d * (mc_B.AB / P_T.AB) ^ (-chi) * t_A * C(4) - T(2);
    F(12) = (1 - lambda_o - lambda_d) * (mc_M.AB / P_T.AB) ^ (-chi) * t_A * C(4) - T(3);
    F(13) = lambda_o * (mc_A.AM / P_T.AM) ^ (-chi) * t_A * C(7) - T(4);
    F(14) = (1 - lambda_o - lambda_d) * (mc_B.AM / P_T.AM) ^ (-chi) * t_A * C(7) - T(5);
    F(15) = lambda_d * (mc_M.AM / P_T.AM) ^ (-chi) * t_A * C(7) - T(6);
    F(16) = (1 - lambda_o - lambda_d) * (mc_A.BM / P_T.BM) ^ (-chi) * t_B * C(8) - T(7);
    F(17) = lambda_o * (mc_B.BM / P_T.BM) ^ (-chi) * t_B * C(8) - T(8);
    F(18) = lambda_d * (mc_M.BM / P_T.BM) ^ (-chi) * t_B * C(8) - T(9);
    F(19) = lambda_d * (mc_A.BA / P_T.BA) ^ (-chi) * t_B * C(2) - T(10);
    F(20) = lambda_o * (mc_B.BA / P_T.BA) ^ (-chi) * t_B * C(2) - T(11);
    F(21) = (1 - lambda_o - lambda_d) * (mc_M.BA / P_T.BA) ^ (-chi) * t_B * C(2) - T(12);
    F(22) = lambda_d * (mc_A.MA / P_T.MA) ^ (-chi) * t_M * C(3) - T(13);
    F(23) = (1 - lambda_o - lambda_d) * (mc_B.MA / P_T.MA) ^ (-chi) * t_M * C(3) - T(14);
    F(24) = lambda_o * (mc_M.MA / P_T.MA) ^ (-chi) * t_M * C(3) - T(15);
    F(25) = (1 - lambda_o - lambda_d) * (mc_A.MB / P_T.MB) ^ (-chi) * t_M * C(6) - T(16);
    F(26) = lambda_d * (mc_B.MB / P_T.MB) ^ (-chi) * t_M * C(6) - T(17);
    F(27) = lambda_o * (mc_M.MB / P_T.MB) ^ (-chi) * t_M * C(6) - T(18);
    

    % Good market clearing
    F(28) = Z_A * L(1) - (C(1) + C(4) + C(7));
    F(29) = Z_B * L(2) - (C(2) + C(5) + C(8));
    F(30) = Z_M * L(3) - (C(3) + C(6) + C(9));

    % Transportation market clearing
    F(31) = Z_TA * LT(1) - T(1) - T(4) - T(7) - T(10) - T(13) - T(16);
    F(32) = Z_TB * LT(2) - T(2) - T(5) - T(8) - T(11) - T(14) - T(17);
    F(33) = Z_TM * LT(3) - T(3) - T(6) - T(9) - T(12) - T(15) - T(18);

    % Labor market clearing
    F(34) = L_US_total - (L(1) + LT(1)) - (L(2) + LT(2));
    % F(35) = L_M_total - (L(3) + LT(3));

    % Spatial equilibrium
    if sigma == 1
        F(36) = C(1) ^ mu * C(2) ^ ((1 - mu) / 2) * C(3) ^ ((1 - mu) / 2) - C(4) ^ ((1 - mu) / 2) * C(5) ^ mu * C(6) ^ ((1 - mu) / 2) ;
    else
        F(36) = (mu ^ (1 / sigma) * C(1) ^ ((sigma - 1) / sigma) + ((1 - mu) / 2) ^ (1 / sigma) * C(2) ^ ((sigma - 1) / sigma) + ((1 - mu) / 2) ^ (1 / sigma) * C(3) ^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1)) - ...
            ((((1 - mu) / 2) ^ (1 / sigma) * C(4) ^ ((sigma - 1) / sigma) + mu ^ (1 / sigma) * C(5) ^ ((sigma - 1) / sigma) + ((1 - mu) / 2) ^ (1 / sigma) * C(6) ^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1)));
    end    

end



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


%% Plotting utiltities
% Create a tiled layout for two plots side by side
figure;

% Plotting the utilities 
nexttile; % Move to the first tile
hold on;
plot(Z_TM_range, U_A, 'b-', 'DisplayName', 'Utility in the US', 'LineWidth', 1.5);
plot(Z_TM_range, U_M, 'r-', 'DisplayName', 'Utility in Mexico', 'LineWidth', 1.5);

% Add labels and title
xlabel('Z_{TM} (Transportation Productivity in Mexico)', 'FontSize', 12);
ylabel('Utility', 'FontSize', 12);
title(sprintf('\\sigma = %.1f, \\mu = %.1f, \\chi = %.1f, \\lambda_o = %.1f, \\lambda_d = %.1f', sigma, mu, chi, lambda_o, lambda_d), 'FontSize', 14);

% Add legend
legend('Location', 'best', 'FontSize', 10);

% Add grid
grid on;
hold off;

% Save the figure
filename = sprintf('sigma_%.1f_mu_%.1f_chi_%.1f_lambda_o_%.1f_lambda_d_%.1f.png', sigma, mu, chi, lambda_o, lambda_d);

% Save the current figure
%saveas(gcf, filename);