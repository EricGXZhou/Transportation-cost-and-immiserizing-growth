clear
clc

%% NOTE:
% This code extends the two-country equilibrium model to include CES aggregation 
% for transportation services. It analyzes how changes in transportation productivity 
% in Mexico (Z_TM) affect the equilibrium utilities of the US and Mexico, 
% accounting for both countries' contributions to transportation with CES preferences.
% The model solves a system of nonlinear equations for each value of Z_TM 
% and plots the resulting utilities for both countries.

%% Parameters
Z_A = 5;
Z_M = 5; 
Z_TA = 5;
d_AM = 1; % d_AM = d_MA (distance)
d_MA = d_AM;
d_AA = 1;
d_MM = 1;
rho1 = 1; % parameter for distance
tau_A_AM = max(d_AA,d_AM); % commuting costs for A in route AM
tau_M_AM = max(d_MA,d_MM); % commuting costs for M in route AM
tau_M_MA = max(d_MA,d_MM); % commuting costs for M in route MA
tau_A_MA = max(d_AM,d_AA); % commuting costs for A in route MA
rho2 = 1; % parameter for commuting costs
tilde_tau_A_AM = 1; % policy for A in route AM
tilde_tau_M_AM = 1; % policy for M in route AM
tilde_tau_M_MA = 1; % policy for M in route MA
tilde_tau_A_MA = 1; % policy for A in route MA

sigma = 1; % elasticity of substitution for consumption
mu = 0.5; % home bias parameter for consumption
chi = 1; % elasticity of substitution for transportation services
lambda = 1; % origin bias for transportation services

t_A = 0.1;
t_M = 0.1;
L_A_total = 1;
L_M_total = 1;


% The range of Z_TM
Z_TM_range = linspace(1,10,100); % Values of Z_TM from 1 to 10

% Results
solutions = cell(length(Z_TM_range), 1);
utilities = zeros(length(Z_TM_range), 2); % Columns: U_A, U_M

%% Solutions
function F = CES(x_log, params)
    x = exp(x_log);  % Log-to-level transformation
    % Unpack variables
    C = x(1:4); % Consumptions
    % C(1): Consumption of goods from country A in country A
    % C(2): Consumption of goods from country M in country A
    % C(3): Consumption of goods from country A in country M
    % C(4): Consumption of goods from country M in country M

    L = x(5:8); % Labor allocations
    % L(1): Labor allocated to the good sector in country A
    % L(2): Labor allocated to the good sector in country M
    % L(3): Labor allocated to the transportation sector in country A
    % L(4): Labor allocated to the transportation sector in country M

    T = x(9:12); % Demand for transportation services
    % T(1): Demand for A's services in route AM
    % T(2): Demand for M's services in route AM
    % T(3): Demand for A's services in route MA
    % T(4): Demand for M's services in route MA

    W_A = x(13); % Wage in country A

    % Unpack parameters
    Z_A = params.Z_A;
    Z_M = params.Z_M;
    Z_TA = params.Z_TA;
    Z_TM = params.Z_TM;
    t_A = params.t_A;
    t_M = params.t_M;
    d_AM = params.d_AM;
    d_MA = params.d_MA;
    tau_A_AM = params.tau_A_AM;
    tau_A_MA = params.tau_A_MA;
    tilde_tau_A_AM = params.tilde_tau_A_AM;
    tilde_tau_A_MA = params.tilde_tau_A_MA;
    sigma = params.sigma;
    mu = params.mu;
    chi = params.chi;
    lambda = params.lambda;
    L_A_total = params.L_A_total;
    L_M_total = params.L_M_total;
    W_M = params.W_M;
    P_M = params.P_M;
    mc_M_AM = params.mc_M_AM;
    mc_M_MA = params.mc_M_MA;
    P_MM = params.P_MM;

    % Delivered prices
    P_A = W_A / Z_A;
    P_AA = P_A;
    mc_A = W_A / Z_TA;
    mc_A_AM = d_AM * tau_A_AM * tilde_tau_A_AM * mc_A;
    mc_A_MA = d_MA * tau_A_MA * tilde_tau_A_MA * mc_A;
    P_T_AM = (lambda * mc_A_AM ^ chi + (1 - lambda) * mc_M_AM ^ chi) ^ (1 / chi);
    P_T_MA = (lambda * mc_M_MA ^ chi + (1 - lambda) * mc_A_MA ^ chi) ^ (1 / chi);
    P_AM = P_A + t_A * P_T_AM;
    P_MA = P_M + t_M * P_T_MA;


    % Equations
    % Utility maximization
    F(1) = (P_AA + ((1 - mu) / mu) * P_AA ^ sigma * P_MA ^ (1 - sigma)) * C(1) - W_A * L_A_total;
    F(2) = (P_AA + ((1 - mu) / mu) * P_AA ^ sigma * P_MA ^ (1 - sigma)) * C(2) - ((1 - mu) / mu) * (P_AA / P_MA) ^ sigma * W_A * L_A_total;
    F(3) = (P_MM + ((1 - mu) / mu) * P_MM ^ sigma * P_AM ^ (1 - sigma)) * C(3) - ((1 - mu) / mu) * (P_MM / P_AM) ^ sigma * W_M * L_M_total;
    F(4) = (P_MM + ((1 - mu) / mu) * P_MM ^ sigma * P_AM ^ (1 - sigma)) * C(4) - W_M * L_M_total;

    % Demand for transportation services
    F(5) = lambda * (mc_A_AM / P_T_AM) ^ (-chi) * t_A * C(3) - T(1);
    F(6) = (1 - lambda) * (mc_M_AM / P_T_AM) ^ (-chi) * t_M * C(3) - T(2);
    F(7) = (1 - lambda) * (mc_A_MA / P_T_MA) ^ (-chi) * t_A * C(2) - T(3);
    F(8) = lambda * (mc_M_MA / P_T_MA) ^ (-chi) * t_M * C(2) - T(4);

    % Good market clearing
    F(9) = Z_A * L(1) - (C(1) + C(3));
    F(10) = Z_M * L(2) - (C(2) + C(4));

    % Transportation market clearing
    F(11) = Z_TA * L(3) - T(1) - T(3);
    F(12) = Z_TM * L(4) - T(2) - T(4);

    % Labor market clearing
    F(13) = L_A_total - (L(1) + L(3));
    F(14) = L_M_total - (L(2) + L(4));

end


for i = 1:length(Z_TM_range)
    % Normalization
    W_M = 1;
    P_M = 1 / Z_M;
    mc_M = 1 / Z_TM_range(i);
    mc_M_AM = d_AM * tau_M_AM * tilde_tau_M_AM * mc_M;
    mc_M_MA = d_MA * tau_M_MA * tilde_tau_M_MA * mc_M;
    P_MM = P_M;


    % Set parameters
    params.Z_A = Z_A;
    params.Z_M = Z_M;
    params.Z_TA = Z_TA;
    params.Z_TM = Z_TM_range(i);
    params.t_A = t_A;
    params.t_M = t_M;
    params.d_AM = d_AM;
    params.d_MA = d_MA;
    params.tau_A_AM = tau_A_AM;
    params.tau_A_MA = tau_A_MA;
    params.tilde_tau_A_AM = tilde_tau_A_AM; 
    params.tilde_tau_A_MA = tilde_tau_A_MA;
    params.sigma = sigma;
    params.mu = mu;
    params.chi = chi;
    params.lambda = lambda;
    params.L_A_total = L_A_total;
    params.L_M_total = L_M_total;
    params.W_M = W_M;
    params.P_M = P_M;
    params.mc_M_AM = mc_M_AM;
    params.mc_M_MA = mc_M_MA;
    params.P_MM = P_MM;

    % Initial guesses for the variables
    x0 = log([ones(4, 1); % Consumptions
          ones(4, 1); % Labor allocations
          ones(4, 1); % Demand for transportation services
          1]); % A's wage


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
    solutions{i} = [x(1:13); W_M];
    
    % Calculate utilities
    if sigma == 1
        U_A = x(1) ^ mu * x(2) ^ (1 - mu);
        U_M = x(4) ^ mu * x(3) ^ (1 - mu);
    else
        U_A = (mu ^ (1 / sigma) * x(1) ^ ((sigma - 1) / sigma) + (1 - mu) ^ (1 / sigma) * x(2) ^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1));
        U_M = (mu ^ (1 / sigma) * x(4) ^ ((sigma - 1) / sigma) + (1 - mu) ^ (1 / sigma) * x(3) ^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1));
    end

    utilities(i, :) = [U_A, U_M];
end

% Label the solutions
label = ["C_AA";"C_MA";"C_AM";"C_MM";"L_A";"L_M";"L_TA";"L_TM";"T_A_AM";"T_M_AM";"T_A_MA";"T_M_MA";"W_A";"W_M"];
for i = 1:length(Z_TM_range)
    solutions{i} = [label, solutions{i}];
end


% Extract utilities for each country
U_A = utilities(:, 1)'; % Utilities in A
U_M = utilities(:, 2)'; % Utilities in Mexico


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
title(sprintf('\\sigma = %.1f, \\mu = %.1f, \\chi = %.1f, \\lambda = %.1f', sigma, mu, chi, lambda), 'FontSize', 14);

% Add legend
legend('Location', 'best', 'FontSize', 10);

% Add grid
grid on;
hold off;

% Save the figure
filename = sprintf('sigma_%.1f_mu_%.1f_chi_%.1f_lambda_%.1f.png', sigma, mu, chi, lambda);

% Save the current figure
%saveas(gcf, filename);