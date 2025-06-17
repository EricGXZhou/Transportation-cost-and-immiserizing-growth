clear
clc

%% NOTE:
% This code analyzes how changes in transportation productivity in Mexico (Z_TM)
% affect the equilibrium utilities of the US and Mexico under two scenarios:
% 1) Counterfactual world (no border frictions)
% 2) Factual world (with border frictions)
% It solves a system of nonlinear equations for each scenario and plots the resulting utilities.

%% Parameters
Z_A = 5;
Z_M = 5; 
Z_TA = 5;
% Z_TM = 1;
t_A = 0.1;
t_M = 0.1;
sigma = 1; % elasticity of substitution
mu = 0.5; % home bias parameter
L_A_total = 1;
L_M_total = 1;


% The range of Z_TM
Z_TM_range = linspace(1,10,100); % Values of Z_TM from 1 to 10

% Results
base_solutions = cell(length(Z_TM_range), 1);
base_utilities = zeros(length(Z_TM_range), 2); % Columns: U_A, U_M
border_solutions = cell(length(Z_TM_range), 1);
border_utilities = zeros(length(Z_TM_range), 2); % Columns: U_A, U_M

%% Counterfactual World
% Define the function for the system of equations
function F = counterfactual(x, params)
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

    W = x(9); % Wage in country A

    % Unpack parameters
    Z_A = params.Z_A;
    Z_M = params.Z_M;
    Z_TA = params.Z_TA;
    Z_TM = params.Z_TM;
    t_A = params.t_A;
    t_M = params.t_M;
    sigma = params.sigma;
    mu = params.mu;
    L_A_total = params.L_A_total;
    L_M_total = params.L_M_total;
    W_M = params.W_M;
    P_M = params.P_M;
    P_TM = params.P_TM;
    P_MA = params.P_MA;
    P_MM = params.P_MM;

    % Delivered prices
    P_A = W / Z_A;
    P_TA = W / Z_TA;
    P_AA = P_A;
    P_AM = P_A + t_A * P_TA;


    % Equations
    % Utility maximization
    F(1) = (P_AA + ((1 - mu) / mu) * P_AA ^ sigma * P_MA ^ (1 - sigma)) * C(1) - W * L_A_total;
    F(2) = (P_AA + ((1 - mu) / mu) * P_AA ^ sigma * P_MA ^ (1 - sigma)) * C(2) - ((1 - mu) / mu) * (P_AA / P_MA) ^ sigma * W * L_A_total;
    F(3) = (P_MM + ((1 - mu) / mu) * P_MM ^ sigma * P_AM ^ (1 - sigma)) * C(3) - ((1 - mu) / mu) * (P_MM / P_AM) ^ sigma * W_M * L_M_total;
    F(4) = (P_MM + ((1 - mu) / mu) * P_MM ^ sigma * P_AM ^ (1 - sigma)) * C(4) - W_M * L_M_total;

    % Good market clearing
    F(5) = Z_A * L(1) - (C(1) + C(3));
    F(6) = Z_M * L(2) - (C(2) + C(4));

    % Transportation market clearing
    F(7) = Z_TA * L(3) - (t_A * C(3));
    F(8) = Z_TM * L(4) - (t_M * C(2));

    % Labor market clearing
    F(9) = L_A_total - (L(1) + L(3));
    F(10) = L_M_total - (L(2) + L(4));

end




for i = 1:length(Z_TM_range)
    % Normalization
    W_M = 1;
    P_M = 1 / Z_M;
    P_TM = 1 / Z_TM_range(i);

    % Delivered prices
    P_MA = P_M + t_M * P_TM;
    P_MM = P_M;


    % Set parameters
    params.Z_A = Z_A;
    params.Z_M = Z_M;
    params.Z_TA = Z_TA;
    params.Z_TM = Z_TM_range(i);
    params.t_A = t_A;
    params.t_M = t_M;
    params.sigma = sigma;
    params.mu = mu;
    params.L_A_total = L_A_total;
    params.L_M_total = L_M_total;
    params.W_M = W_M;
    params.P_M = P_M;
    params.P_TM = P_TM;
    params.P_MA = P_MA;
    params.P_MM = P_MM;

    % Initial guesses for the variables
    x0 = log([2*ones(4, 1); % Consumptions
          2*ones(4, 1); % Labor allocations
          2]); % Wage


    % Solve the system of equations
    options = optimoptions('fsolve', ...
    'Display', 'iter', ...
    'TolFun', 1e-12, ...
    'TolX', 1e-12, ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 5000);
    x = fsolve(@(x) counterfactual(x, params), x0, options);
    % Store the solution
    base_solutions{i} = [x(1:9); W_M];
    
    % Calculate utilities
    if sigma == 1
        U_A = x(1) ^ mu * x(2) ^ (1 - mu);
        U_M = x(4) ^ mu * x(3) ^ (1 - mu);
    else
        U_A = (mu ^ (1 / sigma) * x(1) ^ ((sigma - 1) / sigma) + (1 - mu) ^ (1 / sigma) * x(2) ^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1));
        U_M = (mu ^ (1 / sigma) * x(4) ^ ((sigma - 1) / sigma) + (1 - mu) ^ (1 / sigma) * x(3) ^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1));
    end

    base_utilities(i, :) = [U_A, U_M];
end

% Label the solutions
label = ["C_AA";"C_MA";"C_AM";"C_MM";"L_A";"L_M";"L_TA";"L_TM";"W_A";"W_M"];
for i = 1:length(Z_TM_range)
    base_solutions{i} = [label, base_solutions{i}];
end


% Extract utilities for each country
base_U_A = base_utilities(:, 1)'; % Utilities in A
base_U_M = base_utilities(:, 2)'; % Utilities in Mexico


%% Factual World: Border
% Define the function for the system of equations
function F = factualborder(x, params)
    % Unpack variables
    C = exp(x(1:4)); % Consumptions
    % C(1): Consumption of goods from country A in country A
    % C(2): Consumption of goods from country M in country A
    % C(3): Consumption of goods from country A in country M
    % C(4): Consumption of goods from country M in country M

    L = exp(x(5:8)); % Labor allocations
    % L(1): Labor allocated to the good sector in country A
    % L(2): Labor allocated to the good sector in country M
    % L(3): Labor allocated to the transportation sector in country A
    % L(4): Labor allocated to the transportation sector in country M

    W = exp(x(9)); % Wage in country A

    % Unpack parameters
    Z_A = params.Z_A;
    Z_M = params.Z_M;
    Z_TA = params.Z_TA;
    Z_TM = params.Z_TM;
    t_A = params.t_A;
    t_M = params.t_M;
    sigma = params.sigma;
    mu = params.mu;
    L_A_total = params.L_A_total;
    L_M_total = params.L_M_total;
    W_M = params.W_M;
    P_M = params.P_M;
    P_TM = params.P_TM;

    % Delivered prices
    P_A = W / Z_A;
    P_TA = W / Z_TA;
    P_AA = P_A;
    P_AM = P_A + (t_A / 2) * P_TA + (t_A / 2) * P_TM;
    P_MA = P_M + (t_M / 2) * P_TM + (t_M / 2) * P_TA;
    P_MM = P_M;


    % Equations
    % Utility maximization
    F(1) = (P_AA + ((1 - mu) / mu) * P_AA ^ sigma * P_MA ^ (1 - sigma)) * C(1) - W * L_A_total;
    F(2) = (P_AA + ((1 - mu) / mu) * P_AA ^ sigma * P_MA ^ (1 - sigma)) * C(2) - ((1 - mu) / mu) * (P_AA / P_MA) ^ sigma * W * L_A_total;
    F(3) = (P_MM + ((1 - mu) / mu) * P_MM ^ sigma * P_AM ^ (1 - sigma)) * C(3) - ((1 - mu) / mu) * (P_MM / P_AM) ^ sigma * W_M * L_M_total;
    F(4) = (P_MM + ((1 - mu) / mu) * P_MM ^ sigma * P_AM ^ (1 - sigma)) * C(4) - W_M * L_M_total;

    % Good market clearing
    F(5) = Z_A * L(1) - (C(1) + C(3));
    F(6) = Z_M * L(2) - (C(2) + C(4));

    % Transportation market clearing
    F(7) = Z_TA * L(3) - (t_A / 2 * C(3) + t_M / 2 * C(2));
    F(8) = Z_TM * L(4) - (t_M / 2 * C(2) + t_A / 2 * C(3));

    % Labor market clearing
    F(9) = L_A_total - (L(1) + L(3));
    F(10) = L_M_total - (L(2) + L(4));

end




for i = 1:length(Z_TM_range)
    % Normalization
    W_M = 1;
    P_M = 1 / Z_M;
    P_TM = 1 / Z_TM_range(i);

    % Delivered prices

    % Set parameters
    params.Z_A = Z_A;
    params.Z_M = Z_M;
    params.Z_TA = Z_TA;
    params.Z_TM = Z_TM_range(i);
    params.t_A = t_A;
    params.t_M = t_M;
    params.sigma = sigma;
    params.mu = mu;
    params.L_A_total = L_A_total;
    params.L_M_total = L_M_total;
    params.W_M = W_M;
    params.P_M = P_M;
    params.P_TM = P_TM;

    % Initial guesses for the variables
    x0 = log([2*ones(4, 1); % Consumptions
          2*ones(4, 1); % Labor allocations
          2]); % Wage


    % Solve the system of equations
    options = optimoptions('fsolve', ...
    'Display', 'iter', ...
    'TolFun', 1e-12, ...
    'TolX', 1e-12, ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 5000);
    x = fsolve(@(x) factualborder(x, params), x0, options);
    % Store the solution
    border_solutions{i} = [exp(x(1:9)); W_M];
    
    % Calculate utilities
    if sigma == 1
        U_A = exp(x(1)) ^ mu * exp(x(2)) ^ (1 - mu);
        U_M = exp(x(4)) ^ mu * exp(x(3)) ^ (1 - mu);
    else
        U_A = (mu ^ (1 / sigma) * exp(x(1)) ^ ((sigma - 1) / sigma) + (1 - mu) ^ (1 / sigma) * exp(x(2)) ^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1));
        U_M = (mu ^ (1 / sigma) * exp(x(4)) ^ ((sigma - 1) / sigma) + (1 - mu) ^ (1 / sigma) * exp(x(3)) ^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1));
    end

    border_utilities(i, :) = [U_A, U_M];
end

% Label the solutions
label = ["C_AA";"C_MA";"C_AM";"C_MM";"L_A";"L_M";"L_TA";"L_TM";"W_A";"W_M"];
for i = 1:length(Z_TM_range)
    border_solutions{i} = [label, border_solutions{i}];
end


% Extract utilities for each country
border_U_A = border_utilities(:, 1)'; % Utilities in A
border_U_M = border_utilities(:, 2)'; % Utilities in Mexico


%% Plotting utiltities by worlds
% Create a tiled layout for two plots side by side
figure;

% Plotting the utilities (counterfactual)
nexttile; % Move to the first tile
hold on;
plot(Z_TM_range, base_U_A, 'b-', 'DisplayName', 'Utility in the US (Base)', 'LineWidth', 1.5);
plot(Z_TM_range, base_U_M, 'r-', 'DisplayName', 'Utility in Mexico (Base)', 'LineWidth', 1.5);
plot(Z_TM_range, border_U_A, 'b--', 'DisplayName', 'Utility in the US (Border)', 'LineWidth', 1.5);
plot(Z_TM_range, border_U_M, 'r--', 'DisplayName', 'Utility in Mexico (Border)', 'LineWidth', 1.5);

% Add labels and title
xlabel('Z_{TM} (Transportation Productivity in Mexico)', 'FontSize', 12);
ylabel('Utility', 'FontSize', 12);
title(sprintf('\\sigma = %.1f and \\mu = %.1f', sigma, mu), 'FontSize', 14);

% Add legend
legend('Location', 'best', 'FontSize', 10);

% Add grid
grid on;
hold off;

% Set the paper size for exporting
set(gcf, 'Position', [100, 100, 800, 600]);

% Save the figure
%saveas(gcf, 'sigma1_mu0.5.png');

