function F = CES(x_log, params)
    x = exp(x_log);  % Log-to-level transformation
    % Unpack variables
    % C: [C_AA, C_BA, C_MA, C_AB, C_BB, C_MB, C_AM, C_BM, C_MM]
    %    C_AA: Consumption in Texas of goods from Texas
    %    C_BA: Consumption in Texas of goods from California
    %    C_MA: Consumption in Texas of goods from Mexico
    %    C_AB: Consumption in California of goods from Texas
    %    C_BB: Consumption in California of goods from California
    %    C_MB: Consumption in California of goods from Mexico
    %    C_AM: Consumption in Mexico of goods from Texas
    %    C_BM: Consumption in Mexico of goods from California
    %    C_MM: Consumption in Mexico of goods from Mexico
    C = x(1:9); % Consumptions
    % L: [L_A, L_B, L_M]
    %    L_A: Labor in goods sector in Texas
    %    L_B: Labor in goods sector in California
    %    L_M: Labor in goods sector in Mexico
    L = x(10:12); % Labor allocations in good sector
    % LT: [LT_A, LT_B, LT_M]
    %    LT_A: Labor in transportation sector in Texas
    %    LT_B: Labor in transportation sector in California
    %    LT_M: Labor in transportation sector in Mexico
    LT = x(13:15); % Labor allocations in transportation sector
    % T: [T_A_AB, T_B_AB, T_M_AB, T_A_AM, T_B_AM, T_M_AM, T_A_BM, T_B_BM, T_M_BM, T_A_BA, T_B_BA, T_M_BA, T_A_MA, T_B_MA, T_M_MA, T_A_MB, T_B_MB, T_M_MB]
    %    T_A_AB: Transport services on AB route provided by Texas
    %    T_B_AB: ... by California
    %    T_M_AB: ... by Mexico
    %    T_A_AM: Transport services on AM route provided by Texas
    %    T_B_AM: ... by California
    %    T_M_AM: ... by Mexico
    %    T_A_BM: Transport services on BM route provided by Texas
    %    T_B_BM: ... by California
    %    T_M_BM: ... by Mexico
    %    T_A_BA: Transport services on BA route provided by Texas
    %    T_B_BA: ... by California
    %    T_M_BA: ... by Mexico
    %    T_A_MA: Transport services on MA route provided by Texas
    %    T_B_MA: ... by California
    %    T_M_MA: ... by Mexico
    %    T_A_MB: Transport services on MB route provided by Texas
    %    T_B_MB: ... by California
    %    T_M_MB: ... by Mexico
    T = x(16:33); % Demand for transportation services
    W_A = x(34); % Wage in Texas
    W_B = x(35); % Wage in California

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
    % F(35) = L_M_total - (L(3) + LT(3)); % This equation is redundant because W_M = 1

    % Spatial equilibrium
    if sigma == 1
        F(36) = C(1) ^ mu * C(2) ^ ((1 - mu) / 2) * C(3) ^ ((1 - mu) / 2) - C(4) ^ ((1 - mu) / 2) * C(5) ^ mu * C(6) ^ ((1 - mu) / 2) ;
    else
        F(36) = (mu ^ (1 / sigma) * C(1) ^ ((sigma - 1) / sigma) + ((1 - mu) / 2) ^ (1 / sigma) * C(2) ^ ((sigma - 1) / sigma) + ((1 - mu) / 2) ^ (1 / sigma) * C(3) ^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1)) - ...
            ((((1 - mu) / 2) ^ (1 / sigma) * C(4) ^ ((sigma - 1) / sigma) + mu ^ (1 / sigma) * C(5) ^ ((sigma - 1) / sigma) + ((1 - mu) / 2) ^ (1 / sigma) * C(6) ^ ((sigma - 1) / sigma)) ^ (sigma / (sigma - 1)));
    end    
end
