% Droplet Numerical Model with Latent Heat
% W04 Droplet Evaporation
%
% Update Log:
%   Saaras Pakanati - 11/30/2025


% Initialization
clc; clear; close all;

%% Experimental Data

% Design of Experiment data for experiments conducted in November. 
%   Using this data, we can tune the model to get an accurate qualitative
%   representation of the actual experiments.

exp_data_raw = -1.*[
    0.43866, 0.56136, 0.66188;
    0.44269, 0.55107, 0.67747;
    0.33779, 0.53126, 0.72715;
    0.50985, 0.93962, 1.10661;
    0.53852, 0.85697, 1.12887;
    0.52839, 0.96345, 1.14877;
    0.93719, 1.79751, 2.33329;
    0.79754, 1.83303, 2.64384;
    0.80075, 2.10863, 2.43470
];

expTemp = [80, 90, 100];
expVolume = [100, 200, 300];

% Computing the exerimental means
expMeans = zeros(3,3);
for i = 1:3
    for j = 1:3
        row_start = (i-1)*3 + 1;
        row_end = i*3;
        expMeans(i,j) = mean(exp_data_raw(row_start:row_end, j));
    end
end

% Box-Cox transformation (lambda = 0.21)
lambda = 0.21;
exp_transformed = expMeans .^ lambda;

fprintf('===============================================\n');
fprintf('EXPERIMENTAL DATA SUMMARY\n');
fprintf('===============================================\n\n');
fprintf('Raw Evaporation Rates [mg/s]:\n');
fprintf('Temp\\Vol   100µL    200µL    300µL\n');
for i = 1:3
    fprintf('%3d°C   %.4f  %.4f  %.4f\n', expTemp(i), expMeans(i,:));
end
fprintf('\nBox-Cox Transformed (lambda=%.2f):\n', lambda);
fprintf('Temp\\Vol   100µL    200µL    300µL\n');
for i = 1:3
    fprintf('%3d°C   %.4f  %.4f  %.4f\n', expTemp(i), exp_transformed(i,:));
end

%% Numerical Model Parameters

% Physical Parameters [Do not Tune this pwease! - Saaras]
h_conv = 1e+2;                                  % [W/m^2-K] Convection Coefficient
alpha  = 1.1e-9;                                % [-] Accommodation Coefficient
beta   = 0.975;                                 % [-] Correction Factor
theta_deg = 10;                                 % [deg] Contact Angle

% Physical Parameters (It is what it is)
k      = 0.6515;                                % [W/m-K] Thermal Conductivity
rho_l  = 983.2;                                 % [kg/m^3] Liquid Density
sigma  = 0.06624;                               % [N/m] Surface Tension
hfg    = 2358e3;                                % [J/kg] Latent Heat of Vaporization

%% Run the Numerical Model

model_results = zeros(3,3);                     % [mg/s]
convergence_status = true(3,3);

for i = 1:3
    for j = 1:3
        [m_total, Ti_avg, converged] = run_droplet_model( expVolume(j), expTemp(i), h_conv, alpha, beta, theta_deg, k, rho_l, sigma, hfg);
        model_results(i,j) = m_total * 1e6;      % [kg/s] Evaporation Rate in SI units
        convergence_status(i,j) = converged;
        
        if ~converged
            fprintf('The Ti and M iterations need to be increased!');
        end
        fprintf('\n');
    end
end

%% Plots

figure('Position', [100 100 1400 800]);

% Plot 1: Raw comparison
subplot(2,3,1);
hold on;
for i = 1:2
    plot(expVolume, expMeans(i,:), 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
         'DisplayName', sprintf('Exp %d°C', expTemp(i)));
end
for i = 1:2
    plot(expVolume, model_results(i,:), 's--', 'LineWidth', 1.5, 'MarkerSize', 6, ...
         'DisplayName', sprintf('Model %d°C', expTemp(i)));
end
xlabel('Droplet Volume [µL]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Evaporation Rate [mg/s]', 'FontSize', 12, 'FontWeight', 'bold');
title('Model vs Experiment (Raw)', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 9);
grid on; box on;

% Plot 2: Parity plot (raw)
subplot(2,3,2);
scatter(expMeans(1:6), model_results(1:6), 100, 'filled'); hold on;
plot([0 max(expMeans(1:6))], [0 max(expMeans(1:6))], 'k--', 'LineWidth', 2);
xlabel('Experimental [mg/s]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Model [mg/s]', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Parity Plot (R=%.3f)', corr_raw), 'FontSize', 14);
grid on; box on; axis equal;
% xlim([0 max(expMeans(:))*1.1]);
% ylim([0 max(expMeans(:))*1.1]);

% Plot 3: Transformed parity plot
subplot(2,3,3);
scatter(exp_transformed(:), model_transformed(:), 100, 'filled'); hold on;
plot([min(exp_transformed(:)) max(exp_transformed(:))], ...
     [min(exp_transformed(:)) max(exp_transformed(:))], 'k--', 'LineWidth', 2);
xlabel('Experimental (transformed)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Model (transformed)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Parity (Box-Cox, R=%.3f)', corr_transformed), 'FontSize', 14);
grid on; box on; axis equal;

% Plot 4: Relative errors heatmap
subplot(2,3,4);
imagesc(rel_errors');
colorbar;
colormap(flipud(hot));
set(gca, 'XTick', 1:3, 'XTickLabel', expTemp, ...
         'YTick', 1:3, 'YTickLabel', expVolume);
xlabel('Temperature [°C]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Volume [µL]', 'FontSize', 12, 'FontWeight', 'bold');
title('Relative Error [%]', 'FontSize', 14);
for ix = 1:3
    for jy = 1:3
        text(ix, jy, sprintf('%.1f%%', rel_errors(ix,jy)), ...
             'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
    end
end
box on;

% Plot 5: Temperature trends
subplot(2,3,5);
hold on;
for j = 1:3
    plot(expTemp, expMeans(:,j), 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
         'DisplayName', sprintf('Exp %dµL', expVolume(j)));
    plot(expTemp, model_results(:,j), 's--', 'LineWidth', 1.5, 'MarkerSize', 6, ...
         'DisplayName', sprintf('Model %dµL', expVolume(j)));
end
xlabel('Temperature [°C]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Evaporation Rate [mg/s]', 'FontSize', 12, 'FontWeight', 'bold');
title('Temperature Trend Comparison', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 9);
grid on; box on;

% Plot 6: Normalized trends
subplot(2,3,6);
hold on;
scatter(exp_norm(:), model_norm(:), 100, 'filled');
plot([0 1], [0 1], 'k--', 'LineWidth', 2);
xlabel('Experimental (normalized)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Model (normalized)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Trend Match (RMSE=%.3f)', trend_rmse), 'FontSize', 14);
grid on; box on; axis equal;
xlim([0 1]); ylim([0 1]);

sgtitle('Initial Model Performance Analysis', 'FontSize', 16, 'FontWeight', 'bold');

%% functions

% Calculate the apex height based on the volume we get.
function [h_apex, theta_c] = volume_to_geometry(V_uL, theta_deg)
    
    theta_c = theta_deg * pi / 180;                                 % [rad]
    V = V_uL * 1e-9;                                                % [m^3]
    
    factor = (2 + cos(theta_c)) / (1 - cos(theta_c));
    h_apex = (3 * V / (pi * factor))^(1/3);                         % [m]
end

% Da actual numerical model mate
function [m_total, Ti_avg, converged] = run_droplet_model( V_uL, T_s, h_conv, alpha, beta, theta_deg, k, rho_l, sigma, hfg)

    % Universal constants
    R_g = 8.314462618153240;                                        % [J/K-mol]
    M   = 18.015e-3;                                                % [kg/mol]
    T_inf = 25;                                                     % [C]
    T_inf_K = T_inf + 273.15;                                       % [K]
    T_s_K = T_s + 273.15;                                           % [K]
    
    % Geometry from volume
    [h_apex, theta_c] = volume_to_geometry(V_uL, theta_deg);
    R = h_apex / (1 - cos(theta_c));                                % [m] Radius of Curvature
    y_c = -R + h_apex;                                              % [m] Circle Center Correction
    a = R * sin(theta_c);                                           % [m] Contact Line radius
    
    % Discretization
    dx = 1e-7;                                                      % [m]
    x = -a : dx : a;                                                % [m]
    h = y_c + sqrt(R^2 - x.^2);                                     % [m] Droplet Height Profile
    h(h < 0) = 0;
    
    % Saturation pressure (Antoine)
    A = 8.07131; B = 1730.63; C = 233.426;
    pv_sat = 10^(A - B/(C + T_s)) * 133.322;                        % [Pa]
    pv_inf = 10^(A - B/(C + T_inf)) * 133.322;                      % [Pa]
    pv = pv_inf;                                                    % [Pa] Ambient Vapor Pressure
    
    % Vapor properties
    Tv = T_s_K;                                                     % [K]
    rho_v = (pv * M) / (R_g * Tv);                                  % [kg/m^3]
    
    % Initial guess for interface temperature (no evaporation)
    Ti = T_s_K - (h_conv .* h .* (T_s_K - T_inf_K)) ./ (k + (h_conv .* h));
    Ti_old = Ti;
    
    % Constant curvature for spherical cap
    kappa = 1 / R;
    
    % Iterative solver parameters
    max_iter = 1E+6;
    tol = 1e-6;
    relax = 0.5;
    converged = false;
    
    for iter = 1:max_iter
        % Evaporation model (Bellur-like formulation)
        W1 = pv_sat ./ pv;
        W2 = (1 - (Tv ./ Ti)) .* ((pv_sat .* hfg) ./ pv);
        W3 = (Tv ./ Ti) .* (rho_v ./ rho_l) .* (sigma .* kappa ./ pv);
        W = W1 + W2 + W3;
        
        bellur_1 = (2 .* alpha) ./ (2 - alpha);
        bellur_2 = sqrt(M ./ (2 .* pi .* R_g .* Tv));
        bellur_3 = pv .* beta .* W .* sqrt(Tv ./ Ti) - 1;
        m = bellur_1 .* bellur_2 .* bellur_3;                       % Local mass flux [kg/m^2-s]
        
        % Energy balance at interface: conduction = convection + evaporation
        Ti_new = (k .* T_s_K ./ h + h_conv .* T_inf_K - m .* hfg) ./ (k ./ h + h_conv);
        
        % Relaxation for stability
        Ti = relax .* Ti_new + (1 - relax) .* Ti_old;
        
        max_change = max(abs(Ti - Ti_old));
        if max_change < tol
            converged = true;
            break;
        end
        
        Ti_old = Ti;
    end
    
    % Integrate total evaporation rate over surface (axisymmetric)
    m_total = 0;  % [kg/s]
    for idx = 2 : length(x)
        dh_dx = (h(idx) - h(idx-1)) / dx;                           % Slope
        dS_factor = sqrt(1 + dh_dx^2);                              % Area Correction
        r = abs(x(idx));                                            % Radial Coordinate
        m_total = m_total + 2 * pi * r * m(idx) * dS_factor * dx;
    end
    
    Ti_avg = mean(Ti);  % average interface temperature [K]
end
