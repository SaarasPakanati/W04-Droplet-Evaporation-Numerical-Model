% W04 Droplet Evaporation
% Developed By: Saaras Pakanati

% Initalization
clc; clear; close all;

%% Chemical Properties

% Acetone
% OSHA: https://www.osha.gov/chemicaldata/476
acetone_MW  = 58.08;                                    % [g/mol] Molecular Weight
acetone_rho = 0.784;                                    % [g/mL] Density at 20 C
acetone_PEL = 1000;                                     % [ppm] Permissible Exposure Limit - 8-hour TWA (ppm)
acetone_LEL = 25000;                                    % [ppm] Lower explosive limit (LEL) means the minimum concentration of vapor in air below which propagation of a flame does not occur in the presence of an ignition source.    
acetone_UEL = 128000;                                   % [ppm] Upper explosive limit (UEL) means the maximum concentration of flammable vapor in air above which propagation of flame does not occur on contact with a source of ignition.   
acetone_name = 'Acetone';

% Isopropyl Alcohol 
% OSHA: https://www.osha.gov/chemicaldata/475
IPA_MW  = 60.10;                                        % [g/mol] Molecular Weight
IPA_rho = 0.786;                                        % [g/mL] Density at 20 C
IPA_PEL = 400;                                          % [ppm] Permissible Exposure Limit - 8-hour TWA (ppm)
IPA_LEL = 20000;                                        % [ppm] Lower explosive limit (LEL) means the minimum concentration of vapor in air below which propagation of a flame does not occur in the presence of an ignition source.    
IPA_UEL = 127000;                                       % [ppm] Upper explosive limit (UEL) means the maximum concentration of flammable vapor in air above which propagation of flame does not occur on contact with a source of ignition.   
IPA_name = 'IPA';

% n-Pentane 
% OSHA: https://www.osha.gov/chemicaldata/98
pentane_MW  = 72.15;                                    % [g/mol] Molecular Weight
pentane_rho = 0.626;                                    % [g/mL] Density at 20 C
pentane_PEL = 1000;                                     % [ppm] Permissible Exposure Limit - 8-hour TWA (ppm)    
pentane_LEL = 15000;                                    % [ppm] Lower explosive limit (LEL) means the minimum concentration of vapor in air below which propagation of a flame does not occur in the presence of an ignition source.    
pentane_UEL = 78000;                                    % [ppm] Upper explosive limit (UEL) means the maximum concentration of flammable vapor in air above which propagation of flame does not occur on contact with a source of ignition.   
pentane_name = 'Pentane';

%% Calculation and Plotting

close all;

% Evaporated Liquid Volumes (Based on syringe volumes)
V_liq_mL = [1, 3, 5, 7, 10];                            % [mL] Liquid volumes to simulate

% Initializing Chemical Properties
chem_MW   = [acetone_MW,   IPA_MW,   pentane_MW];       % [g/mol] Molecular Weight
chem_rho  = [acetone_rho,  IPA_rho,  pentane_rho];      % [g/mL] Density
chem_name = {acetone_name, IPA_name, pentane_name};     % [str] Chemical Names

% Choose Limit
chem_limit  = [acetone_PEL,  IPA_PEL,  pentane_PEL];    % [ppm] Limit
% chem_limit  = [acetone_LEL,  IPA_LEL,  pentane_LEL];    % [ppm] Limit
% chem_limit  = [acetone_UEL,  IPA_UEL,  pentane_UEL];    % [ppm] Limit

% Plot parameters for each chemical
maxVolCont = [5000, 5000, 5000];                        % [mL] Max container volume for plot
% maxVolCont = [5, 5, 5];                        % [mL] Max container volume for plot
yMaxPlot   = [25000, 25000, 25000];                        % [ppm] Max y-axis for plot

% Plot Colors
colors = [0.00 0.45 0.74;
          0.85 0.33 0.10;
          0.93 0.69 0.01;
          0.49 0.18 0.56;
          0.00 0.67 0.30];

% Plot Loop
for chem_idx = 1:3

    % Extract Chemical Properties
    MW        = chem_MW(chem_idx);                      % [g/mol] Molecular Weight
    rho       = chem_rho(chem_idx);                     % [g/mL] Density
    PEL       = chem_limit(chem_idx);                   % [ppm] Exposure Limit
    name      = chem_name{chem_idx};                    % [str] Chemical Names
    v_con_max = maxVolCont(chem_idx);                   % [mL] Max container volume
    y_max     = yMaxPlot(chem_idx);                     % [ppm] Max y-axis
    
    % Calculating the Concentration
    for i = 1:length(V_liq_mL)
        v_liq = V_liq_mL(i);                            % [mL] Liquid Volume
        
        % Moles of Liquid
        n = v_liq * rho / MW;                           % [mol]

        % Vapor Volume at STP (NIST: 1 mol = 24.465 L)
        V_vap = n * 24.465;                             % [L]
        
        % Initialize container Volume Array
        V_con = linspace(v_liq/1e3 + 0.05, v_con_max, 2500);        % [L]
        
        % Volume remaining after removing liquid volume (In hindsight, an overkill...)
        V_rem = V_con - v_liq/1e3;                      % [L]
        
        % Vapor Concentration
        C = (V_vap ./ V_rem) * 1e6;                   % [ppm] Vapor concentration
        
        % Saving Data
        C_data{i} = C;
        V_con_data{i} = V_con;
    end
    
    % Plotting
    figure();
    
    % Safe zone
    fill([0 v_con_max v_con_max 0], [0 0 PEL PEL], [0.85 0.95 0.85], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    
    % NOT Safe Zone
    fill([0 v_con_max v_con_max 0], [PEL PEL y_max y_max], [0.95 0.85 0.85], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    hold on;
    
    % Threshold Line
    yline(PEL, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    
    % Data Plot
    for i = 1:length(V_liq_mL)
        C_plot = C_data{i};
        C_plot(C_plot > y_max) = NaN;
        plot(V_con_data{i}, C_plot, 'Color', colors(i,:), 'LineWidth', 2, ...
            'DisplayName', sprintf('%d mL', V_liq_mL(i)));
    end
    
    % Labels
    text(v_con_max*0.75, PEL*0.35, '$\mathrm{Safe}$', ...
        'Interpreter','latex', 'FontName','Serif', 'FontSize', 14, ...
        'Color', [0.0 0.5 0.0], 'FontWeight', 'bold');
    text(v_con_max*0.25, y_max*0.88, '$\mathrm{Exceeds \ Limit}$', ...
        'Interpreter','latex', 'FontName','Serif', 'FontSize', 14, ...
        'Color', [0.7 0.0 0.0], 'FontWeight', 'bold');
    hold off;
    ylabel("$\textbf{Vapor Concentration } \mathbf{(ppm)}$", "Interpreter","latex", 'FontName', 'Serif');
    xlabel("$\textbf{Container Volume } \mathbf{(L)}$",  "Interpreter","latex", 'FontName', 'Serif');
    title(sprintf("%s", name),"Interpreter","latex", 'FontName', 'Serif');
    set(gca, 'FontWeight','bold', 'fontsize', 14, 'FontName', 'Serif');
    legend('Location', 'best');
    axis padded;
    grid minor;
    ylim([0 y_max]);
    xlim([0 v_con_max]);

end
