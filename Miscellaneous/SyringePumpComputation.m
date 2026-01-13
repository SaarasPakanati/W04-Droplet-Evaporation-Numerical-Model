% Syringe Pump Specifications Calculator
% W04 Droplet Evaporation
% Date: Jan 13 2026

% Initializing
clc; clear; close all;

% Inputs
thrRod = "M6";                                  % [-] Standard Metric Size
syrID = 8.66;                                   % [mm] Syringe's Inner Diameter
neeID = 1.3;                                    % [mm] Needle's Inner Diameter
neeLen = 16;                                    % [mm] Needle's Length
tubID = 4.0;                                    % [mm] Tube's Inner Diameter
tubLen = 480;                                   % [mm] Tube's Length
sysFre = 0.5;                                   % [1/s] Software Frequency
microSteps = 16;                                % [#] Sub steps in a  motor's complete step

% Choice Evaluation
if thrRod == "M4"
    thrRodPitch = 0.7;
elseif thrRod == "M5"
    thrRodPitch = 0.8;
elseif thrRod == "M6" || thrRod == "M7"
    thrRodPitch = 1.0;
else
    error("Incorrect threaded rod choice, mate");
end

% Calculations
volStep = (pi*(syrID^2)/4) * thrRodPitch;       % [mm^3] Volume displaced by one pitch length movement
volSignal = volStep / microSteps;               % [mm^3] Volume displaced by one signal
Q = volSignal * sysFre;                         % [mm^3/s] Volumetric Flow Rate

velTube = Q / (pi*(tubID^2)/4);                 % [mm/s] Velocity in the tube's orifce
velNeedle = Q / (pi*(neeID^2)/4);               % [mm/s] Velocity at the needle orfice

% Results
disp("----------------------------------------------------------------------------")
disp("Inputs:")
fprintf("The Inner Diameter of the Syringe        = %.2f mm \n", syrID);
fprintf("The Inner Diameter of the Needle         = %.2f mm \n", neeID);
fprintf("The Inner Diameter of the Tubing         = %.2f mm \n", tubID);
fprintf("The Pitch Length of the Treaded Rod      = %.2f mm \n", thrRodPitch);
disp("----------------------------------------------------------------------------")
disp("Results:")
fprintf("The volumetric flow rate                 = %.2f mm^3/s \n", Q);
fprintf("The velocity of the fluid in the tube    = %.2f mm/s \n", velTube);
fprintf("The velocity of the fluid in the needle  = %.2f mm/s \n", velNeedle);
disp("----------------------------------------------------------------------------")
