%% getConductance finds pipe conductance in the molecular flow regime.
% This formula assumes the following units:
% m    kg
% T    ^oK
% L    cm
% I    cm^3
%
% Conductance is returned with units:
% C    Liters/sec
%
% use 4.65e-26 kg for molecular Nitrogen for example.
% use 293 K for 20 C room temperature.
%
function C = getConductance(m,T,L,I)
    v = 100 * sqrt(1.3801e-23 * T / ( 2*pi*m ) ); % 100 converts to cm/s
    C = 1/2 * I/L * v / 1000; % divide by 1000 to co from cm^3 to L
end