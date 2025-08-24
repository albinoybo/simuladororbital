% --- Constantes Físicas Fundamentais ---
% Este script define as constantes que serão usadas em todo o simulador.

clear; clc; % Limpa o ambiente para garantir consistência

global mu R_terra; % Define como variáveis globais para fácil acesso

% Constante Gravitacional Geocêntrica (km^3/s^2)
mu = 3.986004418e5; % 

% Raio Equatorial da Terra (km)
R_terra = 6378.137; % 

fprintf('Constantes físicas carregadas no ambiente.\n');
fprintf('  mu = %.4e km^3/s^2\n', mu);
fprintf('  R_terra = %.1f km\n\n', R_terra);