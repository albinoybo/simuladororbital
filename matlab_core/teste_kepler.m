% Script para validar a função resolver_kepler

constantes; % Carrega mu e R_terra

% Dados do satélite AQUA 
e_aqua = 0.0002177;
M_aqua_graus = 290.1139;
M_aqua_rad = deg2rad(M_aqua_graus); % Converter para radianos

% Resolver a equação
tolerancia = 1e-8;
E_aqua_rad = resolver_kepler(M_aqua_rad, e_aqua, tolerancia);
E_aqua_graus = rad2deg(E_aqua_rad);

% Exibir resultado no terminal
fprintf('--- Validação do Solucionador de Kepler (AQUA) ---\n');
fprintf('Anomalia Média (M) = %.4f graus\n', M_aqua_graus);
fprintf('Anomalia Excêntrica (E) = %.4f graus\n\n', E_aqua_graus);