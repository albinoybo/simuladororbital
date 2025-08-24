% Script para validar a função kepler_estado

constantes; % Carrega mu e R_terra

% Dados do satélite TERRA 
per_terra = 690;
apo_terra = 694;
a_terra = (per_terra + apo_terra) / 2 + R_terra; % Calcula o semi-eixo maior

e_terra = 0.0002861;
i_terra_rad = deg2rad(98.0190);
raan_terra_rad = deg2rad(92.1562);
aop_terra_rad = deg2rad(41.5034);
M_terra_rad = deg2rad(96.0644);

% Chamar a função de conversão
[r, v] = kepler_estado(a_terra, e_terra, i_terra_rad, raan_terra_rad, aop_terra_rad, M_terra_rad);

% Exibir os vetores no terminal
fprintf('--- Validação da Conversão Kepler -> Vetor de Estado (TERRA) ---\n');
fprintf('Vetor Posição (r) [km]:\n');
fprintf('  X: %.4f\n  Y: %.4f\n  Z: %.4f\n', r(1), r(2), r(3));
fprintf('Magnitude |r|: %.4f km\n\n', norm(r));

fprintf('Vetor Velocidade (v) [km/s]:\n');
fprintf('  VX: %.4f\n  VY: %.4f\n  VZ: %.4f\n', v(1), v(2), v(3));
fprintf('Magnitude |v|: %.4f km/s\n\n', norm(v));