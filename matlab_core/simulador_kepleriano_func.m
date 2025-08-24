function resultados = simulador_kepleriano_func(tempo_total_sim, passo_tempo)
    % --- Simulador Orbital Kepleriano Básico (versão função) ---
    global mu R_terra;
    definir_constantes; % Garante que as globais estão carregadas

    % --- Dados do Satélite (LEO 1) ---
    % ... (mantém os mesmos dados do satélite) ...
    per = 423; apo = 939; a = (per + apo) / 2 + R_terra; e = 0.0365331;
    i = deg2rad(99.2936); raan0 = deg2rad(237.4020); aop0 = deg2rad(179.9420);
    M0 = deg2rad(180.1886); n_rad_s = 14.636810 * 2 * pi / 86400;

    % --- Loop de Simulação ---
    tempos = 0:passo_tempo:tempo_total_sim;
    num_passos = length(tempos);
    
    % Inicializa as estruturas de armazenamento
    posicoes = zeros(3, num_passos);
    velocidades = zeros(3, num_passos);
    altitudes = zeros(1, num_passos);
    excentricidades = zeros(1, num_passos);
    perigeus = zeros(1, num_passos);

    for idx = 1:num_passos
        t = tempos(idx);
        M_atual = mod(M0 + n_rad_s * t, 2*pi);
        
        [r_vec, v_vec] = kepler_para_vetor_estado(a, e, i, raan0, aop0, M_atual);
        posicoes(:, idx) = r_vec;
        velocidades(:, idx) = v_vec;
        
        % ALTERAÇÃO: Calcula os elementos orbitais a cada passo
        elementos_atuais = kepler_para_vetor_estado(r_vec, v_vec);
        
        % Armazena os dados para plotagem
        altitudes(idx) = norm(r_vec) - R_terra;
        excentricidades(idx) = elementos_atuais.e;
        perigeus(idx) = elementos_atuais.altitude_perigeu;
    end
    
    % Retorna uma estrutura com todos os dados
    resultados.tempo = tempos;
    resultados.posicao = posicoes;
    resultados.velocidade = velocidades;
    resultados.altitude = altitudes;
    resultados.excentricidade = excentricidades;
    resultados.perigeu = perigeus;
end