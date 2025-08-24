function elementos = kepler_para_vetor_estado(r_vec, v_vec)
    % Converte vetores de estado (r, v) para elementos Keplerianos clássicos.
    % Baseado no Capítulo 6 de "Introdução à Mecânica Orbital".
    
    global mu R_terra;
    
    r = norm(r_vec);
    v = norm(v_vec);
    
    % Vetor momento angular específico
    h_vec = cross(r_vec, v_vec);
    h = norm(h_vec);
    
    % Inclinação (i)
    i = acos(h_vec(3) / h);
    
    % Vetor do nodo ascendente
    n_vec = cross([0; 0; 1], h_vec);
    n = norm(n_vec);
    
    % Ascensão reta do nodo ascendente (RAAN)
    raan = acos(n_vec(1) / n);
    if n_vec(2) < 0
        raan = 2*pi - raan;
    end
    
    % Vetor de excentricidade
    e_vec = (1/mu) * (cross(v_vec, h_vec) - mu * r_vec / r);
    e = norm(e_vec);
    
    % Argumento do perigeu (aop)
    aop = acos(dot(n_vec, e_vec) / (n * e));
    if e_vec(3) < 0
        aop = 2*pi - aop;
    end
    
    % Anomalia verdadeira (nu)
    nu = acos(dot(e_vec, r_vec) / (e * r));
    if dot(r_vec, v_vec) < 0
        nu = 2*pi - nu;
    end
    
    % Semi-eixo maior (a) e Perigeu/Apogeu
    a = 1 / (2/r - v^2/mu);
    altitude_perigeu = a * (1 - e) - R_terra;
    
    elementos.a = a;
    elementos.e = e;
    elementos.i = i;
    elementos.raan = raan;
    elementos.aop = aop;
    elementos.nu = nu;
    elementos.altitude_perigeu = altitude_perigeu;
end