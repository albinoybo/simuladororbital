function E = resolver_kepler(M, e, tol)
    % Resolve a Equação de Kepler M = E - e*sin(E) usando o método de Newton-Raphson.
    % M: Anomalia Média (em radianos)
    % e: Excentricidade
    % tol: Tolerância para a convergência
    % E: Anomalia Excêntrica (em radianos)

    % Chute inicial
    if M < pi
        E = M + e/2;
    else
        E = M - e/2;
    end

    % Iteração de Newton-Raphson
    delta_E = 1;
    while abs(delta_E) > tol
        f_E = E - e * sin(E) - M; % Equação f(E) = 0
        df_dE = 1 - e * cos(E);   % Derivada f'(E)
        delta_E = f_E / df_dE;
        E = E - delta_E;
    end
end