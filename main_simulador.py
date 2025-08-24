import numpy as np
import matplotlib.pyplot as plt

# Bibliotecas de Astrodinâmica e Unidades
from astropy import units as u
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from tqdm import tqdm

# --- DADOS DOS SATÉLITES ---
SATELITE_SELECIONADO = "LEO 1"

DADOS_SATELITES = {
    "LEO 1": {
        "NOME": "LEO 1",
        "PERIGEU_KM": 423,
        "APOGEU_KM": 939,
        "ECC": 0.0365331,
        "INC_GRAUS": 99.2936,
        "RAAN_GRAUS": 237.4020,
        "AOP_GRAUS": 179.9420,
        "MA_GRAUS": 180.1886,
    },
    # Outros dados de satélite...
}

DADOS_SATELITE = DADOS_SATELITES[SATELITE_SELECIONADO]

# --- CONSTANTES FÍSICAS E DE SIMULAÇÃO ---
RAIO_TERRA_KM = Earth.R.to(u.km).value
MU_TERRA_KM3_S2 = Earth.k.to(u.km**3 / u.s**2).value

TEMPO_TOTAL_SIM_DIAS = 0.2
PASSO_TEMPO_SEGUNDOS = 30.0

# --- NÚCLEO DA SIMULAÇÃO ---
def simulador_kepleriano_python(dados_sat, t_span_dias, passo_s):
    """
    Executa a simulação orbital inteiramente em Python usando a biblioteca Poliastro.
    """
    print("Executando o núcleo da simulação em Python (modo iterativo)...")

    r_p = dados_sat["PERIGEU_KM"] * u.km + Earth.R
    r_a = dados_sat["APOGEU_KM"] * u.km + Earth.R
    a = (r_p + r_a) / 2
    
    ecc = dados_sat["ECC"] * u.one
    inc = dados_sat["INC_GRAUS"] * u.deg
    raan = dados_sat["RAAN_GRAUS"] * u.deg
    argp = dados_sat["AOP_GRAUS"] * u.deg
    M = dados_sat["MA_GRAUS"] * u.deg

    orbita_inicial = Orbit.from_classical(Earth, a, ecc, inc, raan, argp, M)

    tempo_total_s = t_span_dias * 24 * 3600
    tempos_propagacao = np.arange(0, tempo_total_s, passo_s) * u.s
    
    print("Propagando a órbita passo a passo...")
    
    num_pontos = len(tempos_propagacao)
    posicao_hist = np.zeros((3, num_pontos))
    velocidade_hist = np.zeros((3, num_pontos))
    altitude_hist = np.zeros(num_pontos)
    excentricidade_hist = np.zeros(num_pontos)
    perigeu_hist = np.zeros(num_pontos)

    for i, t_voo in enumerate(tqdm(tempos_propagacao, desc="Simulando")):
        orbita_no_tempo_t = orbita_inicial.propagate(t_voo)
        
        posicao_hist[:, i] = orbita_no_tempo_t.r.to(u.km).value
        velocidade_hist[:, i] = orbita_no_tempo_t.v.to(u.km / u.s).value
        altitude_hist[i] = (np.linalg.norm(posicao_hist[:, i]) - RAIO_TERRA_KM)
        excentricidade_hist[i] = orbita_no_tempo_t.ecc.value
        perigeu_hist[i] = (orbita_no_tempo_t.r_p - Earth.R).to(u.km).value

    print("Propagação da órbita concluída.")
    
    resultados = {
        'tempo': tempos_propagacao.to(u.s).value,
        'posicao': posicao_hist,
        'velocidade': velocidade_hist,
        'altitude': altitude_hist,
        'excentricidade': excentricidade_hist,
        'perigeu': perigeu_hist
    }
    
    return resultados

# --- FUNÇÕES AUXILIARES E DE PLOTAGEM ---

def exibir_dados_keplerianos():
    print("\n--- Parâmetros Keplerianos Iniciais ---")
    for chave, valor in DADOS_SATELITE.items():
        print(f"{chave:<15}: {valor}")
    print("-------------------------------------\n")

def plotar_orbita(x, y, z):
    """Cria os gráficos 2D da trajetória orbital com as novas customizações."""
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    fig.canvas.manager.set_window_title('Visualização da Órbita')
    fig.suptitle(f'Trajetória Orbital - {DADOS_SATELITE["NOME"]}', fontsize=16)

    limite_max = max(abs(np.concatenate([x, y, z]))) * 1.1
    limites = [-limite_max, limite_max]
    
    # MODIFICAÇÃO: Estilo da Terra como uma estrela amarela
    estilo_terra = {'marker':'*', 'color':'gold', 'markersize':15, 'zorder': 10}
    
    # Plot 1: Plano XY
    ax1.plot(x, y, color='red') # MODIFICAÇÃO: Cor da órbita
    ax1.plot(0, 0, **estilo_terra)
    ax1.set_title('Plano XY'); ax1.set_xlabel('X (km)'); ax1.set_ylabel('Y (km)')
    ax1.set_xlim(limites); ax1.set_ylim(limites)
    ax1.set_aspect('equal', adjustable='box'); ax1.grid(True, linestyle='--', alpha=0.6)
    
    # Plot 2: Plano XZ
    ax2.plot(x, z, color='red') # MODIFICAÇÃO: Cor da órbita
    ax2.plot(0, 0, **estilo_terra)
    ax2.set_title('Plano XZ'); ax2.set_xlabel('X (km)'); ax2.set_ylabel('Z (km)')
    ax2.set_xlim(limites); ax2.set_ylim(limites)
    ax2.set_aspect('equal', adjustable='box'); ax2.grid(True, linestyle='--', alpha=0.6)

    # Plot 3: Plano YZ
    ax3.plot(y, z, color='red') # MODIFICAÇÃO: Cor da órbita
    ax3.plot(0, 0, **estilo_terra)
    ax3.set_title('Plano YZ'); ax3.set_xlabel('Y (km)'); ax3.set_ylabel('Z (km)')
    ax3.set_xlim(limites); ax3.set_ylim(limites)
    ax3.set_aspect('equal', adjustable='box'); ax3.grid(True, linestyle='--', alpha=0.6)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

def plotar_dados_temporais(tempo, altitude, perigeu, excentricidade, velocidade_vec):
    velocidade_mag = np.linalg.norm(velocidade_vec, axis=0)
    tempo_horas = tempo / 3600.0
    
    fig, axs = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
    fig.canvas.manager.set_window_title('Análise de Parâmetros Orbitais')
    fig.suptitle(f'Parâmetros Orbitais vs. Tempo - {DADOS_SATELITE["NOME"]}', fontsize=16)
    
    axs[0].plot(tempo_horas, altitude)
    axs[0].set_ylabel('Altitude (km)')
    axs[0].grid(True)
    
    axs[1].plot(tempo_horas, perigeu)
    axs[1].set_ylabel('Altitude Perigeu (km)')
    axs[1].grid(True)
    
    axs[2].plot(tempo_horas, excentricidade)
    axs[2].set_ylabel('Excentricidade')
    axs[2].grid(True)
    
    axs[3].plot(tempo_horas, velocidade_mag)
    axs[3].set_ylabel('Velocidade (km/s)')
    axs[3].set_xlabel('Tempo (horas)')
    axs[3].grid(True)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

# --- FUNÇÃO PRINCIPAL ---
def executar_simulacao():
    print("Iniciando o simulador orbital em Python...")
    print(f"Executando simulação para: {DADOS_SATELITE['NOME']}...")
    
    resultados = simulador_kepleriano_python(DADOS_SATELITE, TEMPO_TOTAL_SIM_DIAS, PASSO_TEMPO_SEGUNDOS)
    
    tempo = resultados['tempo']
    posicao = resultados['posicao']
    velocidade = resultados['velocidade']
    altitude = resultados['altitude']
    excentricidade = resultados['excentricidade']
    perigeu = resultados['perigeu']
    
    exibir_dados_keplerianos()
    
    plotar_orbita(posicao[0,:], posicao[1,:], posicao[2,:])
    plotar_dados_temporais(tempo, altitude, perigeu, excentricidade, velocidade)

# --- PONTO DE ENTRADA DO PROGRAMA ---
if __name__ == "__main__":
    executar_simulacao()