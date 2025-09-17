# -*- coding: utf-8 -*-
# Este script Python plota uma matriz de dados de um arquivo de texto.

import numpy as np
import matplotlib.pyplot as plt

def plot_ising_matrix(filename, title):
    """
    Lê uma matriz de spins de um arquivo de texto e a plota.
    
    Args:
        filename (str): O nome do arquivo de dados (ex: '5.000.dat').
        title (str): O título do gráfico.
    """
    try:
        # Lê a matriz de spins do arquivo de texto.
        # np.loadtxt é muito eficiente para este formato de dados.
        matrix = np.loadtxt(filename)

        # Configura o gráfico
        plt.figure(figsize=(8, 8))
        
        # Usa imshow para criar um heatmap.
        # 'interpolation=nearest' garante que os pixels sejam nítidos e quadrados.
        # 'cmap' define a paleta de cores. 'binary' é uma ótima opção.
        plt.imshow(matrix, cmap='binary', interpolation='nearest')
        
        # Adiciona o título e remove os eixos para uma visualização limpa
        plt.title(title, fontsize=16)
        plt.xticks([])
        plt.yticks([])
        
        # Salva o gráfico em um arquivo PNG
        plt.savefig("{}.png".format(filename.replace('.dat', '')))
        print("Plotagem de '{}' concluída. Imagem salva como '{}'".format(filename, filename.replace('.dat', '.png')))

    except FileNotFoundError:
        print("Erro: O arquivo '{}' não foi encontrado.".format(filename))
    except Exception as e:
        print("Ocorreu um erro ao plotar a matriz: {}".format(e))

if __name__ == "__main__":
    # Range de temperaturas a serem plotadas, de 5.0 a 0.5 com passo de 0.1
    # Usamos inteiros para evitar problemas de ponto flutuante no loop
    temperatures = [t / 10.0 for t in range(50, 4, -1)]

    # Loop para plotar cada arquivo de temperatura
    for temp in temperatures:
        # Formata o nome do arquivo com 3 casas decimais
        filename = "{:.3f}.dat".format(temp)
        
        # Define um título dinâmico para cada gráfico
        title = "Matriz de Spins para T = {}".format(temp)
        
        # Chama a função para plotar
        plot_ising_matrix(filename, title)

