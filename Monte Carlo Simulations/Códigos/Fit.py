import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def curve_function(x, a, b, c):
    return a * x**2 + b * x + c

def fit_and_plot_data(file_path):
    
    try:
        df = pd.read_csv(file_path)
        
        
        x_data = df.iloc[:, 0].values
        y_data = df.iloc[:, 5].values
        
        finite_indices = np.isfinite(x_data) & np.isfinite(y_data)
        x_data = x_data[finite_indices]
        y_data = y_data[finite_indices]
        
      
        if len(x_data) < 3:
            print("Dados insuficientes para o ajuste de curva.")
            return

        
        popt, pcov = curve_fit(curve_function, x_data, y_data)
        
        
        x_fit = np.linspace(min(x_data), max(x_data), 100)
        y_fit = curve_function(x_fit, *popt)
        

        plt.figure(figsize=(10, 6))
        plt.plot(x_data, y_data, 'o', label='Dados Originais')
        plt.plot(x_fit, y_fit, '-', label='Curva de Ajuste')
        
        
        equation = "y = {:.2f}x**2 + {:.2f}x + {:.2f}".format(popt[0], popt[1], popt[2])
        plt.text(0.05, 0.95, equation, transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')

        plt.title('Ajuste de Curva a Dados Esparsos')
        plt.xlabel('Eixo X')
        plt.ylabel('Eixo Y')
        plt.legend()
        plt.grid(True)
        plt.show()

        print("Parametros da curva ajustada:")


fit_and_plot_data('fort.8')
# Exemplo de como usar a funcao
# Substitua 'seus_dados.csv' pelo nome do seu arquivo

