#Author: Jailson39
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
#import sympy as sym
import scipy as sp
from pylab import legend, grid
import pylab as pl

###### DEFININDO O SETUP: ######

m = 1                                               # Massa do Corpo(Elemento de Inércia).
k = 8                                               # Constante Elástica da mola(Meio Elástico)
xt0 = 10                                            # Posição Inicial do bloco  ----> x(t0) = xt0 = 10cm
vt0 = 0                                             # Velocidade Inicial do bloco  ----> V(t0) = X'(t0) = X0' = V0 e [V0] = [1cm]/[1s]
gamma = 1.05                                        # Coeficiente de Arrasto do bloco cúbico no ar.
at0 = (k * xt0) / m                                 # Aceleração do bloco na posição em que é imediatamente solto do repouso.
passos = 850                                        # Número de Iterações.
tempomax = 100                                      # Tempo máximo que o eixo t vai correr com o gráfico(t=100 é um bom Coeficiente).
deltat = tempomax / passos                          # Intervalo de tempo é o "passo de iteração"
v = np.empty(passos + 1)                            # Cria um vetor que vai receber os valores de velocidade.
x = np.empty(passos + 1)                            # Cria um vetor que vai receber os valores de posição.
a = np.empty(passos + 1)                            # Cria um vetor que vai receber os valores de aceleração.
v[0] = vt0                                          # Define a velocidade Inicial
x[0] = xt0                                          # Define a Posição Inicial
a[0] = at0                                          # Define a Aceleração Inicial

# Aplicando o Algoritmo de Euler Para Evolução Dinâmica:
for i in range(0, passos):                                                     # PARA O I-ÉSIMO DADO NO ALCANÇE QUE VAI DE ZERO ATÉ O NÚMERO DE PASSOS FAÇA:
    a[i + 1] = (x[i] * k - gamma * v[i]) / m                                   # Implementando a Aceleração Para o Problema.
    v[i + 1] = (-k / m) * deltat * x[i] + v[i] - v[i] * deltat * gamma / m     # Implementando a velocidade do i-ésimo+1 dado de v, de acordo com Euler.
    x[i + 1] = deltat * v[i] + x[i]                                            # Implementando a posição do i-ésimo dado de x, de acordo com Euler.

iterations = 850
tempomax = 100
deltat = tempomax/iterations
t = np.linspace(0, iterations, iterations+1)
x1t0 = 10
v1t0 = 0
momento1t0 = 0
momento2t0 = 0
momento3t0 = 0
momento4t0 = 0
v1 = np.empty(iterations + 1)
x1 = np.empty(iterations + 1)
momento1 = np.empty(iterations + 1)
momento2 = np.empty(iterations + 1)
momento3 = np.empty(iterations + 1)
momento4 = np.empty(iterations + 1)
v1[0] = v1t0
x1[0] = x1t0
momento1[0] = momento1t0
momento2[0] = momento2t0
momento3[0] = momento3t0
momento4[0] = momento4t0


def EspacoDeFase1(tmax, iterations, A, B):
    global tempomax
    tempomax = tmax
    global deltat
    deltat = tempomax/iterations
    GAMMA = [0.1, 1.05, 2.7, 6]
    for i in range(0, 850):
        v1[i + 1] = (-k / m) * deltat * x1[i] + v1[i] - v[i]*deltat * GAMMA[0] / m
        x1[i + 1] = deltat * v1[i] + x1[i]
        momento1[i + 1] = m * v1[i]
    x2, y = np.meshgrid(np.linspace(-15, 15, 30), np.linspace(-30, 30, 30))
    U = A * y                 # Representa o campo do espaco de fase baseado na divisao da EDO de 2ª ordem em um sist. de duas EDOS de 1ª ordem.
    V = B * x2                # Representa o campo do espaco de fase baseado na divisao da EDO de 2ª ordem em um sist. de duas EDOS de 1ª ordem.
    plt.figure(figsize=(12, 6))
    plt.ylim(-30, 30)
    plt.xlim(-15, 15)
    plt.title('$Espaço\:de\:Fases\:Para\:Diferentes\:Coeficientes\:de\:Arrasto\:-$ $({\gamma = 0.1)}$', FontSize = 14)
    plt.quiver(x2, y, U, V, color = 'black')
    plt.plot(x1, momento1, label = '$\gamma = 0.1$', markersize=25, linewidth=1.8, color = 'blue')
    plt.xlabel('$X(t)$', FontSize = 14)
    plt.ylabel('$V(t)$', FontSize = 14)
    legend(fancybox=True, loc='upper left')
    grid(True)
    plt.savefig('EspacoDeFase1.png', transparent=True)
    plt.show()

EF1 = EspacoDeFase1(900, 100000, 1/2, -1)
#ef1 = EF1.transpose()
#plt.plot(x,v)


#print(len(ef1))
#print(len(EF1))

x_2t0 = 10
v_2t0 = 0
v_2 = np.empty(iterations + 1)
x_2 = np.empty(iterations + 1)
v_2[0] = v_2t0
x_2[0] = x_2t0

def EspacoDeFase2(tmax, iterations, A, B):
    global tempomax
    tempomax = tmax
    global deltat
    deltat = tempomax/iterations
    GAMMA = [0, 1.05, 2.7, 6]
    for i in range(0, 850):
        v_2[i + 1] = (-k / m) * deltat * x_2[i] + v_2[i] - v_2[i]*deltat * GAMMA[1] / m
        x_2[i + 1] = deltat * v_2[i] + x_2[i]
        momento2[i + 1] = m * v_2[i]
    x2, y = np.meshgrid(np.linspace(-15, 15, 30), np.linspace(-25, 20, 30))
    U = A*y                 # Representa o campo do espaco de fase baseado na divisao da EDO de 2ª ordem em um sist. de duas EDOS de 1ª ordem.
    V = B*x2                 # Representa o campo do espaco de fase baseado na divisao da EDO de 2ª ordem em um sist. de duas EDOS de 1ª ordem.
    plt.figure(figsize=(12, 6))
    plt.ylim(-25, 20)
    plt.xlim(-15, 15)
    plt.title('$Espaço\:de\:Fases\:Para\:Diferentes\:Coeficientes\:de\:Arrasto\:-\:(\gamma = 1.05)$', FontSize = 14)
    plt.quiver(x2, y, U, V, color = 'black')
    plt.plot(x_2, momento2, label = '$\gamma = 1.05$', linewidth=1.8, color = 'red')
    plt.xlabel('$X(t)$', FontSize = 14)
    plt.ylabel('$V(t)$', FontSize = 14)
    legend(fancybox=True, loc='upper left')
    grid(True)
    plt.savefig('EspacoDeFase2.png', transparent=True)
    plt.show()

EF2 = EspacoDeFase2(999, 100000, 1/2, -1)
#ef2 = EF2.transpose()

x3t0 = 10
v3t0 = 0
v3 = np.empty(iterations + 1)
x3 = np.empty(iterations + 1)
v3[0] = v3t0
x3[0] = x3t0

def EspacoDeFase3(tmax, iterations, A, B):
    global tempomax
    tempomax = tmax
    global deltat
    deltat = tempomax/iterations
    GAMMA = [0, 1.05, 2.7, 6]
    for i in range(0, 850):
        v3[i + 1] = (-k / m) * deltat * x3[i] + v3[i] - v3[i]*deltat * GAMMA[2] / m
        x3[i + 1] = deltat * v3[i] + x3[i]
        momento3[i + 1] = m * v3[i]
    x2, y = np.meshgrid(np.linspace(-15, 15, 30), np.linspace(-20, 20, 30))
    U = A * y                   # Representa o campo associado ao espaco de fases baseado na divisao da EDO de 2ª ordem em um sist. de duas EDOS de 1ª ordem.
    V = B * x2                  # Representa o campo associado ao espaco de fases baseado na divisao da EDO de 2ª ordem em um sist. de duas EDOS de 1ª ordem.
    plt.figure(figsize=(12, 6))
    plt.ylim(-20, 20)
    plt.xlim(-15, 15)
    plt.title('$Espaço\:de\:Fases\:Para\:Diferentes\:Coeficientes\:de\:Arrasto\:-\:(\gamma = 2.7)$', FontSize = 14)
    plt.quiver(x2, y, U, V, color = 'black')
    plt.plot(x3, momento3, label = '$\gamma = 2.7$', linewidth=1.8, color = 'green')
    plt.xlabel('$X(t)$', FontSize = 14)
    plt.ylabel('$V(t)$', FontSize = 14)
    legend(fancybox=True, loc='upper left')
    grid(True)
    plt.savefig('EspacoDeFase3.png', transparent=True)
    plt.show()

EF3 = EspacoDeFase3(999, 100000, 1/2, -1)

x4t0 = 10
v4t0 = 0
v4 = np.empty(iterations + 1)
x4 = np.empty(iterations + 1)
v4[0] = v4t0
x4[0] = x4t0

def EspacoDeFase4(tmax, iterations, A, B):
    global tempomax
    tempomax = tmax
    global deltat
    deltat = tempomax/iterations
    GAMMA = [0, 1.05, 2.7, 6]
    for i in range(0, 850):
        v4[i + 1] = (-k / m) * deltat * x4[i] + v4[i] - v4[i] * deltat * GAMMA[3] / m
        x4[i + 1] = deltat * v4[i] + x4[i]
        momento4[i + 1] = m * v4[i]
    x2, y = np.meshgrid(np.linspace(-15, 15, 30), np.linspace(-20, 20, 30))
    U = A * y                   # Representa o campo associado ao espaco de fases baseado na divisao da EDO de 2ª ordem em um sist. de duas EDOS de 1ª ordem.
    V = B * x2                  # Representa o campo associado ao espaco de fases baseado na divisao da EDO de 2ª ordem em um sist. de duas EDOS de 1ª ordem.
    plt.figure(figsize=(12, 6))
    plt.ylim(-20, 20)
    plt.xlim(-15, 15)
    plt.title('$Espaço\:de\:Fases\:Para\:Diferentes\:Coeficientes\:de\:Arrasto\:-\:(\gamma = 6)$', FontSize = 14)
    plt.quiver(x2, y, U, V, color = 'black')
    plt.plot(x4, momento4, label = '$\gamma = 6$', linewidth=1.8, color = 'purple')
    plt.xlabel('$X(t)$', FontSize = 14)
    plt.ylabel('$V(t)$', FontSize = 14)
    legend(fancybox=True, loc='upper left')
    grid(True)
    #plt.tick_params(axis='y', direction='in', length=10)
    #plt.tick_params(axis='x', direction='in', length=10)
    plt.savefig('EspacoDeFase4.png', transparent=True)
    plt.show()

EF4 = EspacoDeFase4(999, 100000, 1/2, -1)        # A = 1/2, B = -0.06

#print(len(ef1))
#print(len(x))
#print(len(t))
#print(len(v))

#MULTIPLOT COM QUATRO FIGURAS PARA O ESPAÇO DE FASES CORRESPONDENTE A QUATRO VALORES DISTINTOS DE GAMMA:

A = 1/2
B = -1
fig = plt.figure(figsize=(9, 6))
x2, y = np.meshgrid(np.linspace(-15.6, 15.6, 30), np.linspace(-33, 33, 30))
U = A * y                   # Representa o campo associado ao espaco de fases baseado na divisao da EDO de 2ª ordem em um sist. de duas EDOS de 1ª ordem.
V = B * x2                  # Representa o campo associado ao espaco de fases baseado na divisao da EDO de 2ª ordem em um sist. de duas EDOS de 1ª ordem.
ax1 = fig.add_subplot(1, 1, 1)
plt.plot(x1, v1, label = '$\gamma = 0.1$', linewidth=1.8, color='blue')
ax1.plot(x_2, v_2, label = '$\gamma = 1.05$', linewidth=1.8, color='red')
ax1.plot(x3, v3, label = '$\gamma = 2.7$', linewidth=1.8, color='#22DD33')
ax1.plot(x4, v4, label = '$\gamma = 6$', linewidth=1.8, color='purple')
plt.title('$Espaço\:de\:Fases\:Para\:Diferentes\:Valores\:de\: \gamma$', FontSize = 14)
#plt.quiver(x2, y, U, V, color = 'black')
plt.xlim(-15.6, 15.6)
plt.ylim(-33, 33)
plt.xlabel('$X(t)$', FontSize = 14)
plt.ylabel('$V(t)$', FontSize = 14)
legend(fancybox=True, loc='upper left')
grid(True)
#plt.tick_params(axis='y', direction='in', length=10)
#plt.tick_params(axis='x', direction='in', length=10)
plt.savefig('EspacoDeFase4.png', transparent=True)
plt.show()
