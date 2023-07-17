#Trabalho Numérico - TCM (Turma 1) - 2023.1
#Evelyn Soares Pereira - 170102785

import numpy as np
import matplotlib.pyplot as plt

#Problema 1 -> usarei o exemplo do roteiro no código, pois vou fazer algumas funções baseadas em algumas partes do código do roteiro
N = 10
L = 1.0
delta_x = L/N
delta_t = 0.2*delta_x**2
n_final = 10
t_final = 100
t = 0
x = np.linspace(0.0, L, N+1) 

Temp = np.ones(N+1, float) #preencho o array com 1's 
Temp[0] = 0.0 #condições de contorno
Temp[N] = 0.0
Temp_nova = np.copy(Temp)

def edd(i): #função para a equação de diferenças finitas
    Temp_nova[i] = Temp[i]+(delta_t/delta_x**2)*(Temp[i+1]-2.0*Temp[i]+Temp[i-1])
    return Temp_nova[i]

def somatorio(x, t, n_final): #função para o calculo da solução exata do problema 1
    resultado = np.zeros_like(x)
    for n in range(1, n_final):
        resultado += (4*np.sin((2*n-1)*np.pi*x)*np.exp(-((2*n-1)**2)*np.pi**2*t))/((2*n-1)*np.pi)
    return resultado

while t < t_final:
    for i in range(1, N):
        edd(i)
    Temp = np.copy(Temp_nova)
    t += delta_t

resultado = somatorio(x, t, n_final)

fig = plt.figure() #visualização do gráfico do problema 1
ax = fig.add_subplot()
fig.suptitle('t = %.3f' %t, fontsize = 18, fontweight = 'bold')
ax.set_ylabel('T (°C)', fontsize = 18)
ax.set_xlabel('x (m)', fontsize = 18)
plt.plot(x, Temp, '-r', lw = 4) #solução numérica (vermelho)
plt.plot(x, resultado, '-g', lw=4) #solução exata (verde)
plt.savefig('figura_p1.pdf', format = 'pdf', dpi = 1200, bbox_inches = 'tight')
plt.show()

#Problema 2 (so irei redefinir o que vai mudar)
t = 0
Temp = np.zeros(N + 1, float) #preencho o array com 0's 
Temp[0] = 1.0 #nova condição de contorno
Temp_nova = np.copy(Temp)

def edd(i): #função para a equação de diferenças finitas
    Temp_nova[i] = Temp[i]+(delta_t/delta_x**2)*(Temp[i+1]-2.0*Temp[i]+Temp[i-1])
    return Temp_nova[i]

def somatorio(x, t, n_final): #função para o calculo da solução exata do problema 2
    resultado = np.zeros_like(x)
    for n in range(1, n_final):
        resultado += (2*np.sin(n*np.pi*x)*np.exp(-n**2*np.pi**2*t))/((2*n)*np.pi)
    return 1-x-resultado

while t < t_final:
    for i in range(1, N):
        edd(i)
    Temp = np.copy(Temp_nova)
    t += delta_t

resultado = somatorio(x, t, n_final)

#visualização do gráfico do problema 2
fig = plt.figure()
ax = fig.add_subplot()
fig.suptitle('t = %.3f' %t, fontsize = 18, fontweight = 'bold')
ax.set_ylabel('T (°C)', fontsize = 18)
ax.set_xlabel('x (m)', fontsize = 18)
plt.plot(x, Temp, '-r', lw = 4) #soluçao numérica (vermelho)
plt.plot(x, resultado, '-g', lw=4) #solução exata (verde), os resultados devem se sobrepor
plt.savefig('figura_p2.pdf', format = 'pdf', dpi = 1200, bbox_inches = 'tight')
plt.show()

#Problema 3 (so irei redefinir o que vai mudar)
L = 2.0 #x vai de 0 a 2
t = 0
delta_x = L/N
delta_t = 0.2*delta_x**2
x = np.linspace(0.0, L, N+1)

Temp = np.sin((np.pi/2)*x) #nova função de temp
Temp[0] = 0.0 #novas condições de contorno 
Temp_nova = np.copy(Temp)

def edd(i): #função para a equação de diferenças finitas
    Temp_nova[i] = Temp[i]+(delta_t/delta_x**2)*(Temp[i+1]-2.0*Temp[i]+Temp[i-1])
    return Temp_nova[i]

while t < t_final:
    for i in range(1, N):
        edd(i)
    Temp = np.copy(Temp_nova)
    t += delta_t

resultado = np.sin((np.pi/2)*x)*np.exp(-np.pi*np.pi*t/4)

#visualização do gráfico do problema 3
fig = plt.figure()
ax = fig.add_subplot()
fig.suptitle('t = %.3f' %t, fontsize = 18, fontweight = 'bold')
ax.set_ylabel('T (°C)', fontsize = 18)
ax.set_xlabel('x (m)', fontsize = 18)
plt.plot(x, Temp, '-r', lw = 4)  #soluçao numérica (vermelho)
plt.plot(x, resultado, '-g', lw=4)  #soluçao exata (verde)
plt.savefig('figura_p3.pdf', format = 'pdf', dpi = 1200, bbox_inches = 'tight')
plt.show()