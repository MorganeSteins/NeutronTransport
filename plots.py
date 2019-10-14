import matplotlib.pyplot as plt
import numpy as np

def density_no_scattering_homog_point(x,mu,sigmaT) :
    if mu == 0 and x==0:
        return 1
    elif mu==0 :
        return 0
    return (1/mu)*(np.exp(-sigmaT*x/mu))

def density_no_scattering_homog_unif(x,mu,sigmaT):
    if mu==0 :
        return 1
    elif mu<0: 
        return (1./ sigmaT)*(1-np.exp(-sigmaT*(x-1)/mu))
    return (1./ sigmaT)*(1-np.exp(-sigmaT*x/mu))


# Plots de la solution q3 en fonction de x pour différents mu
sigmaT = 1
x = np.linspace(0,1,25)
# liste_x = np.linspace(0,1,100)
# plt.plot(liste_x,[density_no_scattering_homog_point(liste_x[i],1,sigmaT) for i in range(len(liste_x))], label='mu=1')
# plt.plot(liste_x,[density_no_scattering_homog_point(liste_x[i],0.5,sigmaT) for i in range(len(liste_x))],label='mu=0.5')
# plt.plot(liste_x,[density_no_scattering_homog_point(liste_x[i],0.1,sigmaT) for i in range(len(liste_x))],label='mu=0.1')
# plt.legend()
# plt.title("Solution pour différents mu positifs")
# plt.show()

# Plot des erreurs de q4
# text_file = open("erreurs_q4.txt", "r")
# lines = text_file.read().split(',')
# lines.pop()
# print(type(lines[0]))
# x = [1,2,10,20,100,200,1000,2000,10000,20000,100000,200000,1000000,2000000]#,10000000,20000000]
# plt.scatter(np.log10(x),np.log10(np.asarray(lines).astype(np.float)))
# plt.plot(np.log10(x),-0.5*np.log10(x))
# plt.xlabel('Nombre de points en échelle log10')
# plt.ylabel("Erreur en échelle log")
# plt.title("Erreur en fonction du nombre de trajectoires")
# plt.show()

# Intervalles de confiance de q4


# Plot des points de q4
n=25
dx=1/n
sigmaT = 1
N = 10000
text_file = open("Data/phi_q4_"+str(N)+".txt", "r")
lines = text_file.read().split(',')
lines.pop()
lines = np.asarray(lines).astype(np.float)

# freq = np.array([3895,3932,3541,3360,3325,3279,3098,2965,2863,2759,2570,2625,2391,2346,2266,2213,2069,1889,1917,1784,1702,1676,1643,1562,1493])/N
# print(np.sum(freq))
# freq = np.array([392404,376353,361598,347797,333472,321263,308380,297309,284905,273324,262108,251730,241830,233300,224364,215049,205609,198692,191177,182523,176532,169264,162662,155895,150368])/N
x = np.linspace(0,1-1/n,n)
# x = np.concatenate((np.array([0]),np.repeat(np.linspace(1/n,1-1/n,n-1),2),np.array([1])),axis=0)

#plt.plot(np.log10(x),-0.5*np.log10(x))
phi_exact = np.array([density_no_scattering_homog_point(x[i]+dx/2,1,sigmaT) for i in range(len(x))])
#var = np.sqrt(freq*dx*sigmaT*(1-freq*dx*sigmaT)/(N*dx))
var = np.sqrt(phi_exact*dx*sigmaT*(1-phi_exact*dx*sigmaT)/(N*dx))
print(x)
print(len(lines))
# erreur = (lines-phi_exact)**2
# print(np.sqrt(np.sum(erreur)))
# print(x)
# print(np.repeat(lines,2))
# print(len(x),len(np.repeat(lines,2)))
#erreur = np.linalg.norm(lines-phi_exact)
#plt.plot(x,np.repeat(lines,2), label='simulations pour N=$10^7$')
plt.scatter(x+dx/2,lines, linewidth=0.1)
# print(lines)
#plt.plot(x+dx/2,phi_exact, label='solution exacte')
plt.plot(x+dx/2,phi_exact-1.96*var,'r')
plt.plot(x+dx/2,phi_exact+1.96*var,'r')
plt.legend()
plt.title("Simulation sur "+str(n)+" intervalles avec $N=10^4$ particules")
plt.show()

# N = [1,10,100,1000,10000,100000,1000000]
# erreur_L2 = np.zeros(len(N))
# erreur_inf = np.zeros(len(N))
# for i in range(len(N)) :
#     text_file = open("Data/points_q4_"+str(N[i])+".txt", "r")
#     lines = text_file.read().split(',')
#     lines.pop()
#     lines = np.asarray(lines).astype(np.float)
#     phi_exact = [density_no_scattering_homog_point(x[i],1,sigmaT) for i in range(len(x))]
#     erreur_L2[i] = 0.04*np.linalg.norm(lines-phi_exact)
#     erreur_inf[i] = max(abs(lines-phi_exact))

# plt.scatter(np.log10(N),np.log10(erreur_L2), label="Erreur $L^2$")
# plt.scatter(np.log10(N),np.log10(erreur_inf),label="Erreur $L^\infty$")
# plt.plot(np.log10(N),-0.25*np.log10(N)-1, label= 'O($N^{-1/4}$)')
# plt.xlabel('Nombre de points en échelle log10')
# plt.ylabel("Erreur en échelle log")
# plt.legend()
# plt.title("Erreur en fonction du nombre de trajectoires - Source ponctuelle")
# plt.show()

# mu=0.1
# N = [1,10,100,1000,10000,100000,1000000]
# erreur_L2 = np.zeros(len(N))
# erreur_inf = np.zeros(len(N))
# for i in range(len(N)) :
#     text_file = open("Data/points_mu01_q5_"+str(N[i])+".txt", "r")
#     lines = text_file.read().split(',')
#     lines.pop()
#     lines = np.asarray(lines).astype(np.float)
#     phi_exact = [density_no_scattering_homog_unif(x[i],mu,sigmaT) for i in range(len(x))]
#     erreur_L2[i] = 0.04*np.linalg.norm(lines-phi_exact)
#     erreur_inf[i] = max(abs(lines-phi_exact))

# plt.scatter(np.log10(N),np.log10(erreur_L2), label="Erreur $L^2$")
# plt.scatter(np.log10(N),np.log10(erreur_inf),label="Erreur $L^\infty$")
# plt.plot(np.log10(N),-0.25*np.log10(N)-1, label= 'O($N^{-1/4}$)')
# plt.xlabel('Nombre de points en échelle log10')
# plt.ylabel("Erreur en échelle log")
# plt.legend()
# plt.title("Erreur en fonction du nombre de trajectoires - Source uniforme")
# plt.show()

# liste_x = np.linspace(0,1,100)
# plt.plot(liste_x,[density_no_scattering_homog_unif(liste_x[i],1,sigmaT) for i in range(len(liste_x))], label='mu=1')
# plt.plot(liste_x,[density_no_scattering_homog_unif(liste_x[i],0.5,sigmaT) for i in range(len(liste_x))],label='mu=0.5')
# plt.plot(liste_x,[density_no_scattering_homog_unif(liste_x[i],0.1,sigmaT) for i in range(len(liste_x))],label='mu=0.1')
# plt.legend()
# plt.title("Solution pour différents mu positifs")
# plt.show()

# # Question 8 plot de la solution pour N=10000000
# text_file = open("Data/points_q8_"+str(10000000)+".txt", "r")
# sol = text_file.read().split(',')
# sol.pop()
# sol = np.array(sol).astype(float)
# plt.scatter(x,np.array(lines).astype(float),label="Simulations")
# plt.legend()
# plt.title("Calcul Monte-Carlo pour $N=10^7$")
# plt.show()
# print(np.sum(0.04*np.array(sol).astype(float)))

# N = [1,10,100,1000,10000,100000,1000000]
# erreur_L2 = np.zeros(len(N))
# erreur_inf = np.zeros(len(N))
# for i in range(len(N)) :
#     text_file = open("Data/points_q8_"+str(N[i])+".txt", "r")
#     lines = text_file.read().split(',')
#     lines.pop()
#     lines = np.asarray(lines).astype(np.float)
#     erreur_L2[i] = 0.04*np.linalg.norm(lines-sol)
#     erreur_inf[i] = max(abs(lines-sol))
# plt.scatter(np.log10(N),np.log10(erreur_L2), label="Erreur $L^2$")
# plt.scatter(np.log10(N),np.log10(erreur_inf),label="Erreur $L^\infty$")
# plt.plot(np.log10(N),-0.25*np.log10(N)-1, label= 'O($N^{-1/4}$)')
# plt.xlabel('Nombre de points en échelle log10')
# plt.ylabel("Erreur en échelle log p/r à $N=10^7$")
# plt.legend()
# plt.title("Erreur en fonction du nombre de trajectoires - Cas diffusant")
# plt.show()