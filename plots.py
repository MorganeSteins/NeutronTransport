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

# text_file = open("Data/points_q4.txt", "r")
# lines = text_file.read().split(',')
# lines.pop()
# lines = np.asarray(lines).astype(np.float)
# x = np.linspace(0,1,25)
# #plt.scatter(np.log10(x),np.log10(np.asarray(lines).astype(np.float)))
# #plt.plot(np.log10(x),-0.5*np.log10(x))
# phi_exact = [density_no_scattering_homog_point(x[i],1,sigmaT) for i in range(len(x))]
# # erreur = (lines-phi_exact)**2
# # print(np.sqrt(np.sum(erreur)))
# erreur = np.linalg.norm(lines-phi_exact)
# plt.scatter(x,np.asarray(lines).astype(np.float), label='simulations')
# plt.plot(x,phi_exact, label='solution exacte')
# plt.legend()
# plt.title("Tirages moyennés sur 5 simulations")
# plt.show()

# N = [1,10,100,1000,10000,100000,1000000]
# erreur_L2 = np.zeros(len(N))
# erreur_inf = np.zeros(len(N))
# for i in range(len(N)) :
#     text_file = open("Data/points_q4_"+str(N[i])+".txt", "r")
#     lines = text_file.read().split(',')
#     lines.pop()
#     lines = np.asarray(lines).astype(np.float)
#     phi_exact = [density_no_scattering_homog_point(x[i],1,sigmaT) for i in range(len(x))]
#     erreur_L2[i] = np.linalg.norm(lines-phi_exact)
#     erreur_inf[i] = max(abs(lines-phi_exact))

# plt.scatter(np.log10(N),np.log10(erreur_L2), label="Erreur L2")
# plt.scatter(np.log10(N),np.log10(erreur_inf),label="Erreur L infini")
# plt.plot(np.log10(N),-0.25*np.log10(N)+0.2, label= 'O(N^{-1/4})')
# plt.xlabel('Nombre de points en échelle log10')
# plt.ylabel("Erreur en échelle log")
# plt.legend()
# plt.title("Erreur en fonction du nombre de trajectoires - Source ponctuelle")
# plt.show()

mu=1.
N = [1,10,100,1000,10000,100000,1000000]
erreur_L2 = np.zeros(len(N))
erreur_inf = np.zeros(len(N))
for i in range(len(N)) :
    text_file = open("Data/points_mu01_q5_"+str(N[i])+".txt", "r")
    lines = text_file.read().split(',')
    lines.pop()
    lines = np.asarray(lines).astype(np.float)
    phi_exact = [density_no_scattering_homog_unif(x[i],mu,sigmaT) for i in range(len(x))]
    erreur_L2[i] = np.linalg.norm(lines-phi_exact)
    erreur_inf[i] = max(abs(lines-phi_exact))

plt.scatter(np.log10(N),np.log10(erreur_L2), label="Erreur L2")
plt.scatter(np.log10(N),np.log10(erreur_inf),label="Erreur L infini")
plt.plot(np.log10(N),-0.25*np.log10(N)+0.2, label= 'O($N^{-1/4}$)')
plt.xlabel('Nombre de points en échelle log10')
plt.ylabel("Erreur en échelle log")
plt.legend()
plt.title("Erreur en fonction du nombre de trajectoires - Source uniforme")
plt.show()