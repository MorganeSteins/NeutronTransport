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

def density_no_scattering_non_homog(x,mu, sigmaT):
    if (x<0.3):
        return (1/mu)*(np.exp(-x/mu))
    elif (x<0.7):
        return (1/mu)*(np.exp(-3*(x-0.3)/mu))*(np.exp(-0.3/mu))
    else :
        return (1/mu)*np.exp(-(x-0.7)/mu)*(np.exp(-3*(0.7-0.3)/mu))*(np.exp(-0.3/mu))

def density_limite_diffusion(x) :
    return 1.5*x*(1-x)


sigmaT = 1
n_deter=100
n_MC=100
N = 10000000
Nmu=100
dx = 1/n_deter
mu=1

# MC = open("Data/phi_q5_"+str(N)+"_"+str(n_MC)+".txt", "r")
# MC = open("Data/phi_q7_"+str(N)+"_"+str(n_MC)+".txt", "r")
# MC = MC.read().split(',')
# MC.pop()
# MC = np.asarray(MC).astype(np.float)

deter = open("Data/phi_q15_"+str(n_deter)+"_"+str(Nmu)+"_0.100000.txt", "r")
# deter = open("Data/phi_q11_"+str(n_deter)+".txt", "r")
# deter = open("Data/phi_13_"+str(n_deter)+"_1.000000.txt", "r")
deter = deter.read().split(',')
deter.pop()
deter = np.asarray(deter).astype(np.float)

deter_ref = open("Data/phi_q15_"+str(n_deter)+"_"+str(Nmu)+"_0.010000.txt", "r")
# deter = open("Data/phi_q11_"+str(n_deter)+".txt", "r")
# deter = open("Data/phi_13_"+str(n_deter)+"_1.000000.txt", "r")
deter_ref = deter_ref.read().split(',')
deter_ref.pop()
deter_ref = np.asarray(deter_ref).astype(np.float)

# deter_ref2 = open("Data/phi_fast_"+str(n_deter)+"_"+str(Nmu)+"_epsilon0.000100.txt", "r")
# deter_ref2 = deter_ref2.readlines()#.split(',')
# deter_ref2 = np.asarray(deter_ref2).astype(np.float)



# x = np.linspace(0,1,n_deter+1)
# x_cst = np.concatenate((np.array([0]),np.repeat(np.linspace(1/n_deter,1-1/n_deter,n_deter-1),2),np.array([1])),axis=0)

# x = np.linspace(0,1,n_deter+1)
# print(x)
# print(len(x[0:len(x)-1]))
# print(len(deter))
x = np.linspace(0,1-1/n_deter,n_deter)

x_phi = np.concatenate((np.array([x[0]]),x+dx/2,np.array([x[-1]+dx])), axis=0)
# phi_exact = np.array([density_no_scattering_non_homog(x_phi[i],1,sigmaT) for i in range(len(x_phi))])
phi_exact = density_limite_diffusion(x)
# var = np.sqrt(MC*dx*sigmaT*(1-MC*dx*sigmaT)/(N*dx))

# var = np.sqrt(phi_exact*dx*sigmaT*(1-phi_exact*dx*sigmaT)/(N))
# plt.plot(x_cst,np.repeat(deter,2), label='déterministe '+str(n_deter)+' intervalles', linewidth=1)
# plt.scatter(x_cst,np.repeat(deter,2), label='déterministe '+str(n_deter)+' intervalles', linewidth=0.01)
# plt.plot(x,phi_exact, label="solution theorique")
# print(len(x))
# print(len(deter))
# print(np.max(deter))

# print(len(x))
# print(len(deter))
# # plt.scatter(x+dx/2,MC, marker = '*',label='Monte-Carlo pour N='+str(N), linewidth=2.5)
x_double = np.concatenate((np.array([0]),np.repeat(np.linspace(1/n_deter,1-1/n_deter,n_deter-1),2),np.array([1])),axis=0)
# plt.scatter(x,MC, label='Monte-Carlo pour $N=10^4$')
# plt.scatter(x[0:len(x)-1]+dx/2,deter, linewidths=0.01, label='$\mu=1$')
# plt.subplot(1,2,1)
plt.plot(x_double,np.repeat(deter,2), label='$\epsilon=0.1$ : erreur L2=0.13')
plt.plot(x_double,np.repeat(deter_ref,2), label='$\epsilon=0.01$ : erreur L2=0.013',linewidth=2)
# plt.plot(x_double,np.repeat(deter_ref2,2),'k', label='$\epsilon=0.0001$', linewidth=2)
plt.plot(x,phi_exact, 'r--', label='Limite de diffusion')
# plt.legend()
# plt.ylabel("$\phi$")
# plt.xlabel("$x$")
# plt.title("$\phi(x)$ moyenné : deterministe")
# plt.subplot(1,2,2)
# plt.plot(x_double,np.repeat(MC,2), label='Calcul Monte-Carlo')
# plt.plot(x_phi,phi_exact, label="Solution exacte")
# plt.ylim([np.min(deter)-0.1*abs(np.min(deter)), np.max(deter)+0.1*abs(np.max(deter))])
# # print(lines)
# plt.plot(x_phi,phi_exact,'r', label='Solution exacte')
# # plt.plot(x,phi_exact-1.96*var/dx,'r', label='Intervalle de fluctuation')
# # plt.plot(x,phi_exact+1.96*var/dx,'r')
# plt.plot(x,phi_exact, label='solution exacte')
# plt.plot(x,np.array([density_no_scattering_non_homog(x[i],1) for i in range(len(x))]),label="$\mu=1$")
# plt.plot(x,np.array([density_no_scattering_non_homog(x[i],0.7) for i in range(len(x))]),label="$\mu=0.7$")
# plt.plot([0.3,0.3],[0,1.4],'k--', linewidth=1)
# plt.plot([0.7,0.7],[0,1.4],'k--', linewidth=1)
# plt.plot(x,phi_exact-1.96*var/dx,'r',label='intervalle de fluctuation')
# plt.plot(x,phi_exact+1.96*var/dx,'r')
plt.legend()
plt.ylabel("$\phi$")
plt.xlabel("$x$")
plt.title("$\phi(x)$ moyenné : limite de diffusion")
# # # plt.title("Simulation sur "+str(n_deter)+" intervalles avec N="+str(N)+"particules")
# # plt.title("Simulation iterative sans scattering, source ponctuelle")
# # plt.title("$\psi(x)$ moyenné en angle pour "+str(n_deter)+" intervalles")
plt.show()

# print(deter)
# x_phi_t = np.linspace(0,1,len(deter)+1)
# phi_exact = np.array([density_no_scattering_non_homog(x_phi_t[i]+dx/2,1,sigmaT) for i in range(len(x_phi_t)-1)])
# phi_t_exact = np.array([(1/(2*mu))*(phi_exact[i]+phi_exact[i+1]) for i in range(len(x_phi_t)-1)])
# plt.scatter(x_phi_t[0:len(x_phi_t)-1]+dx/2,phi_t_exact)
# plt.scatter(x+dx/2,deter, marker='*')
# plt.show()

errl2 = np.sum(((deter-phi_exact)**2)*dx)
print(np.sqrt(errl2)/np.sqrt(np.sum(phi_exact**phi_exact*dx)))
errl2 = np.sum(((phi_exact-deter_ref)**2)*dx)
print(np.sqrt(errl2)/np.sqrt(np.sum(phi_exact**phi_exact*dx)))

Ns = np.array([2,10,100,500,1000])
errs = np.array([0.1798, 0.006607,0.00017956, 1.65*10**-5, 6.55*10**-6])

# Ns = np.array([10,20,50,100,200,500])
# errs = [0.002997,0.0007519,0.0001205,3.007*10**-5,7.51*10**-6,1.2609*10**-6]
# # errs = [0.00304,0.000759,0.0001898,3.0358*10**-5,7.59*10**-6,1.873*10**-6,4.43*10**-7]
# # errs = np.sqrt([2.382*10**-6,1.498*10**-7,9.37*10**-9, 2.389*10**-10,1.522*10**-11,9.358*10**-13])#,1.106*10**-13,9.097*10**-14])
# plt.loglog(Ns,errs, 'o', label='Erreur $\mathcal{L}^2$')
# plt.loglog(Ns[0:6],0.5/Ns[0:6]**1.5,'--', label='$O(N^{-1.5})$')
# plt.legend(fontsize=15)
# plt.xlabel('$N_\mu$')
# plt.ylabel('Erreur')
# plt.title("Erreur en fonction de $N_\mu$", fontsize=15)
# # # plt.plot(x[0:len(x)-1],phi_t_exact)
# plt.show()