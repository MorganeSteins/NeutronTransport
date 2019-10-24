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

def density_no_scattering_non_homog(x,mu):
    if (x<0.3):
        return (1/mu)*(np.exp(-x/mu))
    elif (x<0.7):
        return (1/mu**2)*(np.exp(-3*(x-0.3)/mu))*(np.exp(-0.3/mu))
    else :
        return  (1/mu**3)*np.exp(-(x-0.7)/mu)*(np.exp(-3*(0.7-0.3)/mu))*(np.exp(-0.3/mu))

sigmaT = 1
n_deter=50
n_MC=25
sigmaT = 1
N = 1000000
dx = 1/n_MC
mu=1

# MC = open("Data/phi_q5_"+str(N)+".txt", "r")
# MC = MC.read().split(',')
# MC.pop()
# MC = np.asarray(MC).astype(np.float)

deter = open("Data/phi_t_q11_"+str(n_deter)+"_1.000000.txt", "r")
# deter = open("Data/phi_q10_"+str(n_deter)+".txt", "r")
deter = deter.read().split(',')
deter.pop()
deter = np.asarray(deter).astype(np.float)

#file_2 = open("Data/phi_q9_"+str(N)+".txt", "r")

# x = np.linspace(0,1-1/n_deter,n_deter)
x_cst = np.concatenate((np.array([0]),np.repeat(np.linspace(1/n_deter,1-1/n_deter,n_deter-1),2),np.array([1])),axis=0)

x = np.linspace(0,1,n_deter+1)
phi_exact = np.array([density_no_scattering_non_homog(x[i],1) for i in range(len(x))])
# var = np.sqrt(lines*dx*sigmaT*(1-lines*dx*sigmaT)/(N*dx))
#var = np.sqrt(phi_exact*dx*sigmaT*(1-phi_exact*dx*sigmaT)/(N*dx))
plt.plot(x_cst,np.repeat(deter,2), label='déterministe '+str(n_deter)+' intervalles', linewidth=1)
# plt.scatter(x_cst,np.repeat(deter,2), label='déterministe '+str(n_deter)+' intervalles', linewidth=0.01)
plt.plot(x,phi_exact, label="solution theorique")

# x = np.linspace(0,1-1/n_MC,n_MC)
# plt.scatter(x+dx/2,MC, marker = '*',label='Monte-Carlo pour N='+str(N), linewidth=2.5)
# x = np.concatenate((np.array([0]),np.repeat(np.linspace(1/n_MC,1-1/n_MC,n_MC-1),2),np.array([1])),axis=0)
# plt.plot(x,np.repeat(MC,2), label='Monte-Carlo pour N='+str(N))
#plt.scatter(x+dx/2,lines, linewidth=0.1)
# print(lines)
#plt.plot(x+dx/2,phi_exact, label='solution exacte')
# plt.plot(x+dx/2,phi_exact-1.96*var,'r')
# plt.plot(x+dx/2,phi_exact+1.96*var,'r')
# plt.plot(x+dx/2,lines-1.96*var,'r',label='intervalle de confiance')
# plt.plot(x+dx/2,lines+1.96*var,'r')
plt.legend()
# plt.title("Simulation sur "+str(n_deter)+" intervalles avec N="+str(N)+"particules")
plt.title("Solution sans scattering et courant entrant unitaire")
plt.show()
