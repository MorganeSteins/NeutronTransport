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


sigmaT = 1
n=100
dx=1/n
sigmaT = 1
N = 50
dx=1/N
text_file = open("Data/phi_q9_"+str(N)+".txt", "r")
lines = text_file.read().split(',')
lines.pop()
lines = np.asarray(lines).astype(np.float)

#file_2 = open("Data/phi_q9_"+str(N)+".txt", "r")

# x = np.linspace(0,1-1/n,n)
# x = np.concatenate((np.array([0]),np.repeat(np.linspace(1/n,1-1/n,n-1),2),np.array([1])),axis=0)

x = np.linspace(0,1,N+1)
phi_exact = np.array([density_no_scattering_homog_point(x[i],1,sigmaT) for i in range(len(x))])
# var = np.sqrt(lines*dx*sigmaT*(1-lines*dx*sigmaT)/(N*dx))
#var = np.sqrt(phi_exact*dx*sigmaT*(1-phi_exact*dx*sigmaT)/(N*dx))
print(len(x))
print(len(lines))
plt.scatter(x,lines, label='simulations pour '+str(N)+' intervalles', linewidth=0.01)
plt.plot(x,phi_exact, label="solution theorique")

# plt.plot(x,np.repeat(lines,2), label='simulations pour N='+str(N))
#plt.scatter(x+dx/2,lines, linewidth=0.1)
# print(lines)
#plt.plot(x+dx/2,phi_exact, label='solution exacte')
# plt.plot(x+dx/2,phi_exact-1.96*var,'r')
# plt.plot(x+dx/2,phi_exact+1.96*var,'r')
# plt.plot(x+dx/2,lines-1.96*var,'r',label='intervalle de confiance')
# plt.plot(x+dx/2,lines+1.96*var,'r')
plt.legend()
# plt.title("Simulation sur "+str(n)+" intervalles avec N="+str(N)+"particules")
plt.title("Solution sans scattering et courant entrant unitaire")
plt.show()
