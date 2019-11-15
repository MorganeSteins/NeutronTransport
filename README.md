## AMS302

Commande de compilation

# Monte carlo

g++ main_monte_carlo.cpp points.cpp MC.cpp aleat_tools.cpp -o monte_carlo.out

# Determininiste

g++ main_deterministe.cpp deterministe.cpp -o deterministe.out -I PATH/TO/EIGEN

# Deterministe DSA

g++ main_deterministe_fast.cpp deterministe.cpp -o deterministe_fast.out -I PATH/TO/EIGEN
