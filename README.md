# AMS302

## Commande de compilation

### Monte carlo

g++ main_monte_carlo.cpp src/utils/points.cpp src/mc/MC.cpp src/utils/aleat_tools.cpp -o monte_carlo.out

### Determininiste

g++ main_deterministe.cpp src/deterministe/deterministe.cpp -o deterministe.out -I PATH/TO/EIGEN

### Deterministe DSA

g++ main_deterministe_fast.cpp src/deterministe/deterministe.cpp -o deterministe_fast.out -I PATH/TO/EIGEN



## How it works
See report (french) for explanation of methods.
- Monte-Carlo for stochastic methods
- Source iteration for deterministic methods

Different sources are usable.