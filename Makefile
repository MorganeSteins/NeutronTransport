x=1
mu=1
type='point'
N = 100000
Nx = 5
Nmu = 2
max_iter = 50
nb_points = 25


all: execute_all

execute_all: 
	@mkdir -p Data
	@echo "Compiling mc, deterministe and dsa in this order ... "
	@g++ main_MC.cpp src/utils/*.cpp src/mc/*.cpp -o mc.out
	@g++ main_deterministe.cpp src/utils/*.cpp src/deterministe/deterministe.cpp -o deter.out
	@g++ main_deterministe_fast.cpp src/utils/*.cpp src/deterministe/deterministe.cpp src/deterministe/fast_IS.cpp -o deter_fast.out
	@echo "Compiled all files"
	@time ./mc.out 
	@time ./deter.out $(Nx) $(Nmu)
	@time ./deter_fast.out $(Nx) $(Nmu)
	@echo "\n Executed all files, done."

mc:
	@mkdir -p Data
	g++ main_MC.cpp src/utils/*.cpp src/mc/*.cpp -o mc.out
	time ./mc.out $(type) $(N) $(max_iter) $(x) $(mu)

deterministe : 
	@mkdir -p Data
	g++ main_deterministe.cpp src/utils/*.cpp src/deterministe/deterministe.cpp -o deter.out
	time ./deter.out $(Nx) $(Nmu)

# WIP
deterministe_dsa :
	@mkdir -p Data
	g++ main_deterministe_fast.cpp src/utils/*.cpp src/deterministe/deterministe.cpp src/deterministe/fast_IS.cpp -o deter_fast.out
	time ./deter_fast.out $(Nx) $(Nmu)

clean:
	rm *.out