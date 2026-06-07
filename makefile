CXX=g++
CXXFLAGS=-std=c++17 -O3 -fopenmp -Wall
GEN_SRC=gen.cpp
GEN_EXE=gen
MESH=ball_mesh.txt
GEN_ARGS=--radius 1.0 --points 820 --min-distance 0.058 --seed 3 --output $(MESH)
EXP_DEMO=explicit.py
IMP_DEMO=implicit.py
EXPERIMENT=experiment.py

.PHONY:clean tet explicit

clean:
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	find . -type f \( -name "*.pyc" -o -name "*.pyo" \) -delete
	rm -f $(GEN_EXE) imgui.ini data/experiment.csv

clear:
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	find . -type f \( -name "*.pyc" -o -name "*.pyo" \) -delete
	rm -f $(GEN_EXE) $(MESH) imgui.ini data/experiment.csv

tet:
	$(CXX) $(CXXFLAGS) "$(GEN_SRC)" -o "$(GEN_EXE)"
	"./$(GEN_EXE)" $(GEN_ARGS)

explicit:
	/home/jerry/miniconda3/envs/graphic/bin/python $(EXP_DEMO)

implicit:
	/home/jerry/miniconda3/envs/graphic/bin/python $(IMP_DEMO)

experiment:
	/home/jerry/miniconda3/envs/graphic/bin/python $(EXPERIMENT)