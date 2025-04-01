# Makefile for sim.cpp

# Compiler and flags
CC = g++
CFLAGS = -std=c++11 -Wall -O2

# Target executable
TARGET = sim

# Default target: build the executable
all: $(TARGET)

$(TARGET): sim.cpp
	$(CC) $(CFLAGS) -o $(TARGET) sim.cpp

# Clean target: remove the executable and any generated CSV files
clean:
	rm -f $(TARGET) *.csv

# Experiment target: run the simulation for five scenarios
experiment: $(TARGET)
	@echo "Running Scenario 1: Infinite buffer, 5 nodes, audio reference"
	./$(TARGET) 5 -1 100 audio 0.1 0.9 0.1 results_scenario1.csv
	@echo "Results saved in results_scenario1.csv"
	@echo "-----------------------------------------------"
	@echo "Running Scenario 2: Finite buffer (K = 100), 5 nodes, audio reference"
	./$(TARGET) 5 100 100 audio 0.1 0.9 0.1 results_scenario2.csv
	@echo "Results saved in results_scenario2.csv"
	@echo "-----------------------------------------------"
	@echo "Running Scenario 3: Finite buffer (K = 100), 10 nodes, audio reference"
	./$(TARGET) 10 100 100 audio 0.1 0.9 0.1 results_scenario3.csv
	@echo "Results saved in results_scenario3.csv"
	@echo "-----------------------------------------------"
	@echo "Running Scenario 4: Finite buffer (K = 100), 5 nodes, video reference"
	./$(TARGET) 5 100 100 video 0.1 0.9 0.1 results_scenario4.csv
	@echo "Results saved in results_scenario4.csv"
	@echo "-----------------------------------------------"
	@echo "Running Scenario 5: Finite buffer (K = 100), 5 nodes, data reference"
	./$(TARGET) 5 100 100 data 0.1 0.9 0.1 results_scenario5.csv
	@echo "Results saved in results_scenario5.csv"
	@echo "Experiment complete. Check the generated CSV files for the results."

.PHONY: all clean experiment
