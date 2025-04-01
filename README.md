# Packet Scheduling and Traffic Flow Analysis

This repository contains a simulation program (`sim.cpp`) designed to study packet scheduling algorithms and traffic flow behavior under different network load scenarios. The project generates CSV output files that capture performance metrics, which can then be analyzed or visualized (e.g., in a Jupyter notebook).

## Table of Contents

- [Packet Scheduling and Traffic Flow Analysis](#packet-scheduling-and-traffic-flow-analysis)
  - [Table of Contents](#table-of-contents)
  - [Overview](#overview)
  - [Directory Structure](#directory-structure)
  - [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Building the Simulation](#building-the-simulation)
  - [Running Experiments](#running-experiments)
  - [Results](#results)
  - [Visualization](#visualization)
  - [License](#license)
  - [Acknowledgments](#acknowledgments)

---

## Overview

The primary goal of this project is to simulate how different packet scheduling strategies affect network performance. By varying input traffic types (e.g., audio, video, data) and packet sizes, we can observe how queues behave under various offered loads. The simulation outputs performance metrics such as:

- End-to-end delay  
- Packet blocking ratio  
- Throughput  
- Queue utilization  

These metrics are stored in CSV files for further analysis.

---

## Directory Structure

```
.
├── Makefile
├── results_scenario1.csv
├── results_scenario2.csv
├── results_scenario3.csv
├── results_scenario4.csv
├── results_scenario5.csv
├── sim.cpp
├── visual.ipynb
└── README.md
```

- **Makefile**: Contains rules to build and run the simulation.  
- **sim.cpp**: Main C++ simulation code.  
- **results_scenarioX.csv**: Output CSV files containing experiment results for each scenario.  
- **visual.ipynb**: Jupyter notebook for analyzing and visualizing the results.  
- **README.md**: This documentation file.

---

## Getting Started

### Prerequisites

- A C++ compiler (e.g., `g++`, `clang++`) supporting C++11 or later.
- [GNU Make](https://www.gnu.org/software/make/).
- (Optional) [Python 3](https://www.python.org/) and [Jupyter](https://jupyter.org/) if you want to use `visual.ipynb` for visualization.

### Building the Simulation

1. **Clone or download** this repository.
2. **Navigate** to the project root directory in your terminal.
3. **Run** `make experiment` (or just `make`) to build and execute the simulation.

   ```bash
   make experiment
   ```

   This will:
   - Compile `sim.cpp` into an executable (named `sim` by default).
   - Run the simulation for each scenario.
   - Generate CSV files in the project directory.

---

## Running Experiments

- **Default Experiments**: The command `make experiment` executes predefined scenarios (e.g., scenario1, scenario2, etc.) that vary traffic types and parameters. The results are stored in the respective `results_scenarioX.csv` files.
- **Custom Experiments**: To modify or add scenarios, update the `Makefile` and/or the parameters in `sim.cpp` (e.g., packet sizes, offered load range, scheduling algorithms, etc.) and then rerun `make experiment`.

---

## Results

After running the simulation, you will find CSV files in the project directory, such as:

- `results_scenario1.csv`
- `results_scenario2.csv`
- `results_scenario3.csv`
- `results_scenario4.csv`
- `results_scenario5.csv`

Each file typically contains columns like:
- `Offered_Load`
- `AvgDelay`
- `BlockingRatio`
- `AvgBacklog`
- `N_audio`
- `N_video`
- `N_data`


---

## Visualization

Use the included Jupyter notebook (`visual.ipynb`) to analyze and plot the data:

1. **Install** the required Python libraries (e.g., `pandas`, `matplotlib`, `numpy`):
   ```bash
   pip install pandas matplotlib numpy
   ```
2. **Open** the notebook:
   ```bash
   jupyter notebook visual.ipynb
   ```
3. **Run** the cells to generate plots and statistics from the CSV files.


