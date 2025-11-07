# Creutz Demon — 2D Ising Model Simulation

**C++ implementation of the Creutz demon algorithm applied to the 2D Ising spin lattice.**  
The project studies the thermodynamic behavior of a 2D Ising system using a **microcanonical Monte Carlo** simulation and analyzes energy histograms, magnetization, and temperature dependence.

📄 **Note:** The full report (*sprawozdanie*) is written in **Polish** and available in the file [`Sprawozdanie.pdf`](./Sprawozdanie.pdf).

---

## 🧠 Project Overview

This project implements the **Creutz demon algorithm** for the **two-dimensional Ising model**.  
It allows simulation of spin interactions on a lattice and analysis of phase transition behavior near the critical temperature.

The algorithm introduces an auxiliary variable — the **demon** — that exchanges energy with the system while keeping total energy constant, enabling a **microcanonical ensemble** simulation.

---

## ⚙️ Implementation Details

- **Language:** C++17  
- **Environment:** CLion 2024.3.4 (JetBrains)  
- **Compiler:** GCC  
- **Key libraries:** `<vector>`, `<map>`, `<random>`, `<cmath>`, `<fstream>`, `<iostream>`, `<filesystem>`  

The program reads input parameters from `stdin` or a file and generates:
- `histogram/` – energy distribution of the demon (`E_demon`, `N(E)`, `ln(N)`)  
- `magnetization/` – time evolution of magnetization during Monte Carlo sweeps  
- `mT.txt` – summary file containing demon energy, temperature, average magnetization ⟨m⟩, and slope of the energy distribution  

---

## 🧩 Input Parameters

Example input format:
Lx Ly numSweeps numEnergies
E_demon_1 E_demon_2 ... E_demon_N

37 37 10000 10
56 128 256 512 768 960 1112 1280 1400 1520

---

## 📊 Output Data

After the simulation, the program generates:

| File                                   | Description                                              |
|----------------------------------------|----------------------------------------------------------|
| `mT.txt`                               | Table of results: demon energy, temperature, ⟨m⟩, slope   |
| `histogram/histogram_E=...txt`         | Histogram of demon energy                                |
| `magnetization/magnetization_E=...txt` | Magnetization as a function of Monte Carlo steps         |

---

## 🧪 Results Summary

- The simulation reproduces the **expected phase transition** in the 2D Ising model.  
- The computed **critical temperature** was  T_c^(sim) ≈ 2.3237 J/kB
- The magnetization ⟨m⟩ decreases sharply near the critical point, showing a clear transition from the ferromagnetic to the paramagnetic phase.

---

## 🧾 Report (in Polish)

The complete analysis, theoretical background, methodology, and discussion of results are presented in the Polish-language report:

📘 **`DemonCreutz.pdf`**

---

## 📂 Repository Structure
CreutzIsing2D/
├── main.cpp # Source code
├── Sprawozdanie.pdf # Full Polish report
├── mT.txt # Summary results file
├── histogram/ # Energy distribution data
├── magnetization/ # Magnetization over MC sweeps
└── README.md # Project documentation

---

## 🧑‍💻 Author
Created for educational purposes by Avui.

