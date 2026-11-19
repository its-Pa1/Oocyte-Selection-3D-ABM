// ****************************************************************************
//
//              Ovarian model GRN (Simuscale-based project)
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package
// Author: Pawan Kumar
// This file was implemented by Pawan Kumar and used to produce results
// reported in the article cited in the main README.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// ****************************************************************************

// ============================================================================
//                                   Includes
// ============================================================================



#ifndef SIMUSCALE_GRN_H__
#define SIMUSCALE_GRN_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <random>
#include <cstring>
#include <memory>
#include <array>

/**
 * @file GRN.h
 * @brief Header file for the Gene Regulatory Network (GRN) simulation using the PDMP framework.
 * 
 * ## Overview of the Gene Regulatory Network (GRN) Simulation
 * This file defines the `GRN` class, which models the dynamics of a gene regulatory network (GRN).
 * The GRN class simulates the behavior of mRNA and protein concentrations in a cell under the 
 * influence of internal interactions (via a gene interaction matrix) and an external signal. 
 * The simulation is based on the **Piecewise Deterministic Markov Process (PDMP)** framework, 
 * which combines deterministic dynamics with stochastic events.
 * 
 * ## Overview of the PDMP Framework
 * The PDMP framework is used to model systems where deterministic dynamics (e.g., ODEs) are 
 * interspersed with stochastic events (e.g., jumps or bursts). This hybrid approach is particularly 
 * suited for simulating gene regulatory networks, where:
 * - Deterministic dynamics describe the continuous evolution of mRNA and protein concentrations.
 * - Stochastic events represent discrete bursts in mRNA production.
 * 
 * ### Key Components of the PDMP Framework in This Implementation:
 * - **Deterministic Evolution**:
 *   - The `Intracellular_ExactEvol` function solves the coupled ODEs for mRNA and protein concentrations.
 *   - These equations describe the decay of mRNA and proteins over time, as well as the production of proteins from mRNA.
 * - **Stochastic Bursts**:
 *   - The `advance` function models stochastic jumps in mRNA levels, where the size of each burst is drawn from an exponential distribution.
 *   - The time between stochastic events is sampled from an exponential distribution with a rate equal to the total propensity (`Tau`).
 * - **Adaptive Thinning**:
 *   - Used to efficiently simulate stochastic events by dynamically adjusting upper bounds (`PPmax`) for propensities.
 * - **Gene Interactions**:
 *   - The `GenesInteractionsMatrix_` defines how genes interact with each other, influencing their activation or inhibition.
 * - **External Signal**:
 *   - An external signal can influence specific genes (e.g., PTN in this case), modifying their activation levels.
 * 
 * ## Key Features of the GRN Class
 * - **mRNA and Protein Dynamics**:
 *   - Simulates mRNA and protein concentrations for each gene in the network.
 *   - Includes gene-specific decay rates for mRNA (`d0`) and proteins (`d1`).
 * - **Gene Interactions**:
 *   - Models interactions between genes using a gene interaction matrix (`GenesInteractionsMatrix_`).
 *   - Supports both activation and inhibition effects.
 * - **Stochastic Bursting**:
 *   - Models stochastic bursts in mRNA production, where the burst size is drawn from an exponential distribution.
 * - **External Signal**:
 *   - Allows an external signal to influence specific genes, modifying their activation levels dynamically.
 * - **Time Evolution**:
 *   - Supports both deterministic evolution (via ODEs) and stochastic evolution (via bursts).
 * - **Scalability**:
 *   - Designed for large-scale simulations of cellular populations, where each cell can have its own GRN instance.
 * 
 * ## Applications
 * This framework is particularly suited for:
 * - Simulating gene regulatory networks in individual cells.
 * - Modeling large-scale cellular populations in the context of SimuScale.
 * - Studying the effects of gene interactions and external signals on cellular behavior.
 */
/**
 * @class GRN
 * @brief Class for simulating a Gene Regulatory Network (GRN).
 * 
 * The `GRN` class models the dynamics of a gene regulatory network, including:
 * - mRNA and protein concentrations for each gene.
 * - Gene-specific decay rates for mRNA and proteins.
 * - Interactions between genes through a gene interaction matrix.
 * - The effect of an external signal on specific genes.
 * - Stochastic bursting events for mRNA production.
 * 
 * The class provides methods for initializing the GRN, simulating its dynamics over time, and
 * handling stochastic events such as mRNA bursts. It is designed to be flexible and efficient,
 * allowing it to be applied to individual cells in a large-scale simulation framework.
 */

 struct GeneConstants {
    double K0;  // Precomputed K0 for the gene
    double K1;  // Precomputed K1 for the gene
    double d0;  // Decay rate (d0) for the gene
};

class GRN {
public:
    // Constructors and destructor declarations
    /**
     * @brief Constructor for the GRN class.
     * 
     * Initializes the GRN with the given configuration, kinetics, and interactions files.
     * 
     * @param config_file Path to the configuration file.
     * @param kinetics_file Path to the kinetics file.
     * @param interactions_file Path to the interactions file.
     */

    GRN(const std::string& config_file, const std::string& kinetics_file, const std::string& interactions_file);
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Constructor for the GRN class with an external signal.
     * 
     * Initializes the GRN with the given configuration, kinetics, and interactions files,
     * and sets an initial external signal.
     * 
     * @param config_file Path to the configuration file.
     * @param kinetics_file Path to the kinetics file.
     * @param interactions_file Path to the interactions file.
     * @param external_signal Initial value of the external signal.
     */

    GRN(const std::string& config_file, const std::string& kinetics_file, const std::string& interactions_file, double external_signal);
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Destructor for the GRN class.
     * 
     * Frees dynamically allocated memory used by the GRN.
     */
    virtual ~GRN();
    //----------------------------------------------------------------------------------------------------------------------------------
    // Public method declarations

    // Interface methods to access GRN state
    /**
     * @brief Get the number of genes in the GRN.
     * 
     * @return int Number of genes in the GRN.
     */
    int getNumberOfGenes() const { return Number_Of_Genes_; }
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Get the number of parameters in the GRN.
     * 
     * @return int Number of parameters in the GRN.
     */
    double getProteinLevel(int i) const { return Protein_array_[i]; }
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Get the mRNA level for a specific gene.
     * 
     * @param i Index of the gene.
     * @return double mRNA level for the gene.
     */
    double getMRNALevel(int i) const { return mRNA_array_[i]; }
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Set the mRNA level for a specific gene.
     * 
     * @param i Index of the gene.
     * @param value New mRNA level for the gene.
     */
    void setMRNALevel(int i, double value) { mRNA_array_[i] = value; }
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Set the protein level for a specific gene.
     * 
     * @param i Index of the gene.
     * @param value New protein level for the gene.
     */
    void setProteinLevel(int i, double value) { Protein_array_[i] = value; }

    //----------------------------------------------------------------------------------------------------------------------------------

    /**
     * @brief Get the physical protein concentration for a specific gene.
     *
     * @param gene_index Index of the gene.
     * @return double Physical protein concentration for the gene.
     */
    double getPhysicalProteinConcentration(int gene_index) const {
    if (gene_index < 0 || gene_index >= Number_Of_Genes_) return 0.0;
    return (cell_volume_ > 0) ? (Protein_array_[gene_index] / cell_volume_) : 0.0;}
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Get the physical mRNA concentration for a specific gene.
     *
     * @param gene_index Index of the gene.
     * @return double Physical mRNA concentration for the gene.
     */
    double getPhysicalMRNAConcentration(int gene_index) const {
    if (gene_index < 0 || gene_index >= Number_Of_Genes_) return 0.0;
    return (cell_volume_ > 0) ? (mRNA_array_[gene_index] / cell_volume_) : 0.0;}
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Set the external signal value.
     * 
     * @param external_signal New value of the external signal.
     */
    void setExternalSignal(double external_signal) { external_signal_ = external_signal; }
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Get the external signal value.
     * 
     * @return double Current value of the external signal.
     */
    double getExternalSignal() const { return external_signal_; }
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Get the time of the next jump event.
     * 
     * @return double Time of the next jump event.
     */
    double getNextJumpTime() const { return Time_NextJump_;}
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Set the time of the next jump event.
     * 
     * @param time New time for the next jump event.
     */
    void setNextJumpTime(double time) { Time_NextJump_ = time;}
    //----------------------------------------------------------------------------------------------------------------------------------
    // time marching method
    /**
     * @brief Advances the GRN state by a specified time step using the PDMP framework.
     * 
     * ## PDMP Workflow
     * 1. **Deterministic Evolution**:
     *    - The system evolves deterministically up to the next stochastic event or the end of the time step.
     *    - This is handled by the `Intracellular_ExactEvol` function, which solves the ODEs for mRNA and protein concentrations.
     * 
     * 2. **Calculate Propensities**:
     *    - The total propensity (`Tau`) is calculated as the sum of individual propensities for all genes.
     *    - Adaptive thinning is used to calculate an upper bound (`PPmax`) for protein levels.
     * 
     * 3. **Sample Time to Next Jump**:
     *    - The time to the next stochastic event (`DeltaT`) is sampled from an exponential distribution with rate `1 / Tau`.
     * 
     * 4. **Stochastic Jump**:
     *    - A gene is selected for the stochastic jump based on the probabilities of individual propensities.
     *    - The mRNA level for the selected gene is updated with a burst size drawn from an exponential distribution.
     * 
     * 5. **Repeat Until End of Time Step**:
     *    - The process repeats until the total time step (`dt`) is completed.
     * 
     * ## Parameters
     * @param dt Time step for the simulation.
     * @param global_time Current global time in the simulation.
     * @param volume_ratio Ratio of the new cell volume to the old cell volume.
     */
    void advance(const double& dt, const double& global_time, double volume_ratio = 1.0);
        //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Set the cell volume for concentration calculations.
     * 
     * This method updates the cell volume, which is used to convert between
     * molecular amounts and concentrations.
     * 
     * @param volume New cell volume (in µm³ / unit³).
     */
    void setCellVolume(double volume) { 
        cell_volume_ = (volume > 0) ? volume : 1.0;  // Ensure volume is always positive
}
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Get the current cell volume.
     * 
     * @return double Current cell volume (in µm³ / unit³).
     */
    double getCellVolume() const { return cell_volume_; }
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Get molecular amount for mRNA from concentration.
     * 
     * Converts the stored mRNA concentration to absolute amount by multiplying
     * with current cell volume.
     * 
     * Formula: amount = concentration × volume
     * 
     * @param gene_index Index of the gene (0, 1, ..., Number_Of_Genes_-1).
     * @return double Absolute mRNA amount (concentration × volume).
     */
    double getMRNAAmount(int gene_index) const {
        // if (gene_index < 0 || gene_index >= Number_Of_Genes_) {
        //     std::cerr << "Error: Invalid gene index " << gene_index << std::endl;
        //     return 0.0;
        // }
        return mRNA_array_[gene_index] * cell_volume_;
    }
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Get molecular amount for protein from concentration.
     * 
     * Converts the stored protein concentration to absolute amount by multiplying
     * with current cell volume.
     * 
     * Formula: amount = concentration × volume
     * 
     * @param gene_index Index of the gene (0, 1, ..., Number_Of_Genes_-1).
     * @return double Absolute protein amount (concentration × volume).
     */
    double getProteinAmount(int gene_index) const {
        // if (gene_index < 0 || gene_index >= Number_Of_Genes_) {
        //     std::cerr << "Error: Invalid gene index " << gene_index << std::endl;
        //     return 0.0;
        // }
        return Protein_array_[gene_index] * cell_volume_;
    }
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Set mRNA amount (converts to concentration and stores).
     * 
     * Converts the given absolute mRNA amount to concentration by dividing
     * by current cell volume, then stores in GRN.
     * 
     * Formula: concentration = amount / volume
     * 
     * @param gene_index Index of the gene.
     * @param amount Absolute mRNA amount to set.
     */
    void setMRNAAmount(int gene_index, double amount) {
        // if (gene_index < 0 || gene_index >= Number_Of_Genes_) {
        //     std::cerr << "Error: Invalid gene index " << gene_index << std::endl;
        //     return;
        // }
        double concentration = (cell_volume_ > 0) ? (amount / cell_volume_) : 0.0;
        mRNA_array_[gene_index] = concentration;
    }
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Set protein amount (converts to concentration and stores).
     * 
     * Converts the given absolute protein amount to concentration by dividing
     * by current cell volume, then stores in GRN.
     * 
     * Formula: concentration = amount / volume
     * 
     * @param gene_index Index of the gene.
     * @param amount Absolute protein amount to set.
     */
    void setProteinAmount(int gene_index, double amount) {
        // if (gene_index < 0 || gene_index >= Number_Of_Genes_) {
        //     std::cerr << "Error: Invalid gene index " << gene_index << std::endl;
        //     return;
        // }
        double concentration = (cell_volume_ > 0) ? (amount / cell_volume_) : 0.0;
        Protein_array_[gene_index] = concentration;
    }
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Update GRN concentrations after a pure volume change (no molecular transfer).
     * 
     * When cell volume changes but molecular amounts stay constant (e.g., Strategy 1
     * in mass transfer), concentrations must be adjusted to maintain conservation.
     * 
     * Conservation principle:
     *   amount = old_conc × old_volume = new_conc × new_volume
     * 
     * Therefore:
     *   new_conc = old_conc × (old_volume / new_volume)
     * 
     * This method updates all mRNA and protein concentrations in the GRN.
     * 
     * @param old_volume Cell volume before the change (µm³).
     * @param new_volume Cell volume after the change (µm³).
     */
    void updateConcentrationsForVolumeChange(double old_volume, double new_volume) ;
    //----------------------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------------------
    // ============================================================================
    //                      Declaration of variables and methods
    // ============================================================================
    // Protected member variables
    int Number_Of_Genes_; /**< Number of genes in the GRN. */
    int Number_Of_Parameters_; /**< Number of parameters in the GRN. */
    double External_Signal_Coeff_; /**< Coefficient for the external signal in the sigma equation. */
    double Basal_activity_Phi_; /**< Basal activity (Phi) for the sigma equation. */
    std::vector<double> KinParam_; /**< Array of kinetic parameters. */
    std::vector<double> GenesInteractionsMatrix_; /**< Matrix of gene interactions. */
    std::vector<double> mRNA_array_; /**< Array of mRNA concentrations. */
    std::vector<double> Protein_array_; /**< Array of protein concentrations. */
    std::mt19937 gen_; /**< Random number generator. */
    std::uniform_real_distribution<> dis_; /**< Uniform distribution for random numbers. */
    int ithGene_; /**< Index of the current gene being processed. */
    double Time_NextJump_; /**< Time of the next jump event. */
    std::vector<int> TrueJumpCounts_array_; /**< Array of true jump counts for each gene. */
    int PhantomJumpCounts_; /**< Count of phantom jumps. */
    double external_signal_; /**< Current value of the external signal. */
    double cell_volume_;  /**< Current cell volume for concentration calculations. */

    std::vector<GeneConstants> gene_constants_; /**< Precomputed constants for each gene. */
    double d1_; /**< Common decay rate for proteins (d1). */
    double a2_; /**< Shared parameter (a2). it is k_off/S_{0,1} as in the Herbach et al. 2021 */
    //----------------------------------------------------------------------------------------------------------------------------------
    //Protected method declarations

    /**
     * @brief Read the configuration file.
     * 
     * Reads the configuration file to determine the number of genes and parameters.
     * 
     * ### Configuration File Format (`config_file`)
     * - Specifies the number of genes, number of parameters, and the external signal coefficient.
     * - Example:
     *   ```
     *   Number_Of_Genes: 2
     *   Number_Of_Parameters: 6
     *   External_Signal_Coeff: 1.0
     *   ```
     * 
     * @param config_file Path to the configuration file.
     */
    void readConfig(const std::string& config_file);
    //----------------------------------------------------------------------------------------------------------------------------------
    
    /**
     * @brief Reads the kinetics and interactions files and initializes parameters.
     * 
     * ### Kinetics File Format (`kinetics_file`)
     * - Defines the kinetic parameters for mRNA and protein dynamics.
     * - Example:
     *   ```
     *   0.1  # d0_1: mRNA decay rate for gene 1
     *   0.05 # d1: Protein decay rate (shared by all genes)
     *   0.0  # a0: Basal transcription rate
     *   10.0 # a1: Transcriptional activation rate
     *   5.0  # a2: Burst size
     *   0.11 # d0_2: mRNA decay rate for gene 2
     *   ```
     * 
     * ### Interactions File Format (`interactions_file`)
     * - Defines the interaction matrix between genes.
     * - Example:
     *   ```
     *   0.0  1.0
     *   -1.0 0.0
     *   ```
     *   - Gene 1 activates Gene 0, no self-activation of Gene 0.
     *   - Gene 0 inhibits Gene 1, no self-activation of Gene 1.
     * 
     * @param kinetics_file Path to the kinetics file.
     * @param interactions_file Path to the interactions file.
     */
    void get_GeneParams(const std::string& kinetics_file, const std::string& interactions_file);
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Initializes the arrays for mRNA, protein, and jump counts.
     */
    void initializeArrays();   
    //----------------------------------------------------------------------------------------------------------------------------------

    /**
     * @brief Initializes the gene constants based on the kinetics parameters.
     */
    void initializeGeneConstants();
    //----------------------------------------------------------------------------------------------------------------------------------
   /**
     * @brief Calculates the sigma value for a specific gene.
     * 
     * Sigma represents the activation level of a gene based on:
     * - Interactions with other genes (from `GenesInteractionsMatrix_`).
     * - The influence of an external signal (if applicable).
     * 
     * ## Formula
     * \f[
     * \text{Sigma} = \frac{1}{1 + e^{-\text{initval}}}
     * \f]
     * where:
     * - `initval` is the sum of contributions from other genes and the external signal.
     * 
     * ## Parameters
     * @param i Index of the gene.
     * @param P Protein array.
     * @param AcceptNegative Flag to accept negative interaction values.
     * @param external_signal External signal value.
     * @param volume_ratio Volume ratio for the gene.
     * @return double Sigma value for the gene.
     */
    double GetSigma_i(int i, std::vector<double>& P, bool AcceptNegative = false, double external_signal = 0.0, double volume_ratio = 1.0);
    //----------------------------------------------------------------------------------------------------------------------------------
    /**
     * @brief Exact evolution of the intracellular system for mRNA and protein concentrations.
     * 
     * This function solves the coupled ODEs for mRNA and protein concentrations using explicit solutions (eq (6)-(7) of in Herbach et al. 2021):
     * 
     * \f[
     * M(t) = M(0) \cdot e^{-d_0 \cdot t}
     * \f]
     * 
     * \f[
     * P(t) = \frac{S_1}{d_0 - d_1} \cdot M(0) \cdot (e^{-d_1 \cdot t} - e^{-d_0 \cdot t}) + P(0) \cdot e^{-d_1 \cdot t}
     * \f]
     *  
     * ## Parameters
     * @param t Time step for the evolution.
     * @param P Protein array (updated in place).
     * @param M mRNA array (updated in place).
     */
    void Intracellular_ExactEvol(double thetime, std::vector<double>& P, std::vector<double>& M); // or a ODE solver e.g. gls ones 
    //----------------------------------------------------------------------------------------------------------------------------------

private:
    // Private helper method declarations

    /**
     * @brief Calculates the value of S1 for a specific gene.
     * 
     * S1 represents the rate of protein production based on the mRNA decay rate (`d0`), 
     * protein decay rate (`d1`), and the burst size parameter (`a2`).
     * 
     * ## Formula
     * \[
     * S1 = \frac{d_0 \cdot d_1 \cdot a_2}{K_1}
     * \]
     * 
     * ## Parameters
     * @param gene_index Index of the gene.
     * @return double Value of S1.
     */
    double calculateS1(int gene_index);
    //----------------------------------------------------------------------------------------------------------------------------------
};

// ============================================================================
//                                Constructors  
// ============================================================================
GRN::GRN(const std::string& config_file, const std::string& kinetics_file, const std::string& interactions_file)
    : gen_(std::random_device()()), dis_(0.01, 1.0), ithGene_(-1), Time_NextJump_(0.0), external_signal_(0.0), cell_volume_(1.0) {
    readConfig(config_file);
    get_GeneParams(kinetics_file, interactions_file);
    initializeArrays();
    initializeGeneConstants();
}

GRN::GRN(const std::string& config_file, const std::string& kinetics_file, const std::string& interactions_file, double external_signal)
    : gen_(std::random_device()()), dis_(0.01, 1.0), ithGene_(-1), Time_NextJump_(0.0), external_signal_(external_signal), cell_volume_(1.0) {
    readConfig(config_file);
    get_GeneParams(kinetics_file, interactions_file);
    initializeArrays();
    initializeGeneConstants();
}
// ============================================================================
//                                Destructor
// ============================================================================
GRN::~GRN() {
}
// ============================================================================
//                                  Methods
// ============================================================================



// Read the configuration file
void GRN::readConfig(const std::string& config_file) {
    std::ifstream file(config_file);
    if (!file) {
        std::cerr << "Error: file " << config_file << " could not be opened" << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key;
        if (std::getline(iss, key, ':')) {
            std::string value;
            if (std::getline(iss, value)) {
                if (key == "Number_Of_Genes") {
                    Number_Of_Genes_ = std::stoi(value);
                } else if (key == "Number_Of_Parameters") {
                    Number_Of_Parameters_ = std::stoi(value);
                } else if (key == "External_Signal_Coeff") {
                    External_Signal_Coeff_ = std::stod(value);
                } else if (key == "Basal_activity_Phi"){
                    Basal_activity_Phi_ = std::stod(value);
                }
            }
        }
    }

    file.close();
}
//----------------------------------------------------------------------------------------------------------------------------------
// Initialize gene constants based on the kinetics parameters
// This function precomputes mRNA decay rates (d0), K0 and K1 for each gene based on the given kinetics parameters
void GRN::initializeGeneConstants() {
    d1_ = KinParam_[1];
    a2_ = KinParam_[4];

    for (int i = 0; i < Number_Of_Genes_; i++) {
        double d0 = (i % 2 == 0) ? KinParam_[0] : KinParam_[5]; // Use d0_1 for even indices and d0_2 for odd indices, first is d0_1 for KDCH3 and d0_2 for PTN
        double a0 = KinParam_[2];
        double a1 = KinParam_[3];

        GeneConstants constants;
        constants.d0 = d0;
        constants.K0 = a0 * d0;
        constants.K1 = a1 * d0;

        gene_constants_.push_back(constants);
    }
}
//----------------------------------------------------------------------------------------------------------------------------------
void GRN::initializeArrays() {
    mRNA_array_.resize(Number_Of_Genes_, 0.0);
    Protein_array_.resize(Number_Of_Genes_, 0.0);
    TrueJumpCounts_array_.resize(Number_Of_Genes_, 0);

    // fallback initialization:
    double init_mrna_0 = 2.0;
    double init_protein_0 = 0.3;
    double init_mrna_1 = 2.0;
    double init_protein_1 = 0.3;

    const std::string init_file = "full_trajectories_signal_0.0_5000h_pro_oocyte.txt";
    std::ifstream init_stream(init_file);
    if (init_stream) {
        std::string line;
        if (std::getline(init_stream, line)) {
            std::vector<std::array<double, 5>> candidates;
            while (std::getline(init_stream, line)) {
                if (line.empty() || line[0] == '#') {
                    continue;
                }
                std::istringstream iss(line);
                double time = 0.0;
                double m0_mean = 0.0;
                double m0_std = 0.0;
                double p0_mean = 0.0;
                double p0_std = 0.0;
                double m1_mean = 0.0;
                double m1_std = 0.0;
                double p1_mean = 0.0;
                double p1_std = 0.0;
                if (!(iss >> time >> m0_mean >> m0_std >> p0_mean >> p0_std >> m1_mean >> m1_std >> p1_mean >> p1_std)) {
                    continue;
                }
                if (time >= 1000.0) {
                    candidates.push_back({time, m0_mean, p0_mean, m1_mean, p1_mean});
                }
            }

            if (!candidates.empty()) {
                std::uniform_int_distribution<size_t> pick(0, candidates.size() - 1);
                const auto& chosen = candidates[pick(gen_)];
                init_mrna_0 = chosen[1];
                init_protein_0 = chosen[2];
                init_mrna_1 = chosen[3];
                init_protein_1 = chosen[4];
            }
        }
    }

    for (int j = 0; j < Number_Of_Genes_; j++) {
        if (j == 0) {
            // Initialize at the pro-oocyte gene.
            Protein_array_[j] = init_protein_0;
            mRNA_array_[j] = init_mrna_0;
        } else if (j == 1) {
            // Initialize at the pro nurse gene.
            Protein_array_[j] = init_protein_1;
            mRNA_array_[j] = init_mrna_1;
        } else {
            Protein_array_[j] = 1.0 * dis_(gen_);
            mRNA_array_[j] = 1.0 * dis_(gen_);
        }
        TrueJumpCounts_array_[j] = 0;
    }
    PhantomJumpCounts_ = 0;
    Time_NextJump_ = 0.0;
    ithGene_ = -1;
}
//----------------------------------------------------------------------------------------------------------------------------------
void GRN::get_GeneParams(const std::string& kinetics_file, const std::string& interactions_file) {
    // Validate parameters
    if (Number_Of_Parameters_ <= 0 || Number_Of_Genes_ <= 0) {
        std::cerr << "Error: Invalid parameter values" << std::endl;
        exit(1);
    }

    // Open kinetics file
    std::ifstream kin_file(kinetics_file);
    if (!kin_file) {
        std::cerr << "Error: file " << kinetics_file << " could not be opened" << std::endl;
        exit(1);
    }

    // Resize KinParam_ vector
    KinParam_.resize(Number_Of_Parameters_, 0.0);

    std::string line;
    int param_index = 0;
    while (std::getline(kin_file, line) && param_index < Number_Of_Parameters_) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip empty lines and comments
        }
        std::istringstream iss(line);
        while (param_index < Number_Of_Parameters_ && iss >> KinParam_[param_index]) {
            param_index++;
        }
    }
    kin_file.close();

    // Open interactions file
    std::ifstream inter_file(interactions_file);
    if (!inter_file) {
        std::cerr << "Error: file " << interactions_file << " could not be opened" << std::endl;
        exit(1);
    }

    // Resize GenesInteractionsMatrix_ vector
    GenesInteractionsMatrix_.resize(Number_Of_Genes_ * Number_Of_Genes_, 0.0);

    int row = 0, col = 0;
    while (std::getline(inter_file, line)) {
        if (line.empty() || line[0] == '#') continue; // Skip comments & empty lines

        std::istringstream iss(line);
        while (iss >> GenesInteractionsMatrix_[row * Number_Of_Genes_ + col]) {
            col++;
            if (col == Number_Of_Genes_) { // Move to next row
                col = 0;
                row++;
                if (row == Number_Of_Genes_) break; // Stop when matrix is filled
            }
        }
    }
    inter_file.close();

    // Check if the matrix was fully read
    if (row < Number_Of_Genes_ || col != 0) {
        std::cerr << "Error: Unexpected end of file or invalid data in " << interactions_file << std::endl;
        exit(1);
    }
}
//----------------------------------------------------------------------------------------------------------------------------------
// Get the sigma value for gene i
double GRN::GetSigma_i(int i, std::vector<double>& P, bool AcceptNegative, double external_signal, double volume_ratio) {
    double initval = Basal_activity_Phi_;  // Basal activity for other genes -0.5
    for (int j = 0; j < Number_Of_Genes_; j++) {
        double interIJ = GenesInteractionsMatrix_[i * Number_Of_Genes_ + j]; // J acts on I
        if ((interIJ > 0) || AcceptNegative)  {
            initval += P[j] * volume_ratio * interIJ;
        }
    }
    // Apply external signal only to gene index 1 (i.e. second gene PTN)
    if (i == 1) {
        initval += external_signal * External_Signal_Coeff_;  // plus external signal
    }
    return (1 / (1 + exp(-initval))); //   check it..............................ß

   
}
//----------------------------------------------------------------------------------------------------------------------------------
// Exact evolution of the intracellular system, think better for a complex system of ODEs
void GRN::Intracellular_ExactEvol(double t, std::vector<double>& P, std::vector<double>& M) {
    // Create copies of the current mRNA and protein arrays
    std::vector<double> M_Old = M; // Copy current mRNA array
    std::vector<double> P_Old = P; // Copy current protein array
    
    // Explicit ODEs solutions for each gene
    for (int i = 0; i < Number_Of_Genes_; i++) {
        // Use gene_constants_ for gene-specific decay rate (d0)
        double d0 = gene_constants_[i].d0;
        
        // Calculate S1 using the helper method (which already uses gene_constants_)
        double S1 = calculateS1(i);

        // Update mRNA - simple exponential decay
        M[i] = M_Old[i] * exp(-t * d0); 
        
        // Update protein - analytical solution of the coupled ODEs
        P[i] = ((S1 / (d0 - d1_)) * M_Old[i] * (exp(-t * d1_) - exp(-t * d0))
                + P_Old[i] * exp(-t * d1_));
    }
}
//----------------------------------------------------------------------------------------------------------------------------------
// Advance the GRN state by a specified time step
void GRN::advance(const double& dt, const double& global_time, double volume_ratio) {
    // Calculate the end time for this time step
    double end_time = global_time + dt;
    
    // Use modern C++ containers for memory management
    std::vector<double> PPmax(Number_Of_Genes_, 0.0);
    std::vector<double> ProbaArray(Number_Of_Genes_ + 1, 0.0);
    std::vector<double> KKon(Number_Of_Genes_, 0.0);
    std::vector<double> Sigma(Number_Of_Genes_, 0.0);
    
    // Initialize Time_NextJump_ if necessary (first call or we've reached the previous jump)
    if (Time_NextJump_ <= global_time) {
        // Calculate initial Tau for determining the next jump time
        double Tau = 0.0;
        for (int i = 0; i < Number_Of_Genes_; i++) {
            // Calculate S1 for the current gene
            double S1 = calculateS1(i);

            // Calculate PPmax for the current gene [Adaptive thinning]
            double thetime = log(gene_constants_[i].d0 / d1_) / (gene_constants_[i].d0 - d1_ + 1e-10);
            PPmax[i] = Protein_array_[i] + (S1 / (gene_constants_[i].d0 - d1_ + 1e-10)) * 
                       mRNA_array_[i] * (exp(-thetime * d1_) - exp(-thetime * gene_constants_[i].d0));
            
            // Calculate Sigma for the upper bound (AcceptNegative=false)
            Sigma[i] = GetSigma_i(i, PPmax, true, external_signal_, volume_ratio);
            
            // Calculate KKon for the current gene
            KKon[i] = (1 - Sigma[i]) * gene_constants_[i].K0 + Sigma[i] * gene_constants_[i].K1 + exp(-10 * log(10));
            Tau += KKon[i];
        }
        
        // Calculate time to next jump
        double DeltaT = std::exponential_distribution<double>(Tau)(gen_);
        Time_NextJump_ = global_time + DeltaT;
    }
    
    // Current position in time (starts at global_time)
    double current_time = global_time;
    
    // Core simulation loop - process jumps until beyond the end time
    while (Time_NextJump_ < end_time) {
        // Evolve deterministically up to the next jump time
        double time_to_evolve = Time_NextJump_ - current_time;
        if (time_to_evolve > 0) {
            Intracellular_ExactEvol(time_to_evolve, Protein_array_, mRNA_array_);
        }
        
        // Update current time to the jump time
        current_time = Time_NextJump_;
        
        // Calculate PPmax and Tau for determining the next jump
        double Tau = 0.0;
        for (int i = 0; i < Number_Of_Genes_; i++) {
            // Calculate S1 for the current gene
            double S1 = calculateS1(i);

            // Calculate PPmax for the current gene [Adaptive thinning]
            double thetime = log(gene_constants_[i].d0 / d1_) / (gene_constants_[i].d0 - d1_ + 1e-10);
            PPmax[i] = Protein_array_[i] + (S1 / (gene_constants_[i].d0 - d1_ + 1e-10)) * 
                       mRNA_array_[i] * (exp(-thetime * d1_) - exp(-thetime * gene_constants_[i].d0));
            
            // Calculate Sigma for the Kon upper bound (AcceptNegative=false)
            Sigma[i] = GetSigma_i(i, PPmax, true, external_signal_, volume_ratio);
            
            // Calculate KKon for the current gene
            KKon[i] = (1 - Sigma[i]) * gene_constants_[i].K0 + Sigma[i] * gene_constants_[i].K1 + exp(-10 * log(10));
            Tau += KKon[i];
        }
        
        // Calculate probabilities for each gene
        double Proba_no_jump = 1.0;
        for (int i = 0; i < Number_Of_Genes_; i++) {
            // Use actual Sigma values (AcceptNegative=true) for dynamics
            Sigma[i] = GetSigma_i(i, Protein_array_, true, external_signal_, volume_ratio);
            ProbaArray[i] = ((1 - Sigma[i]) * gene_constants_[i].K0 + Sigma[i] * gene_constants_[i].K1) / Tau;
            Proba_no_jump -= ProbaArray[i];
        }
        
        // Add the probability of no jump
        ProbaArray[Number_Of_Genes_] = std::max(0.0, Proba_no_jump);
        
        // Normalize ProbaArray to ensure it sums to 1
        double sum = 0.0;
        for (int i = 0; i < Number_Of_Genes_ + 1; i++) {
            sum += ProbaArray[i];
        }
        if (sum > 0) { // Avoid division by zero
            for (int i = 0; i < Number_Of_Genes_ + 1; i++) {
                ProbaArray[i] /= sum;
            }
        }
        
        // Select the next gene to jump
        std::discrete_distribution<int> DiscreteJumpDistrib(ProbaArray.begin(), ProbaArray.end());
        ithGene_ = DiscreteJumpDistrib(gen_);
        
        // Apply the stochastic burst if a true jump occurs
        if (ithGene_ < Number_Of_Genes_ && ithGene_ >= 0) {
            TrueJumpCounts_array_[ithGene_]++;
            double burst_amount_ = std::exponential_distribution<double>(a2_)(gen_);
            mRNA_array_[ithGene_] += burst_amount_/volume_ratio; // Adjust burst amount based on volume change
        } else {
            PhantomJumpCounts_++;
        }
        
        // Calculate time to next jump
        double DeltaT = std::exponential_distribution<double>(Tau)(gen_);
        
        // Update the time of next jump
        Time_NextJump_ += DeltaT;
    }
    
    // Evolve deterministically up to the end time if needed
    if (current_time < end_time) {
        Intracellular_ExactEvol(end_time - current_time, Protein_array_, mRNA_array_);
    }
}
//----------------------------------------------------------------------------------------------------------------------------------
double GRN::calculateS1(int gene_index) {
    const GeneConstants& constants = gene_constants_[gene_index];
    return constants.d0 * d1_ * a2_ / constants.K1;
}
//----------------------------------------------------------------------------------------------------------------------------------
void GRN::updateConcentrationsForVolumeChange(double old_volume, double new_volume) {
    if (old_volume <= 0 || new_volume <= 0) {
        std::cerr << "Error: Volumes must be positive values." << std::endl;
        return;
    }

    double volume_ratio = old_volume / new_volume;

    // Update mRNA concentrations
    for (int i = 0; i < Number_Of_Genes_; i++) {
        mRNA_array_[i] *= volume_ratio;
    }

    // Update protein concentrations
    for (int i = 0; i < Number_Of_Genes_; i++) {
        Protein_array_[i] *= volume_ratio;
    }
}

#endif // SIMUSCALE_GRN_H__