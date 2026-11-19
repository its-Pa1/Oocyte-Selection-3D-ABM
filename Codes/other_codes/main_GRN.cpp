/**
 * @file main.cpp
 * @brief Main file for simulating the Gene Regulatory Network (GRN).
 * 
 * This file contains the main function that simulates the GRN using the `advance` function
 * and outputs the results to a file. It handles the introduction of an external signal
 * at a specified time and synchronizes the state of the GRN.
 */

#include "GRN.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <sstream>

// Structure to store data points
struct DataPoint {
    double time;
    std::vector<double> values; // mRNA and protein values
};

/**
 * @brief Read signal strength from config file
 * 
 * @param config_file Path to configuration file
 * @return double Signal strength value, or 0.0 if not found
 */
double readSignalStrength(const std::string& config_file) {
    std::ifstream file(config_file);
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open config file " << config_file << std::endl;
        return 0.0;
    }
    
    std::string line;
    while (std::getline(file, line)) {
        // Remove whitespace and comments
        size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }
        
        // Look for SIGNAL_STRENGTH parameter
        if (line.find("SIGNAL_STRENGTH") != std::string::npos) {
            size_t equals_pos = line.find('=');
            if (equals_pos != std::string::npos) {
                std::string value_str = line.substr(equals_pos + 1);
                // Remove whitespace
                value_str.erase(0, value_str.find_first_not_of(" \t"));
                value_str.erase(value_str.find_last_not_of(" \t") + 1);
                
                try {
                    double signal_strength = std::stod(value_str);
                    std::cout << "Read signal strength: " << signal_strength << std::endl;
                    return signal_strength;
                } catch (const std::exception& e) {
                    std::cerr << "Error parsing signal strength: " << e.what() << std::endl;
                    return 0.0;
                }
            }
        }
    }
    
    std::cerr << "Warning: SIGNAL_STRENGTH not found in config file, using default 0.0" << std::endl;
    return 0.0;
}

/**
 * @brief Main function to simulate the Gene Regulatory Network (GRN).
 * 
 * This function initializes the GRN, introduces an external signal at a specified time,
 * and outputs the mRNA and protein concentrations to a file.
 * 
 * @return int Returns 0 on success, or 1 if an error occurs.
 */
int main() {
    // Example configuration file, kinetics file, and interactions file
    std::string config_file = "config.txt";
    std::string kinetics_file = "kinetics.txt";
    std::string interactions_file = "interactions.txt";
    std::string output_file = "concentrations.txt";
    std::string random_output_file = "random_before_signal.txt";

    // Create an instance of GRN
    try {
        double total_time = 240.0;
        double dt = 0.5;
        double after_signal_time = 80.0;
        double timeToStartSignal = 60.0;
        
        // Read signal strength from config file
        double signal_strength = readSignalStrength(config_file);
        
        // Vector to store data points after signal introduction
        std::vector<DataPoint> before_signal_data;

        GRN grn(config_file, kinetics_file, interactions_file);
        GRN grn_with_signal(config_file, kinetics_file, interactions_file, signal_strength);

        std::cout << "Simulating GRN with signal strength: " << signal_strength << std::endl;
        
        // Using advance function and saving results to concentrations.txt
        std::ofstream out(output_file);
        if (!out) {
            std::cerr << "Error: could not open output file " << output_file << std::endl;
            return 1;
        }
        
        double current_time = 0.0;
        bool signal_started = false;
        while (current_time <= total_time) {
            if (current_time < timeToStartSignal) {
                grn.advance(dt, current_time);
            } else {
                if (!signal_started) {
                    // Synchronize the state of grn_with_signal with grn only once
                    for (int i = 0; i < grn.getNumberOfGenes(); i++) {
                        grn_with_signal.setMRNALevel(i, grn.getMRNALevel(i));
                        grn_with_signal.setProteinLevel(i, grn.getProteinLevel(i));
                    }
                    // Also synchronize the Time_NextJump_ to maintain stochastic continuity
                    grn_with_signal.setNextJumpTime(grn.getNextJumpTime());
                    signal_started = true;
                }
                grn_with_signal.advance(dt, current_time);
            }
            
            // Output the current state
            out << current_time;
            
            // Create a data point for after-signal storage
            DataPoint data_point;
            data_point.time = current_time;
            
            for (int i = 0; i < grn.getNumberOfGenes(); i++) {
                if (current_time < timeToStartSignal) {
                    double mrna = grn.getMRNALevel(i);
                    double protein = grn.getProteinLevel(i);
                    out << " " << mrna << " " << protein;
                    data_point.values.push_back(mrna);
                    data_point.values.push_back(protein);
                } else {
                    double mrna = grn_with_signal.getMRNALevel(i);
                    double protein = grn_with_signal.getProteinLevel(i);
                    out << " " << mrna << " " << protein;
                    
                    // Store the values if after signal
                    data_point.values.push_back(mrna);
                    data_point.values.push_back(protein);
                }
            }
            out << std::endl;
            
            // Store the data point if before signal
            if (current_time < timeToStartSignal && current_time > after_signal_time) {
                before_signal_data.push_back(data_point);
            }
            
            current_time += dt;
        }

        out.close();
        
        // Process and export 32 random data points before signal
        if (before_signal_data.size() >= 32) {
            // Random number generator
            std::random_device rd;
            std::mt19937 gen(rd());
            
            // Randomly select 32 indices
            std::vector<size_t> indices(before_signal_data.size());
            for (size_t i = 0; i < indices.size(); ++i) {
                indices[i] = i;
            }
            std::shuffle(indices.begin(), indices.end(), gen);
            
            // Open output file for random data
            std::ofstream random_out(random_output_file);
            if (!random_out) {
                std::cerr << "Error: could not open random output file " << random_output_file << std::endl;
                return 1;
            }
            
            // Write header
            random_out << "Time";
            int num_genes = grn.getNumberOfGenes();
            for (int i = 0; i < num_genes; ++i) {
                random_out << " mRNA_" << i << " Protein_" << i;
            }
            random_out << std::endl;
            
            // Write the 32 randomly selected data points
            for (size_t i = 0; i < 32; ++i) {
                const DataPoint& data = before_signal_data[indices[i]];
                
                random_out << data.time;
                for (const double& value : data.values) {
                    random_out << " " << value;
                }
                random_out << std::endl;
            }
            
            random_out.close();
            std::cout << "Successfully exported 32 random data points before signal to " << random_output_file << std::endl;
        } else {
            std::cerr << "Error: Not enough data points after signal introduction (need 32, have " 
                      << before_signal_data.size() << ")" << std::endl;
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}