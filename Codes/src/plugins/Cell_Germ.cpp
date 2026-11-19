// ****************************************************************************
//
//              SiMuScale - Multi-scale simulation framework
//
// ****************************************************************************
//
// Copyright: See the AUTHORS file provided with the package
// E-mail: simuscale-contact@lists.gforge.inria.fr
// Original Authors : Samuel Bernard, Carole Knibbe, David Parsons
//
// Additional contributions: Pawan Kumar
// GRN, Cell_Germ, Cell_Somatic, updated view files, Python scripts,
// and universal CMake / cross-platform support.
// These implementations were used to produce results reported in the
// article cited in the main README.
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
#include "Cell_Germ.h"
//#include "CellSize.h"
#include "Simulation.h"
#include "params/SimulationParams.h"
#include "Population.h"
//#include "GRN.h"
#include <cmath> // Include for mathematical functions
#include <fstream> // Include for file operations
#include <iostream> // Include for input/output operations
#include <unistd.h>  // For getcwd
#include <cstdlib>   // For system
#include <algorithm>
#include <random>
#include <limits>
#include <unordered_map>
/**
 * @file Cell_Germ.cpp
 * @brief Implementation of the Cell_Germ class, a specific cell formalism plugin within the Simuscale framework.
 * This class models germ cells with specific behaviors and interactions, including a gene regulatory network (GRN) and volume transfer dynamics.
 * The Cell_Germ class inherits from the base Cell class and implements its own version of the Update, Save, Load, and other relevant methods to simulate germ cell behavior in a multi-scale simulation context.
 * The Cell_Germ class includes attributes for the GRN, internal state variables, and methods for handling volume transfer, death, and interactions with neighboring cells.
 */

// ============================================================================
//                       Definition of static attributes
// ============================================================================
constexpr CellFormalism Cell_Germ::classId_;
constexpr char Cell_Germ::classKW_[];
std::map<int, double> Cell_Germ::transfer_rates_map_;
bool Cell_Germ::transfer_rates_assigned_ = false;
double Cell_Germ::reference_volume_ = 0.0;
double Cell_Germ::reference_external_radius_ = 0.0;
double Cell_Germ::reference_internal_radius_ = 0.0;
bool Cell_Germ::reference_set_ = false;
double Cell_Germ::transfer_rate_min_ = 36.0;
double Cell_Germ::transfer_rate_max_ = 60.0;
double Cell_Germ::transfer_rate_mean_ = 48.0;     
double Cell_Germ::transfer_rate_std_ = 8.0;      
std::mt19937 Cell_Germ::rng_;
bool Cell_Germ::rng_initialized_ = false;
double Cell_Germ::reference_nuclei_factor_ = 0.25; // Default to 25% of external radius
// ============================================================================
//                                Constructors
// ============================================================================
/**
 * Create a new object with the provided constitutive elements and values
 */
Cell_Germ::Cell_Germ(CellType cellType,
                   const MoveBehaviour& move_behaviour,
                   const Coordinates<double>& pos,
                   double initial_volume,
                   double volume_min,
                   double doubling_time) :
    Cell(cellType, move_behaviour,
         pos, CellSize(initial_volume, volume_min),
         doubling_time),
         donor_cell_(0),
          activation_start_time_(-1.0) {
          internal_state_.resize(2, 0.0);
          read_additional_params();
          grn_ = new GRN("config.txt", "kinetics.txt", "interactions.txt");
          grn_->setCellVolume(initial_volume);
          setRefVolumeRadius();


}
Cell_Germ::Cell_Germ(const Cell_Germ &model):
  Cell(model), donor_cell_(model.donor_cell_), 
    activation_start_time_(model.activation_start_time_) {
    internal_state_ = model.internal_state_;
    Germ_mRNA_array_ = model.Germ_mRNA_array_;
    Germ_Protein_array_ = model.Germ_Protein_array_;
    donor_target_hours_ = model.donor_target_hours_;
    donor_rate_decimal_ = model.donor_rate_decimal_;
    donor_rate_initialized_ = model.donor_rate_initialized_;
    grn_ = new GRN(*model.grn_); // a copy constructor
}

Cell_Germ::Cell_Germ(gzFile backup_file) :
    Cell(backup_file) {
    internal_state_.resize(2, 0.0);
    grn_ = new GRN("config.txt", "kinetics.txt", "interactions.txt");
    Load(backup_file);
}
// ============================================================================
//                                 Destructor
// ============================================================================
// Cell_Germ::~Cell_Germ() = default;
Cell_Germ::~Cell_Germ() {
  delete grn_;
}

// ============================================================================
//                                   Methods
// ============================================================================
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
void Cell_Germ::InternalUpdate(const double& dt) {
  if (!concentrations_initialized_) {
    load_initial_concentrations();
    concentrations_initialized_ = true;
  }
  AssignTransferRates();
  setRefVolumeRadius();
  double global_time = Simulation::sim_time();
  double neighbourhood_score = get_output(InterCellSignal::S_nbd); // or 
  grn_->setExternalSignal(neighbourhood_score);

  // Advance the GRN simulation
  // receiver_volume_ratio = std::max(1.0, outer_volume() / reference_volume_);
  grn_->advance(dt, global_time);

  // Update cell state based on GRN output
  updateCellDiffState();

  // Record current time
  static thread_local double last_recorded_time = -1.0;
  if (global_time != last_recorded_time) {
    recordNurseCellActivation(global_time);
    last_recorded_time = global_time;
  }
  
  // Add final time for any simulation duration
  static thread_local bool final_recorded = false;
  static thread_local double cached_final_time = -1.0;
  
  if (!final_recorded) {
    // Calculate final time once
    if (cached_final_time < 0) {
      cached_final_time = Simulation::instance().max_timestep() * dt;
    }
    
    // Check if we're at the last normal timestep
    if (global_time + 1e-10 >= (cached_final_time - dt)) {
      recordNurseCellActivation(cached_final_time);
      final_recorded = true;
    }
  }
    AssignDonorReceiverRoles();
    // Execute volume transfer if this cell is a donor (only after transfer phase has started)
    if (donor_cell_ == 1) {
      auto cell_list = Simulation::instance().pop().cell_list();
      std::vector<Cell*> TotalCells(cell_list.begin(), cell_list.end());
      TransferVolumeToClosest(TotalCells, dt);
    }

    // Update death status
    updateDeathStatus();

    // After potential transfers and volume changes, refresh the local mirrored arrays
    // so they reflect the current GRN state and concentrations following any volume change.
    Germ_mRNA_array_[0] = grn_->getMRNALevel(0);
    Germ_mRNA_array_[1] = grn_->getMRNALevel(1);
    Germ_Protein_array_[0] = grn_->getProteinLevel(0);
    Germ_Protein_array_[1] = grn_->getProteinLevel(1);
  
}

//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
NeighbourCounts Cell_Germ::get_neighbour_cells_population() const {
  auto TotalCells = Simulation::instance().pop().cell_list();
  NeighbourCounts counts = {0, 0, 0, 0, 0, 0};

  // Calculate the reference neighborhood radius
  double refence_nbd_radius = nbd_factor_ *reference_external_radius_; // external_radius(); // nbd_factor_ * current_external_radius(); reference_external_radius_
  double sphere_nbd_radius = refence_nbd_radius;
  // std::cout << "[Debug] Sphere radius: %f\n", sphere_nbd_radius);

  for (auto cell : TotalCells) {
    Coordinates<double> diff = cell->pos() - pos();
    double distance = diff.norm();
    double neighbor_radius = cell->external_radius();

    if (cell->cell_type() == SOMATIC) {
      // Check for spherical region
      if (distance <= sphere_nbd_radius + neighbor_radius) {
        counts.nbh_somatic_sphere++;
      }
    } 
    else if (cell->cell_type() == GERM) {
      const Cell_Germ* germ_cell = dynamic_cast<const Cell_Germ*>(cell);
      if (germ_cell && germ_cell->internal_state_[NurseCellActivation] == 0) {
        // Check for spherical region
        if (distance <= sphere_nbd_radius + neighbor_radius) {
          counts.nbh_germ_sphere++;
        }
    }
    }
  }
  return counts;
}
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
double Cell_Germ::get_output(InterCellSignal signal) const {
  // Initialize concentrations if not already done, see if 
  if (!concentrations_initialized_) {
    Cell_Germ* non_const_this = const_cast<Cell_Germ*>(this);
    non_const_this->load_initial_concentrations();
    non_const_this->concentrations_initialized_ = true;
    
    // update the cell state based on these initial concentrations
    non_const_this->updateCellDiffState();
    
    // // Debug 
    // std::cout << "[DEBUG] Cell " << this->id() << " initialized in get_output()" << std::endl;
  }
  // Return the requested output based on the signal type
  if (signal == InterCellSignal::GERM)
    return germ_indentity;
  else if (signal == InterCellSignal::DEATH)
    return -1.0;
  else if (signal == InterCellSignal::CYCLE)
    return static_cast<double>(isCycling_);
  else if (signal == InterCellSignal::S_nbd) {
      NeighbourCounts counts = get_neighbour_cells_population();
      
      int total_neighbors = counts.nbh_germ_sphere + counts.nbh_somatic_sphere;
      if (total_neighbors == 0 || alpha_ + beta_ == 0.0 ) { //|| counts.nbh_somatic_sphere == 0.0
          return 0.0;  // Return 0 if no neighbors or invalid parameters
      }
      
      return (alpha_ * counts.nbh_germ_sphere + beta_ * counts.nbh_somatic_sphere) /
            ((alpha_ + beta_) * total_neighbors);
  }
  else if (signal == InterCellSignal::NURSE_CELL_ACTIVATION) {
    return internal_state_[NurseCellActivation];
  } 
  else if (signal == InterCellSignal::TOTAL_GERM_NBD_COUNT) {
    NeighbourCounts counts = get_neighbour_cells_population();
    return static_cast<double>(counts.nbh_germ_sphere);
  } 
  else if (signal == InterCellSignal::TOTAL_SOMATIC_NBD_COUNT) {
    NeighbourCounts counts = get_neighbour_cells_population();
    return static_cast<double>(counts.nbh_somatic_sphere);
  }
  else if (signal == InterCellSignal::GERM_P1) {
    return grn_->getProteinLevel(0);
  } 
  else if (signal == InterCellSignal::GERM_P2) {
    return grn_->getProteinLevel(1);
  } 
  else if (signal == InterCellSignal::GERM_mRNA1) {
    return grn_->getMRNALevel(0);
  } 
  else if (signal == InterCellSignal::GERM_mRNA2) {
    return grn_->getMRNALevel(1);
  } 
  else if (signal == InterCellSignal::DONOR_CELL) {
    return static_cast<double>(donor_cell_); 
  } 
  else {
    return 0.0; // if signal is unknown, return 0.0, not a failure
  }
}

//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------

void Cell_Germ::read_additional_params() {
  std::ifstream infile("additional_param.txt");
  if (!infile) {
    std::cerr << "ERROR: Could not open additional_param.txt" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  while (std::getline(infile, line)) {
    // Remove comments (anything after '#')
    size_t comment_pos = line.find('#');
    if (comment_pos != std::string::npos) {
      line = line.substr(0, comment_pos);
    }

    // Skip empty lines
    if (line.empty()) {
      continue;
    }

    // Parse the parameter and its value
    std::istringstream iss(line);
    std::string param;
    if (iss >> param) {
      if (param == "nbdFactor" || param == "alpha" || param == "beta" || param == "PTN_threshold") {
        double value;
        if (iss >> value) {
          if (param == "nbdFactor") {
            nbd_factor_ = value;
          } else if (param == "alpha") {
            alpha_ = value;
          } else if (param == "beta") {
            beta_ = value;
          } else if (param == "PTN_threshold") {
            PTN_threshold_ = value;
          }
        } 
        else {
          std::cerr << "ERROR: Invalid value for parameter " << param << " in additional_param.txt: " << line << std::endl;
          exit(EXIT_FAILURE);
        }
      } 
      else if (param == "transfer_rate_min") {       
        double value;
        if (iss >> value) {
          transfer_rate_min_ = value;
        } else {
          std::cerr << "ERROR: Invalid value for parameter " << param << " in additional_param.txt: " << line << std::endl;
          exit(EXIT_FAILURE);
        }
      } 
      else if (param == "transfer_rate_max") {   
        double value;
        if (iss >> value) {
          transfer_rate_max_ = value;
        } else {
          std::cerr << "ERROR: Invalid value for parameter " << param << " in additional_param.txt: " << line << std::endl;
          exit(EXIT_FAILURE);
        }
      } 
      else if (param == "transfer_rate_mean") {       
        double value;
        if (iss >> value) {
          transfer_rate_mean_ = value;
        } else {
          std::cerr << "ERROR: Invalid value for parameter " << param << " in additional_param.txt: " << line << std::endl;
          exit(EXIT_FAILURE);
        }
      } 
      else if (param == "transfer_rate_std") {   
        double value;
        if (iss >> value) {
          transfer_rate_std_ = value;
        } else {
          std::cerr << "ERROR: Invalid value for parameter " << param << " in additional_param.txt: " << line << std::endl;
          exit(EXIT_FAILURE);
        }
      }       
      else if (param == "nurse_cell_activation_file") {
        std::string value_str;
        if (iss >> value_str) {
          nurse_cell_activation_file_ = value_str; // Store the filename
        } else {
          std::cerr << "ERROR: Invalid value for parameter " << param << " in additional_param.txt: " << line << std::endl;
          exit(EXIT_FAILURE);
        }
      } 
      else if (param == "reference_nuclei_factor") {
        double value;
        if (iss >> value) {
          reference_nuclei_factor_ = value;
        } else {
          std::cerr << "ERROR: Invalid value for parameter " << param << " in additional_param.txt: " << line << std::endl;
          exit(EXIT_FAILURE);
        }
      }
      else {
        std::cerr << "ERROR: Unknown parameter in additional_param.txt: " << param << std::endl;
        exit(EXIT_FAILURE);
      }
    } 
    else {
      std::cerr << "ERROR: Invalid line in additional_param.txt: " << line << std::endl;
      exit(EXIT_FAILURE);
    }
  }
}
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
void Cell_Germ::updateCellDiffState() {
  int state = internal_state_[NurseCellActivation];
  // Thresholds for differentiation state
    if (!reference_set_ || reference_volume_ <= 0.0) {
      // Reference not yet initialized, skip threshold check
      internal_state_[NurseCellActivation] = state;
      return;
  }

  if (state == 10) {
    // If already activated, maintain the state
    internal_state_[NurseCellActivation] = state;
    return;
  }
  // Get current cell state
  double current_volume = outer_volume();
  double volume_ratio = current_volume / reference_volume_;
  // std::cout << "[DEBUG] Cell " << id() << " volume ratio: " << volume_ratio << std::endl;


  double current_PTN_concentration = grn_->getProteinLevel(1);
  double effective_threshold_conc = PTN_threshold_ * volume_ratio;

  if (current_PTN_concentration > (effective_threshold_conc)) {
    state = 10; // Activated
    // std::cout << "[ACTIVATE] Cell " << id() << " becomes donor!" << std::endl;
  } 

  internal_state_[NurseCellActivation] = state;
}
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
void Cell_Germ::load_initial_concentrations() {
  std::ifstream infile("random_before_signal.txt");
  if (!infile.is_open()) {
    std::cerr << "ERROR: Could not open random_before_signal.txt" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Get the global ID of this cell
  int global_id = this->id();

  // Get the list of all cells in the simulation
  auto TotalCells = Simulation::instance().pop().cell_list();

  // Create a mapping from global ID to local germ cell index
  std::vector<int> germ_cell_ids;
  for (const auto& cell : TotalCells) {
    if (cell->cell_type() == GERM) {
      germ_cell_ids.push_back(cell->id());
    }
  }

  // Sort the germ cell IDs to ensure consistent ordering
  std::sort(germ_cell_ids.begin(), germ_cell_ids.end());

  // Find the position of this cell's ID in the sorted list
  int local_index = -1;
  for (size_t i = 0; i < germ_cell_ids.size(); i++) {
    if (germ_cell_ids[i] == global_id) {
      local_index = i;
      break;
    }
  }

  if (local_index == -1) {
    std::cerr << "ERROR: Cell ID " << global_id << " not found in germ cell list" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Skip the header line
  std::string line;
  if (!std::getline(infile, line)) {
    std::cerr << "ERROR: File random_before_signal.txt is empty or improperly formatted." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Skip to the line corresponding to this cell's local index
  int current_line = 0;
  while (current_line <= local_index) {
    if (!std::getline(infile, line)) {
      std::cerr << "ERROR: Not enough lines in random_before_signal.txt for cell index " << local_index << std::endl;
      exit(EXIT_FAILURE);
    }

    // Skip comments and empty lines
    if (line.empty() || line[0] == '#') {
      continue;
    }

    current_line++;
  }

  // Parse the line to extract mRNA and protein values
  std::istringstream iss(line);
  double time, mRNA_0, protein_0, mRNA_1, protein_1;
  if (!(iss >> time >> mRNA_0 >> protein_0 >> mRNA_1 >> protein_1)) {
    std::cerr << "ERROR: Invalid format in random_before_signal.txt for cell index " << local_index << std::endl;
    exit(EXIT_FAILURE);
  }

  // Assign values to the member arrays
  Germ_mRNA_array_ = {mRNA_0, mRNA_1};
  Germ_Protein_array_ = {protein_0, protein_1};

  // Set these values in the GRN
  grn_->setMRNALevel(0, mRNA_0);
  grn_->setMRNALevel(1, mRNA_1);
  grn_->setProteinLevel(0, protein_0);
  grn_->setProteinLevel(1, protein_1);
  grn_->setCellVolume(outer_volume());

}
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
void Cell_Germ::recordNurseCellActivation(double global_time) {
  static thread_local bool file_initialized = false;
  static thread_local double last_recorded_time = -1.0;

  // // Debug: Print when function is called
  // std::cout << "[DEBUG] recordNurseCellActivation called at time: " << global_time << std::endl;

  // Create directory structure if needed
  if (!file_initialized) {
    std::string dir_path = nurse_cell_activation_file_;
    size_t last_slash = dir_path.find_last_of("/\\");
    
    if (last_slash != std::string::npos) {
      dir_path = dir_path.substr(0, last_slash);
      std::string mkdir_cmd = "mkdir -p \"" + dir_path + "\"";
      int result = system(mkdir_cmd.c_str());
      if (result != 0) {
        std::cerr << "ERROR: Failed to create directory: " << dir_path << std::endl;
      }
    }
    
    std::ofstream outfile(nurse_cell_activation_file_);
    if (!outfile.is_open()) {
      std::cerr << "ERROR: Could not create file: " << nurse_cell_activation_file_ << std::endl;
      return;
    }
    outfile << "TimeStep\tNurseCellCount\tTotalGermCells\tTotalProOocytes" << std::endl;
    outfile.close();
    // std::cout << "[DEBUG] Created file: " << nurse_cell_activation_file_ << std::endl;
    file_initialized = true;
  }

  // Skip if same time step
  if (global_time == last_recorded_time) {
    // std::cout << "[DEBUG]: Skipping duplicate time: " << global_time << std::endl;
    return;
  }

  // Count nurse cells with detailed debugging
  int nurse_cell_count = 0;
  int total_germ_cells = 0;
  int total_pro_oocytes = 0;
  auto TotalCells = Simulation::instance().pop().cell_list();

  // std::cout << "[DEBUG]: Total cells in simulation: " << TotalCells.size() << std::endl;

  for (const auto& cell : TotalCells) {
    if (cell->cell_type() == GERM) {
      total_germ_cells++;
      const Cell_Germ* germ_cell = dynamic_cast<const Cell_Germ*>(cell);
      if (germ_cell) {
        // Debug: Print activation state
        auto activation_value = germ_cell->internal_state_[NurseCellActivation];
        // std::cout << "[DEBUG]: Germ cell activation value: " << activation_value << std::endl;

        if (activation_value == 10) {
          nurse_cell_count++;
          // std::cout << "[DEBUG]: Found activated nurse cell!" << std::endl;
        } else if (activation_value == 0) {
          total_pro_oocytes++;
        }
      }
      //  else {
      //     std::cout << "[DEBUG]: Failed to cast to Cell_Germ" << std::endl;
      // }
    }
  }

  // std::cout << "[DEBUG]: Total germ cells: " << total_germ_cells 
  //           << ", Activated nurse cells: " << nurse_cell_count << std::endl;

  // Append data to file
  std::ofstream outfile(nurse_cell_activation_file_, std::ios::app);
  if (!outfile.is_open()) {
    std::cerr << "ERROR: Could not open file for appending: " << nurse_cell_activation_file_ << std::endl;
    return;
  }
  outfile << global_time << "\t" << nurse_cell_count << "\t" << total_germ_cells << "\t" << total_pro_oocytes << std::endl;
  outfile.close();

  // std::cout << "[DEBUG]: Recorded time: " << global_time << ", count: " << nurse_cell_count << std::endl;

  // Update last recorded time
  last_recorded_time = global_time;
}
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------

void Cell_Germ::AssignDonorReceiverRoles() {
    // By default all Germ cells are receiver cells
    receiver_cell_ids_.clear();
    receiver_cell_ids_.push_back(id());

    double current_time = Simulation::sim_time();

    // Update this cell's donor status
    if (internal_state_[NurseCellActivation] == 10) {
        if (activation_start_time_ < 0) {
            activation_start_time_ = current_time;
        }
        
        double activation_duration = current_time - activation_start_time_;
        if (activation_duration >= DONOR_ACTIVATION_TIME) {
            donor_cell_ = 1;
        }
    } else {
        activation_start_time_ = -1.0;
        donor_cell_ = 0;
      donor_target_hours_ = 0.0;
      donor_rate_decimal_ = 0.0;
      donor_rate_initialized_ = false;
    }

    // Global donor list management - collect all donors from all germ cells
    static std::vector<int> global_donor_list;
    global_donor_list.clear();
    
    auto all_cells = Simulation::instance().pop().cell_list();
    for (const auto& cell : all_cells) {
        if (cell->cell_type() == GERM) {
            const Cell_Germ* germ_cell = dynamic_cast<const Cell_Germ*>(cell);
            if (germ_cell && germ_cell->donor_cell_ == 1) {
                global_donor_list.push_back(germ_cell->id());
            }
        }
    }
    
    // Update donor list for this cell
    donor_cell_ids_ = global_donor_list;
  
}
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------

void Cell_Germ::AssignTransferRates(){
  if (!transfer_rates_assigned_) {
    // Clear the transfer rates map
    if (!rng_initialized_) {
      rng_.seed(std::random_device{}());
      rng_initialized_ = true;
    }
    transfer_rates_map_.clear();
    // Create normal distribution
    std::normal_distribution<double> normal_dist(transfer_rate_mean_, transfer_rate_std_);
 

    // Assign transfer rates for each germ cell using cell ID as key
    for (const auto& cell : Simulation::instance().pop().cell_list()) {
      if (cell->cell_type() == GERM) {
        const Cell_Germ* germ_cell = dynamic_cast<const Cell_Germ*>(cell);
        if (germ_cell) {
          int cell_id = germ_cell->id();
          double transfer_rate;
          int attempts = 0;
          const int max_attempts = 1000;
          bool sampled_successfully = false;

          // Try to sample from truncated Gaussian
          while (attempts < max_attempts) {
            transfer_rate = normal_dist(rng_);
            attempts++;
            
            // Check if sample is within bounds
            if (transfer_rate >= transfer_rate_min_ && 
                transfer_rate <= transfer_rate_max_) {
              sampled_successfully = true;
              break;  // Success! Exit the loop
            }
          }
          if (!sampled_successfully) {
            // Fallback to uniform distribution if Gaussian sampling fails, this should be rare
            std::uniform_real_distribution<double> uniform_dist(transfer_rate_min_, transfer_rate_max_);
            transfer_rate = uniform_dist(rng_);
            // std::cerr << "Warning: Truncated Gaussian sampling failed for cell " 
            //           << cell_id << " after " << max_attempts 
            //           << " attempts. Using uniform fallback: " 
            //           << transfer_rate << "%" << std::endl;
          }
          transfer_rates_map_[cell_id] = transfer_rate;

          


          // Debug output
          // std::cout << "[DEBUG] Assigned transfer rate " << transfer_rate 
          //           << " to cell ID " << cell_id << std::endl;
        }
      }
    }
    transfer_rates_assigned_ = true;

    // std::cout << "[DEBUG] Transfer rates assigned for " << transfer_rates_map_.size() 
    //           << " germ cells" << std::endl;
  }
}
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------

double Cell_Germ::getTransferRate(int cell_id) const {
    auto it = transfer_rates_map_.find(cell_id);
    if (it != transfer_rates_map_.end()) {
        return it->second;
    } else {
        std::cerr << "Warning: Transfer rate not found for cell ID " << cell_id << std::endl;
        return 0.5; // Default transfer rate
    }
}
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------

// Get transfer rate for this cell
double Cell_Germ::getMyTransferRate() const {
    return getTransferRate(this->id());
}
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
Cell_Germ* Cell_Germ::FindClosestCell(const std::vector<Cell*>& all_cells) {
    Cell_Germ* closest_cell = nullptr;
    double min_distance = std::numeric_limits<double>::max();
    double search_radius = nbd_factor_ * external_radius(); // same as neighbourhood
    
    Coordinates<double> my_pos = pos();
    
    for (Cell* cell : all_cells) {
        if (cell == this) continue;  // Skip self

        Cell_Germ* vt_cell = dynamic_cast<Cell_Germ*>(cell);
        if (vt_cell == nullptr) continue;  // Skip non-Germ cells

        // Only transfer to receiver cells (skip donors)
        if (vt_cell->donor_cell_) continue;

        // if a cell is set to die, let it die, you cant change someone's destiny
        if (vt_cell->isDying()) continue;
        
        double distance = (my_pos - cell->pos()).norm();
        double max_allowed_distance = search_radius + vt_cell->external_radius(); // same as neighbourhood check
        if (distance < min_distance && distance <= max_allowed_distance) {
            min_distance = distance;
            closest_cell = vt_cell;
        }
    }
    
    return closest_cell;
}
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
void Cell_Germ::setRefVolumeRadius() {
  if (!reference_set_) {
      reference_volume_ = this->outer_volume();
      reference_external_radius_ = this->external_radius();
      double reference_internal_volume = reference_nuclei_factor_ * reference_volume_;
      reference_internal_radius_ = std::cbrt((3.0 * reference_internal_volume) / (4.0 * M_PI));
      reference_set_ = true;
      //std::cout << "[DEBUG] Cell " << id() << " reference volume set to " << reference_volume_ << std::endl;
  }
}
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
void Cell_Germ::ReceiveVolume(double volume, double donor_mRNA_0, double donor_mRNA_1, double donor_protein_0, double donor_protein_1) {
    double old_volume = outer_volume();
    double min_vol = size_.volume_min();
    double new_volume = old_volume + volume;

  
    // std::cout << "[DEBUG] RECEIVE: Cell " << id()
    //           << " | Current: " << old_volume
    //           << " | Min: " << min_vol
    //           << " | Target: " << new_volume << std::endl;
    
    // Initialize static variables if needed (this will work) didnt :P
    // CellSize::set_volume_max_min_ratio(2.0);

    // Dynamically increase ratio if needed
    double required_ratio = (new_volume / min_vol);
    if (required_ratio > 2.2) {
        // std::cout << "[DEBUG] Expanding ratio to " << required_ratio << " to accommodate volume " << new_volume << std::endl;
        CellSize::set_volume_max_min_ratio(required_ratio);
    }
    // Get current molecular amounts before volume change
    double old_amount_mRNA_0 = grn_->getMRNAAmount(0);
    double old_amount_mRNA_1 = grn_->getMRNAAmount(1);
    double old_amount_protein_0 = grn_->getProteinAmount(0);
    double old_amount_protein_1 = grn_->getProteinAmount(1);

       // Calculate old concentration for logging
    double old_conc_PTN = (old_volume > 0) ? (old_amount_protein_1 / old_volume) : 0;
    
        
    // Add incoming molecular amounts from donor
    double new_amount_mRNA_0 = old_amount_mRNA_0 + donor_mRNA_0;
    double new_amount_mRNA_1 = old_amount_mRNA_1 + donor_mRNA_1;
    double new_amount_protein_0 = old_amount_protein_0 + donor_protein_0;
    double new_amount_protein_1 = old_amount_protein_1 + donor_protein_1;
    
    // Update volume
    size_.set_volume(new_volume);
    double actual_new_volume = outer_volume();

    // Update GRN volume
    grn_->setCellVolume(actual_new_volume);
    
    // Update GRN with new total amounts
    grn_->setMRNAAmount(0, new_amount_mRNA_0);
    grn_->setMRNAAmount(1, new_amount_mRNA_1);
    grn_->setProteinAmount(0, new_amount_protein_0);
    grn_->setProteinAmount(1, new_amount_protein_1);
    
    // Calculate new concentration
    double new_conc_PTN = grn_->getProteinLevel(1);
    double conc_change_percent = (old_conc_PTN > 0) ? 
                                  ((new_conc_PTN / old_conc_PTN) - 1.0) * 100.0 : 0;

    
    // // DEBUG LOGGING
    // std::cout << "[RECEIVER-S2] Cell " << id() << " Complete Transfer" << std::endl;
    // std::cout << "  Volume: " << old_volume << " -> " << actual_new_volume 
    //           << " (+" << ((actual_new_volume - old_volume) / old_volume * 100.0) << "%)" << std::endl;
    // std::cout << "  PTN: " << old_protein_1 << " + " << donor_protein_1 
    //           << " = " << new_protein_1_total << std::endl;
    // std::cout << "  PTN conc: " << old_conc_PTN << " -> " << new_conc_PTN
    //           << " (" << std::showpos << conc_change_percent << std::noshowpos << "%)" 
    //           << std::endl;        


  }
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------

void Cell_Germ::TransferVolumeToClosest(const std::vector<Cell*>& all_cells, double dt) {

    // Early exit if there is no donor cell
    if (donor_cell_!=1) {
        return;
    }
    
    // Find the closest receiver cell
    Cell_Germ* closest_cell = FindClosestCell(all_cells);
    if (closest_cell == nullptr) {
        return;  // No receiver cells found
    }
    
    // Get current cell dimensions for minimum size protection
    double current_volume = outer_volume();
    double nucleus_volume = reference_nuclei_factor_ * reference_volume_;
    
    // Initialize when cell first becomes donor
    if (!donor_rate_initialized_) {
        // Get timeline from truncated Gaussian (mean 48h)
        donor_target_hours_ = getMyTransferRate();
        
        // Calculate rate to reach FIXED nucleus from CURRENT volume
        // This adapts to whatever size the cell has when it becomes donor
        double volume_ratio = nucleus_volume / current_volume;
        if (donor_target_hours_ <= 0.0 || volume_ratio <= 0.0) {
          return;
        }

        donor_rate_decimal_ = 1.0 - std::pow(volume_ratio, 1.0 / donor_target_hours_);
        donor_rate_initialized_ = true;
        
    }
    
    double rate_decimal = donor_rate_decimal_;
    
    // Safety check
    double max_transferable = current_volume - nucleus_volume;
    if (max_transferable <= 1e-6) {
        return;
    }
    
    // Exponential decay transfer
    double remaining_fraction = std::pow(1.0 - rate_decimal, dt);
    double transfer_fraction = 1.0 - remaining_fraction;
    double desired_transfer = current_volume * transfer_fraction;
    
    double volume_to_transfer = std::min(desired_transfer, max_transferable);
    
    if (volume_to_transfer < 1e-6) {
        return;
    }

    double transfer_ratio = volume_to_transfer / current_volume;
    
    // Get current donor molecular AMOUNTS using GRN methods
    double donor_amount_mRNA_0 = grn_->getMRNAAmount(0);
    double donor_amount_mRNA_1 = grn_->getMRNAAmount(1);
    double donor_amount_protein_0 = grn_->getProteinAmount(0);
    double donor_amount_protein_1 = grn_->getProteinAmount(1);

    double transferred_mRNA_0 = 0.0;
    double transferred_mRNA_1 = 0.0;
    double transferred_protein_0 = 0.0;
    double transferred_protein_1 = 0.0;

    double old_donor_volume = current_volume;
    double new_donor_volume = current_volume - volume_to_transfer;

    // Get old donor PTN concentration
    double old_donor_conc_PTN = grn_->getProteinLevel(1);
      
    transferred_mRNA_0 = donor_amount_mRNA_0 * transfer_ratio;
    transferred_mRNA_1 = donor_amount_mRNA_1 * transfer_ratio;
    transferred_protein_0 = donor_amount_protein_0 * transfer_ratio;
    transferred_protein_1 = donor_amount_protein_1 * transfer_ratio;
    
    // Calculate new donor amounts
    double new_donor_amount_mRNA_0 = donor_amount_mRNA_0 - transferred_mRNA_0;
    double new_donor_amount_mRNA_1 = donor_amount_mRNA_1 - transferred_mRNA_1;
    double new_donor_amount_protein_0 = donor_amount_protein_0 - transferred_protein_0;
    double new_donor_amount_protein_1 = donor_amount_protein_1 - transferred_protein_1;
    
    // Update donor volume
    this->size_.set_volume(new_donor_volume);
    
    // Update GRN volume
    grn_->setCellVolume(new_donor_volume);
    
    // Set new amounts in GRN (automatically converts to concentrations)
    grn_->setMRNAAmount(0, new_donor_amount_mRNA_0);
    grn_->setMRNAAmount(1, new_donor_amount_mRNA_1);
    grn_->setProteinAmount(0, new_donor_amount_protein_0);
    grn_->setProteinAmount(1, new_donor_amount_protein_1);
    
    // Get new donor PTN concentration (should be ~same as old)
    double new_donor_conc_PTN = grn_->getProteinLevel(1);
    double donor_conc_change = (old_donor_conc_PTN > 0) ? 
                                ((new_donor_conc_PTN / old_donor_conc_PTN) - 1.0) * 100.0 : 0;
    
    // // Logging (uncomment for debugging)
    // std::cout << "[DONOR-S2] Cell " << id() << " Complete Transfer" << std::endl;
    // std::cout << "  Volume: " << old_donor_volume << " → " << new_donor_volume
    //           << " (-" << (transfer_ratio * 100.0) << "%)" << std::endl;
    // std::cout << "  PTN amount: " << donor_amount_protein_1 << " - " << transferred_protein_1
    //           << " = " << new_donor_amount_protein_1 << std::endl;
    // std::cout << "  PTN conc: " << old_donor_conc_PTN << " → " << new_donor_conc_PTN
    //           << " (" << std::showpos << donor_conc_change << std::noshowpos 
    //           << "%) [Expected: ~0%]" << std::endl;
          
    // Transfer to receiver
    closest_cell->ReceiveVolume(volume_to_transfer, 
                                 transferred_mRNA_0, transferred_mRNA_1,
                                 transferred_protein_0, transferred_protein_1);

    // Detailed logging
    // std::cout << "[DEBUG: TRANSFER] Time: " << Simulation::sim_time() 
    //           << " | Donor " << id() << " -> Receiver " << closest_cell->id()
    //           << " | Rate: " << rate_percent << "% per hour"
    //           << " | Transferred: " << volume_to_transfer << " µm³"
    //           << " | Donor: " << current_volume << " -> " << outer_volume()
    //           << " | Receiver: " << (closest_cell->outer_volume() - volume_to_transfer) 
    //           << " -> " << closest_cell->outer_volume() << std::endl;


}
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
void Cell_Germ::updateDeathStatus() {
    if (!isDying() && this->donor_cell_==1) {
        // Set death threshold
      double current_volume = outer_volume();
      double nucleus_volume = reference_nuclei_factor_ * reference_volume_;
      double dying_volume_threshold = 1.5 * nucleus_volume;
      double death_volume_threshold = nucleus_volume;

        if (current_volume < dying_volume_threshold && current_volume > death_volume_threshold) {
                // std::cout << "[DEBUG] Cell " << this->id() << " is DYING with radius " << this->external_radius() << std::endl;
            set_to_dying();
        }
    }
    
    if (isDying()) {
      double death_volume_threshold = reference_nuclei_factor_ * reference_volume_;
        if (outer_volume() <= death_volume_threshold) {
                // std::cout << "[DEBUG] Cell " << this->id() << " is DEAD with radius " << this->external_radius() << std::endl;
            set_to_dead();
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
bool Cell_Germ::isDying() const {
  return is_dying_;
}
//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
bool Cell_Germ::isDead() const {
  return is_dead_;
}
//-----------------------------------------------------------------------------------------------------------------------------

Coordinates<double> Cell_Germ::MotileDisplacement(const double& dt) {
  // <TODO> Motile cell autonomous displacement behaviour
  Coordinates<double> displ = dt*sigma_*orientation();
  // Coordinates<double> displ { sqrt(dt)*sigma_t*Alea::gaussian_random(),
  //                             sqrt(dt)*sigma_t*Alea::gaussian_random(),
  //                             - dt*0.01 };
  return displ;
}

Cell* Cell_Germ::Divide(void) {
  // <TODO> Determine how the cell should divide </TODO>
  Cell_Germ* newCell = new Cell_Germ(*this);

  SeparateDividingCells(this, newCell);

  return newCell;
}

void Cell_Germ::Save(gzFile backup_file) const {
  // Write my classId
  gzwrite(backup_file, &classId_, sizeof(classId_));

  // Write the generic cell stuff
  Cell::Save(backup_file);

  // Write the specifics of this class
  // <TODO> Save your specific attributes </TODO>
}

void Cell_Germ::Load(gzFile backup_file) {
  // <TODO> Load your specific attributes </TODO>
}

double Cell_Germ::gaussian_field_weight(InterCellSignal signal) {
  if (signal == InterCellSignal::GERM)
    return germ_indentity;
  else
    return 0.0;
}

bool Cell_Germ::isCycling() const {
  // <TODO> Is the cell currently cycling ? </TODO>
  return isCycling_;
}

bool Cell_Germ::isDividing() const {
  // <TODO> Should the cell divide now ? </TODO>
  if ( isCycling() && outer_volume() > size_.volume_min() ) 
    return true;
  else
    return false;
}
void Cell_Germ::UpdateCyclingStatus() {

  if ( nn_ ) {
    isCycling_ = true;
  }
  else {
    isCycling_ = false;
  }
}

bool Cell_Germ::StopCycling() {
  // <TODO> Should the cell stop cycling now ? </TODO>
  return true; // cell is not cycling.........
}

bool Cell_Germ::StartCycling() {
  // <TODO> Should the cell start cycling now ? </TODO>
  return false;
}

// Register this class in Cell
bool Cell_Germ::registered_ =
    Cell::RegisterClass(classId_, classKW_,
                        [](CellType type,
                           const MoveBehaviour& move_behaviour,
                           const Coordinates<double>& pos,
                           double initial_volume,
                           double volume_min,
                           double doubling_time){
                          return static_cast<Cell*>(
                              new Cell_Germ(type, move_behaviour,
                                           pos, initial_volume, volume_min,
                                           doubling_time));
                        },
                        [](gzFile backup_file){
                          return static_cast<Cell*>(new Cell_Germ(backup_file));
                        }
    );