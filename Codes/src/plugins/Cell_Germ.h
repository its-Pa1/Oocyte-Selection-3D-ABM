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

// <TODO> Modify the include guard </TODO>
#ifndef SIMUSCALE_CELL_GERM_H__
#define SIMUSCALE_CELL_GERM_H__


// ============================================================================
//                                   Includes
// ============================================================================
#include "Cell.h"
#include "GRN.h"
struct NeighbourCounts {
    double nbh_somatic_cylinder;
    double nbh_somatic_sphere;
    double nbh_somatic_ellipsoid;
    double nbh_germ_cylinder;
    double nbh_germ_sphere;
    double nbh_germ_ellipsoid;
};

/**
 * @file Cell_Germ.h
 * @brief Declaration of the Cell_Germ class, a specific cell formalism plugin within the Simuscale framework.
 * This class models germ cells with specific behaviors and interactions, including a gene regulatory network (GRN) and volume transfer dynamics.
 * The Cell_Germ class inherits from the base Cell class and implements its own version of the Update, Save, Load, and other relevant methods to simulate germ cell behavior in a multi-scale simulation context.
 * The Cell_Germ class includes attributes for the GRN, internal state variables, and methods for handling volume transfer, death, and interactions with neighboring cells.
 * 
 */


/**
 * @class  Cell_Germ
 * @brief  A specific cell formalism plugin representing germ cells in the Simuscale framework
 * I
 * 
 * 
 */
class Cell_Germ : public Cell {
  private:
    GRN* grn_; // A pointer to the gene regulatory network

 public :
  typedef enum {
    NurseCellActivation, NeutralState, DifferentiationState // -1, 0, 1 flags
  } StateVariable;
  // ==========================================================================
  //                               Constructors
  // ==========================================================================
  Cell_Germ() = default; //< Default ctor
  Cell_Germ(const Cell_Germ &); //< Copy ctor
  Cell_Germ(Cell_Germ &&) = delete; //< Move ctor
  Cell_Germ(CellType cellType,
                   const MoveBehaviour& move_behaviour,
                   const Coordinates<double>& pos,
                   double initial_volume,
                   double volume_min,
                   double doubling_time);
  Cell_Germ(gzFile backup_file);

  // ==========================================================================
  //                                Destructor
  // ==========================================================================
  virtual ~Cell_Germ();

  // ==========================================================================
  //                                Operators
  // ==========================================================================
  /// Copy assignment
  Cell_Germ &operator=(const Cell_Germ &other) = delete;

  /// Move assignment
  Cell_Germ &operator=(Cell_Germ &&other) = delete;

  // ==========================================================================
  //                              Public Methods
  // ==========================================================================
  void InternalUpdate(const double& dt) override;
  Cell* Divide(void) override;

  void Save(gzFile backup_file) const override;
  void Load(gzFile backup_file) override;
  NeighbourCounts get_neighbour_cells_population() const;
  Cell_Germ* FindClosestCell(const std::vector<Cell*>& all_cells);
  void TransferVolumeToClosest(const std::vector<Cell*>& all_cells, double dt);
  void ReceiveVolume(double volume, double donor_mRNA_0 = 0.0, double donor_mRNA_1 = 0.0,
                      double donor_protein_0 = 0.0, double donor_protein_1 = 0.0);   
  // ==========================================================================
  //                                Accessors
  // ==========================================================================
  double get_output(InterCellSignal signal) const override;

  bool isCycling() const override;
  bool isDividing() const override;
  bool isDying() const override;
  bool isDead() const override;

  // Fast Gaussian Transform -----------------------------------------
  double gaussian_field_weight(InterCellSignal signal) override;
  // END Fast Gaussian Transform -------------------------------------

 protected :
  // ==========================================================================
  //                            Protected Methods
  // ==========================================================================
  void UpdateCyclingStatus() override;
  bool StopCycling() override;
  bool StartCycling() override;

  Coordinates<double> MotileDisplacement(const double &dt) override;
  void updateDeathStatus();

  // ==========================================================================
  //                               Attributes
  // ==========================================================================
public:
  // <TODO> Modify the class id and keyword </TODO>
  // Class ID, uniquely identifies this class
  static constexpr CellFormalism classId_ = 0x97fb477e; // CRC32 Hash 97fb477e
  // Keyword to be used in param files
  static constexpr char classKW_[] = "GERM";
  // void InitializeConcentrations();

 protected:
  // <TODO> Add your formalism-specific attributes here </TODO>

  double germ_indentity = 1.0;
  double sigma_ = 1.5;
  int32_t nn_ = 0;
  std::vector<double> internal_state_;
  std::vector<double> Germ_mRNA_array_; /**< Array of mRNA concentrations. */
  std::vector<double> Germ_Protein_array_; /**< Array of protein concentrations. */


 private:
  /** dummy attribute - allows to register class in Simuscale statically */
  static bool registered_;
  double alpha_ = 1.0;
  double beta_ = 1.0;
  double nbd_factor_ = 1.0;
  double PTN_threshold_ = 10.0;
  void read_additional_params();
  void updateCellDiffState();
  void load_initial_concentrations(); // Method to load initial mrna and protein concentrations (random ones from a file)
  void recordNurseCellActivation(double global_time);
  void AssignDonorReceiverRoles();
  void AssignTransferRates();
  double getTransferRate(int cell_id) const;
  double getMyTransferRate() const;
  void setRefVolumeRadius();
  std::string nurse_cell_activation_file_;
  //int GetDifferentiationState();
  //void assign_update_GRN();
  bool concentrations_initialized_ = false;
  int donor_cell_;  // -1: somatic, 0: germ (receiver), 1: donor
  // bool donor_cell_ = false; /**< Flag indicating if the cell is a donor. */
  std::vector<int> receiver_cell_ids_; /**< List of receiver cell IDs. */  // once a receiver stays a receiver
  std::vector<int> donor_cell_ids_; /**< List of donor cell IDs. */ // these will die
  static constexpr double DONOR_ACTIVATION_TIME = 0.0; //  hibernation period
  double activation_start_time_ = -1.0;  // Time when NurseCellActivation started
  static std::map<int, double> transfer_rates_map_;  // Map of cell_id -> transfer_rate
  static bool transfer_rates_assigned_;
  static double reference_volume_;
  static double reference_external_radius_;
  static double reference_internal_radius_;
  static bool reference_set_;
  static double transfer_rate_mean_; 
  static double transfer_rate_std_; 
  static double transfer_rate_min_;
  static double transfer_rate_max_;   
  // bool pro_oocyte_;
  static std::mt19937 rng_;
  static bool rng_initialized_;
  static double reference_nuclei_factor_;
  double donor_target_hours_ = 0.0;
  double donor_rate_decimal_ = 0.0;
  bool donor_rate_initialized_ = false;

};

#endif // SIMUSCALE_CELL_GERM_H__
