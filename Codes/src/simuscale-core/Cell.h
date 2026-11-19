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

#ifndef SIMUSCALE_CELL_H__
#define SIMUSCALE_CELL_H__


// =================================================================
//                              Includes
// =================================================================
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cfloat>

#include <list>

#include <zlib.h>

#include "CellType.h"
#include "Coordinates.h"
#include "Observable.h"
#include "Alea.h"
#include "InterCellSignal.h"
#include "InterCellSignals.h"
#include "CellSize.h"
#include "movement/MoveBehaviour.h"

// ============================================================================
//                         Class declarations, Using etc
// ============================================================================
class Mobile;

using real_type = InterCellSignals::real_type_;

/*!
  \brief An abstraction for any kind of cells.
*/
class Cell : public Observable
{
  friend class Mobile;
  friend class Motile;

 public :


  using CellFactory = Cell* (*)(CellType type,
                                const MoveBehaviour& move_behaviour,
                                const Coordinates<double>& pos,
                                double initial_volume,
                                double volume_min,
                                double doubling_time);
  using CellFactories = std::map<CellFormalism, CellFactory>;
  using CellLoader = Cell* (*)(gzFile backup_file);
  using CellLoaders = std::map<CellFormalism, CellLoader>;
  static bool RegisterClass(CellFormalism formalism,
                            const char param_file_keyword[],
                            CellFactory factory,
                            CellLoader loader);

  // =================================================================
  //                             Constructors
  // =================================================================
  Cell() = default;
  Cell(const Cell &model); //< Copy ctor
  Cell(CellType cell_type, 
       const MoveBehaviour& move_behaviour,
       Coordinates<double> pos,
       CellSize&& size,
       double doubling_time);
  Cell(gzFile backup_file);


  // =================================================================
  //                             Destructor
  // =================================================================
  virtual ~Cell() = default;


  // =================================================================
  //                            Public Methods
  // =================================================================
  void Update(const double& dt);
  double Distance(const Cell* other) const;
  void AddNeighbour(Cell * neighbour);
  void ResetNeighbours();

  /**
   * Make the cell divide
   *
   * The cell will continue to exist as one of the daughter cells
   * @return the second daughter cell
   */
  virtual Cell* Divide() = 0;

  void ComputeInteractions();
  void ComputeGaussianFields();
  void AddGaussianField(InterCellSignal signal, std::vector<real_type>& value);
  void ResetInteractions();
  void Move(const double& dt);

  void UpdateMinMaxSignal() const;
    
  /**
   * Save the cell into a backup
   *
   * @param backup_file a write-opened gz file
   */
  virtual void Save(gzFile backup_file) const;

  /**
   * Load a cell's generic data from a backup
   *
   * @param backup_file a read-opened gz file
   * @see LoadCell(gzFile backup_file)
   */
  virtual void Load(gzFile backup_file);

  /**
   * Load a cell of any known formalism from a backup (Polymorphic factory)
   *
   * Return a heap allocated Cell instance of the same derived type as the
   * original cell reading from a gzFile
   *
   * \throws std::out_of_range if the derived type is unknown
   */
  static Cell* MakeCell(CellFormalism formalism,  
                        const MoveBehaviour& move_behaviour,
                        CellType type,
                        Coordinates<double> pos,
                        double initial_volume,
                        double volume_min,
                        double doubling_time);

  static Cell* LoadCell(gzFile backup_file);
  static void SaveStatic(gzFile backup_file);
  static void LoadStatic(gzFile backup_file);

  // =================================================================
  //                              Accessors
  // =================================================================
  int32_t id() const {return id_;};
  const Coordinates<double>& pos() const {return pos_;};
  double pos_x() const {return pos_.x;};
  double pos_y() const {return pos_.y;};
  double pos_z() const {return pos_.z;};
  const std::vector<Cell*>& neighbours() const {return neighbours_;};

  const Coordinates<double>& orientation() const { return orientation_; };
  double orientation_x() const { return orientation_.x; };
  double orientation_y() const { return orientation_.y; };
  double orientation_z() const { return orientation_.z; };

  const CellType& cell_type() const {return cell_type_;};
  double outer_volume() const {return size_.volume();};
  double internal_radius() const {return size_.internal_radius();};
  double external_radius() const {return size_.external_radius();};

  virtual bool isCycling() const { return false; }; // isCycling: deprecated, make false by default
  virtual bool isDividing() const = 0;
  virtual bool isDying() const { return is_dying_; };
  virtual bool isDead() const;

  virtual double get_output(InterCellSignal signal) const = 0; // deprecated, use local_signal

  virtual double local_signal(InterCellSignal signal) const;

  static void set_volume_max_min_ratio(double volume_max_min_ratio) {
    CellSize::set_volume_max_min_ratio(volume_max_min_ratio);
  }
  static void set_radii_ratio(double radii_ratio) {
    CellSize::set_radii_ratio(radii_ratio);
  }

  // FGT Fast Gaussian Transform -----------------------------------------
  virtual double                            gaussian_field_weight(InterCellSignal signal);
  virtual const Coordinates<double>&        gaussian_field_source(InterCellSignal signal); 
  virtual std::vector<Coordinates<double>>  gaussian_field_targets(InterCellSignal signal);
  // END Fast Gaussian Transform -------------------------------------

  virtual const char* cell_formalism() const { return NULL;};

 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================
  virtual void InternalUpdate(const double& dt) = 0;
  static void SeparateDividingCells(Cell* c1, Cell* c2);
  virtual void Grow(const double& dt); // Cell::Grow has a default behaviour
  void Move(Coordinates<double> dpos);
  virtual Coordinates<double> MotileDisplacement(const double& dt); // deprecated -> motility
  virtual Coordinates<double> motility(const double& dt);           // replaces MotileDisplacement

  virtual void UpdateCyclingStatus() { /* do nothing */ };          // deprecated, no replacement
  virtual bool StopCycling() { return false; };                     // deprecated, no replacement
  virtual bool StartCycling() { return false; };                    // deprecated, no replacement 

  virtual void AddMechForce(const Cell* other);
  void AddChemCom(const Cell* other);
  void ResetMechForce();
  void ResetInSignals();

  inline void set_to_dying() { is_dying_ = true; };
  inline void set_to_dead()  { is_dead_  = true; }; 

  double getInSignal(InterCellSignal inSignal) const;               // deprecated, -> get_local_signal

  double get_local_signal(InterCellSignal inSignal) const; 

  std::vector<real_type> get_diffusive_signal(InterCellSignal inSignal) const;


 private:
  // =================================================================
  //                           Private Methods
  // =================================================================

 protected:
  // =================================================================
  //                              Attributes
  // =================================================================
  /** Maximum ID that has been attributed so far */
  static uint32_t maxID_;

  /** Cell ID */
  int32_t id_;

  /** Coordinates in space */
  Coordinates<double> pos_;

  /** Orientation in space: normalized Coordinates */
  Coordinates<double> orientation_;

  /** Size */
  CellSize size_;

  /** Cell Type */
  CellType cell_type_;

  /** Time needed for the volume of the cell to be doubled */
  double doubling_time_;

  /**
   * Factor by which the volume is multiplied at each time step during growth
   *
   * It is a simple transformation of the (exponential) growth rate when both
   * the time_step and doubling time are known.
   * It represents the common ratio of the geometric progression of the volume
   *
   * This attribute exists only for the sake of performance since it can be
   * computed once and for all
   */
  double growth_factor_;

  /**
   * Defines how the cell moves (strategy pattern)
   */
  const MoveBehaviour* move_behaviour_ = nullptr;

  /**
   * Combination of all the mechanical forces exerted on the cell
   */
  Coordinates<double> mechanical_force_;

  /**
   * Whether the cell is cycling or not
   */
  bool isCycling_ = false;  // deprecated, do not use anymore

  /**
   * Whether the cell is alive or not
   *
   * set is_dying_ to flag the cell as dying. This will
   * mark the cell as 'D' in the output trajectory file
   * but keep the cell in the cell list.
   *
   * getter: isDying() 
   * setter: set_to_dying() 
   */
  bool is_dying_ = false;
  
  /** 
   * set is_dead_ to flag the cell as dead. Cell will be removed
   * from the cell population by the core simulator, 
   * and will not be listed in the trajectory file.
   *
   * getters: isDead()
   * setters: set_to_dead().
   */
  bool is_dead_ = false;

  /**
   * Cells that are in physical contact with the object
   */
  std::vector<Cell*> neighbours_;

  /**
   * Chemical inputs perceived by the cell
   *
   * This is basically the sum of the outputs emitted by neighbouring cells
   */
  InterCellSignals intrinsic_inputs_;

  /**
   * Container for derived class factories
   *
   * This is a method for technical reasons only, it returns a reference to
   * the actual "attribute" that can be used normaly, including as an lvalue
   */
  static CellFactories& factories();

  /**
   * Container for derived class loader factories
   *
   * This is a method for technical reasons only, it returns a reference to
   * the actual "attribute" that can be used normaly, including as an lvalue
   */
  static CellLoaders& loaders();
};

#endif // SIMUSCALE_CELL_H__
