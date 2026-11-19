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

#ifndef SIMUSCALE_CELL_SOMATIC_H__
#define SIMUSCALE_CELL_SOMATIC_H__


// =================================================================
//                              Includes
// =================================================================
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

#include "Cell.h"
#include "CellType.h"


/*!
  This class implements a specific model based on a predefined set of ODEs

  This model is described in the following article:
  Bernard S, Gonze D, Čajavec B, Herzel H, Kramer A (2007)
  Synchronization-Induced Rhythmicity of Circadian Oscillators in the
  Suprachiasmatic Nucleus.
  PLOS Computational Biology 3(4): e68. doi: 10.1371/journal.pcbi.0030068
*/
class Cell_Somatic : public Cell {
 public :

  // =================================================================
  //                             Constructors
  // =================================================================
  // Ctors are protected because we either encapsulate them into factories
  // or call them directly from within the class
 protected:
  Cell_Somatic() = delete;
  Cell_Somatic(const Cell_Somatic &model);

  Cell_Somatic(CellType cellType,
           const MoveBehaviour& move_behaviour,
           const Coordinates<double>& pos,
           double initial_volume,
           double volume_min,
           double doubling_time);
  Cell_Somatic(gzFile backup_file);
 public:

  // =================================================================
  //                             Destructor
  // =================================================================
  virtual ~Cell_Somatic(void) noexcept;


  // =================================================================
  //                            Public Methods
  // =================================================================
  void InternalUpdate(const double& dt) override;
  Cell* Divide(void) override;

  void Save(gzFile backup_file) const override;
  void Load(gzFile backup_file) override;


  // =================================================================
  //                              Accessors
  // =================================================================
  double get_output(InterCellSignal signal) const override;

  bool isCycling() const override;
  bool isDividing() const override;
  bool isDying() const override;

 protected :
  // =================================================================
  //                           Protected Methods
  // =================================================================


  // =================================================================
  //                               Attributes
  // =================================================================

 public:
  // <TODO> Modify the class id and keyword </TODO>
  // Class ID, uniquely identifies this class
  static constexpr CellFormalism classId_ = 0x5ae37d71; // CRC32 Hash
  // Keyword to be used in param files
  static constexpr char classKW_[] = "SOMATIC";

 protected:
  double sigma_ = 0.0;

 private:
  /** dummy attribute - allows to register class in Simuscale statically */
  static bool registered_;
};

#endif // SIMUSCALE_CELL_SOMATIC_H__
