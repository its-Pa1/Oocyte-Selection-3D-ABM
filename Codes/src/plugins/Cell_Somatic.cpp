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

// =================================================================
//                              Includes
// =================================================================
#include "Cell_Somatic.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "Simulation.h"
#include "Alea.h"
#include "Coordinates.h"

// =================================================================
//                    Definition of static attributes
// =================================================================
constexpr CellFormalism Cell_Somatic::classId_;
constexpr char Cell_Somatic::classKW_[];

// =================================================================
//                             Constructors
// =================================================================
Cell_Somatic::Cell_Somatic(const Cell_Somatic &model):
    Cell(model) {
}
/**
 * \brief Initialization ctor
 *
 * Create a new Cell_ODE object with randomly initialized internal
 * state and per_scal_ parameter
 */
Cell_Somatic::Cell_Somatic(CellType cellType,
                   const MoveBehaviour& move_behaviour,
                   const Coordinates<double>& pos,
                   double initial_volume,
                   double volume_min,
                   double doubling_time) :
    Cell(cellType, move_behaviour,
         pos, CellSize(initial_volume, volume_min),
         doubling_time) {}

Cell_Somatic::Cell_Somatic(gzFile backup_file) :
    Cell(backup_file) {
  Load(backup_file);
}

// =================================================================
//                             Destructor
// =================================================================
Cell_Somatic::~Cell_Somatic(void) noexcept {
}

// =================================================================
//                             Accessors
// =================================================================
double Cell_Somatic::get_output(InterCellSignal signal) const {
  switch (signal) {
    case InterCellSignal::S_nbd:
      return NAN; // Just to visualize it  
    default:
      return 0.0;
      // printf("%s:%d: error: case not implemented.\n", __FILE__, __LINE__);
      //exit(EXIT_FAILURE);
  }
}

// =================================================================
//                            Public Methods
// =================================================================
// WARNING : intrinsic inputs MUST have been computed before !
void Cell_Somatic::InternalUpdate(const double& dt) {
  // Update the intracellular state of the cell
  
}

Cell* Cell_Somatic::Divide(void) {
  Cell_Somatic* newCell = new Cell_Somatic(*this);
  SeparateDividingCells(this, newCell);
  return newCell;
}

bool Cell_Somatic::isCycling() const {
  return false;
}

bool Cell_Somatic::isDividing() const {
  return false;
}

bool Cell_Somatic::isDying() const {
  return false;
}



void Cell_Somatic::Save(gzFile backup_file) const {
  // Write my classId
  gzwrite(backup_file, &classId_, sizeof(classId_));

  // Write the generic cell stuff
  Cell::Save(backup_file);

  // This reduced SOMATIC formalism currently has no plugin-specific state.
}

void Cell_Somatic::Load(gzFile backup_file) {
  // This reduced SOMATIC formalism currently has no plugin-specific state.
}



// =================================================================
//                           Protected Methods
// =================================================================


// =================================================================
//                           Private Methods
// =================================================================
bool Cell_Somatic::registered_ =
    Cell::RegisterClass(classId_, classKW_,
                        [](CellType type,
                           const MoveBehaviour& move_behaviour,
                           const Coordinates<double>& pos,
                           double initial_volume,
                           double volume_min,
                           double doubling_time){
                          return static_cast<Cell*>(
                              new Cell_Somatic(type, move_behaviour,
                                           pos, initial_volume, volume_min,
                                           doubling_time));
                        },
                        [](gzFile backup_file){
                          return static_cast<Cell*>(new Cell_Somatic(backup_file));
                        }
    );
    
