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

// ============================================================================
//                                   Includes
// ============================================================================
#include "Cell_Boilerplate.h"

// ============================================================================
//                       Definition of static attributes
// ============================================================================
constexpr CellFormalism Cell_Boilerplate::classId_;
constexpr char Cell_Boilerplate::classKW_[];

// ============================================================================
//                                Constructors
// ============================================================================
/**
 * Create a new object with the provided constitutive elements and values
 */
Cell_Boilerplate::Cell_Boilerplate(CellType cellType,
                   const MoveBehaviour& move_behaviour,
                   const Coordinates<double>& pos,
                   double initial_volume,
                   double volume_min,
                   double doubling_time) :
    Cell(cellType, move_behaviour,
         pos, CellSize(initial_volume, volume_min),
         doubling_time) {
}

/**
 * Restore an object that was backed-up in backup_file
 */
Cell_Boilerplate::Cell_Boilerplate(gzFile backup_file) :
    Cell(backup_file) {
  Load(backup_file);
}

// ============================================================================
//                                 Destructor
// ============================================================================
Cell_Boilerplate::~Cell_Boilerplate() = default;

// ============================================================================
//                                   Methods
// ============================================================================
void Cell_Boilerplate::InternalUpdate(const double& dt) {
  // <TODO>
  // Describe what happens inside the cell when taking a time step
  // </TODO>

  (void)dt; // to avoid compiler warning: unused parameter 'dt'

  RBP1_ = getInSignal(InterCellSignal::BP1);
  std::vector<real_type> v = get_diffusive_signal(InterCellSignal::BP1);
  RBP1_ += v[0];

  // Default dying behaviour
  // If is_dying_ is set, then cell becomes dead
  if ( is_dying_ ) {
    set_to_dead();
    return; // it should die immediately
  }

  // Default dying behaviour
  // If isDying evals to true, then set is_dying_
  // isDying should be called only once.
  if ( isDying() ) {
    set_to_dying();
    return; // if dying, return immediately
  }
}

Cell* Cell_Boilerplate::Divide(void) {
  // <TODO> Determine how the cell should divide </TODO>
  return nullptr;
}

void Cell_Boilerplate::Save(gzFile backup_file) const {
  // Write my classId
  gzwrite(backup_file, &classId_, sizeof(classId_));

  // Write the generic cell stuff
  Cell::Save(backup_file);

  // Write the specifics of this class
  // <TODO> Save your specific attributes </TODO>
}

void Cell_Boilerplate::Load(gzFile backup_file) {
  // <TODO> Load your specific attributes </TODO>
  (void)backup_file; // to avoid compiler warning: unused parameter 'backup_file'
}

double Cell_Boilerplate::get_output(InterCellSignal signal) const {
  // <TODO> Determine the amount of <signal> the cell is secreting </TODO>
  // example: cell expresses BP1
  if ( signal == InterCellSignal::BP1 ) {
    return .1;
  }
  else if ( signal == InterCellSignal::RBP1 ) {
    return RBP1_;
  }
  else {
    return 0.0;
  }
}

bool Cell_Boilerplate::isCycling() const {
  // <TODO> Is the cell currently cycling ? </TODO>
  return false;
}

bool Cell_Boilerplate::isDividing() const {
  // <TODO> Should the cell divide now ? </TODO>
  return false;
}

bool Cell_Boilerplate::isDying() const {
  // <TODO> Should the cell die now ? </TODO>
  return false;
}

bool Cell_Boilerplate::isDead() const {
  // Should the cell die ?
  return is_dead_;
}

void Cell_Boilerplate::UpdateCyclingStatus() {
}

bool Cell_Boilerplate::StopCycling() {
  // <TODO> Should the cell stop cycling now ? </TODO>
  return false;
}

bool Cell_Boilerplate::StartCycling() {
  // <TODO> Should the cell start cycling now ? </TODO>
  return false;
}


Coordinates<double> Cell_Boilerplate::MotileDisplacement(const double& dt) {
  // <TODO> Motile cell autonomous displacement behaviour
  Coordinates<double> displ { 0.0, 0.0, 0.0 };
  (void)dt; // to avoid compiler warning: unused parameter 'dt'
  return displ;
}


double Cell_Boilerplate::gaussian_field_weight(InterCellSignal signal) {
  // <TODO> set the strength of the diffusive signals
  (void)signal; 
  // example: cell expresses diffusive signal BP1
  if ( signal == InterCellSignal::BP1 ) {
    return 1.5;
  }
  else {
    return 0.0;
  }
}

const Coordinates<double>& Cell_Boilerplate::gaussian_field_source(InterCellSignal signal) {
  // <TODO> set the position of the source of the diffusive signals
  (void)signal; 
  return pos();  // default: cell position

}

std::vector<Coordinates<double>> Cell_Boilerplate::gaussian_field_targets(InterCellSignal signal) {
  // <TODO> set the vector of positions for sensing diffusive signals 
  (void)signal;
  std::vector<Coordinates<double>> targets;
  targets.push_back(pos());                   // default: cell position and
  targets.push_back(pos() + orientation());   // a position along the cell orientation
  return targets;
}


// Register this class in Cell
bool Cell_Boilerplate::registered_ =
    Cell::RegisterClass(classId_, classKW_,
                        [](CellType type,
                           const MoveBehaviour& move_behaviour,
                           const Coordinates<double>& pos,
                           double initial_volume,
                           double volume_min,
                           double doubling_time){
                          return static_cast<Cell*>(
                              new Cell_Boilerplate(type, move_behaviour,
                                           pos, initial_volume, volume_min,
                                           doubling_time));
                        },
                        [](gzFile backup_file){
                          return static_cast<Cell*>(new Cell_Boilerplate(backup_file));
                        }
    );
