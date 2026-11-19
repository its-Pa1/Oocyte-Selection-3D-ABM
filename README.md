# A Multiscale Agent-Based Model for Oocyte Selection during Ovarian Development

This repository provides the computational framework, source code, and supplementary materials associated with the manuscript: 
**"A Multiscale Agent-Based Model for Oocyte Selection during Ovarian Development"** (Link to be added upon preprint/publication).

---

## Project Overview
The model integrates a stochastic gene regulatory network (GRN) with germ–germ and germ–somatic cell interactions to investigate how local cellular environments and cyst geometry influence germ cell fate determination.

It also accounts for cytoplasmic transfer, nurse-cell death, and the resulting changes in germ cell volumes and spatial organisation.

## Directory Structure
- **`Codes/`**
  - C++ source code of the multiscale model.
  - Cross-platform build instructions (CMake/Makefile).
  - Python scripts and wrappers for running simulations
    and reproducing the manuscript figures.
- **`Supplementary-Data/`**: 
  - Spatiotemporal evolution animation of cysts in 3D.
  - `Robustness_Analysis_Table.txt`: Results for $\pm 10$% perturbations in baseline parameters.

## Framework and Repositories

The model is developed using the [Simuscale framework](https://gitlab.inria.fr/bernard1/simuscale).

The project is also maintained at the following locations:

- **Inria GitLab:**
  [Simuscale – pawan branch](https://gitlab.inria.fr/bernard1/simuscale/-/tree/pawan)

- **GitHub:**
  [Ovarian Development Model](https://github.com/its-Pa1/Ovarian-Development-Model)

## Installation and Usage

1. Clone the repository:

   ```bash
   git clone https://gitlab.inria.fr/pawkumar/oocyte-selection-3d-abm.git
   ```

2. Navigate to the `Codes/` directory.

3. Follow the platform-specific installation and execution instructions provided in `Codes/README.md`.


## License and Attribution

The computational model builds upon the [Simuscale framework](https://gitlab.inria.fr/bernard1/simuscale), licensed under the GNU Affero General Public License v3.0
(AGPLv3).

The combined software is distributed under the GNU AGPLv3. See the [LICENSE](LICENSE) file for details. Original Simuscale copyright and licensing notices are
retained in `Codes/`.

Supplementary materials (simulation videos and numerical results) are distributed under the Creative Commons Attribution 4.0 International License (CC BY 4.0).
See [Supplementary-Data/LICENSE](Supplementary-Data/LICENSE).

## Citation

If you use this model or the associated simulation results, please cite:

**A Multiscale Agent-Based Model for Oocyte Selection during Ovarian Development**

*Preprint DOI coming soon.*

## Contact

For questions regarding the model, installation, or simulation data, please open an issue in this repository or contact the authors.
