## Supplementary Videos — Cyst Evolution

The following animations show the full temporal evolution of the germline cyst
from t = 0 to t = 240 h for three different cyst geometries.

| File | Geometry |
|---|---|
| `sphere_cyst.mp4` | Spherical Cyst |
| `cylindrical_cyst.mp4` | Cylindrical Cyst |
| `ellipsoidal_cyst.mp4` | Ellipsoidal Cyst |
| `all_cysts.mp4` | All these three cysts together |

### Color Code

The cells are colored using the **cividis colormap**:

| Hex | Color swatch | Cell type |
|---|---|---|
| `#FDE737` | ![#FDE737](https://placehold.co/15x15/FDE737/FDE737.png) | Germ cell  |
| `#4B73B0` | ![#4B73B0](https://placehold.co/15x15/4B73B0/4B73B0.png) | Somatic cell |
| `#7B3F5E` | ![#7B3F5E](https://placehold.co/15x15/7B3F5E/7B3F5E.png) | Dying cell | 

### Reading the Cell Counter

Each frame displays the simulation time and the **total cell count**
(germ + somatic). The somatic cell number is fixed at **128** throughout
the simulation. The number of surviving germ cells is therefore:


For example, a counter showing **135** means **7 germ cells** survived
alongside the 128 somatic cells. At t = 240 h, these surviving germ
cells are the oocytes.



## `Robustness_Analysis_Table.txt` — Robustness Analysis Results

This file contains the results of a one-at-a-time (OAT) robustness analysis, where each of the **17 model parameters** is independently perturbed by ±10% from its baseline value while all other parameters are held fixed.

### Experimental Design

- **Total parameter sets:** 35 — 1 baseline + 2 perturbations (−10%, +10%) × 17 parameters
- **Simulations per parameter set:** 35 independent stochastic runs
- **Readout time point:** end of simulation at **t = 240 h**

### Parameters Swept

| Symbol | Description | Baseline |
| :--- | :--- | :--- |
| $[P_n]_0$ | $P_n$ initial activation threshold | 0.400 |
| $\zeta$ | External signal coefficient | 50.000 |
| $\alpha$ | Germ cell weight | 0.030 |
| $\beta$ | Somatic cell weight | 0.200 |
| $\lambda$ | Sensing radius coefficient | 1.000 |
| $V_\text{nucleus}$ | Nucleus volume (death threshold) | 0.200 |
| $d_{0,o}$ | mRNA degradation rate (oocyte gene) | 0.100 |
| $d_1$ | Protein degradation rate | 0.020 |
| $a_0$ | Coefficient of $k_0$ | 1.000 |
| $a_1$ | Coefficient of $k_1$ | 10.000 |
| $a_2$ | Burst amplitude factor ($1/\text{mean burst size}$) | 2.500 |
| $d_{0,n}$ | mRNA degradation rate (nurse gene) | 0.100 |
| $\Phi$ | Basal activation level | −1.000 |
| $\theta_{oo}$ | Self-interaction, gene o | 5.000 |
| $\theta_{on}$ | Gene n → gene o interaction | −4.000 |
| $\theta_{no}$ | Gene o → gene n interaction | −8.000 |
| $\theta_{nn}$ | Self-interaction, gene n | 5.000 |

> **Note on composite parameters:** In the code, the effective rates governing the system dynamics are defined as $k_0 = a_0 \cdot d_{0,o}$, $k_1 = a_1 \cdot d_{0,o}$, and mean burst size $\gamma^{-1} = 1 / a_2$. 
>
> Because these effective rates share underlying base parameters, a standard ±10% OAT perturbation on the base parameters inherently subjects the system to coupled variations in the effective dynamics:
> * Perturbing $a_0$, $a_2$, or $d_1$ scales the respective effective rate independently by exactly ±10%.
> * Perturbing $d_{0,o}$ by ±10% shifts **both** $k_0$ and $k_1$ simultaneously by the exact same amount.
> * Perturbing $a_1$ forces an inverse coupling on protein production ($s_1 \propto 1/a_1$): increasing $a_1$ by 10% increases $k_1$ by 10% but decreases $s_1$ by ~9.1%. Decreasing $a_1$ by 10% decreases $k_1$ by 10% but increases $s_1$ by ~11.1%.
> 
> Therefore, this OAT design successfully tests not only isolated parameter shifts, but also these inherent simultaneous and inverse couplings within the system's effective rates.

### Output Columns

| Column | Description |
|---|---|
| `Mean_Nurse` | Mean number of nurse cells at t = 240 h |
| `CV_Nurse` | Coefficient of variation of nurse cell count |
| `Mean_Germ` | Mean total germ cell count at t = 240 h |
| `CV_Germ` | Coefficient of variation of total germ cells |
| `Mean_Oocyte` | Mean oocyte count at t = 240 h |
| `CV_Oocyte` | Coefficient of variation of oocyte count |
| `N/G_Ratio` | Nurse-to-total-germ cell ratio |

### Notes

- Each row corresponds to one unique parameter set (one perturbed parameter, rest at baseline)
- Row 1 is the **baseline** with all parameters at their nominal values
- The file is tab-separated with fixed-width padding for human readability