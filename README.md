##  Water and Hydrate by TIP5P Molecular Dynamics Simulation ##

We are studying the microwave heating of water and ice by means of the molecular dynamics TIP5P model (with if_xyz2=.true.), and that of the ice and methane hydrate (with if_xyz1=.true.) by the same code.

This is the second generation code against the first one of 3-point SPC/E shake/rattle algorithm 
to study the microwave heating of water and ice. It was shown clearly in the first generation code that 
the ice below 273 K becomes a crystal and is not melted against the strong electric field of microwaves 
(JCP, 2007, Ref. 1).
The new results by the 5-body rotation coordinate system and TIP5P-Ewald model are shown below 
and the PDF in the arXiv Library: https://arxiv.org/abs/2311.01182 (2023, 2024).

The five-point molecules for water are used known as the TIP5P-Ewald model. 
A four-water molecule is specified to calculate positive two hydrogens q_H= 0.241e 
and negative two hydrogens q_L= -0.241e with e the electron charge. 
The fifth oxygen atom called the dummy site q_O=0 is to correlate with adjacent molecules 
using the Lennard-Jones potential Psi(r)= eps_A/r^12 -eps_B /r^6 (Ref. 2).  

The fortran code with MPI is given in the file @p3mtip5p07a.f03 with additional files as param_tip5D07a.h, TIP507_config.start0, and initial coordinates 1cx666a.exyz and quaternions 1cx666a.q. Its description of the code is shown at README.md and also PDF documents of this repository. The latter documents are "Water_TIP5P_Simulation.pdf" for numerical coding, and for the physics simulation run by "Water_and_hydrate_molecules_by_TIP5P_code.pdf". 
The freezing ice state by microwaves, which is our theory discovery in JCP 2007 mentioned 
above, remains basically the same due to the structure of six-membered ice ! 

Methane hydrate is simulated by switching to if_xyz1=.true. of the TIP5P code, 
like with the SPC/E code (Ref. 3). We need the initial coordinate file mh3.exyz and quaternions mh3.q. 
It runs up to t= 1.7 x 10^(-8) s with E_x = 3. x 10^7 V/cm, which is terminated suddenly due to 
collapse of methane hydrate. Energies and scatter plots of molecules, the distribution of cosine's 
in the x-direction are shown in the PDF file "MDTip5pWater-8.pdf", the new version at https://arxiv.org/abs/2311.01182v8 (2024).

As natural resources, however, the methanes are 30 times environmentally more hazardous 
materials than carbon dioxide. The methane hydrate that will be mined and then burnt in air 
should be confined back to the deep interior of the earth.

### Procedure of Water Molecules by Five-Points Model ###

The important points of the five-body molecules are summarized.
(LaTeX is used for this and the next section.)

A. Five sites are one oxygen (O), hydrogen of H1 and H2 (H), and 
negative hydrogen virtual sites of L1 and L2 (L). 
Their charges are $0$, $0.241 e$,  $0.241 e$, $-0.241 e$ and $-0.241 e$, respectively. The L1 and L2 are called the dummy sites.

B. Separate the translational part $\bm{R}_{j}$, $\bm{V}_{j}$ for molecules ($j= 1,N/5$), where the rotation part $\bm{r}_{i}= (x_{i}, y_{i}, z_{i})$  ($i= 1,N$) for atoms is used  for the five sites.
The separation is done at the starting step only;
once determined at $t=0$, they become constant in time.

C. The half time step for the molecules is first executed for a predictor step, and
the full time step is made for advance of time.

D. Before the end of the step, the forces are calculated by atom positions, where the dummy L sites are calculated by algebraic vector operation.

E. After correction of quaternions, the kinetic and Coulombic energies are calculated, 
and go to the beginning of the cycle. The leap-frog method is used for the plasmas and water. 

### Each Step of Molecular Dynamics Simulation Code ###

Each step of the cycle corresponds to (i) translational motion (Step 1), 
(ii) rotational motion (Step 2-4), and (iii) addition of the fields (Step 5-8).

0) Read positions $(x,y,z)$ and quaternions from files by 
$ read(17) x_{i}, y_{i}, z_{i} (i=1,N), where the dummy sites are obtained 
by algebraic operation, and $ read(30) e_{0j}, e_{1j}, e_{2j}, e_{3j} (j=1,N/5).
The initial configuration of water and ice "GenIce" is constracted by
Dr. M. Matsumote (Ref. 4). This step is executed only at the first time step. 

2) The position $\bm{R}_{j}$ and velocity $\bm{V}_{j} $ of each molecule ($j=1,N/5$) are advanced by summation over five sites of forces $\bm{F}_{k}$ for the translational motion
$(k=1,N)$,
\begin{equation}
d\bm{V}_{j}/dt= (1/m_{j})\sum_{k=1}^{5} \bm{F}_{k}, \ \ 
d\bm{R}_{j}/dt= \bm{V}_{j}.
\end{equation}

3) Next 2)-5) steps are made for half a time steps $\Delta t_{1}=\Delta t/2$ by prediction, 
and then for a full time step $\Delta t_{2}= \Delta t$ by correction.
The angular momentum of rotational motion is calculated by summation over the torque of five sites at a time step $\Delta t_{1}$ or $\Delta t_{2}$,
\begin{equation}
 d\bm{L}_{j}/dt_{n}= \sum_{k=1}^{5}(y_{k}F_{k}^{z} -z_{k}F_{k}^{y},
z_{k}F_{k}^{x} -x_{k}F_{k}^{z}, x_{k}F_{k}^{y} -y_{k}F_{k}^{x})
\end{equation}

4) The angular frequency $\omega_{j,\alpha}$ is connected to the angular mementum 
and inertia moment $Im_{j,\alpha}$ with $\alpha= x,y,z$ and the matrix 
$A_{\alpha,\beta}$ by,
\begin{equation} 
\omega_{j,\alpha}= (A_{\alpha 1}L_{x} +A_{\alpha 2}L_{y} +A_{\alpha 3}L_{z})/
Im_{j,\alpha}
\end{equation}

\begin{equation}
\begin{matrix}
A_{11}= e_{0}^2 +e_{1}^2 -e_{2}^2 -e_{3}^2, &
A_{12}= 2(e_{1}e_{2} +e_{0}e_{3}),&
A_{13}= 2(e_{1}e_{3} -e_{0}e_{2}),\\
%
A_{21}= 2(e_{1}e_{2} -e_{0}e_{3}),&
A_{22}= e_{0}^2 -e_{1}^2 +e_{2}^2 -e_{3}^2, &
A_{23}= 2(e_{2}e_{3} +e_{0}e_{1}),\\
%
A_{31}= 2(e_{1}e_{3} +e_{0}e_{2}),&
A_{32}= 2(e_{2}e_{3} -e_{0}e_{1}),&
A_{33}= e_{0}^2 -e_{1}^2 -e_{2}^2 +e_{3}^2, 
\end{matrix}
\end{equation}

4) The time derivative of quaternion $\bm{q}_{j}$ is given by the angular frequency by,
\begin{equation} 
d\bm{q}_{j}/dt_{n}= (1/2) \Delta t_{n} 
\begin{pmatrix}
-e_{1} \omega_x -e_{2} \omega_y -e_{3} \omega_z \\
e_{0} \omega_x -e_{3} \omega_y +e_{2} \omega_z \\
e_{3} \omega_x +e_{0} \omega_y -e_{1} \omega_z \\
-e_{2} \omega_x +e_{1} \omega_y +e_{0} \omega_z
\end{pmatrix}
\end{equation}

5) Get a rotational matrix $ A_{\alpha \beta}(e_{0},e_{1},e_{2},e_{3})$
in half a time steps for prediction and go back to Step 2.
In the correction step it is made for a full time step and go to Step 6.

6) The three sites $\bm{r}_{i}$ and the position $\bm{R}_{j}$ are connected by,
\begin{equation}
  \bm{r}_{i}= \bm{R}_{j} +
     \begin{pmatrix}
     A_{11}  & A_{21}  & A_{31} \\
     A_{12}  & A_{22}  & A_{32} \\
     A_{13}  & A_{23}  & A_{33} % \nonumber
     \end{pmatrix}
     \begin{pmatrix}    
     x_{i} \\ y_{i} \\ z_{i} % \nonumber
     \end{pmatrix}
\end{equation}
The position of dummy sites are calulated from known three sites by algebraic operation.

7) The forces of positions are calculated from Coulombic and Lennard-Jones potentials using the five sites. 

8) Correction to a normalization of quaternions is made at every 10 steps. Go to the new time step Step 1.  

Note that a time step is important and it will be $\Delta t= 0.025 \tau$ or less by the time advancing scheme; otherwise the code will be inaccurate or go to overflow. 

The equations of A_{ij} and e0(i),e1(i),e2(i),e3(i) are written in the PDF file, "Water_TIP5P_Simulation.pdf". Checks of equations of the TIP5P code are also shown.

### The Lennard-Jones Potential ###

With the Coulombic interactions, the Lennard-Jones 12-6 potential is adopted,
A=3.85x10^(-8) erg Ang^(12), B=4.36x10^(-11) erg Ang^(6) for TIP5P-Ewald sum's case, 
and A=4.17x10^(-8) erg Ang^(12), B=4.24x10^(-11) erg Ang^(6) for TIP4P.

Other parameters are: r(OH)=0.9572 Ang, Delta(HOH)=104.52 Ang. r(OM)=0.15 Ang is used 
for the TIP4P case, where the equipartition line of the virtual M sites is on the plane 
that equally separates the HOH angle. 

### In Order to Start a Run ###

To start a simulation of water cluster with the TIP5P code, the adjacent 4x4 hydrogen pairs of molecules 
are summed electrostatically, and also the oxygen pairs are coupled by TIP5P Lennard-Jones potentials (Ref. 5).

1. To get an initial state, we make the size of at least a 6x6x6 water cluster for numerical stability.
Short-range Coulombic and long-range forces are best separated for interactions, where the short-range 
forces for r>r_LJ ranges are made to be spatially dumped.

2. Around a given temperature, a dryrun is executed at least for 5 periods that is 10^(-9) seconds.
The long dryrun of a simulation is very important !!

3. Then, we may be ready to apply the electric field E_x= E_0*sin(\omega *time) in the x-direction to excite the 
electric dipole interactions of water. For the moment, we give the electric field 10 GHz where 
the electric field E_0 and electric dipole p_0 are of the order of 5x10^(-3) eV.
The E_0 value can be as small as 5x10^(-4) eV for the 6x6x6 water cluster.

### To Obtain the Initial Equilibrium for 298 K ###

To give uniform random noises to molecules, an equilibration takes extremely a long time of periods.
Thus, we try to use salt ions where 6x6x6 water clusters are fully strengthend with 64 Na(+) and 64 Cl(-) ions 
at t=1000, and waited by t=1700. Then, the salt ions are gradually removed by t=4700, and only water 
molecules are continued at least for 5 periods of t=50,000 before the real simulation run is again started from t=0.

The simulation run is shown in the next section in the PDF document with color pictures, 
"Molecular_Dynamics_Simulation_of_Water_and_Ice_by_TIP5P_Code.pdf".

It is noted that we impose the NVE ensemble simulation so that the box size is constant.
In the water heating of our daily experiment, however, the system is open and an energy goes out 
to the surrounding system. It is possible to have rapid heating of hot bubbling water.

### Figures of TIP5P-Model Simulations ###

Molecular dynamics simulations of water and ice, and/or methane hydrate, are shown by using the TIP5P-Ewald model, 
and the current file is "Water_and_hydrate_by_molecular_dynamics_TIP5P_code.pdf" in this repository (Ref. 5).

Figures in the color PDF file are water and ice and/or methane hydrate that include the energy of molecules, 
3-dimensional scatter plots, pair distribtion functions of oxygen and hydrogen atoms, and 
cosine distribution functions in the x-direction.

Generally speaking, water is heated by microwaves where its heating efficiency is highest at 0 Celsius.
Below 273 K, however, ice is frozen and not heated by the electric field.
Methane hydrate that is excited by microwaves application is heated and eventually collapsed. 

### Necessary Files of the Repository ###

Files for simulation runs: 

!*   1. @p3mtip5p07a.f03 : simulation code                       

!*   2. param_tip5p_D07a.h : parameter file 1, be constants      

!*   3. TIP507_config.start0 : parameter file 2, kstart=0 or 2   

!*      and/or TIP507_config.start1 : kstart=1 or 3              

!*   4. Initial files of water and ice: 1c666a.exyz and 1c666a.q

Post-processing files:                                

!*   * @lplotip507.f03 - the final history of energy             

!*   * @dipol_seqtip507.f03 - dipole Ex field                    

!*   * @iceplotip507.f03  - for 3D plot for x,y,z                

!*   * @wat_radtip507.f03 - radial distribution functions        
                                                         
References of Numerical Technique:

!* 1) M.Tanaka, J.Comput.Phys., vol. 79, 206 (1988).   

!* 2) M.Tanaka, J.Comput.Phys., vol.107, 124 (1993).   

!* 3) M.Tanaka, Comput.Phys.Comm., vol.87, 117 (1995). 

!*  4) M.Tanaka, Comput.Phys.Comm., vol.241, 56 (2019). 


### References ### 
 
1. M. Tanaka and M. Sato, Microwave heating of water, ice and saline solution: Molecular dynamics study, J.Chem.Phys., 126, 034509 (2007).

2. "Classical Mechanics", H. Goldstein, C. Poolee, J. Safko, 3rd Edition, Pearson Education Inc., England (2003); 古典力学，吉岡書店 (2006).

3. M. Tanaka, M. Sato, S. Nakatani, Microwave heating and collapse of methane hydrate by molecular dynamics simulations, https://arxiv.org/1909.01024, Cornell University, USA (2019).

4. M. Matsumoto, T. Yagasaki, and H. Tanaka,"GenIce: Hydrogen-Disordered Ice
Generator", J. Comput. Chem. 39, 61-64 (2017).

5. M. Tanaka, Molecular dynamics simulation of methane hydrate by means of the TIP5P-Ewald model,
https://arxiv.org/2311.01182, Cornell University Library, USA (2023, 2024).


