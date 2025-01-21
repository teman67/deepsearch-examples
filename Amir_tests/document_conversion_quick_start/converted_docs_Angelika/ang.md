## ScienceDirect

Comput. Methods Appl. Mech. Engrg. 354 (2019) 689-705

www.elsevier.com/locate/cma

## Design concepts for generating optimised lattice structures aligned with strain trajectories

Stephen Daynes$^{a}$, Stefanie Feih$^{$_{a,}$∗}$, Wen Feng Lu$^{b}$, Jun Wei $_{a,}$∗

$^{a}$Singapore Institute of Manufacturing Technology, 2 Fusionopolis Way, 138634, Singapore $^{b}$Department of Mechanical Engineering, National University of Singapore, 117975, Singapore

Received 15 January 2019; received in revised form 30 May 2019; accepted 31 May 2019 Available online 14 June 2019

## Highlights

· A novel lattice infill methodology creates optimally aligned lattice structures.

· Lattice structures are automatically generated using principal strain trajectories.

· The lattices conform to arbitrary 3D design space volumes.

· Optimally aligned lattice structures have higher stiffness and strength.

· Additively manufactured lattices validate the optimisation predictions.

## Graphical Abstract

## Abstract

Additively manufactured lattice structures enable the realisation of light-weight, multi-functional, structures. For example, lattices can be used for high stiffness and buckling resistance in sandwich structures or as support material for additive manufacturing. Topology optimisation and additive manufacturing are two technologies that allow the design, optimisation and manufacture of complex lattice designs. In this work, a new lattice optimisation methodology is presented that tailors the

0045-7825/ c ⃝ 2019 Elsevier B.V. All rights reserved.

size, shape and orientation of individual lattice trusses in three-dimensional space by using principal strain fields obtained from topology optimisation. This new method of generating functionally graded lattices is shown both numerically and experimentally to be capable of generating lattice structures with greatly improved stiffness and strength when compared to lattice structures with a uniform lattice infill. Upper and lower relative density thresholds and minimum truss member sizes are included in the optimisation workflow to ensure that the optimised lattice designs are compatible with additive manufacturing process constraints. The functional grading method is also shown to be capable of generating conformal lattice structures in three dimensions, even for complex loading conditions and arbitrary volume boundaries. c ⃝ 2019 Elsevier B.V. All rights reserved.

Keywords: Lattice structures; Functional grading; Infill; Topology optimisation; Additive manufacture

## 1. Introduction

Additive manufacturing enables the generation of complex lattice structures for a variety of applications, often when a light-weight core is required to withstand external loads. It has been demonstrated that cellular structures are particularly well suited to multi-functional applications where cell topology can be tailored to provide sufficient stiffness and strength while possessing other beneficial characteristics [1,2]. Multifunctional applications include buckling resistant sandwich panel cores [3,4], controlling fluid flows within a structure [5], absorbing impact energy [6,7], heat transfer or insulation [8,9], acoustic damping [10], and surface texturing for biomedical implants [11,12]. Such cellular structures can be found throughout the natural world, with examples including bamboo [13], sea sponges [14,15] and the internal cores of bones [16].

The most common cell arrangements are periodic, such as Kagome lattice cores [17] and triply periodic minimal surfaces [18], or stochastic, such as foams where the cellular structure is randomly distributed [19]. Functionally graded versions of these topologies are also possible where there is a variation in relative density and hence tailored properties through the structure. Periodic truss-type lattice structures with a relative density variation through the thickness have been shown to have gradual crushing behaviour with the sequential collapse of the graded layers [20]. Both periodic and stochastic functionally graded structures have also been investigated for orthopaedic hip implants [21,22] where there is a requirement to match the stiffness and strength of bone and also to tailor the surface texture to promote bone growth. Additive manufacturing technologies are particularly well suited to the realisation of such complex cellular structures [23,24]. When compared to traditional manufacturing techniques, the available design space for additively manufactured parts is largely free from manufacturing constraints or limitations on design complexity.

Computational methods have been developed to address the difficult challenge of generating functionally graded lattice structures. Such approaches include a voxel-based method [25] and a 'space warping' method to compress cell size in high stress areas [26]. A dithering approach combined with Voronoi tessellations has also been proposed to generate graded stochastic patterns [27,28]. The potential to generate highly complex parts also enables extensive use of topology optimisation (TO) as part of the design process [29,30]. With TO, an optimal relative density distribution can be converted into a cellular structure with a tailored relative density varying between cells [31,32]. TO approaches employing cellular structures often utilise homogenisation approaches to simplify unit cell representation for modelling of the structural performance at macro-level [33]. A different approach is the multi-scale concurrent design method where cell level details are optimised concurrently with their macroscopic distribution [34,35]. Cell orientation has also been identified as an important factor to account for during the optimisation of lattice structures. Tang and Zhao developed a method to optimise the distribution of lattice orientations in a discrete number of 'sub functional volumes' [36]. Reinhart et al. developed a different approach by arranging lattice trusses to align with 'flux of force' within a solid [37] and showed reduced stress concentrations occurred at truss connections when compared with a uniform lattice design. The study of such optimal truss arrangements can be traced back to Michell, who defined the minimum mass truss structure to sustain a given force [38].

This paper addresses outstanding questions regarding best practice for arranging lattice cells for maximum stiffness (minimum compliance) in three-dimensional space including printing constraints, as well as benchmarking of the resulting mechanical lattice properties against solid topology designs of the same mass. We begin with an outline of the general state-of-the-art lattice optimisation approach in commercial finite element analysis (FEA) software. FEA is capable of generating functionally graded lattice structures, where the diameter of each lattice

member is uniquely sized while maintaining the same cell size. This approach is then compared with our in-house developed method that can additionally tailor the lattice cell size, shape and orientation based on optimal strain distributions in three-dimensional space. Additionally, we highlight the importance of incorporating lower and upper density threshold values to separate void and solid regions from the lattice structure. These density thresholds are related in this work to respective additive manufacturing process constraints. Two engineering case studies are presented to demonstrate that functionally graded lattices with the additional design freedom of lattice arrangement can lead to significantly higher stiffness.

## 2. Numerical and experimental methodology

This section contains five subsections. The first subsection introduces the TO problem description and the procedure to generate optimal relative density distributions. The second subsection introduces lattice manufacturing constraints in the form of lower and upper relative density bounds. This is followed by details of the method for generating lattices aligned with principal strain trajectories. These graded lattices fill regions of intermediate relative density. The fourth subsection gives details of parametric (size and shape) optimisation to 'fine tune' the lattice design. The last subsection describes the design validation in this work as undertaken by Selective Laser Sintering (SLS). It is noted that the presented consideration of the manufacturing limits in Section 2.5 applies to most AM processes.

## 2.1. Optimisation problem description

We use Altair's OptiStruct 2017 in-built TO functions as a basis for our work. The software uses the Solid Isotropic Material Penalisation (SIMP) method [39] whereby the optimum material distribution is parameterised by relative density. This method was chosen because SIMP enables both a discrete solid-void representation of optimal material distribution and also a greyscale representation of the optimal relative density where lattice regions can be introduced. The linear static formulation for the relative density-based TO problem is given in Eq. (1):

subject to :

K (ρ) U = F (ρ)

(1b)

0 ≤ ρ ≤ 1

(1d)

Here c is the objective function, ρ is a vector containing the relative density design variables, the displacement vector is given by U and the global stiffness matrix for the FEA problem is K. The objective is to minimise structural compliance, c, which is measured by elastic strain energy. Compliance can be considered a reciprocal measure of the stiffness. The displacement vector U is found by solving the static elasticity equation (Eq. (1b)). To satisfy equilibrium, the stiffness matrix and the displacement matrix are related to the force vector F. The constraint function g in Eq. (1c) links the material volume fraction in the design space V and the allowable volume fraction upper bound V$_{f}$.

A power law penalisation is applied to the stiffness-relative density relationship. With this power greater than unity, the law will penalise the formation of intermediate relative density regions since the relative stiffness of the constituent solid material will typically be higher than that of a lattice [40]. In this work the relationship between solid and lattice stiffness is based on the following semi-empirical expression [41-43]:

The stiffness of the solid material is given by K$_{s}$. The penalty factor p can be specified prior to optimisation. A penalty of p = 4 is used for solid topology designs, resulting in a mostly binary result of void and solid regions. No penalty can be applied (p = 1), although this idealised state is hard to achieve with a lattice structure, where a proportion of the trusses will not be aligned with the applied loads, which in turn introduces a degree of structural inefficiency. An exponent of p = 1. 8 is characteristic of the open cell lattice structures used in this work [19]. Such a penalty factor typically leads to an optimised design with a mixture of fully dense and void regions with some

Fig. 1. Simply supported floor beam with a uniformly distributed pressure load (dimensions in mm). (For interpretation of the references to colour in this figure legend, the reader is referred to the web version of this article.)

Fig. 2. SIMP results before (left) and after (right) relative density thresholds are applied, using a penalty of (a) p = 1, (b) p = 1. 25 and (c) p = 1. 8. Volumes exceeding upper and lower relative density bounds are shown in black and white respectively. (For interpretation of the references to colour in this figure legend, the reader is referred to the web version of this article.)

intermediate relative density regions. A less conservative penalty of p = 1. 25 is also available in OptiStruct and will be evaluated later in the paper.

The methodology is presented in this section using our first case study of a simply supported floor beam, see Fig. 1, which is subject to a design space volume fraction constraint, V$_{f}$, of 1/3. The beam has a 1 mm thick upper skin (non-design space) that is subject to a uniformly distributed pressure load. The design space (grey region) is meshed with 2.5 mm second order hexahedral elements prior to TO. The upper 1 mm thick floor is modelled with 2.5 mm second order quadrilateral shell elements. A minimum thickness constraint of 10 mm is set during TO, allowing the generation of two complete lattice cells.

The optimal relative density distributions in the beam using the three penalty factors of p = 1, 1.25 and 1.8, are shown on the left hand side of Fig. 2. It can be seen that a larger proportion of intermediate relative density regions is generated when reducing penalty factors.

## 2.2. Lattice relative density bounds

Even with recent advances in additive manufacturing, it is impractical to manufacture a lattice structure that can approach either limits of 0% or 100% relative density. Low densities are difficult to achieve since any additive manufacturing process will have a minimum print resolution, which in turn will determine the minimum feasible truss diameter. Such minimum feature sizes are a result of melt pool and powder characteristics in the case of powder-based printing processes. In this work we set a minimum relative density bound of 10%, which

Fig. 3. Cubic unit cell with (a) 10% relative density showing minimum printable feature size of 0.5 mm radius (1 mm diameter) and (b) 60% relative density showing trusses in close proximity which will result in powder entrapment and larger intersecting volumes of neighbouring trusses.

approximately corresponds to the relative density of a cubic cell with a 5 mm edge length and a minimum truss diameter of 1 mm, see Fig. 3a. It will be shown in the following section that the lattice generation method in this work creates three orthogonal sets of lattice trusses. Therefore, using a cubic unit cell, which is also composed of three orthogonal sets of lattice trusses, is a reasonable first approximation to evaluate relative density bounds for manufacturing feasibility. After the TO step in Section 2.1, this lower relative density bound is applied by removing elements with low densities and re-meshing the remaining volume. After the application of relative density bounds the remaining solid regions are re-meshed to have 2.5 mm second order tetrahedron elements. The TO is then performed for a second time, since the reduction of design volume-when the lower bound is applied-has some influence on the optimised relative density distribution. The lower bounds are shown by the white regions on the right hand side of Fig. 2.

Similarly, high relative densities of lattice regions can also be difficult to manufacture due to print resolution limitations, see Fig. 3b. The open cell nature of truss-like lattice structures have been found to be useful for avoiding powder entrapment, when compared with using wall-like infill designs [44]. However, printing trusses in close proximity can still result in the formation of internal cavities, which can trap supporting powder and/or trusses merging together due to thermal conduction occurring between close neighbouring trusses. For these reasons, an upper relative density bound of 60% is used in this work, as shown by the black regions on the right hand side of Fig. 2. Lattice regions with approximately 60% relative density will have internal pore features less than 2 mm across.

An additional advantage of applying an upper density bound is that the SIMP penalty factor approach as per Equation (2), used for approximating relative lattice stiffness in terms of relative density, becomes less accurate with increasing relative density. This is caused by the intersecting volumes of neighbouring trusses becoming larger as the truss diameters increase relative to the unit cell size [45], also indicated in Fig. 3b. Modelling the stiffness for these higher relative densities with increased accuracy is possible by considering additional terms in the penalty factor approximation [46]; however, this approach is currently not available in OptiStruct.

## 2.3. Lattice generation based on principal strain trajectories

Usually with TO a threshold relative density is applied to the optimised relative density distribution and any elements falling below this value are deleted and any elements above this value are converted to solids, resulting in a binary void-solid solution. However, in this work these intermediate relative densities are represented by lattice trusses in the form of beam elements. This solution can be considered as a greyscale interpretation of the idealised optimal relative density distribution. The conventional approach is to populate these intermediate relative density regions with a repeated uniform array of lattice unit cells. The research presented in this work is an advancement

on this state-of-the-art by generating graded lattice infill arrangements. Lattice truss spacing, length and orientation become additional design variables. This in-house developed approach uses the optimised strain data to generate an interconnected network of trusses aligned with principal strain trajectories [47,48]. After TO, the strain tensor is available for each node in the solid FEA model:

Each row describes three strains acting on a plane that is aligned with one of the global coordinate axes. This data enables the three principal strain magnitudes to be determined (ε $_{j}$, j = 1, 2, 3) in addition to their orientation [49,50]. This is done in Matlab by calculating the eigenvalues and eigenvectors respectively. Three orthogonal, normalised, eigenvectors are calculated (n$_{j}$, j = 1, 2, 3), each one corresponding to a principal strain direction. Matlab is also used to trace the resulting strain trajectories from initial starting seed locations a$_{0}$. This is done using a fourth-order Runge-Kutta scheme [51]. A step size ± ∆ s is used to trace along the principal strain trajectories, set at 2 mm in this work.

In Eq. (4), a$_{i}$$_{+}$$_{1}$ is the next point along the trajectory to be calculated and n$_{j}$ ⏐ a$_{i}$ is a principal strain direction at point a$_{i}$. The seed locations, shown by the red markers on the left hand side of Fig. 4, are equally spaced at 2.5 mm intervals within the design volume. Principal strain trajectories are traced, from each seed point, until the boundary of the lattice region is reached. If the seed point lies outside the intermediate relative density volume then no tracing of trajectories is initiated.

A close-up view of a lattice region for the beam with the penalty factor p = 1 is shown in Fig. 5. The close-up lattice region shown has dimensions of 20 mm × 20 mm × 5 mm, with the shorter dimension in the throughthickness direction. The three strain trajectories emanating from a single seed point are shown in Fig. 5a. At the point of intersection all three trajectories are orthogonal. All of the trajectories from all seed points are shown in Fig. 5b. Principal strain trajectories generated for each model are shown on the left hand side of Fig. 4. Note that one third of the trajectories generated are perpendicular to the page for this plane stress loading condition. Comparing Figs. 2 and 4, it can be seen that the principal strain trajectories have a tendency to move closer in areas of higher relative density.

The next step is determining how to appropriately space and join the strain trajectories to create an interconnected lattice structure. Here, an average spacing criterion is implemented whereby the distances and relative orientation of neighbouring trajectories are compared, see Fig. 5c. If the average distance is less than a prescribed threshold value (5 mm in this work) and the trajectories are aligned in the same direction, then the trajectory in question is deleted. This deletion process is done in order of descending trajectory length. Longer principal strain trajectories are preferable in the final lattice design, since they are more likely to form continuous load paths and well formed, conformal, void regions.

Next, the remaining principal strain trajectories are connected together, see Fig. 5c. The trajectories are colour coded in red, green and blue in Fig. 5 to represent their relative orientation. If two principal strain trajectories of different orientations are in close proximity to one another, then they are connected by an additional line. In this work the threshold distance is set to 2.5 mm, corresponding to half the average trajectory spacing criterion. These additional lines, shown as bold black lines in Fig. 5c, are orthogonal to both strain trajectories that they connect to. The orientations of these connecting lines are determined by taking the cross product of the two strain trajectories at

Fig. 4. Seed locations and all principal strain trajectories generated (left). Remaining trusses after deletion and connection steps (right). Penalty factors of (a) p = 1, (b) p = 1. 25 and (c) p = 1. 8. (For interpretation of the references to colour in this figure legend, the reader is referred to the web version of this article.)

Fig. 5. Principal strain trajectory generation in 3D space under plane stress conditions; (a) three orthogonal trajectories from a single seed point, (b) all trajectories from all seed points and (c) trajectory deletion by enforcing average distance criterion followed by connecting remaining trajectories that are orthogonal and in close proximity. (For interpretation of the references to colour in this figure legend, the reader is referred to the web version of this article.)

their closest points. Finally, size and shape (parametric) design variables are assigned to the lattice nodes in Matlab before being exported back into OptiStruct in text file format (*.fem). It is noted that this process has been fully automated for arbitrary 3D design volumes.

Fig. 6. Size optimisation with shape variables highlighted in yellow (left). Size and shape optimisation (right). Penalty factors of (a) p = 1, (b) p = 1. 25 and (c) p = 1. 8. (For interpretation of the references to colour in this figure legend, the reader is referred to the web version of this article.)

## 2.4. Size and shape (parametric) optimisation

The average spacing approach presented in Section 2.3 is a convenient means to generate a lattice structure that is limited by the printable feature constraints shown Fig. 3. However, average distance spacing does not guarantee an optimal load re-distribution between trusses. Further improvements in performance can be achieved by fine tuning the lattice size and shape parameters. Each lattice member is modelled as a tapered truss element. The diameter at each end of the truss then forms a unique nodal design variable that can be size optimised, ensuring that connecting nodes have the same diameter. During size optimisation, the minimum lattice truss diameter is constrained to 1 mm to ensure the minimum feature dimensions remain printable with high geometric accuracy.

Shape design variables are also assigned to some of the lattice nodes. These nodes are shown in yellow on the left hand side of Fig. 6. The shape variables allow the nodes to move in two orthogonal directions for the case of plane stress problems and in all three directions for 3D stress problems. A node-to-surface freeze contact is implemented to connect the lattice nodes that are in close proximity (< 1 mm) to the solid and shell regions. Nodes that lie on a solid-lattice interface are not assigned shape variables. This ensures that lattice and solid/shell regions remain connected during size and shape optimisation.

During this parametric optimisation the same problem formulation, as described in Eq. (1) is implemented, except now it is the related truss diameters and positions that determine the relative density distribution. The final size and shape optimised designs are shown on the right hand size of Fig. 6. As a typical example, for the beam with p = 1 in Fig. 6 we generate 1644 size variables and 2198 shape variables.

We use one of the in-built optimisation algorithms in OptiStruct, known as the Method of Feasible Directions (MFD), for this parametric optimisation step. An advantage of this algorithm is that it always produces a feasible design if the initial design is itself feasible [52]. The algorithm improves the design at each iteration, so that an

Fig. 7. Optimisation workflow to generate topology optimised structures with lattice regions aligned with strain trajectories.

acceptable solution can be quickly achieved. It can be seen in Fig. 6 that the shape optimisation procedure has the effect of smoothing the final design. Discontinuous lines at lattice-void interfaces tend to become vanishingly small.

The final step in the process is to export the optimised structure for printing. The structure is exported from HyperMesh as a *.fem file with its updated size and shape parameters. The file is converted into *.stl format for printing via Materialise 3-matic software. This software enables the conversion of generation of tapered trusses, with smooth connections at lattice nodes. Lattice and solid regions are finally merged in 3-matic to create a single manifold surface. A summary of the optimisation workflow is shown in Fig. 7.

## 2.5. Experimental details

The lattice structures were printed by SLS on an EOS P395 machine. This is a powder based additive manufacturing process that uses a laser as a localised heat source to sinter polymeric powder particles into solid components. The raw powder material used in this work is polyamide PA 2200 (Nylon 12). Standard EOS parameters were used for the printing process with a layer thickness of 120 µ m [53]. The density and tensile modulus of the laser sintered parts are 0.93 g/cm 3 and 1700 MPa, respectively. The stress-strain relationship of the material is mostly linear up to yield strength. For all experimental validation tests, an Instron 5982 universal test machine was used with a 10kN load cell and at a test speed of 1 mm/min.

Printing accuracy was assessed for the relevant diameter ranges of circular structures via optical microscopy and found to be within +5%-10% of the design diameter. For the standard processing parameters used, lattice trusses tend to print larger than the design diameter, and this increase can be attributed to the partial attachment of sintered powder particles to the surface. For specimens printed for design validation, removal of powder after printing required several iterations of sand blasting the specimens, mechanical scraping and brushing. However, additional loose powder remained present within smaller cavities. These particles add to the surface roughness and some additional mass of the structure, but do not contribute structurally to an increase in stiffness and strength.

## 3. Results and discussion

In this section two case studies are presented that compare state-of-the art topology and lattice optimisation with our in-house developed method for generating functionally graded lattices.

Fig. 8. Benchmark studies. (a) SIMP with p = 4 and (b) the corresponding final design after 50% relative density threshold is applied. (c) Uniform truss diameter lattice aligned with strain trajectories, plotted for a uniform solid volume. (For interpretation of the references to colour in this figure legend, the reader is referred to the web version of this article.)

## 3.1. Floor beam

Referring back to the case study shown in Fig. 1, the lattice designs shown in Fig. 6 are benchmarked against a classical solid-void SIMP optimisation by setting the penalty factor high to p = 4 and re-running this optimisation problem. The minimum feature size is also reduced to 1 mm in this solid-void case, so it is comparable to the minimum allowable lattice diameters in the other designs. The optimal relative density distribution from this SIMP optimisation (p = 4) can be seen in Fig. 8a, showing that there are almost no intermediate relative density regions due to the high penalisation applied. The final design is shown in Fig. 8b after applying a 50% relative density threshold.

For an additional comparison, the result of applying the strain trajectory approach to a volume with a uniform density is shown in Fig. 8c. For this example, the truss diameters are generated of uniform size and scaled so that the resultant beam has the same mass as the other benchmark studies presented. The load path alignments between Figs. 8b and 8c show some similarities; however, the compliance of the latter design is 1.61 times higher. One of the causes for such a high compliance in Fig. 8c is the poor load transfer between the thin lattice members and the simply supported boundary conditions at the lower corners. The design in Fig. 8c will also be difficult to print with a powder bed fusion processes due to the close proximity of neighbouring lattice trusses near the mid-span on the beam, leading to particle fusion and entrapped powder.

Reverting back the beam designs presented in Fig. 6, the compliance results of these three optimised lattice designs with strain trajectories and different penalty factors are compared with the solid design, see Fig. 9, to determine the influence of the penalty factor on the structural response. The results have been normalised with respect to the solid design (p = 4). Of the three penalty factors available for lattice optimisation in OptiStruct, the p = 1. 8 penalty results in the final design with the lowest compliance when compared to the other available penalty factors. It can also be seen in Fig. 9 that the compliance of the final solid topology design using p = 4 results in the highest of all the four penalty factors considered. The cause of this is attributed to the high penalty factor approach leaving larger areas of the upper side of the beam, where the pressure load is applied, poorly supported by long slender members when applying the lower threshold density cut-off of 0.5.

Two designs containing uniform lattice unit cell topologies are also compared with the best performing optimised lattice design (Fig. 10) for further benchmarking purposes. The first of these benchmark lattice structures consists

Fig. 9. Relative compliance from FEA of the final optimised lattice (p = 1, p = 1. 25 and p = 1. 8) designs. Results are normalised to the solid design with p = 4.

Fig. 10. Printed models used for experimental validation; (a) cubic unit cells, (b) tetrahedron unit cells and (c) lattice design aligned with principal strain trajectories (p = 1. 8).

of cubic unit cells with 5 mm side lengths. Lattice trusses are placed along the edges of these cubic unit cells. The second of these two benchmark structures contains tetrahedron unit cells with an average edge length of 9 mm. Edge lengths vary slightly between cells in order for the cells to fit within the rectangular beam profile. Similarly, for the cubic lattice unit cells, trusses are placed along the edges of the tetrahedron cells. A 9 mm edge length is used for the tetrahedron unit cell to achieve a similar lower relative density bound of 10% when a minimum truss diameter constraint of 1 mm is imposed. For both of these uniform unit cell arrangements, each lattice node (truss diameter) is then size-optimised to produce a functionally graded lattice.

For design validation, 4-point bend test fixtures were used as a uniformly distributed load is difficult to achieve experimentally. For testing purposes, the uniform load is approximated using a 4-point bend test with loads applied at 1/3 and 2/3 span. This loading arrangement results in a similar bending moment distribution, maximum deflection and identical maximum bending moments at 1/2 span when the 4-point bending loads are scaled by a factor of 0.75. The load-displacement curves of the three lattice structures tested are shown in Fig. 11. There is close agreement between the initial stiffness estimated by FEA for the four structures. All lattice structures displayed an initial linear-elastic response followed by buckling in the lattice regions. The optimised design results in the highest failure load, in addition to the highest stiffness.

Fig. 11. Four-point bend load-displacement curves (FEA predictions shown with dashed lines).

Table 1 Beam compliance results, relative to the solid design, c / c$_{solid}$.

Table 2 Characteristic truss diameter ranges (FEA).

|                       |   Cubic |   Tetrahedron |   Strain trajectories ( p = 1 . 8) |
|-----------------------|---------|---------------|------------------------------------|
| Distributed load, FEA |    2.2  |          1.56 |                               0.63 |
| 4-point bend, FEA     |    2.65 |          2.11 |                               0.69 |
| 4-point bend, Expt.   |    2.52 |          2.08 |                               0.76 |

Table 1 shows the comparison of compliance for distributed and 4-point bend loading, based on FEA. The 4-point bend test is displacement controlled, which can lead to a non-uniform distribution of load applied along the loading rollers. This non-uniformity is negligible for the optimised design, where load is placed over solid regions in the structure. However, for the two benchmark designs, this discrepancy is more significant, which is considered a contributing factor to their higher compliance values. The compliance of all three beam designs also differs slightly in this 4-point loading configuration, compared to uniformly distributed loading, as a result of localised deformations around the loading points. However, in terms of relative improvements against the respective benchmark designs, the results show consistent trends.

|           |   Cubic |   Tetrahedron |   Strain trajectories ( p = 1 . 8) |
|-----------|---------|---------------|------------------------------------|
| Mean [mm] |    1.48 |          1.61 |                               1.18 |
| S.D [mm]  |    0.61 |          0.66 |                               0.23 |

The experimental compliance results for these three designs are also provided in Table 1. Overall the agreement between experimental testing and numerical analysis is considered good in absolute values and very similar in terms of relative stiffness improvement. It can be seen that the optimised design has the lowest compliance. This highlights that implementing the relative density bounds and arranging lattice cells are both beneficial for reducing compliance. It is also interesting to note that of the two uniform designs, it is the tetrahedron that performs better under distributed loading. The uniform square design has a bending-dominated behaviour whereas the tetrahedron unit cell is stretch-dominated [41,42]. The optimised design does not show bending-dominated behaviour since its arrangement of cells of varying shape, size and orientation prevents it from folding in a mechanism-like manner at its node regions.

Characteristic truss diameter ranges between these three beams are provided in Table 2. These values were extracted from the FEA models after size optimisation. Of the three designs, the trusses aligned with strain trajectories have the smallest average size. This is a consequence of the upper relative density threshold being

Fig. 12. Spider bracket case study: (a) design space, (b) benchmark solid optimisation with p = 4, (c) compression testing set-up. Dimensions in mm. (For interpretation of the references to colour in this figure legend, the reader is referred to the web version of this article.)

Table 3 Comparison of compliance for simulated and tested structures. Results are given relative solid design, c / c$_{solid}$.

|            |   p = 4 (solid) |   p = 1 . 8 (solid-lattice) |
|------------|-----------------|-----------------------------|
| FEA        |            1    |                        1.05 |
| Experiment |            0.96 |                        1.19 |

applied, which limits lattice generation to intermediate relative density regions. It is also interesting to note that the strain trajectory approach results in the lowest standard deviation in truss diameters. This result is also evident by observing Fig. 10, where the diameters of the aligned trusses in Fig. 10c appear more self-similar. This design feature is favourable for reduced stress concentrations at connecting nodes.

## 3.2. Spider bracket

The second case study is for a three-dimensional spider bracket, as shown in Fig. 12. A vertical load is applied to the central pad and the four lower pads are given simply supported boundary conditions. These 10 mm thick pads form the non-design space, shown in red in Fig. 12a. The same meshing, optimisation and lattice generation parameters are used in this second case study. The only difference in the setup of the optimisation problem description is that the design space volume fraction in Eq. (1c) is reduced to 10%.

The SIMP method with a high penalty factor of p = 4 tends to work well for such applications where there are discrete loading points to connect, see Fig. 12b. The SIMP result with a penalty of p = 1. 8 and a 10% lower relative density bound and a 60% upper relative density bound applied is shown in Fig. 13a. Only the penalty factor of p = 1. 8 is presented for this case study since it has already been shown in the previous case study to be the most effective penalty available for modelling lattices aligned to principal strain trajectories (Fig. 9). It can be seen in this p = 1. 8 design that the web regions primarily consist of intermediate densities, see Fig. 13b. The final printed lattice design after parametric size and shape optimisation is presented in Fig. 13c, showing that our in-house developed method is capable of generating a conformal, interconnected, lattice network in complex three dimensional spaces.

The boundary conditions of the experimental compression testing set-up are shown in Fig. 12c. A 6 mm diameter screw is attached to the central pad, which is then clamped to the moving cross-head of the Instron machine. The four arms are free to slide on a lower steel base plate. The relative compliance results are presented in Table 3, normalised with respect to the solid model with p = 4. For this case study the compliance of the solid design is the lowest, but it is noted that the compliance results of the two final FEA models are within 5%.

The solid design in Fig. 12b has large overhang regions that would require support material with metallic powderbed fusion processes [54], such as selective layer melting or electron beam melting. Support structure is removed after printing, resulting in additional waste material and labour-intensive post processing efforts. In contrast, the unsupported length in the lattice design is much smaller, effectively eliminating the need for additional support at

Fig. 13. Optimisation workflow for generating a bracket with lattice infill; (a) remaining intermediate relative density region after application of lower and upper relative density thresholds, (b) principal strain trajectories in remaining volume, (c) final printed structure comprising both solid and aligned lattice regions. (For interpretation of the references to colour in this figure legend, the reader is referred to the web version of this article.)

minimal increase of compliance. Hence, the potential for these optimised lattice designs to have a multi-functional application as a support structure during printing is investigated in the following.

A numerical investigation is carried out using the overhang constraint available in OptiStruct, see Fig. 14. The default overhang constraint 'CONSTR' option is strictly applied with the default larger overhang 'STEP' length activated to promote faster solution convergence. A 45 degree overhang constraint is implemented in both the solid and lattice optimisation workflows. This overhang angle is measured from the build direction, and a larger angle consequently implies more design freedom. As expected, the addition of this overhang constraint increases the compliance of the solid topology solution by 12%, see Fig. 14c. However, the inclusion of the overhang constraint with lattice optimisation reduces the compliance by 13%, hence the lattice design results in a better solution when realistic manufacturing constraints are considered. In Fig. 14d it can be seen that a solid region forms on the overhanging face.

## 4. Conclusions

Topology and lattice optimisation are increasingly being used as part of the design process for generating high performance additively manufactured parts. Methods for the design and optimisation of complex lattice structures have significant potential to realise enhanced or new functionality in a range of applications. We have presented a new lattice generation method based on lattices aligned with principal strain trajectories capable of full 3D design solutions. This method is shown to be an effective means of generating lightweight, high stiffness, lattice designs. The inclusion of upper and lower relative density bounds is an important consideration to ensure that the final lattice design is printable, within the constraints of a given additive manufacturing process. These relative density bounds are included in the optimisation workflow to ensure a printable lattice is generated.

The results presented are for two specific case studies, though the methodology presented is applicable more generally. For the floor beam example, where there is a requirement to support a uniformly distributed pressure load, our lattice optimisation method exceeds the state-of-the-art topology and lattice optimisation results in terms of minimising compliance. For the spider bracket example, where an additional overhang constraint is included, the compliance of our lattice optimisation solution is shown to be lower than the conventional solid TO result.

In addition to potentially further minimising compliance, the lattice generation procedure presented in this paper has other advantages over the current state-of-the art:

(1) Ability to closely achieve the original topology optimised compliance after lattice generation, using a stiffness penalty factor that accurately describes the lattice infill.

(2) Creation of conformal infill lattices for complex 3D volumes. Smooth outer surfaces of lattice regions that promote better strength and aesthetic appearance are created.

(3) Ability to create built-in functional support structure, hence reducing material waste and post-processing time.

Areas of future work include the study of objective functions other than minimising compliance and other constraints, such as buckling or yielding. The lattice structures in the presented case studies create a significantly larger surface area compared to equivalent solid topology designs. This benefits multi-functional applications where high mechanical performance is required in addition to heat transfer considerations.

Fig. 14. Influence of the 45 degree overhang constraint on compliance showing; (a) p = 4, (b) p = 1. 8, (c) p = 4 with overhang constraint and (d) p = 1. 8 with overhang constraint.

## Acknowledgements

The authors would like to acknowledge the support from the Agency for Science, Technology and Research (A*STAR) and Science and Engineering Research Council (SERC) of Singapore through the Additive Manufacturing Centre (AMC) Initiative-SIMTech-led R&D projects (SERC Grant No 142 6800088).

## References

[1] A. Evans, J. Hutchinson, N. Fleck, M. Ashby, H. Wadley, The topological design of multifunctional cellular metals, Prog. Mater. Sci. 46 (3) (2001) 309-327, http://dx.doi.org/10.1016/S0079-6425(00)00016-5.

[2] H. Wadley, Multifunctional periodic cellular metals, Phil. Trans. R. Soc. A 364 (1838) (2006) 31-68, http://dx.doi.org/10.1098/rsta. 2005.1697.

[3] G. Zu, J. Zhai, T. Zeng, Z. Wang, S. Cheng, D. Fang, Response of composite sandwich beams with graded lattice core, Compos. Struct. 119 (2015) 666-676, http://dx.doi.org/10.1016/j.compstruct.2014.09.042.

[4] C. Thomsen, F. Wang, O. Sigmund, Buckling strength topology optimization of 2d periodic materials based on linearized bifurcation analysis, Comput. Methods Appl. Mech. Engrg. 339 (2018) 115-136, http://dx.doi.org/10.1016/j.cma.2018.04.031.

[5] T. Lu, Heat transfer efficiency of metal honeycombs, Int. J. Heat Mass Transfer 42 (11) (1999) 2031-2040, http://dx.doi.org/10.1016/ S0017-9310(98)00306-8.

[6] T. Schaedler, A. Jacobsen, A. Torrents, A. Sorensen, J. Lian, J. Greer, L. Valdevit, W. Carter, Ultralight metallic microlattices, Science 334 (6058) (2011) 962-965, http://dx.doi.org/10.1126/science.1211649.

[7] I. Maskery, N.T. Aboulkhair, A.O. Aremu, C.J. Tuck, I.A. Ashcroft, Compressive failure modes and energy absorption in additively manufactured double gyroid lattices, Addit. Manuf. 16 (2017) 24-29, http://dx.doi.org/10.1016/j.addma.2017.04.003.

[8] J. Tian, T. Kim, T. Lu, H. Hodson, D. Queheillalt, D. Sypeck, H. Wadley, The effects of topology upon fluid-flow and heat-transfer within cellular copper structures, Int. J. Heat Mass Transfer 47 (14) (2004) 3171-3186, http://dx.doi.org/10.1016/j.ijheatmasstransfer. 2004.02.010.

[9] H. Brooks, K. Brigden, Design of conformal cooling layers with self-supporting lattices for additively manufactured tooling, Addit. Manuf. 11 (2016) 16-22, http://dx.doi.org/10.1016/j.addma.2016.03.004.

[10] F. Scarpa, M. Ouisse, M. Collet, K. Saito, Kirigami auxetic pyramidal core: mechanical properties and wave propagation analysis in damped lattice, J. Vib. Acoust. 135 (4) (2013) 041001, http://dx.doi.org/10.1115/1.4024433.

[11] Y. Koizumi, A. Okazaki, A. Chiba, T. Kato, A. Takezawa, Cellular lattices of biomedical co-cr-mo-alloy fabricated by electron beam melting with the aid of shape optimization, Addit. Manuf. 12 (2016) 305-313, http://dx.doi.org/10.1016/j.addma.2016.06.001.

[12] L. Wang, J. Kang, C. Sun, D. Li, Y. Cao, Z. Jin, Mapping porous microstructures to yield desired mechanical properties for application in 3d printed bone scaffolds and orthopaedic implants, Mater. Des. 133 (2017) 62-68, http://dx.doi.org/10.1016/j.matdes.2017.07.021.

[13] U. Wegst, H. Bai, E. Saiz, A. Tomsia, R. Ritchie, Bioinspired structural materials, Nature Mater. 14 (2014) 23-36, http://dx.doi.org/ 10.1038/nmat4089.

[14] D. Jang, L. Meza, F. Greer, J. Greer, Fabrication and deformation of three-dimensional hollow ceramic nanostructures, Nature Mater. 12 (10) (2013) 893-898, http://dx.doi.org/10.1038/nmat3738.

[15] J. Weaver, J. Aizenberg, G. Fantner, D. Kisailus, A. Woesz, P. Allen, K. Fields, M. Porter, F. Zok, P. Hansma, P. Fratzl, D. Morse, Hierarchical assembly of the siliceous skeletal lattice of the hexactinellid sponge euplectella aspergillum, J. Struct. Biol. 158 (1) (2007) 93-106, http://dx.doi.org/10.1016/j.jsb.2006.10.027.

[16] J. Wolff, The classic: on the inner architecture of bones and its importance for bone growth. 1870, Clin. Orthop. Relat. Res. 468 (4) (2010) 1056-1065, http://dx.doi.org/10.1007/s11999-010-1239-2.

[17] I. Ullaha, J. Elambasseril, M. Brandt, S. Feih, Performance of bio-inspired kagome truss core structures under compression and shear loading, Compos. Struct. 118 (2014) 294-302, http://dx.doi.org/10.1016/j.compstruct.2014.07.036.

[18] L. Zhang, S. Feih, S. Daynes, S. Chang, M. Wang, J. Wei, W. Lu, Energy absorption characteristics of metallic triply periodic minimal surface sheet structures under compressive loading, Addit. Manuf. (2018) http://dx.doi.org/10.1016/j.addma.2018.08.007.

[19] L. Gibson, M. Ashby, Cellular Solids: Structure and Properties, Cambridge University Press, 1999.

[20] I. Maskery, A. Hussey, A. Panesar, A. Aremu, C. Tuck, I. Ashcroft, R. Hague, An investigation into reinforced and functionally graded lattice structures, J. Cell. Plast. 53 (2) (2017) 151-165, http://dx.doi.org/10.1177/0021955X16639035.

[21] S. Khanoki, D. Pasini, Multiscale design and multiobjective optimization of orthopedic hip implants with functionally graded cellular material, J. Biomech. Eng. 134 (3) (2012) 031004, http://dx.doi.org/10.1115/1.4006115.

[22] L. Murr, S. Gaytan, F. Medina, H. Lopez, E. Martinez, B. Machado, D. Hernandez, L. Martinez, M. Lopez, R. Wicker, J. Bracke, Next-generation biomedical implants using additive manufacturing of complex, cellular and functional mesh arrays, Phil. Trans. R Soc. A 368 (1917) (2010) 1999-2032, http://dx.doi.org/10.1098/rsta.2010.0010.

[23] I. Maskery, L. Sturm, A. Aremua, A. Panesar, C. Williams, C. Tuck, R. Wildman, I. Ashcroft, R. Hague, Insights into the mechanical properties of several triply periodic minimal surface lattice structures made by polymer additive manufacturing, Polymer (2017) (in press).

[24] L. Zhang, S. Feih, S. Daynes, Y. Wang, M.Y. Wang, J. Wei, W.F. Lu, Buckling optimization of kagome lattice cores with free-form trusses, Mater. Des. (2018) http://dx.doi.org/10.1016/j.matdes.2018.02.026, (in press).

[25] A. Aremu, J. Brennan-Craddock, A. Panesar, I. Ashcroft, R. Hague, R. Wildman, C. Tuck, A voxel-based method of constructing and skinning conformal and functionally graded lattice structures suitable for additive manufacturing, Addit. Manuf. 13 (2017) 1-13, http://dx.doi.org/10.1016/j.addma.2016.10.006.

[26] Y. Chen, 3d texture mapping for rapid manufacturing, Computer-Aided Des. Appl. 4 (6) (2013) 761-771, http://dx.doi.org/10.1080/ 16864360.2007.10738509.

[27] D. Brackett, I. Ashcroft, R. Hague, A dithering based method to generate variable volume lattice cells for additive manufacturing, in: 22nd Annual International Solid Freeform Fabrication Symposium, Austin, TX, 2011.

[28] D. Brackett, I. Ashcroft, R. Wildman, R. Hague, An error diffusion based method to generate functionally graded cellular structures, Comput. Struct. 138 (2014) 102-111, http://dx.doi.org/10.1016/j.compstruc.2014.03.004.

[29] A. Panesar, I. Ashcroft, D. Brackett, R. Wildman, R. Hague, Design framework for multifunctional additive manufacturing: coupled optimization strategy for structures with embedded functional systems, Addit. Manuf. 16 (2017) 98-106, http://dx.doi.org/10.1016/j. addma.2017.05.009.

[30] A. Panesar, M. Abdi, D. Hickman, I. Ashcroft, Strategies for functionally graded lattice structures derived using topology optimisation for additive manufacturing, Addit. Manuf. 19 (2018) 81-94, http://dx.doi.org/10.1016/j.addma.2017.11.008.

[31] Y. Wang, H. Xu, D. Pasini, Multiscale isogeometric topology optimization for lattice materials, Comput. Methods Appl. Mech. Engrg. 316 (2017) 568-585, http://dx.doi.org/10.1016/j.cma.2016.08.015.

[32] J. Wu, A. Clausen, O. Sigmund, Minimum compliance topology optimization of shell-infill composites for additive manufacturing, Comput. Methods Appl. Mech. Engrg. 326 (2017) 358-375, http://dx.doi.org/10.1016/j.cma.2017.08.018.

[33] J. Robbins, S. Owen, B. Clark, T. Voth, An efficient and scalable approach for generating topologically optimized cellular structures for additive manufacturing, Addit. Manuf. 12 (2016) 296-304, http://dx.doi.org/10.1016/j.addma.2016.06.013.

[34] Y. Wang, F. Chen, M. Wang, Concurrent design with connectable graded microstructures, Comput. Methods Appl. Mech. Engrg. 317 (2017) 84-101, http://dx.doi.org/10.1016/j.cma.2016.12.007.

[35] Y. Wang, L. Zhang, S. Daynes, H. Zhang, S. Feih, M.Y. Wang, Design of graded lattice structure with optimized mesostructures for additive manufacturing, Mater. Des. 142 (2018) 114-123, http://dx.doi.org/10.1016/j.matdes.2018.01.011.

[36] Y. Tang, Y. Zhao, Lattice-skin structures design with orientation optimization, in: Proceedings of the Solid Freeform Fabrication Symposium, Austin, TX, 2015.

[37] G. Reinhart, S. Teufelhart, F. Riss, Investigation of the geometry-dependent anisotropic material behavior of filigree struts in alm-produced lattice structures, Physics Procedia 39 (2012) 471-479, http://dx.doi.org/10.1016/j.phpro.2012.10.063.

[38] A. Michell, The limits of economy of material in frame-structures, Lond. Edinb. Dubl. Phil. Mag. J. Sci. 8 (47) (1904) 589-597, http://dx.doi.org/10.1080/14786440409463229.

[39] O. Sigmund, A 99 line topology optimization code written in matlab, Struct. Multidiscip. Optim. 21 (2) (2001) 120-127.

[40] J. Berger, H. Wadley, R. McMeeking, Mechanical metamaterials at the theoretical limit of isotropic elastic stiffness, Nature (2017) http://dx.doi.org/10.1038/nature21075.

[41] V. Deshpande, M. Ashby, N. Fleck, Foam topology: bending versus stretching dominated architectures, Acta Mater. 49 (6) (2001) 1035-1040, http://dx.doi.org/10.1016/S1359-6454(00)00379-7.

[42] M. Ashby, The properties of foams and lattices, Phil. Trans. R. Soc. A 364 (1838) (2006) 15-30, http://dx.doi.org/10.1098/rsta.2005. 1678.

[43] I. Maskery, A. Aremu, M. Simonelli, C. Tuck, R. Wildman, I. Ashcroft, R. Hague, Mechanical properties of ti-6al-4v selectively laser melted parts with body-centred-cubic lattices of varying cell size, Exp. Mech. 55 (7) (2015) 1261-1272, http://dx.doi.org/10.1007/ s11340-015-0021-5.

[44] J. Wu, N. Aage, R. Westermann, O. Sigmund, Infill optimization for additive manufacturing-approaching bone-like porous structures, IEEE Trans. Vis. Comput. Graph. 24 (2) (2018) 1127-1140.

[45] C. Portela, J. Greer, D. Kochmann, Impact of node geometry on the effective stiffness of non-slender three-dimensional truss lattice architectures, Extreme Mech. Lett. 22 (2018) 138-148, http://dx.doi.org/10.1016/j.eml.2018.06.004.

[46] L. Meza, G. Phlipot, C. Portela, A. Maggi, L. Montemayor, A. Comella, D. Kochmann, J. Greer, Reexamining the mechanical property space of three-dimensional lattice architectures, Acta Mater. 140 (2017) 424-432, http://dx.doi.org/10.1016/j.actamat.2017.08.052.

[47] S. Daynes, S. Feih, W. Jun, Method and system of manufacturing a load-bearing structure and a load-bearing structure manufactured thereof. Singapore Patent Application No. PCT/SG2017/050637, 21 2017.

[48] S. Daynes, W. Jun, L. Wen Feng, S. Feih, Optimisation of functionally graded lattice structures using isostatic lines, Mater. Des. 127 (2017) 215-223, http://dx.doi.org/10.1016/j.matdes.2017.04.082.

[49] S.G.J. Timoshenko, Theory of Elasticity, McGraw-Hill, New York, 1951.

[50] P. Frankovský, F. Trebuˇna, J. Kostka, F. Šimˇcák, O. Ostertag, W. Papacz, Characteristic entities in photostress method, Am. J. Mech. Eng. 2 (7) (2014) 239-243, http://dx.doi.org/10.12691/ajme-2-7-13.

[51] D. Kelly, C. Reidsema, A. Bassandeh, G. Pearce, M. Lee, On interpreting load paths and identifying a load bearing topology from finite element analysis, Finite Elem. Anal. Des. 47 (8) (2011) 867-876, http://dx.doi.org/10.1016/j.finel.2011.03.007.

[52] X. Chen, M. Kostreva, Methods of feasible directions: a review, in: Progress in Optimization, Springer, Boston, MA, 2000, pp. 205-219, http://dx.doi.org/10.1007/978-1-4613-0301-5_14.

[53] PA 2200 Material data sheet, EOS GmbH-Electro Optical Systems, Munich, Germany, 2008.

[54] T. Wang, J. Gao, Z. Kang, Level set-based topology optimization with overhang constraint: Towards support-free additive manufacturing, Comput. Methods Appl. Mech. Engrg. 339 (2018) 591-614, http://dx.doi.org/10.1016/j.cma.2018.04.040.

