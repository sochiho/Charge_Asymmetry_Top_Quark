# Charge Asymmetry Top Quark

Durham University Level 4 Theoretical Physics project 

There is a large discrepancy between the measured charge asymmetry at the
Tevatron and the theoretical predictions from the Standard Model for top quark pair
production. This asymmetry arises from a difference in the production rates of top
and antitop quarks in the forward hemisphere. In this work, a parton-level Monte
Carlo event generator is constructed to directly compute the charge asymmetry for
the Tevatron and at the Large Hadron Collider (LHC), which is calculated at first
non-vanishing order (&alpha;<sub>s</sub>). For the Tevatron, an asymmetry of 4-5% is predicted
in the lab frame and ∼ 7% in the tt̄ frame which differs from the experimental
result of 15% measured by CDF. At the LHC, cut-independent and cut-dependent
asymmetries are introduced to study the asymmetry in suitable kinematic regions.
Through these cuts the top quark is found to be predominantly produced at higher
rapidities while the antitop is more centrally produced.

The Monte Carlo (MC) integration method was used due its fast convergence for high dimensional integrals which is particularly relevant for particle physics. The other aspect which makes it appealing is how it can be utilised for constructing event generators due to the randomness of the numerical method being analogous to the random physical process. The benefit of generating events is that multiple observables can be binned simultaneously into histograms without analytically performing Jacobian transforms which can be mathematically arduous.

References to research papers and Lepage's VEGAS integrator, as well as the relevant parton distribution function are provided in the 'Top_Quark_Charge_Asymmetry.pdf'physics report.

In this project, a MC event generator was produced for the simpler leptonic process, e<sup>+</sup> e<sup>-</sup> &rarr; Z/&gamma; &rarr; &mu;<sup>+</sup> &mu;<sup>-</sup>. This process allowed us to build a good foundation to be expanded for top quark production. The files for producing (vegas_ee.py) and histogramming (ee_binning.py) theoretical data can be found in the ee_mumu.

The phase-space slicing method, which requires a soft cutoff energy (E<sub>cut</cut>) was employed to separate the phase-space being integrated into two regions labelled 'hard' and 'soft'. This was neccesary to solve the infared divergences that occur from integration momenta tending to zero in phase space integrals for one-loop diagrams and in real gluon emission. The files for showing that the soft and hard contributions are independent of E<sub>cut</cut> are given in the folder 2_3_body.

The folder DAT_Bin contains the relevant files for producing thoeretical data for calculating the charge asymmetry at the LHC and the Tevatron. It also contains files for binning the data into histograms. The program is capable of producing distributions for the rapidity, invariant mass, psuedo-rapidity in the lab and rest frame. Other observables can also be easily implemented. 

The graphs and tables produced by these programs are included in the folder Graphs_tables.
