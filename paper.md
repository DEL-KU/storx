---
title: 'STORX: An Open-Source Object-Oriented Framework for Shape and Topology Optimization in MATLAB'
tags:
  - MATLAB
  - topology optimization
  - shape optimization
  - finite element analysis
  - structural optimization
  - object-oriented programming
  - engineering education
authors:
  - name: Amir M. Mirzendehdel
    orcid: 0000-0002-4407-1877
    affiliation: 1
  - name: Krishnan Suresh
    orcid: 0000-0002-9688-9697
    affiliation: 2
affiliations:
 - name: University of Kansas, USA
   index: 1
 - name: University of Wisconsin-Madison, USA
   index: 2
date: 07 July 2026
bibliography: paper.bib
---

# Summary

STORX (Shape and Topology Optimization for Research and Experimentation) is an open-source, object-oriented MATLAB framework and, to the best of our knowledge, the first to unify parametric shape, level-set shape, and multiple families of topology optimization within a single consistent architecture. It provides the full design optimization pipeline in a single codebase: boundary representation of general 2D geometry, structured mesh generation, finite element analysis (elasticity, thermal, thermoelasticity, and fluid physics), and multiple optimization paradigms, including parametric shape optimization [@suresh2021design], level-set shape optimization[@van2013level,@allaire2005structural], density-based topology optimization (SIMP, RAMP), level-set topology optimization (standard and modified Hamilton-Jacobi), and topological-sensitivity-driven methods (ESO, BESO, and Pareto-tracing) [@sigmund2013topology,@feijoo2005topological,@suresh2013efficient]. 

Rather than implementing each method as an isolated script, STORX organizes them under shared abstract class hierarchies: a common `simulation2d` interface for all physics solvers, a `functional` interface for objectives and constraints, and a `mfgConstraints` interface for manufacturing constraints such as minimum feature size, retained regions, and symmetry. This separation makes explicit which parts of the pipeline (state solve, functional evaluation, sensitivity analysis, regularization, constrained update) are shared across shape and topology optimization (SO/TO) methods and which are method-specific, and it allows new physics solvers, objectives, and constraints to be added as derived classes without modifying the core code. Optimized designs can be exported to DXF and watertight STL for downstream manufacturing (e.g., 3D printing). The repository is organized as a sequence of numbered chapters with matching worked examples that build up the framework progressively, from geometry and analysis to each family of optimization methods.

# Statement of need

Several well-known MATLAB codes serve as standard entry points to topology optimization, including the 99-line SIMP code [@sigmund2001], its 88-line successor [@andreassen2011], a discrete level-set implementation [@challis2010], and a Pareto-tracing code based on topological sensitivity [@suresh2010]. For a comprehensive review of such educational codes published through 2021, see [@wang2021comprehensive].

These codes are valuable precisely because of their brevity, but the same "script-per-method" style limits their use beyond a single method. Each script hardcodes one physics model, one mesh (typically a rectangular domain), and one optimization algorithm together, so components cannot easily be swapped or methods compared within a common structure. For researchers, this means that prototyping a new algorithm, objective, or constraint typically requires rebuilding the geometry, meshing, and analysis layers from scratch, or invasively modifying a monolithic script. It also makes controlled comparisons between SO/TO methods difficult, since each implementation carries its own physics, discretization, and numerical settings. For students and instructors, the same coupling obscures the design-optimization pipeline that all SO/TO methods share.

STORX addresses this gap by separating representation, meshing, analysis, and optimization into distinct, well-defined classes, and by defining the core software interfaces via abstract base classes. This architecture enables controlled comparisons between SO/TO methods under consistent physical and numerical conditions on general 2D geometries, and provides a readable, extensible foundation for prototyping new SO/TO algorithms: new objective functionals, design/manufacturing constraints, and physics solvers are implemented as derived classes that conform to established interfaces, without rewriting the rest of the pipeline. Key routines are vectorized for efficiency, with optional explicit nested-loop implementations retained for transparency and ease of learning. A companion article provides a longer technical description of the underlying algorithms and software design [@mirzendehdel2026storx].

# Target audience

STORX is intended for researchers who want an extensible codebase for developing and benchmarking new SO/TO algorithms without building the representation, meshing, and FEA layers from scratch, and for senior undergraduate and graduate students in mechanical, aerospace, or civil engineering studying structural, shape, or topology optimization. The chapter-organized examples also give instructors runnable material to accompany lectures.

# Evidence of use

STORX and its chapter-organized examples were used as the primary computational framework in AE700: Shape and Topology Optimization, a graduate-level course at the University of Kansas, in Fall 2025. Students used the framework to work through boundary representation, meshing, and finite element analysis before implementing and comparing parametric shape optimization, density-based topology optimization, level-set topology optimization, and topological-sensitivity-driven methods, following the chapter sequence in the `00-examples` directory.

The framework's modular architecture enabled rapid incorporation of new formulations: for example, following an invited guest lecture introducing local volume fraction constraints [@wu2017], students were able to formulate and deploy this constraint within STORX as part of a homework assignment. Students also used STORX as the computational basis for independent research-oriented final projects spanning antenna shape optimization, fluid-structure interaction for wing design, and topology optimization of auxetic metamaterials, extending the framework with their own physics models, optimization formulations, and example problems.

# Acknowledgements

A.M. gratefully acknowledges the Aerospace Department at the University of Kansas for supporting this work.
A.M. also thanks the students of the graduate-level course AE700, Shape and Topology Optimization (Fall 2025), for their feedback during the development of STORX, and Dr. Jun Wu, Dr. Julian Norato, Dr. Faez Ahmed, Dr. Aaditya Chandrasekhar, and Dr. Krishnan Suresh for their insightful guest lectures, which enriched the course and contributed to the framework's further development.

# References
