---
title: 'AGATE.jl: A probabilistic programming framework for tuneable aquatic ecosystems'
tags:
  - Julia
  - biogeochemistry
  - probabilistic programming
  - process-based ecological modelling
  - species distribution modelling
  - earth science
  - ecology
  - geography
authors:
  - name: Joost DeVries
    orcid: https://orcid.org/0000-0003-3427-6921
    equal-contrib: false
    affiliation: 1 
  - name: Radka Jersakova 
    orcid: https://orcid.org/0000-0001-6846-7158
    equal-contrib: false
    affiliation: 2 
  - name: Christopher Follet
    orcid: https://orcid.org/0000-0002-7164-1660
    equal-contrib: true
    affiliation: 3
  - name: Fanny Monteiro
    orcid: https://orcid.org/0000-0002-7164-1660
    equal-contrib: true
    affiliation: 1 
  - name: Levi John Wolf
    orcid: https://orcid.org/0000-0003-0274-599X
    equal-contrib: true
    affiliation: 1     
affiliations:
 - name: University of Bristol, United Kingdom
   index: 1
 - name: Alan Turing Institute, United Kingdom
   index: 2
 - name: University of Liverpool, United Kingdom
   index: 3
date: 30 March 2025
bibliography: paper.bib
---

# Summary

Aquatic ecosystems are critical in regulating global carbon sequestration, biodiversity and environmental resources. However, their inherent processes are highly complex and uncertain, limiting our understanding of ecosystem functions under threat in a changing climate. To tackle these issues, we have developed a next-generation aquatic ecosystem modelling framework, combining cutting-edge ecosystem theory and probabilistic programming. This new framework supports new thinking about how we represent aquatic ecosystems, making models more accessible, extendable and accurate in quantifying and predicting uncertainties. This package will support new scientific insights into aquatic ecosystems while building the social and technical infrastructure to support the aquatic modelling community.

# Statement of need

Since aquatic ecosystems are vital for many aspects of our lives, it is essential to develop digital twins of these systems that accurately represent their processes. Such twins can help inform climate policies and promote habitat restoration, thereby facilitating our society’s transition to net zero emissions. However, current aquatic ecosystem models in climate change intercomparison projects, such as IPCC, show large disagreements in response to climate change. There are even large inter-model disagreements in historical ocean primary productivity , which poses problems for model validation and constraining uncertainties in carbon dioxide removal policies . Despite large uncertainties in current ecosystem models, they are not captured in model outputs, which rely on large sets of deterministic partial differential equations describing the system state and its changes. The deterministic nature of these models also limits our ability to conduct inferences about the processes within the model itself, for example, comparing the relative effect of predation versus growth on biodiversity. Likewise, their nature prevents us from modelling the uncertainty intrinsic to measurements we make of our natural environment. Finally, current state-of-the-art aquatic ecosystem models such as the `MITgcm-DARWIN` are implemented in FORTRAN, making them hard to access, interpret, extend and couple with other models. As a result, inter-model ecosystem comparison is challenging, with ecosystem representation entangled with other processes inherent to each model, such as physics and biogeochemistry. 

Thus, there is an urgent need to develop aquatic ecosystem models that combine process knowledge based on centuries of ecological theory and advanced probabilistic programming methods to provide flexible, extensible, performant and composable models of the aquatic environment that can capture intrinsic uncertainties about ecosystem processes.

# Models

We developed `AGATE.jl` (Aquatic GCM-Agnostic Tunable Ecosystems) to as a framework of probabilistic programming primitives to support the development and extension of aquatic ecosystem models. In this framework, model parameters representing traits, such as size and photosynthesis, are mechanistically defined and co-vary based on trait trade-offs based on well-established ecological theories . Such trade-offs represent the costs and benefits of investment in different traits. For example, predator defence is a trait that comes at the cost of additional resource requirements. In this vein, `AGATE.jl` implements many of the basic processes required to model aquatic ecosystems as follows.

## Underlying ecosystem model

Basic idea of the BGC model as a biology model (like NPZD) integrated with a geochemical model from Oceananigans/OceanBioMe. 

## bio processes

### Photosynthesis

### predation

### calcification

## geochemical processes

short discussion of how the package integrates w/ other bits

# Example




<!-- # FOR REFERENCE
# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
-->