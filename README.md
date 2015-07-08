---
bibliography: 'tube-dynamics-partial-control.bib'
nocite: |
    @sabuco2010partial, @Sabuco2012c, @Sabuco2012, @Coccolo2013,
    @Zambrano2008
output:
  html_document:
    theme: united
    toc: True
  md_document:
    toc: True
    variant: 'markdown\_github'
  pdf_document:
    highlight: zenburn
    toc: True
title: README for partial control of ship dynamics
...

### Usage

MATLAB files with **run** as prefix are driver script that depend on
other functions to perform certain computation. MATLAB files with
**func** as prefix are functions and can be called from MATLAB script or
command line. Check if it has global variables (not the best practice,
so I use it sparingly) and define accordingly. I have refrained from
changing the nomenclature used by @Sabuco2012c as that is not the goal
of this research.

### References

Papers which developed and implemented the safe set sculpting
algorihtms:
