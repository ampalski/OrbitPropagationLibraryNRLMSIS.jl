# OrbitPropagationLibraryNRLMSIS

[![Build Status](https://github.com/ampalski/OrbitPropagationLibraryNRLMSIS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ampalski/OrbitPropagationLibraryNRLMSIS.jl/actions/workflows/CI.yml?query=branch%3Amain)

Implement's a Julia version of NRL's NRLMSISE-00 atmospheric model. 

# Example Use
```
    using OrbitPropagationLibraryNRLMSIS
    using OrbitPropagationLibrarySOFA
    jd0, _ = datevec2jdate([2020.0, 2, 23, 0, 0, 0])
    lat = 55.0 # deg
    lon = 45.0 # deg
    alt = 100.0 # km

    density = getDensity(jd0, lat, lon, alt) / 1000
```

# Planned Improvements
* Predictive space weather values
