# pySRIMP: Super-Resolution for Icy Moons and Planetary bodies

## Description
pySRIMP is a Python pipeline for building super-resolution (SR) maps of planetary
surfaces from overlapping hyperspectral observations. It takes a set of calibrated
image cubes acquired under different viewing and illumination conditions and
combines them into a single map at finer spatial sampling than the individual
observations.

The pipeline runs in three stages:

1. **Navigation** – computes accurate projections for every cube from the
   spacecraft geometry (SPICE kernels), so that observations can be co-registered
   into global mosaics.
2. **Photometric correction** – fits and applies a global photometric model that
   compensates for the varying incidence, emission and phase angles across
   acquisitions.
3. **Super-resolution** – merges all overlapping, corrected observations of a
   target into the final SR map.

pySRIMP was developed and validated on the Cassini VIMS observations of Enceladus,
where it yields the first global infrared map of the moon at enhanced resolution
(Ntinos et al., 2026). The methodology is not specific to that data set: with
appropriate kernels and calibrated cubes it can be applied to other icy moons and
instruments, including future hyperspectral data from missions such as JUICE.

Mapping products at this level of detail matter because they link surface
morphology and composition to internal activity, which is central to constraining
the structure and potential habitability of the outer solar system's icy moons.

## Pipeline
In the following schematic the user can find the SR pipeline workflow as implemented in the _main.py module. The pipeline consists of three blocks corresponding to the three main steps of the methodology: Navigation, Photometric Correction and Super-Resolution. An additional block to the right, lists auxiliary modules used for constants, filtering of the data and visualisation.

![Pipeline flowchart](sr_pipeline_flowchart.svg)

## Usage
The SR pipeline has been tested on the Enceladus VIMS data set (but we strongly encourage broad implementation on varying instruments and targets). The names of the calibrated cubes used can be found in this repository in the Enceladus_cubes_ds4.txt file (Enceladus_cubes_flybys_ds4.txt contains the corresponding flyby names) and can be downloaded from the [Cassini VIMS Data Portal by Observatoire des Sciences de l'Univers Nantes Atlantique (Osuna)](https://vims.univ-nantes.fr/). The list of the SPICE kernels used can be found in the following metakernels: cas_2004_v26.tm, cas_2005_v27.tm, cas_2006_v26.tm, cas_2007_v24.tm, cas_2008_v22.tm, cas_2009_v19.tm, cas_2011_v17.tm, cas_2012_v14.tm, cas_2013_v12.tm, cas_2014_v08.tm, cas_2015_v08.tm, cas_2016_v06.tm, cas_2017_v03.tm (most updated at the time of the development) and can be downloaded from the [Navigation and Ancillary Information Facility (NAIF)](https://naif.jpl.nasa.gov/naif/).

## Contributing
Development takes place on the IPGP GitLab at <https://pss-gitlab.math.univ-paris-diderot.fr/ntinos/super-resolution-for-icy-moons>.
The GitHub repository at <https://github.com/cntinos/Super-Resolution-for-Icy-Moons> is an automatically updated mirror.

In case the institutional GitLab requires an approved account, please report bugs, ask questions and request features by opening an issue on **GitHub**. Pull requests on GitHub are welcome; they are merged on GitLab and mirrored back.

For other enquiries contact C. Ntinos at ntinos@ipgp.fr or christos.ntinos@community.isunet.edu.

## Authors and acknowledgment
Accompanying paper:

[Christos Ntinos](https://orcid.org/0000-0003-4709-3665), [Sebastien Rodriguez](https://orcid.org/0000-0003-1219-0641), [Stephane Le Mouelic](https://orcid.org/0000-0001-5260-1367), [Rozenn Robidel](https://orcid.org/0000-0002-5544-8940), [Benoit Seignovert](https://orcid.org/0000-0001-6533-275X), [Thomas Cornet](https://orcid.org/0000-0001-5971-0056), [Nicolas Altobelli](https://orcid.org/0000-0003-4244-6302), [Gabriel Tobie](https://orcid.org/0000-0002-4494-0294), and [Claire Vallat](https://orcid.org/0000-0001-6921-868X) (2026). Super-resolution maps of icy moons: Application to the Cassini-VIMS observations of Enceladus with updated navigation and photometric correction. *Planetary Research*, **1**(2), 349.
[https://doi.org/10.53480/arz2-gn59](https://doi.org/10.53480/arz2-gn59)

## License
TBD

## Project status
Project is under development.
