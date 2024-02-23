# bling_baltic

This repository contains a version of the Biology Light Iron Nutrients and Gases (BLING) model (Galbraith et al. 2010, https://doi.org/10.5194/bg-7-1043-2010), originally developed as a component of the Modular Ocean Model (MOM, https://github.com/mom-ocean/MOM5). 

The BLING code in this repository is compatible with NEMO version 4.2.1 (https://forge.nemo-ocean.eu/nemo/nemo/-/releases/4.2.1).

This code was adapted from previous NEMO 3.6-compatible BLING code developed at the **University of Alberta**: this is stored at **https://borealisdata.ca/dataset.xhtml?persistentId=doi:10.5683/SP3/DMGYXI** and is described in detail by **Laura Castro de la Guardia in her PhD thesis (https://doi.org/10.7939/R31G0J98H)**. 

Between versions 3.6 and 4.2.1 of NEMO there were substantial changes to the TOP module of NEMO, which this new version of NEMO-BLING accounts for.  Other modifications include a basic representation of nitrogen limitation, inspired by the MITgcm version of BLING (https://github.com/MITgcm/MITgcm/tree/master/pkg/bling).

All code necessary to compile BLING with NEMO 4.2.1 is available in the MY_SRC directory.  The directory 'EXP_BLING' contains some example namelists and xml files for a low resolution test case covering the Baltic Sea.  Input and boundary files are not included.

-------------------------------------

To compile:

1) mkdir cfgs/BLING_TEST
2) mkdir cfgs/BLING_TEST/MY_SRC
3) Copy all MY_SRC files from this repository to cfgs/BLING_TEST/MY_SRC
4) Copy cpp_BLING_BALTIC.fcm from this repository to cfgs/BLING_TEST

*./makenemo -m [insert_your_architecture file] -r ORCA2_ICE_PISCES -n BLING_TEST*

-------------------------------------

References:

*BLING with MOM*:
  - Galbraith, Eric D., Anand Gnanadesikan, John P. Dunne, and Michael R. Hiscock. "Regional impacts of iron-light colimitation in a global biogeochemical model." Biogeosciences 7, no. 3 (2010): 1043-1064.
  - Galbraith, Eric D., John P. Dunne, Anand Gnanadesikan, Richard D. Slater, Jorge L. Sarmiento, Carolina O. Dufour, Gregory F. De Souza et al. "Complex functionality with minimal computation: Promise and pitfalls of reduced‐tracer ocean biogeochemistry models." Journal of Advances in Modeling Earth Systems 7, no. 4 (2015): 2012-2028.

***BLING with NEMO***:
  - **Castro de la Guardia, Laura. 2018. Modelling the response of Arctic and SubArctic marine systems to climate warming. Doctor of Philosophy, University of Alberta.**
  - **Castro de la Guardia, Laura, Yarisbel Garcia‐Quintana, M. Claret, Xianmin Hu, E. D. Galbraith, and Paul G. Myers. "Assessing the role of high‐frequency winds and sea ice loss on Arctic phytoplankton blooms in an ice‐ocean‐biogeochemical model." Journal of Geophysical Research: Biogeosciences 124, no. 9 (2019): 2728-2750.**
  - **Castro de la Guardia, Laura, Tania Hernández Fariñas, Christian Marchese, Martí Amargant-Arumí, Paul G. Myers, Simon Bélanger, Philipp Assmy, Rolf Gradinger, and Pedro Duarte. "Assessing net primary production in the northwestern Barents Sea using in situ, remote sensing and modelling approaches." Progress in Oceanography 219 (2023): 103160.**

*BLING with MITgcm*:
  - Verdy, Ariane, and Matthew R. Mazloff. "A data assimilating model for estimating Southern Ocean biogeochemistry." Journal of Geophysical Research: Oceans 122, no. 9 (2017): 6968-6988.
  - Twelves, Andrew G., Daniel N. Goldberg, Sian Frances Henley, Matthew R. Mazloff, and Daniel C. Jones. "Self‐shading and meltwater spreading control the transition from light to iron limitation in an Antarctic coastal polynya." Journal of Geophysical Research: Oceans 126, no. 2 (2021): e2020JC016636.
