# reliable-network-inference

Data and code for "Reliable Network Inference From Unreliable Data: A Tutorial on Latent Network Modeling Using STRAND", published in Psychological Methods

----------------------------
# The manuscript can be found at:

- Psychological Methods: https://psycnet.apa.org/record/2023-51200-001
- psyArxiv (preprint): https://psyarxiv.com/mkp2y/

# Requirements for analyses:

- R: https://cran.r-project.org
- STRAND: https://github.com/ctross/STRAND
- cmdstanr: https://mc-stan.org/cmdstanr/

# Packages used for data processing and visualisation:

- Rethinking: https://xcelab.net/rm/statistical-rethinking/
- iGraph: https://igraph.org/r/
- Rcolorbrewer: https://cran.r-project.org/web/packages/RColorBrewer/index.html
- tidyverse: https://www.tidyverse.org

To reproduce the model validation parameter sweeps and example analyses using empirical data, please go into the 'Code/' folder of the repository. It would probably be best to first review the different scripts for the many parameter sweeps. Then, if you would like to reproduce the results reported in the publication, call the "run_all.R" file that can be found in "Code/". This can be done by simply calling:

``````````
source("Code/run_all.R")
``````````

However, please note that you will need to specify the correct working directory in the file. It is important to note that all analyses were run on the high performance computational cluster available in the Department of Human Behavior, Ecology and Culture at the Max Planck Institute for Evolutionary Anthropology--and they take a long time to finish. 

The project is maintained by Cody Rosee (cody_ross@eva.mpg.de) and Daniel Redhead (daniel_redhead@eva.mpg.de) and STRAND is hosted at https://github.com/ctross/STRAND.
