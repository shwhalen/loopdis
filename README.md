# LD vs. Hi-C Pipeline

This repository contains code for reproducing analyses from the paper "[Most regulatory interactions are not in linkage disequilibrium](https://doi.org/10.1101/272245)", currently in review.  The resulting files are 622 gigabytes and so cannot be bundled with the source code.  In addition, the pipeline relies on integrating large amounts of data from other published studies, organized in a particular hierarchy, that also cannot be feasibly bundled.  Finally, the pipeline requires up to 80 gigabytes of RAM and weeks of compute time to complete from start to finish.

Thus, the aim of this repository is to allow auditing our analyses, rather than provide a tool or end-to-end-pipeline.  We are happy to arrange a transfer of our pre-generated results to interested parties.

These analyses utilize several command line tools, as well as R and Python libraries.  A full list of tools and the versions we used are provided in the Methods section of the paper.  Prior to posting our first draft, all analyses were run from scratch using the latest version of all required tools and libraries.  Python 3.6.4 was provided by the [Miniconda](https://conda.io/miniconda.html) distribution; R 3.4.3 was compiled from source using gcc 7.2.1.
