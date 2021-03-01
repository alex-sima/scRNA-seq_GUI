INTRODUCTION
------------

An R Shiny App for visualizing scRNA sequencing data. Find and visualize differentially expressed or user-inputted genes under certain treatments.

Files
 * pairWiseTests.rds: R object with the coefficient and p values of tests on genes that are significant after cytokine treatment (not publicly available)
 * treatments.so.rds: R Seurat object with QC, Normalization, and PCA performed on it. (not publicly available) Cells have been clustered and non-linear dimension reductions (UMAP/tSNE) have been run. See https://github.com/satijalab/seurat for more info about Seurat.
 * my_data.rds: Scaled R Seurat object for HeatMap visualization. (not publicly available)
 * MarKDOWn.png: Image of documentation containing detailed instructions on how to use the app (can also be found on the "Documentation" tab of the application)
 * app.R: R Shiny application source code.

REQUIREMENTS
------------

This application requires:

 * R version 4.0+ (https://www.r-project.org)
 * RStudio (https://rstudio.com/products/rstudio/download/)
   - RStudio is recommended, but you can also launch the Shiny app in Terminal
 * The following R packages (see installation for more details):
   - shiny, shinyEventLogger, Seurat, ggplot2, dplyr, grid, shinycssloaders, shinythemes, DT 

Recommended general system requirements/specs:
 * OS: Windows 10 (64-bit) or Mac OS X 10.13+
 * Processor (CPU): Intel Core i5 (sixth generation or newer) or equivalent
 * Memory: 8 GB RAM

INSTALLATION (not publicly available yet)
------------
 
 * Download the entire scRNA sequencing Shiny App box folder
 * Open RStudio (or Terminal) 
 * If not installed yet, install the required R packages (run the following commands in the RStudio console, or launch R and run the commands in terminal):
   - Shiny: run "install.packages("shiny")"
   - devtools: run "install.packages("devtools")"
   - shinyEventLogger: run "devtools::install_github("kalimu/shinyEventLogger")"
     * shinyEventLogger has to be installed from GitHub because it's no longer hosted on CRAN
   - Seurat: run "install.packages("Seurat")"
   - ggplot2: run "install.packages("ggplot2")"
   - dplyr: run "install.packages("dplyr")"
   - grid: run "install.packages("grid")"
   - shinycssloaders: run "install.packages("shinycssloaders")"
   - shinythemes: run "install.packages("shinythemes")"
   - DT: run "install.packages("DT")"
 * Navigate to the "app.R" file in the "scRNA_shiny" sub-folder. Open the file. (In terminal, navigate to the Shiny App box folder)
 * In the top right corner of the code screen, click the "Run App" button. Alternatively, you can also run "runApp("scRNA_shiny/")" in the console or in terminal
 * An event in the console should say, "Hello World!" after the logging statements. This indicates that the app has finished launching successfully. The app should open in a new window.

CONTACT
-------

 * For bug reports, suggestions, or any further questions, contact alex.sima@ucsf.edu

  
