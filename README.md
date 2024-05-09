# Differential Gene Expression and Preliminary ORA

Welcome to the repository for Assignment 2 of the Bioinformatics course at BCB420. This project involves analyzing normalized gene expression data from a previous assignment, ranking genes by differential expression, and conducting thresholded over-representation analysis (ORA) to identify key biological themes.

[Link to original assignment guidelines](https://github.com/bcb420-2020/General_Course_Info/wiki/Assignment%232)

[Link to original jounral entry](https://github.com/bcb420-2022/Sabbir_Hossain/wiki/Journal-Entry-Assignment-%232:--Differential-Gene-expression-and-Preliminary-ORA)

[Link to rendered HTML document for this project](https://html-preview.github.io/?url=https://github.com/itsSabbir/DifferentialGeneexpressionAndPreliminaryORA/blob/main/A2_Sabbir_Hossain.html)

## Objective

The main objective of this assignment is to build on the normalized expression data from Assignment 1 by performing differential gene expression analysis and ORA. This standalone project emphasizes extracting meaningful biological insights from the data and presenting them comprehensively.

## Project Timeline

- **Estimated Time**: 12 hours
- **Time Used**: 18 hours
- **Project Start Date**: 2022/03/10
- **Project Completion Date**: 2022/04/16

## Repository Overview

This repository contains the R Notebook and all associated outputs for the analysis conducted. The project is documented through detailed journal entries and the notebook itself, providing insights into the decision-making process and technical details.

### Key Features

- Differential expression analysis using normalized data from Assignment 1
- Ranking of genes by p-values and adjustment using multiple hypothesis correction methods
- Thresholded ORA to highlight significant biological themes
- Visualization of results through various plots and heatmaps

## Differential Gene Expression

### Model Design
The analysis began by deciding on a model design for calculating differential expression. The design was based on an MDS plot which indicated a strong clustering, suggesting a significant relationship between genes and their expression under the conditions studied.

### Analysis Steps
1. **Fitting a Linear Model**: A linear model was fitted to the data to compute differential expression.
2. **P-value Calculation and Correction**: Genes were ranked by p-values, which were then adjusted using the Benjamini-Hochberg correction method to control the false discovery rate.
3. **Visualization**: Differentially expressed genes were visualized using volcano plots and MA plots. A heatmap was also generated to visualize the top hits and examine the clustering of conditions.

## Thresholded Over-Representation Analysis (ORA)

### ORA Execution
The significantly up-regulated and down-regulated genes were analyzed separately using gene set enrichment analysis tools. The analysis helped in identifying the dominant biological pathways affected by the conditions under study.

### Methodology
- **Gene Set Enrichment Tool**: g:Profiler was used due to its comprehensive database and user-friendly output.
- **Annotation Data**: Latest annotation data was employed to ensure accuracy in the pathway analysis.

## Results Interpretation

The results of the ORA were compared against the original study to check for consistency and relevance. The analysis supported the conclusions from the original study, providing a deeper understanding of the biological mechanisms involved in the conditions being studied.

## How to Use This Repository

1. **Clone the Repository**: Obtain all the files by cloning this repository to your local machine.
2. **Explore the R Notebook**: The notebook contains detailed code along with explanations for each step of the analysis.
3. **Review the Results**: Results are documented within the notebook and can also be viewed directly in the generated HTML output.

## Repository Contents

- **R Notebook (.Rmd and .html formats)**: Complete analysis workflow with detailed annotations.
- **Data/**: Contains processed datasets used in the analysis.
- **Figures/**: All plots and visual outputs generated during the analysis.
- **docs/**: Additional documentation and reference materials.

## References and Further Reading

This project utilizes a wide range of resources and literature to ensure a robust analysis framework. Key references include:

1. **GABAergic Function Markers Study**: [Link to study](https://www.frontiersin.org/articles/10.3389/fpsyt.2022.827972/full)
2. **R Programming for Data Science**: [Link to resource](https://bookdown.org/rdpeng/rprogdatascience/)
3. **Bioinformatics Tools and Packages**: Detailed in the R Notebook and cited appropriately in the analysis.
4. **Geistlinger et al., 2020**: Discusses benchmarks for gene set enrichment analysis. [PubMed](https://www.ncbi.nlm.nih.gov/pubmed/32026945)

## Conclusion

This README provides a comprehensive overview of the differential gene expression analysis and preliminary ORA conducted in this project. The analysis not only reinforced the findings from prior studies but also introduced new insights into the data, showcasing the power of bioinformatics in understanding complex biological datasets.

