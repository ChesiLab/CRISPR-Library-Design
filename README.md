File Information:

To perform library design, first run preprocessing scripts to generate guide pools:
- preprocessCRISPick.Rmd
    - Input SNP locations to create SNP regions for CRISPick design input.
- preprocessGDO.Rmd
    - Takes raw CRISPick output file, extracts relevant columns, adds guide locations and row index. 

For the final library design, run:
- finalDesign.Rmd

To compare different library designs, go through:
- pickComparison.Rmd

libraryDesignFunctions.R contains functions to design the library and also perform some simple quality control of library designs.

Final library design files located in 'output' folder:
- finalPromoterGDOa30.xlsx: contains promoter targeting positive control guide designs.
- finalTargetGDOi30.xlsx: contains guides targeting candidate CREs.
