These scripts take single base pair genomic locations and design an sgRNA library targeting these locations using CRISPick. Compatible with human hg19/GRCh37. See sgRNA Library Design pdf for underlying motivations for library design and QC.

---

To perform library design:

1) Run **preprocessCRISPick.Rmd**: preprocessing to generate regions for sgRNA design
    - From 'input' folder: input SNP locations to create SNP regions for CRISPick design input.
    - To create SNP "regions" for sgRNA design:
        - First, expands from SNP location by 200, 300, 400, 500, and 600 bp. 
        - Second, if any two or more regions are overlapping, they are combined into one larger region. i.e. if one or more SNPs are very closeby each other.
    - Exports regions to /preprocessed folder as regions###.xlsx for human readability and as loc###.fasta for CRISPick input.
2) Input regions into **CRISPick**:
    - Go to: https://portals.broadinstitute.org/gppx/crispick/public
    - Select Reference Genome: Human GRCh37. Pick your relevant mechanism, enzyme, and On Target Scorer.
    - Upload each regiion###.xlsx file of interest.
    - Set CRISPick Quota - typically 50 is comprehensive.
    - Save CRISPick outputs to the /crispick_output folder.
3) Run **preprocessGDO.Rmd**: converts to friendly format for library design.
    - From 'crispick_output' folder: takes raw CRISPick output file, extracts relevant columns, adds guide locations and row index. 
    - Saves processed files to the /preprocessed folder.
4) For the final library design, run: **targetDesign.Rmd**, contains five sections:
    - Set Directory Information + Assemble Library: performs basic library assembly
    - Remaining Sections: examples of adding or removing individual guides

libraryDesignFunctions.R contains functions to design the library and perform simple quality control of library designs.


Final library design files located in 'output' folder:
- finalPromoterGDOa30.xlsx: contains promoter targeting positive control guide designs.
- finalTargetGDOi30.xlsx: contains guides targeting candidate CREs.

---

To compare different library parameters, go through /analyses folder:
- pickComparison.Rmd: Can compare different libraries directly with simple QC.
- overlapLeniency.Rmd: Compares impact of different guide overlap leniencies. i.e. if we allow 2 vs 4bp of overlap between the guides in the library.

---

Additional scripts:
- adtlROIDesign.Rmd: Designs guides that tile along a large region rather than a 1bp SNP location.
- promoterDesign.Rmd: Incorporates list of promoter-targeting guides from CRISPick into your final library. Performs some QC to ensure no guides are overlapping with existing guides in the library. Picks top performing guides.
- postprocessToOrder.Rmd: Postprocessess library of interest into format for ordering library from Cellecta.