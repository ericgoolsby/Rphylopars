## Test environments
* Local Windows 8 install, R 3.1.2
* win-builder (devel and release)

## R CMD check results
* There were no ERRORs or WARNINGs. 
* There were 2 NOTES on the local Windows check
    * checking CRAN incoming feasibility ... NOTE
        * Maintainer: 'Eric W. Goolsby <eric.goolsby.evolution@gmail.com>'
          New submission
            * EXPLANATION: This is a new submission
    * checking package dependencies ... NOTE
        * No repository set, so cyclic dependency check skipped
            * EXPLANATION: OK?
* There was 1 NOTE on the win-builder (devel and release) check
    * checking CRAN incoming feasibility ... NOTE
        * Maintainer: 'Eric W. Goolsby <eric.goolsby.evolution@gmail.com>'
        New submission
            * EXPLANATION: This is a new submission
        * Possibly mis-spelled words in DESCRIPTION:
        Phylogenetic (3:8)
        datasets (7:72)
        intraspecific (7:126)
        pPCA (7:390)
        phylogenetic (7:35, 7:276, 7:301, 7:346)
            * EXPLANATION: This spelling is correct