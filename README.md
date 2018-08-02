# Linkage and QTL mapping of GBS data using Stacks2
## Analysis Approach
1. Assess <font face='consolas'>Stacks2</font> _de novo_ parameters for maximum coverage and minimum error rates (inferred from replicate samples)  
2. Combine replicate read files to one file (to increase coverage for them)
3. Call variants using <font face='consolas'>Stacks2</font> _de novo_ pipeline using optimal parameters  
4. Map catalog back to a reference genome (_Lens culinaris_ v1.2; UoS)
5. Integrate aligned tags into the catalog and match chromosomal position to each loci  
6. Filter output files (remove samples and loci with too much missing information)  
7. Perform linkage mapping and order markers inside each chromosome/linkage group using <font face='consolas'>ASMap</font> or <font face='consolas'>OneMap</font> R packages
8. QTL analysis and visualisation
9. Annotation of SNPs under the QTLs
