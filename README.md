### 🌎 Population Genomics &nbsp; &nbsp; &nbsp; 🔍 Single Nucleotide Variants &nbsp; &nbsp; &nbsp; 📈 Introgression

Hybridization — the interbreeding of individuals from populations that differ by one or more heritable characters — can be studied to better understand the evolutionary processes associated with speciation, as well as the way species are responding to environmental change. Black-capped and Carolina Chickadees hybridize in a narrow hybrid zone that extends from New Jersey to Kansas, and that is moving northwards at a rate of ~1 km per year (Taylor et al. 2014; Wagner et al. 2020). 

![Population Map](https://github.com/kbfeldmann/WGS_chickadee_hybridization/assets/47021794/7da254e1-519b-4704-b11f-3acd3c02c8b0)  
**Figure 1:** Map of the Black-capped (blue) and Carolina (red) Chickadee range distributions. The hybrid zone, represented by the area of range overlap (purple), spans from New Jersey to Kansas. Points represent the populations sampled: allopatric Black-capped Chickadees (DNY, INY and HRP; blue circles), allopatric Carolina Chickadees (ECW and LSU; red squares), and sympatric populations (JSP, LHU and DSU; purple diamonds).

Hybrid chickadees have deficient learning and memory abilities compared to either parental species (McQuillan et al. 2018), suggesting that cognition may be a postzygotic reproductive isolating barrier — where hybrids with deficient cognitive abilities experience negative selection. To determine if cognition is acting as a post-zygotic reproductive isolating barrier in chickadees, we used high-resolution, whole-genome data to examine patterns of introgression across the hybrid zone and determined the biological processes associated with loci experiencing reduced introgression. We also used whole genome data to determine chickadee ancestry across our hybrid zone transect and compared hybrid zone position to previous studies.

**Research Question:** *Does cognitive breakdown in hybrid offspring act as a postzygotic, reproductive isolating barrier in a rapidly moving hybrid zone?*

We found 1) that chickadee ancestry and geographic cline analysis indicated continued northward movement of the chickadee hybrid zone in Pennsylvania (Taylor et al. 2014; Wagner et al. 2020), 2) that genomic cline analysis revealed reduced introgression of loci on the Z chromosome, and 3) that loci experiencing reduced introgression across the hybrid zone are related to cognitive and metabolic function.

To learn more, check out my graduate-level thesis: [click here](https://github.com/kbfeldmann/WGS_chickadee_hybridization/blob/main/Masters_Thesis.pdf)

## Bioinformatics Pipeline

```
# Bioinformatics workflow to identify single-nucleotide variants.
trim_QC.sh   # Trim Reads
align_sort.sh   # Align Reads to Reference Genome
merge_bams.sh   # Merge Bam Files
call_variants.sh   # Call Variants
```
![Bioinformatics Pipeline](https://github.com/kbfeldmann/WGS_chickadee_hybridization/assets/47021794/20600853-ff17-431b-b362-2340e6436e2f)
**Figure 2:** Brief description of the bioinformatics pipeline developed to call single-nucleotide variants. Filenames for each step (excluding Raw Data) can be found in the above code block.

## Whole Genome Analyses

```
# Different programs used to analyse the whole-genome single nucleotide variants.
/HZAR
/GGHYBRID
/ADMIXTURE
/OUTFLANK
/PCA

# Additional programs explored, but not used for thesis.
/MANHATTAN
/HYBRIDDETECTIVE
```

![WGA1](https://github.com/kbfeldmann/WGS_chickadee_hybridization/assets/47021794/9073eac1-3adb-4934-936d-0659fd7fc9f2)
![WGA2](https://github.com/kbfeldmann/WGS_chickadee_hybridization/assets/47021794/24bbc978-2e77-41c1-b0a2-c18e8b6a64d3)
**Figure 3:** Brief description of the programs explored and implemented for whole-genome analyses. Programs shown in blue were used for thesis research. Programs shown in black were explored, but not used in thesis.

Presentation for Faculty-Level Ecologists and Evolutionary Biologists: insert link here
