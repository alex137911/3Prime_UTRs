# 3 Prime UTR Characterization
Sample code from summer internship in the [Computational Rare Disease Group](https://rarediseasegenomics.org/) at the University of Oxford.

### A. preliminaryAnalysis.R

Preliminary analysis of 3' UTR composition in MANE v1.1 transcripts.
- 521 393 transcripts in **MANE v1.1**, 19 735 canonical transcripts (18 872 genes) with "three_prime_UTR" tag
- 894 canonical transcripts (834 genes) with 3' UTR introns

<ins>Colour palette:</ins>

![#BAE4C4](https://placehold.co/15x15/BAE4C4/BAE4C4.png) `#BAE4C4`, ![#B7E3D4](https://placehold.co/15x15/B7E3D4/B7E3D4.png) `#B7E3D4`, ![#A4DCD7](https://placehold.co/15x15/A4DCD7/A4DCD7.png) `#A4DCD7`, ![#89D3D3](https://placehold.co/15x15/89D3D3/89D3D3.png) `#89D3D3`, ![#89C1D3](https://placehold.co/15x15/89C1D3/89C1D3.png) `#89C1D3`, ![#4EB3D3](https://placehold.co/15x15/4EB3D3/4EB3D3.png) `#4EB3D3`, ![#2B8CBE](https://placehold.co/15x15/2B8CBE/2B8CBE.png) `#2B8CBE`, ![#0868AC](https://placehold.co/15x15/0868AC/0868AC.png) `#0868AC`, ![#094081](https://placehold.co/15x15/094081/094081.png) `#094081`, ![#042042](https://placehold.co/15x15/042042/042042.png) `#042042`

<ins>Tasks:</ins>
- [x] Calculate 3' UTR exon count per gene
- [x] Calculate total (cumulative) exon length per gene
- [x] Calculate 3' UTR intron count per gene


### 3' UTR Composition in MANE Transcripts 
![image](https://github.com/alex137911/3Prime_UTRs/assets/78928747/31abcd79-1024-4206-a7f8-6fb99eb91280)


### B. 3primeUTR_LOEUF.R

Investigating the differences in 3' UTR composition between genes that are tolerant/intolerant to loss-of-function mutation.
- Uses **gnomAD v2.1.1**[^1]
- 17 058 genes with 3' UTR and LOEUF score (merge using gene_name and selecting canonical transcripts)

<ins>Tasks:</ins>
- [x] Calculate LOEUF deciles from gnomAD oe_lof_upper column
- [x] Model 3' UTR exon and intron length across LOEUF deciles
- [x] Run linear regression on median length of 3' UTR exons and introns across LOEUF deciles

[^1]: Note that gnomAD v2.1.1 uses hg19 build, while MANE v1.1 uses hg38 genome assembly


### 3' UTR Exon Length Across LOEUF Deciles
![image](https://github.com/alex137911/3Prime_UTRs/assets/78928747/91803b6f-e87e-4da3-b4d4-f676c1f23b01)


### C1. cDNAconverter.R

Modified Ruby's code to convert 3' UTR genomic coordinates to cDNA coordinates. Also determines the presence of stop introns and calculates their coordinates.
- Uses **Ensembl v110**
- A compiled list of all (_n = 31_) stop introns and their coordinates is available in the Data folder of this repository


### C2. 3primeUTR_NMD.R

Determines transcript likelihood to undergo/escape nonsense-mediated mRNA decay (NMD) based on the **exon junction complex** (EJC) **dependent model** for governing NMD; that is, if a premature termination codon (PTC) is located more than 55 bp upstream from the last exon-exon junction, it is predicted to elicit NMD and degrade the mRNA.[^2] Similarly, **it is hypothesized that EJCs within the 3' UTR may cause the original stop codon to be considered as a PTC if located more than 55 bp away downstream**.
- 834 genes with 3' UTR introns (thus, 834 genes which introduce 3' UTR EJCs)

<ins>Colour palette:</ins>

![#EE5496](https://placehold.co/15x15/EE5496/EE5496.png) `#EE5496`, ![#FFC000](https://placehold.co/15x15/FFC000/FFC000.png) `#FFC000`, ![#78DAD5](https://placehold.co/15x15/78DAD5/78DAD5.png) `#78DAD5`, ![#5D576B](https://placehold.co/15x15/5D576B/5D576B.png) `#5D576B`  

<ins>Tasks:</ins>
- [x] Calculate 3' UTR intron length
- [x] Determine EJC position
- [x] Calculate EJC distance from stop codon
- [x] Classify genes as likely to escape NMD (_Category 1_, all EJCs < 55 bp from stop codon), likely to undergo NMD (_Category 2_, all EJC > 55 bp from stop codon), or _Category 3_ (where the gene has 3' UTR EJCs in both _Category 1_ and _2_)
- [x] Cross-reference genes which are impacted by NMD with NHS PanelApp green genes and LOEUF scores

[^2]: M. W. Popp and L. E. Maquat, “Leveraging Rules of Nonsense-Mediated mRNA Decay for Genome Engineering and Personalized Medicine,” Cell, vol. 165, no. 6, pp. 1319–1322, Jun. 2016, doi: 10.1016/j.cell.2016.05.053.


### Calculating EJC Position Relative to CDS End
![image](https://github.com/alex137911/3Prime_UTRs/assets/78928747/0fd099f8-bc8c-41bf-8385-f488b15dffa5)


### D. 3primeUTR_dORF.R

Investigating dORF composition across 3' UTRs.

<ins>Tasks:</ins>
- [x] Determine number of dORFs across NMD categories
- [x] Determine dORF count across LOEUF deciles
- [ ] Calculate distance between dORFs start and end to 3' UTR EJC positions

### E. 3primeUTR_diseaseGenes

Compare 3' UTR composition across different disease gene groups.
- Uses **Gene2Phenotype** developmental disorder dataset
- Uses **COSMIC Cancer Gene census** dataset
- Uses gene dosage sensitivity dataset from: _R. L. Collins et al., “A cross-disorder dosage sensitivity map of the human genome,” Cell, vol. 185, no. 16, pp. 3041-3055.e25, Aug. 2022, doi: 10.1016/j.cell.2022.06.036_.
- Uses **NHS PanelApp** green genes


### Total 3' UTR Exon Length Across Disease Genes
![image](https://github.com/alex137911/3Prime_UTRs/assets/78928747/1d3f40c9-f4d7-4a13-82d7-a311536cc82a)

