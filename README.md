# Jola DNA Research
## Tracing African American Paternal Ancestry to the Jola-Fonyi People of The Gambia

**Researcher:** Kelland Drumgoole, Louisiana, USA  
**YFull ID:** YF133256  
**Haplogroup:** E-FTF75935 / E-Z15174* (downstream of E-CTS246 → E-Z906 → E-Z958 → E-Z15174)  
**Repository:** github.com/plantedlabs/joladnaresearch  

---

## Overview

This repository documents an African American genealogy research project that traced a direct paternal line from Louisiana, USA, back to the **Jola-Fonyi people of the Western Division of The Gambia, West Africa** — a connection severed by the transatlantic slave trade in the late 1700s.

Using a combination of:
- **Big Y-700 DNA testing** (FTDNA)
- **YFull whole-genome haplogroup analysis**
- **Systematic screening of publicly available genomic data** (Gambian Genome Variation Project, PRJEB3252)
- **Documentary genealogy** (census records, slave schedules, plantation records)

...this project identified living Jola-Fonyi men in The Gambia as the closest known genetic relatives of an African American family whose paternal ancestor was enslaved and transported to Virginia approximately 250 years ago.

---

## Key Finding

**On September 5, 2026, this project completed the first systematic Y-DNA haplogroup screen of all male samples in the Gambian Genome Variation Project (GGVP) Jola dataset.**

Result: **E-Z15174 is present in exactly 3 out of 47 Jola males sampled from the Western Division of The Gambia — a frequency of 6.4%.**

Those 3 men are confirmed genetic relatives of Kelland Drumgoole on the direct paternal line. No other man in the entire GGVP Jola dataset shares this ancient lineage.

---

## Confirmed Paternal Lineage

| Generation | Person | Birth | Location |
|------------|--------|-------|----------|
| You | Kelland Drumgoole | — | Louisiana, USA |
| Father | Leroy Drumgoole Jr | — | Louisiana |
| Grandfather | Leroy Drumgoole Sr | Jan 24, 1932 | Pine Bluff, Jefferson Co., Arkansas |
| Great Grandfather | Arthur Drumgoole | 1859 | Catahoula Parish, LA → Jefferson Co., AR |
| 2x Great Grandfather | Woodley Drumgoole | ~1809 | Virginia → Trinity, Catahoula Parish, LA |
| 3x Great Grandfather | Unknown | ~1770s | Brunswick County, Virginia |
| **4x Great Grandfather** | **Jola-Fonyi ancestor** | **~1740s** | **Western Division, The Gambia** |

---

## Haplogroup — E-Z15174*

E-Z15174 is one of the rarest and most geographically specific Y-DNA haplogroups documented in African diaspora research.

**Age:** Approximately 7,000-11,000 years before present  
**Origin:** West Africa, Senegambia region  
**Current distribution:** Almost exclusively found among Jola-Fonyi (dyo) and Mandinka (mnk) speakers in the Western Division of The Gambia

The branch structure:
```
E-M96
  └── E-P147
      └── E-M132
          └── E-Z965
              └── E-CTS246
                  └── E-Z906
                      └── E-Z958
                          └── E-Z15174 (formed ~11,000 ybp, TMRCA ~7,000 ybp)
                              ├── E-Z15174* ← KELLAND DRUMGOOLE (YF133256)
                              │   └── E-FTF75935 (private terminal SNP)
                              └── E-Z5991 (formed ~7,000 ybp, TMRCA ~2,300 ybp)
                                  ├── E-FT402211 (TMRCA ~375 ybp)
                                  │   ├── SC_GMJOL5309851 🇬🇲 Jola-Fonyi
                                  │   └── SC_GMJOL5309804 🇬🇲 Jola-Fonyi (dyo)
                                  └── E-Z15196 (TMRCA ~1,000 ybp)
                                      ├── HG03027 🇬🇲 Mandinka (mnk)
                                      ├── HG03024 🇬🇲 Mandinka (mnk)
                                      └── SC_GMJOL5309829 🇬🇲 Jola-Fonyi (dyo)
```

**Key insight:** The shared common ancestor between Kelland's branch (E-Z15174*) and the Jola trio branches (E-Z5991) lived approximately **7,000 years ago** — making this one of the oldest documented paternal connections between an African American and a specific West African ethnic community.

---

## Confirmed Jola Y-DNA Matches

Three Jola-Fonyi men from the Gambian Genome Variation Project are confirmed Y-DNA relatives:

| Sample ID | Population | Language | Haplogroup | Trio Role |
|-----------|-----------|----------|------------|-----------|
| SC_GMJOL5309804 | GWJ (Gambian Jola) | dyo (Jola-Fonyi) | E-Z15196 | Father, Trio 1 |
| SC_GMJOL5309829 | GWJ (Gambian Jola) | dyo (Jola-Fonyi) | E-Z15196 | Father, Trio 2 |
| SC_GMJOL5309851 | GWJ (Gambian Jola) | — | E-FT402211 | Father, Trio 3 |

All three samples are from the GGVP study:  
**ENA Study:** PRJEB3252 / ERP001781  
**Title:** Low coverage sequencing of the Jola from Gambia Western Division (GWJ)

---

## Original Research — GGVP Y-DNA Haplogroup Screen

### What Was Done

On September 5, 2026 a systematic mpileup analysis was conducted on all 47 male samples in the GGVP Jola dataset to screen for haplogroup E-Z15174.

**Method:**
- Accessed CRAM alignment files via EBI FTP server
- Used samtools mpileup with GRCh38 reference genome
- Checked 4 haplogroup-defining SNP positions per sample:

| SNP | GRCh38 Position | Ref Allele | Derived Allele |
|-----|-----------------|------------|----------------|
| Z15174 | chrY:4,202,358 | G | T |
| Z15234 | chrY:8,414,659 | C | T |
| Z15231 | chrY:14,719,978 | G | A |
| Z15232 | chrY:14,765,245 | T | C |

**Result:** Only 3 of 47 Jola males carry E-Z15174 (6.4%). These are the same 3 samples already confirmed as YFull matches. No additional E-Z15174 carriers were found in the dataset.

**Scientific significance:** E-Z15174 is not a broadly distributed Jola haplogroup — it is specific to a small subset of Jola families in the Western Division. This suggests the three matched families likely come from the same village or extended clan.

### How To Replicate

**Requirements:**
- Linux/WSL2 environment
- samtools installed (`sudo apt-get install samtools`)
- Internet access to ftp.1000genomes.ebi.ac.uk

**Run the screen:**
```bash
# Download chrY reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrY.fa.gz
gunzip chrY.fa.gz
samtools faidx chrY.fa

# Run the full scan
chmod +x run_all.sh
./run_all.sh
```

Results are saved to `all_results.txt`.

---

## Documentary Evidence

### Key Records

- **1870 Census:** "W Brungood," Ward 8, Trinity, Catahoula Parish. Age 61, born Virginia. House Carpenter. *(Woodley Drumgoole)*
- **1880 Census:** Woodley Drumgoole, age 76, 1st Ward Catahoula Parish. Farmer.
- **Dromgoole Plantation:** "Canaan" plantation, Brunswick County, Virginia. 1854 slave inventory held at Duke University Libraries.

### Virginia Connection

The Dromgoole family of Brunswick County, Virginia enslaved people in the late 1700s and early 1800s. Woodley Drumgoole (born ~1809 Virginia) is documented in Louisiana records from 1870. His father — the Jola ancestor — was transported to Virginia between approximately 1770 and 1800.

---

## Researcher Outreach

| Researcher | Institution | Status |
|-----------|-------------|--------|
| Mamadou Jallow | MRC Unit The Gambia (former) / AQA Medical Diagnostics | In contact via WhatsApp |
| Fatou Sissay-Joof | MRC Unit The Gambia | Emailed |
| Prof. Umberto d'Alessandro | MRC Unit The Gambia (former Director) | Responded — referred to Data Management |
| Prof. Ed Clarke | MRC Unit The Gambia (Director) | CC'd by d'Alessandro |
| Dr. Ellen Leffler | University of Utah | Responded — actively investigating |
| Sanger Data Sharing | Wellcome Sanger Institute | Responded — referred to HG programme |
| Laura Allen | Sanger Human Genetics | Responded — added P&M team |
| Vikki | Sanger Parasites & Microbes | Responded — confirmed MRC are data owners |
| MalariaGEN | MalariaGEN Consortium | Emailed |

---

## Scientific Publications Using This Data

- **Leffler et al. (2017)** — Resistance to malaria through structural variation of red blood cell invasion receptors. *Science*, 356: eaam6393. [PMID: 28522690](https://www.ncbi.nlm.nih.gov/pubmed/28522690)
- **Band et al. (2019)** — New insights into malaria susceptibility from genomes of 17,000 individuals from Africa, Asia, and Oceania. *Nature Communications*.
- **Provyn, H. (2024)** — My E-CTS246/Z906 Client's Rare Gambian-related 7,000-10,000 Year Old Haplogroup. PhyloGeographer/Mygrations. [Link](https://phylogeographer.com/my-e-cts246-z906-clients-rare-gambian-related-7000-10000-year-old-haplogroup-based-on-clade-finder-and-str-analysis/)

---

## Data Sources

All genomic data used in this project is publicly available:

- **GGVP CRAM files:** http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/gambian_genome_variation_project/data/JOLA/
- **GGVP VCF calls:** http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/gambian_genome_variation_project/release/20200217_biallelic_SNV/
- **ENA Study:** https://www.ebi.ac.uk/ena/browser/view/PRJEB3252
- **IGSR Sample Portal:** https://www.internationalgenome.org/data-portal/sample/SC_GMJOL5309829
- **YFull Tree:** https://www.yfull.com/tree/E-Z15174/

---

## Acknowledgments

This research would not be possible without:

- The Jola families of the Western Division of The Gambia who volunteered blood samples for the GGVP study
- Mamadou Jallow, Fatou Sissay-Joof, and Umberto d'Alessandro (MRC Unit The Gambia) for sample collection
- The Wellcome Sanger Institute for sequencing and data release
- MalariaGEN and the IGSR for data hosting
- Hunter Provyn (PhyloGeographer) for professional haplogroup analysis
- Dr. Ellen Leffler (University of Utah) for ongoing assistance

*"The authors would like to thank Muminatou Jallow, Fatou Sissay-Joof, Umberto d'Alessandro from the Gambia Genome Variation Project (MRC-LSHTM Unit in Fajara, The Gambia); Jim Stalker, Katja Kivinen, Eleanor Drury, Kirk Rockett, Dominic Kwiatkowski plus members of the DNA pipelines Teams from The Wellcome Sanger Institute; Anna Jeffreys, Kate Rowland, Christina Hubbart, Christopher Spencer, Gavin Band, Quang Si Le, Ellen Leffler, Kirk Rockett and Dominic Kwiatkowski from the MRC Centre for Global Genomics and Health, University of Oxford; and all the participants who volunteered blood samples for this project."*

---

## Contact

**Kelland Drumgoole**  
Louisiana, USA  
YFull ID: YF133256  
GitHub: github.com/plantedlabs/joladnaresearch

---

## License

This methodology and code is released under MIT License for use by other African American genealogy researchers conducting similar ancestral reconnection research.

The genomic data referenced belongs to its respective owners and is subject to their terms of use.
