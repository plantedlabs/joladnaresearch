JOLA Y-DNA HAPLOGROUP SCANNER
Drumgoole Family Research Project
===================================

WHAT THIS DOES
--------------
This script checks 40 Jola males from the Gambian Genome Variation Project
whose Y-DNA haplogroup has never been examined. It looks for the specific
SNPs that define Kelland Drumgoole's Jola paternal lineage (E-Z15174).

Any positive hits would be new members of the Jola family — potentially
even on the exact E-Z15174* branch which currently has only Kelland as
a known member.

REQUIREMENTS
------------
1. Windows 10 or 11 with WSL2 (Ubuntu) installed
   OR a Linux/Mac computer
2. Internet connection (accesses EBI FTP server)
3. samtools installed

SETUP INSTRUCTIONS
------------------

Step 1 — Enable WSL2 (Windows only)
Open PowerShell as Administrator and run:
  wsl --install
Restart your computer.

Step 2 — Open Ubuntu terminal
Click Start → search "Ubuntu" → open it

Step 3 — Install samtools
  sudo apt-get update && sudo apt-get install -y samtools python3

Step 4 — Copy this script to Ubuntu
Either drag and drop the file, or run:
  cp /mnt/c/Users/[YourName]/Downloads/check_jola_males.py ~/

Step 5 — Run the script

TEST FIRST (3 samples, quick):
  python3 check_jola_males.py --test

CHECK UNCHECKED FATHERS ONLY (3 samples):
  python3 check_jola_males.py --fathers

RUN FULL SCAN (40 samples — takes 1-2 hours):
  python3 check_jola_males.py

RUN EVERYTHING (43 samples):
  python3 check_jola_males.py --all

READING THE RESULTS
-------------------
The script will print results as it runs and save a report file.

POSSIBLE_DERIVED = This man may carry your haplogroup SNP — IMPORTANT
ANCESTRAL       = Reference allele only — not your haplogroup
LOW_COVERAGE    = Not enough reads at this position to call
TIMEOUT         = Network issue — try again
ERROR           = File access problem

A POSSIBLE_DERIVED result at Z15174 would mean that man is potentially
on the E-Z15174 branch — your Jola lineage. Multiple POSSIBLE_DERIVED
results across several SNPs strongly suggests a true match.

WHAT TO DO WITH RESULTS
-----------------------
Share any POSSIBLE_DERIVED results with:
- Your YFull page (YF133256) — to see if they cluster with your branch
- Ellen Leffler (leffler@genetics.utah.edu)
- Mamadou Jallow (WhatsApp)

This would be genuine new scientific discovery — no researcher has
done this analysis on this dataset before.

SNP POSITIONS CHECKED (GRCh38)
-------------------------------
Z15174  chrY:4,202,358   (defines E-Z15174 branch — your Jola connection)
Z15234  chrY:8,414,659   (confirms E-Z15174 placement)
Z15231  chrY:14,719,978  (downstream E-Z15174 marker)
Z15232  chrY:14,765,245  (downstream E-Z15174 marker)

CONFIRMED MATCHES (already on YFull tree, not re-checked)
----------------------------------------------------------
SC_GMJOL5309804  E-Z15196 / Jola-Fonyi (dyo)
SC_GMJOL5309829  E-Z15196 / Jola-Fonyi (dyo)
SC_GMJOL5309851  E-FT402211

CONTACT
-------
Research log maintained with Claude (Anthropic)
YFull ID: YF133256
Haplogroup: E-FTF75935 / E-Z15174*
