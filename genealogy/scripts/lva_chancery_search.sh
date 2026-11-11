#!/usr/bin/env bash
# Query the Library of Virginia Chancery Records Index (public, no login) and list case-detail links.
# Usage: bash lva_chancery_search.sh SURNAME [LOCALITY] [enslaved]
#   SURNAME   value for the index "Surname" field (enslaved people are indexed by given name, e.g. Woodley)
#   LOCALITY  e.g. "BRUNSWICK CO" (default All)
#   enslaved  pass the word "enslaved" to restrict to cases indexed with enslaved individuals
# Case images: https://old.lva.virginia.gov/chancery/case_images.asp?CFN=<CFN> lists per-page PDFs under /ChanceryImg/.
set -euo pipefail
SN="$1"; LOC="${2:-All}"; ENS=(); [[ "${3:-}" == enslaved ]] && ENS=(-d chkEnslaved=on)
curl -sS -L --max-time 90 -e https://old.lva.virginia.gov/chancery/ --data-urlencode "LOC=$LOC" \
  -d plselect=Exact -d pfselect=Exact -d dlselect=Exact -d dfselect=Exact -d sn1select=Exact \
  --data-urlencode "surname1=$SN" -d sn2select=Exact "${ENS[@]}" -d intPage=1 -d Submit=Search \
  "https://old.lva.virginia.gov/chancery/default.asp" | grep -o 'case_detail.asp?CFN=[0-9-]*' | sort -u \
  | sed 's|^|https://old.lva.virginia.gov/chancery/|'
