#!/usr/bin/env bash
# Fast SNP->gene annotator using bedtools + GENCODE v19 (GRCh37/hg19)
#
# Usage:
#   ./Step1_annotate_genes_bedtools.sh \
#       --input /path/to/GWAS_file.txt.gz \
#       [--gtf /path/to/gencode.v19.annotation.gtf.gz] \
#       [--extend 10000]
#
# Requirements (in GenePrioritiser_env):
#   bedtools, sort (coreutils), gzip/zcat, python (with pandas installed)

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)

# Defaults
GTF_DEFAULT="$REPO_ROOT/utils/gencode.v19.annotation.gtf.gz"
EXTEND=10000

print_usage() {
  echo "Usage: $0 --input <gwas.txt.gz> [--gtf <gencode.gtf.gz>] [--extend <bp>]"
  echo "Produces: Annotated_GWAS_<phenotype>.csv and variant_data_<phenotype>.csv"
  echo "          in config.variant_output_directory"
}

# ----------------------------
# Parse CLI arguments
# ----------------------------
INPUT=""
GTF_PATH=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input|-i)
      INPUT="$2"; shift 2;;
    --gtf)
      GTF_PATH="$2"; shift 2;;
    --extend)
      EXTEND="$2"; shift 2;;
    --help|-h)
      print_usage; exit 0;;
    *)
      echo "Unknown arg: $1" >&2
      print_usage
      exit 1;;
  esac
done

# ----------------------------
# Resolve GWAS inputs
# ----------------------------
if [[ -z "$INPUT" ]]; then
  CONFIG_SH="$REPO_ROOT/config/config.sh"
  if [[ ! -f "$CONFIG_SH" ]]; then
    echo "Error: no --input provided and $CONFIG_SH not found" >&2
    print_usage
    exit 2
  fi
  # shellcheck disable=SC1090
  source "$CONFIG_SH"
  if [[ -z "${GWAS_PATHS+x}" ]]; then
    echo "Error: GWAS_PATHS not defined in $CONFIG_SH" >&2
    exit 2
  fi
  INPUTS=("${GWAS_PATHS[@]}")
else
  INPUTS=("$INPUT")
fi

# ----------------------------
# Resolve GTF path
# ----------------------------
if [[ -z "$GTF_PATH" ]]; then
  GTF_PATH="$GTF_DEFAULT"
fi

if ! command -v bedtools >/dev/null 2>&1; then
  echo "Error: bedtools not found in PATH. Install via 'conda install -c bioconda bedtools'" >&2
  exit 3
fi

if [[ ! -f "$GTF_PATH" ]]; then
  echo "Error: GTF not found at $GTF_PATH" >&2
  echo "Please download GENCODE v19 (GRCh37) GTF to $GTF_PATH or pass --gtf <path>" >&2
  exit 4
fi

# ----------------------------
# Determine output directory from Python config
# ----------------------------
OUT_DIR=$(PYTHONPATH="$REPO_ROOT/config" python - <<'PY'
import os
try:
    import config
    print(getattr(config, 'variant_output_directory', os.getcwd()))
except Exception:
    print(os.getcwd())
PY
)

mkdir -p "$OUT_DIR"

# ----------------------------
# Build / cache gene BED from GENCODE
# ----------------------------
GENE_BED_CACHE="$REPO_ROOT/utils/gencode_v19.genes.ext${EXTEND}.sorted.bed"

if [[ ! -f "$GENE_BED_CACHE" ]]; then
  echo "Building gene BED cache from $GTF_PATH (extend ${EXTEND}bp)..."
  # Set environment variable EXT for perl ($ENV{EXT})
  EXT="$EXTEND" gzip -dc "$GTF_PATH" \
    | perl -ne '
        next if /^#/;
        @f = split("\t");
        next unless $f[2] eq "gene";
        my ($chr,$start,$end,$attr) = ($f[0], $f[3], $f[4], $f[8]);
        # NORMALISE: strip leading "chr" so chromosomes match GWAS (1..22,X,Y)
        $chr =~ s/^chr//;
        my ($gn) = ($attr =~ /gene_name "([^"]+)"/);
        $gn = defined $gn ? $gn : ".";
        my $ext = $ENV{EXT} // 0;
        my $s = $start - 1 - $ext;
        $s = 0 if $s < 0;
        my $e = $end + $ext;
        print join("\t", $chr, $s, $e, $gn) . "\n";
      ' \
    | sort -k1,1 -k2,2n > "$GENE_BED_CACHE"
  echo "Wrote gene bed: $GENE_BED_CACHE"
else
  echo "Using cached gene bed: $GENE_BED_CACHE"
fi

# ----------------------------
# Process each GWAS file
# ----------------------------
for INP in "${INPUTS[@]}"; do
  # Resolve file path relative to repo root if not absolute
  if [[ "$INP" = /* ]]; then
    FILEPATH="$INP"
  else
    FILEPATH="$REPO_ROOT/$INP"
  fi

  if [[ ! -f "$FILEPATH" ]]; then
    echo "Warning: GWAS file not found: $FILEPATH; skipping" >&2
    continue
  fi

  WORKDIR="$(pwd)/tmp_annotation_$(basename "$FILEPATH" .gz)_$$"
  mkdir -p "$WORKDIR"
  echo "Processing $FILEPATH in $WORKDIR"

  SNPS_BED="$WORKDIR/snps.bed"
  ORIG_MAP="$WORKDIR/orig_rows.tsv"
  HDR_FILE="$WORKDIR/header.txt"

  echo "Parsing GWAS input and extracting SNP positions (Python)..."
  export PY_FILEPATH="$FILEPATH"
  export PY_SNPS="$SNPS_BED"
  export PY_ORIG="$ORIG_MAP"
  export PY_HDR="$HDR_FILE"

  python - <<'PY'
import os, gzip

fp = os.environ['PY_FILEPATH']
snps = os.environ['PY_SNPS']
orig = os.environ['PY_ORIG']
hdr_out = os.environ['PY_HDR']

def find_idx(hdr, names):
    for n in names:
        if n in hdr:
            return hdr.index(n)
    return None

def split_line(line: str):
    """
    Split a line on tabs if present, otherwise on any whitespace.
    This makes the parser robust to both TSV and space-separated files.
    """
    line = line.rstrip('\n')
    if '\t' in line:
        return line.split('\t')
    else:
        return line.split()

with gzip.open(fp, 'rt') as fh, open(snps, 'w') as sf, open(orig, 'w') as of, open(hdr_out, 'w') as hf:
    header = fh.readline().rstrip('\n')
    hf.write(header + '\n')

    # robust header parsing
    if '\t' in header:
        cols = header.split('\t')
    else:
        cols = header.split()

    # more generous detection of chromosome / position columns
    chrom_idx  = find_idx(cols, ['CHROM', 'chr', 'Chromosome', 'CHR', 'chrom'])
    pos_idx    = find_idx(cols, ['GENPOS', 'pos', 'position', 'POS', 'BP', 'bp'])
    # marker / ID column (for 10:100000625:SNP style)
    marker_idx = find_idx(cols, ['MarkerName', 'SNP', 'ID'])

    n_snps = 0
    for lineno, line in enumerate(fh, start=2):
        parts = split_line(line)
        chrom = ''
        pos = ''

        # Use explicit CHROM / POS if present
        if chrom_idx is not None and pos_idx is not None \
           and chrom_idx < len(parts) and pos_idx < len(parts):
            if parts[chrom_idx] and parts[pos_idx]:
                chrom = parts[chrom_idx]
                pos = parts[pos_idx]

        # Fallback: parse from MarkerName / ID (e.g. "10:100000625:SNP")
        if (not chrom or not pos) and marker_idx is not None and marker_idx < len(parts):
            marker = parts[marker_idx]
            if ':' in marker:
                mparts = marker.split(':')
                if len(mparts) >= 2:
                    chrom = mparts[0]
                    pos   = mparts[1]

        if chrom and pos:
            chrom = str(chrom).lstrip('chr')
            try:
                p = int(float(pos))
            except Exception:
                continue
            sf.write(f"{chrom}\t{p-1}\t{p}\t{lineno}\n")
            of.write(f"{lineno}\t{line}")
            n_snps += 1

print(f"Extracted {n_snps} SNP positions from {fp}")
PY

  if [[ ! -s "$SNPS_BED" ]]; then
    echo "No SNP positions extracted from $FILEPATH; check file format" >&2
    rm -rf "$WORKDIR"
    continue
  fi

  SNPS_SORTED="$WORKDIR/snps.sorted.bed"
  sort -k1,1 -k2,2n "$SNPS_BED" > "$SNPS_SORTED"

  echo "Running bedtools closest (this is the fast step, but can be I/O heavy on large GWAS)..."
  CLOSEST_OUT="$WORKDIR/snps.closest.tsv"
  bedtools closest -a "$SNPS_SORTED" -b "$GENE_BED_CACHE" -D a -t first -sorted > "$CLOSEST_OUT"

  echo "Building gene map (line -> gene, distance)..."
  awk -F"\t" 'BEGIN{OFS="\t"} { line=$4; gene=$8; if(gene=="." || gene=="") gene="NoGene"; dist=$NF; print line,gene,dist }' \
    "$CLOSEST_OUT" \
    | sort -k1,1n > "$WORKDIR/gene_map_sorted.tsv"

  sort -k1,1n "$ORIG_MAP" > "$WORKDIR/orig_rows_sorted.tsv"

  # Derive phenotype from input filename (suffix after last underscore, without .txt/.gz)
  base=$(basename "$FILEPATH")
  base=${base%.gz}
  base=${base%.txt}
  phenotype=${base##*_}
  OUT_CSV="$OUT_DIR/Annotated_GWAS_${phenotype}.csv"
  VAR_CSV="$OUT_DIR/variant_data_${phenotype}.csv"

  echo "Joining original rows with gene map and writing outputs via Python..."
  export MERGE_WORK="$WORKDIR"
  export MERGE_INPUT="$FILEPATH"
  export MERGE_OUT="$OUT_CSV"
  export MERGE_VAR="$VAR_CSV"

  python - <<'PY'
import os, gzip

work       = os.environ['MERGE_WORK']
input_path = os.environ['MERGE_INPUT']
out_csv    = os.environ['MERGE_OUT']
var_csv    = os.environ['MERGE_VAR']

gene_map_file = os.path.join(work, 'gene_map_sorted.tsv')

# Load gene map: original file line number -> gene
gene_map = {}
with open(gene_map_file) as f:
    for ln in f:
        parts = ln.rstrip('\n').split('\t')
        try:
            idx = int(parts[0])
        except Exception:
            continue
        gene = parts[1]
        gene_map[idx] = gene

# Open input GWAS (gz or plain)
if input_path.endswith('.gz'):
    fh = gzip.open(input_path, 'rt')
else:
    fh = open(input_path, 'rt')

header = fh.readline().rstrip('\n')
cols = header.split('\t') if '\t' in header else header.split()

# Determine effect column index if present
effect_candidates = ['Effect', 'BETA', 'beta', 'Beta']
effect_idx = None
for i, c in enumerate(cols):
    if c in effect_candidates:
        effect_idx = i
        break

with open(out_csv, 'w') as outf, open(var_csv, 'w') as vf:
    # Main annotated GWAS output
    outf.write(','.join(cols) + ',Gene\n')
    # Variant-level output
    if effect_idx is not None:
        vf.write('Gene,Effect\n')

    line_no = 2  # header is line 1
    kept = 0
    for line in fh:
        line = line.rstrip('\n')
        gene = gene_map.get(line_no, 'NoGene')
        if gene != 'NoGene':
            fields = line.split('\t') if '\t' in line else line.split()
            outf.write(','.join(fields) + ',' + gene + '\n')
            if effect_idx is not None and effect_idx < len(fields):
                vf.write(gene + ',' + fields[effect_idx] + '\n')
            kept += 1
        line_no += 1

print(f'Wrote annotated {kept} SNPs to {out_csv}')
if effect_idx is None:
    print('No Effect column found; variant_data file will be empty (no Effect column)')
else:
    print(f'Wrote variant data to {var_csv}')

fh.close()
PY

  echo "Cleaning up temporary files: $WORKDIR"
  rm -rf "$WORKDIR"
  echo "Done for $FILEPATH"
done
