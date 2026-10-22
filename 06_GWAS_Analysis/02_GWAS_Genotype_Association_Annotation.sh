#!/usr/bin/env bash
set -euo pipefail

# Project folders and software
working_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
r_input_dir="$working_dir/r_input"
r_output_dir="$working_dir/r_output"
annotation_resources="$working_dir/annotation_resources"
phenotype_dir="$working_dir/bash_output_phenotypes"
vcf_dir="$working_dir/bash_output_vcf"
genotype_dir="$working_dir/bash_output_genotype"
incidence_dir="$working_dir/bash_output_incidence"
intensity_dir="$working_dir/bash_output_intensity"
annotation_dir="$working_dir/bash_output_annotation"
mkdir -p "$phenotype_dir" "$vcf_dir" "$genotype_dir" "$incidence_dir" "$intensity_dir" "$annotation_dir"

PLINK="${PLINK:-plink}"
BCFTOOLS="${BCFTOOLS:-bcftools}"
ANNOVAR="${ANNOVAR:-$annotation_resources/table_annovar.pl}"
ANNOVAR_DB="${ANNOVAR_DB:-$annotation_resources/annovar_db}"
GFF3="${GFF3:-$annotation_resources/Merged_Annotation_FSJgenome_July2024_FINAL.braker.interpro.gff3}"

# Convert phenotype CSVs to PLINK FID/IID/phenotype files
python3 - "$r_output_dir" "$phenotype_dir" <<'PY'
import csv
import pathlib
import sys

input_dir = pathlib.Path(sys.argv[1])
output_dir = pathlib.Path(sys.argv[2])
for name, column in (("incidence", "Incidence"), ("intensity", "logIntensity")):
    with (input_dir / f"mhc_{name}.csv").open(newline="") as source, \
         (output_dir / f"{name}_pheno.txt").open("w") as target:
        for row in csv.DictReader(source):
            sample = row["SampleID"].strip()
            value = row[column].strip()
            if sample and value and value.upper() not in {"NA", "NAN"}:
                target.write(f"{sample}\t{sample}\t{value}\n")
PY

# Extract shared individuals and summarize the VCF
"$BCFTOOLS" view -S "$r_output_dir/shared_sampleIDs.txt" -Oz \
  -o "$vcf_dir/MHC_snps.vcf.gz" "$r_input_dir/snps_afterLD_DP3_unrelated.vcf.gz"
"$BCFTOOLS" index "$vcf_dir/MHC_snps.vcf.gz"
"$BCFTOOLS" stats "$vcf_dir/MHC_snps.vcf.gz" > "$vcf_dir/MHC_snps.stats"

# Convert shared-sample VCF to PLINK binary format
"$PLINK" --vcf "$vcf_dir/MHC_snps.vcf.gz" --make-bed --allow-extra-chr \
  --out "$genotype_dir/genotype_base"

# Incidence: SNP filtering, PCA, and PC covariates
"$PLINK" --bfile "$genotype_dir/genotype_base" --geno 0.20 --mac 5 \
  --make-bed --allow-extra-chr --out "$incidence_dir/incidence_qc"
"$PLINK" --bfile "$incidence_dir/incidence_qc" --pca 10 --allow-extra-chr \
  --out "$incidence_dir/incidence_pca"
awk 'BEGIN{OFS="\t"; print "FID","IID","PC1","PC2"} {print $1,$2,$3,$4}' \
  "$incidence_dir/incidence_pca.eigenvec" > "$incidence_dir/incidence_pcs.txt"

# Incidence: additive logistic GWAS and FDR-significant SNPs
"$PLINK" --bfile "$incidence_dir/incidence_qc" \
  --pheno "$phenotype_dir/incidence_pheno.txt" --1 \
  --covar "$incidence_dir/incidence_pcs.txt" --covar-name PC1,PC2 \
  --logistic hide-covar --adjust --allow-no-sex --allow-extra-chr \
  --out "$incidence_dir/incidence_gwas"
awk 'NR==1 || ($9!="NA" && $9<0.05)' \
  "$incidence_dir/incidence_gwas.assoc.logistic.adjusted" > "$incidence_dir/incidence_FDR_significant.txt"

# Intensity: infected-bird subset, SNP filtering, PCA, and PC covariates
awk 'BEGIN{OFS="\t"} {print $1,$2}' "$phenotype_dir/intensity_pheno.txt" \
  > "$intensity_dir/intensity_keep.txt"
"$PLINK" --bfile "$genotype_dir/genotype_base" --keep "$intensity_dir/intensity_keep.txt" \
  --geno 0.20 --mac 5 --make-bed --allow-extra-chr --out "$intensity_dir/intensity_qc"
"$PLINK" --bfile "$intensity_dir/intensity_qc" --pca 10 --allow-extra-chr \
  --out "$intensity_dir/intensity_pca"
awk 'BEGIN{OFS="\t"; print "FID","IID","PC1","PC2"} {print $1,$2,$3,$4}' \
  "$intensity_dir/intensity_pca.eigenvec" > "$intensity_dir/intensity_pcs.txt"

# Intensity: additive linear GWAS and FDR-significant SNPs
"$PLINK" --bfile "$intensity_dir/intensity_qc" \
  --pheno "$phenotype_dir/intensity_pheno.txt" \
  --covar "$intensity_dir/intensity_pcs.txt" --covar-name PC1,PC2 \
  --linear hide-covar --adjust --allow-no-sex --allow-extra-chr \
  --out "$intensity_dir/intensity_gwas"
awk 'NR==1 || ($9!="NA" && $9<0.05)' \
  "$intensity_dir/intensity_gwas.assoc.linear.adjusted" > "$intensity_dir/intensity_FDR_significant.txt"

# Collect significant SNP IDs, beta coefficients, and adjusted p-values
awk 'NR>1 && $9!="NA" && $9<0.05 {print $2}' \
  "$intensity_dir/intensity_gwas.assoc.linear.adjusted" > "$annotation_dir/significant_intensity_snps.txt"
awk 'NR==FNR {wanted[$1]=1; next} FNR==1 || ($2 in wanted)' \
  "$annotation_dir/significant_intensity_snps.txt" "$intensity_dir/intensity_gwas.assoc.linear" \
  > "$annotation_dir/significant_intensity_effects.txt"
awk 'NR==FNR {wanted[$1]=1; next} FNR==1 || ($2 in wanted)' \
  "$annotation_dir/significant_intensity_snps.txt" "$intensity_dir/intensity_gwas.assoc.linear.adjusted" \
  > "$annotation_dir/significant_intensity_adjusted.txt"

# Prepare VCF, ANNOVAR input, and SNP ID mapping
"$BCFTOOLS" view -i "ID=@$annotation_dir/significant_intensity_snps.txt" \
  "$vcf_dir/MHC_snps.vcf.gz" -Ov -o "$annotation_dir/significant_intensity_snps.vcf"
awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$2,$2,$4,$5}' \
  "$annotation_dir/significant_intensity_snps.vcf" > "$annotation_dir/significant_intensity_snps.avinput"
awk 'BEGIN{OFS="\t"} !/^#/ {print $3,$1,$2,$4,$5}' \
  "$annotation_dir/significant_intensity_snps.vcf" > "$annotation_dir/significant_intensity_snps_IDmap.txt"

# Annotate intensity-associated SNPs and extract overlapping genes
if [[ -s "$annotation_dir/significant_intensity_snps.avinput" ]]; then
  perl "$ANNOVAR" "$annotation_dir/significant_intensity_snps.avinput" "$ANNOVAR_DB" \
    -buildver fsj -out "$annotation_dir/significant_intensity_snps_annotated" -remove \
    -protocol refGene -operation g -nastring .

  python3 - "$annotation_dir" <<'PY'
import csv
import pathlib
import re
import sys

root = pathlib.Path(sys.argv[1])
with (root / 'significant_intensity_snps_annotated.fsj_multianno.txt').open(newline='') as handle:
    genes = set()
    for row in csv.DictReader(handle, delimiter='\t'):
        genes.update(re.findall(r'g\d+', row.get('Gene.refGene', '') or ''))
with (root / 'significant_gene_ids.txt').open('w') as handle:
    for gene in sorted(genes):
        handle.write(gene + '\n')
PY

  if [[ -s "$annotation_dir/significant_gene_ids.txt" ]]; then
    grep -F -f "$annotation_dir/significant_gene_ids.txt" "$GFF3" \
      > "$annotation_dir/significant_genes_full.gff3" || true
    awk '$3=="gene"' "$annotation_dir/significant_genes_full.gff3" \
      > "$annotation_dir/significant_genes_only.gff3"
  else
    : > "$annotation_dir/significant_genes_full.gff3"
    : > "$annotation_dir/significant_genes_only.gff3"
  fi
else
  : > "$annotation_dir/significant_gene_ids.txt"
  : > "$annotation_dir/significant_genes_full.gff3"
  : > "$annotation_dir/significant_genes_only.gff3"
fi

# Export detailed and gene-level functional annotation tables
python3 - "$annotation_dir/significant_genes_full.gff3" "$annotation_dir" <<'PYTHON'
import csv
import pathlib
import re
import sys

source = pathlib.Path(sys.argv[1])
output_dir = pathlib.Path(sys.argv[2])
full_columns = ["Gene_ID", "Chromosome", "Feature", "Start", "End", "Strand",
                "Name", "Gene_name", "Product", "GO_terms", "InterPro_IDs",
                "UniProt", "Dbxref", "All_attributes"]
summary_columns = ["Gene_ID", "Chromosome", "Gene_name", "Product",
                   "GO_terms", "InterPro_IDs", "UniProt", "Dbxref"]

records = []
with source.open() as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        chrom, _, feature, start, end, _, strand, _, attributes = parts[:9]
        attrs = {}
        for field in attributes.split(";"):
            if "=" in field:
                key, value = field.split("=", 1)
                attrs[key] = value
        identifiers = ",".join((attrs.get("ID", ""), attrs.get("Parent", "")))
        found = re.search(r"(?<![A-Za-z0-9])g\d+", identifiers)
        if not found:
            continue
        gene = found.group()
        dbxref = attrs.get("Dbxref", "")
        row = {
            "Gene_ID": gene, "Chromosome": chrom, "Feature": feature,
            "Start": start, "End": end, "Strand": strand,
            "Name": attrs.get("Name", ""),
            "Gene_name": attrs.get("gene", attrs.get("gene_name", "")),
            "Product": attrs.get("product", ""),
            "GO_terms": ";".join(sorted(set(re.findall(r"GO:\d+", attrs.get("Ontology_term", ""))))),
            "InterPro_IDs": ";".join(sorted(set(re.findall(r"IPR\d+", dbxref)))),
            "UniProt": attrs.get("uniprot_id", ""),
            "Dbxref": dbxref, "All_attributes": attributes,
        }
        records.append(row)

with (output_dir / "significant_genes_functional_annotation_FULL.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=full_columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(records)

by_gene = {}
for row in records:
    gene = row["Gene_ID"]
    if gene not in by_gene:
        by_gene[gene] = {key: "" for key in summary_columns}
        by_gene[gene]["Gene_ID"] = gene
        by_gene[gene]["Chromosome"] = row["Chromosome"]
    summary = by_gene[gene]
    for field in summary_columns[2:]:
        if not summary[field] and row[field]:
            summary[field] = row[field]

with (output_dir / "significant_genes_functional_annotation_SUMMARY.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=summary_columns, delimiter="\t")
    writer.writeheader()
    writer.writerows(by_gene[gene] for gene in sorted(by_gene))
PYTHON
