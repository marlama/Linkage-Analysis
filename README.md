# Linkage-Disequilibrium
Pipeline for quality Control and auxiliary scripts to run Linkage analysis in families using Merlin software (http://csg.sph.umich.edu/abecasis/Merlin/tour/input_files.html). 

# Control Quality

## Filter just SNPs common with  Reference data
### For Positions
awk '{print "chr"$1"\t"$4}' FileReference.bim > File_UpdateID_POS.txt
 
for chr in (1..22) ; do  vcftools --gzvcf TargetFile_chr${chr}.vcf.gz --positions File_UpdateID_POS.txt --recode --stdout | gzip -c > TargetFile_SNPsReference_chr${chr}.vcf.gz ; done

## Merged chr and filter Biallelic Only

bcftools concat TargetFile_SNPsReference_chr${chr}.vcf.gz  | bcftools view --max-alleles 2 -Oz -o TargetFile_ReferenceSNPs_BiallelicOnly.vcf.gz

## Transform to bfile
plink2 --vcf TargetFile_ReferenceSNPs_BiallelicOnly.vcf.gz --make-bed --out TargetFile_ReferenceSNPs_BiallelicOnly

## Run SmartQC
https://github.com/ldgh/Smart-Cleaning
- Remove chr 0
- Remove duplicate data
- Remove missing data
- Infer individual sex
- Remove A|T and C|G variants
- Remove 100% heterozigotes variants
- dbSNPname
- liftOver

## Maybe will need add Case/Control information and update rs IDs before merge

## QC specific
### Filter Family Individuals

for fam in 1 2 3 ; do plink --bfile Data_Autossomic_SmartQC_ReferenceSNPs_BiallelicOnly --keep Family${fam}.txt --make-bed --out Family${fam} ; done

### Remove monomorphic markers
- freqNonMonomorphic = 1/((N*2)*10) #N= Number of family members, 10x to have sure that just remove monomorphic
- freqNonMonomorphic =1/18 = 0.05

plink --bfile Family1 --maf 0.05 --make-bed --out Family1_NoMonomorphic

### Remove Markers Within Regions That are Known for Long LD
https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)

plink --bfile Family1_NoMonomorphic --exclude SNPs2Remove_LD.txt --make-bed --out Family1_NoMonomorphic_NoLDRegions

### Make a list of SNPs in LD to be removed for the reference pop
for pop in African EastAsia European NativeAmerican SouthAsia ; do plink --bfile Reference --indep-pairwise 100 10 0.1 --out ${pop}_Reference_LD0.1 ; done

### Prun-in using the list for the correspondent ancestry and missing data
for fam in 1 2 3 ; do plink --bfile Family${fam}_NoMonomorphic_NoLDRegions --extract European_Reference_LD0.1.prune.in --geno 0.01 --mind 0.01 --make-bed --out Family${fam}_NoMonomorphic_NoLDRegions_LD0.1_GenoMind_Eur ; done

### Add Centimorgan Positions

for chr in (1..22) ; do awk  -v chr="$chr" '{if($1=="chr"chr){print $2"\t"$3"\t"$4}}' /MapaGenetico/recomb-hg38/genetic_map_GRCh38_merged.tab > /MapaGenetico/genetic_map_GRCh38_chr${chr}_combined.txt ; done

for chr in (1..22) ; do cat /MapaGenetico/Header.txt /MapaGenetico/genetic_map_GRCh38_chr${chr}_combined.txt > /MapaGenetico/genetic_map_GRCh38_chr${chr}_combined_Header.txt ; done

for fam in 1 2 3 ; do for chr in (1..22) ; do plink --bfile Family${fam}_NoMonomorphic_NoLDRegions_LD0.1_GenoMind_Eur --chr ${chr} --make-bed --out Family${fam}_NoMonomorphic_NoLDRegions_LD0.1_GenoMind_Eur_chr${chr} ; done ; done

for fam in 1 2 3 ; do for chr in (1..22) ; do plink --bfile Family${fam}_NoMonomorphic_NoLDRegions_LD0.1_GenoMind_Eur_chr${chr} --cm-map /MapaGenetico/genetic_map_GRCh38_chr${chr}_combined_Header.txt ${chr} --make-bed --out Family${fam}_NoMonomorphic_NoLDRegions_LD0.1_GenoMind_Eur_chr${chr}_map ; done ; done

# Linkage Analysis
