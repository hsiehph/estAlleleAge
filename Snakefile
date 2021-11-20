import os, re
import pandas as pd
import numpy as np
import glob

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

shell.executable("/bin/bash")
shell.prefix("source %s/env.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)

if config == {}:
	configfile: "%s/config.yaml" % SNAKEMAKE_DIR

# parameters for input preparation
beagle=config["beagle"]

SUB_CHR_NAME = config["SUB_CHR_NAME"]
FILE_SAMPLES = config["FILE_SAMPLES"]
PATH2MASK = config["PATH2MASK"]
PATH2ANCESTRAL = config["PATH2ANCESTRAL"]

dict_chrom_refVCF = {}
for f in glob.glob(config["PATH2refVCF"] + "*.gz"):
	if f.endswith("vcf.gz"):
		chrom = re.search(".+?(chr[0-9]+).+?",f).group(1)
		print (chrom)
#		chrom = f.split(".")[1].split("_")[0]
		dict_chrom_refVCF[chrom] = f

dict_chrom_geneMap = {}
for f in glob.glob(config["PATH2geneMap"] + "*.map"):
	if f.endswith("map"):
		chrom = f.split(".")[1]
		dict_chrom_geneMap[chrom] = f

dict_chrom_RelateGeneMap = {}
for f in glob.glob(config["PATH2GENEMAP"] + "/*.txt"):
	if f.endswith("txt"):
		chrom = f.split("/")[-1].split("_")[7].split(".")[0]
		dict_chrom_RelateGeneMap[chrom] = f

# parameters for Relate 
PATH_RELATE = config["PATH_RELATE"]
PATH2RelateOut = config["PATH2RelateOut"]
Ne = config["Ne"]
mu = config["mu"]

if not os.path.exists("log"):
	os.makedirs("log")


dict_chrom_mask = {}
for f in glob.glob(config["PATH2MASK"] + "/*.fa"):
	if f.endswith("fa"):
		chrom = f.split("/")[-1].split("_")[-1].split(".")[0]
#		chrom = f.split("/")[-1].split("_")[0][3:]
		dict_chrom_mask[chrom] = f

dict_chrom_ancestral = {}
for f in glob.glob(config["PATH2ANCESTRAL"] + "/*.fa"):
	if f.endswith("fa"):
		chrom = f.split("/")[-1].split("_")[-1].split(".")[0]
#		chrom = f.split("/")[-1].split("_")[0][3:]
		dict_chrom_ancestral[chrom] = f


def _get_unphasedVCF(wildcards):
	return config["unphasedVCF"]["chr"+wildcards.CHROM][wildcards.NAME]

DATA = []
for chrom in config["unphasedVCF"].keys():
	for name in config["unphasedVCF"][chrom].keys():
		DATA.append(name)
	break

rangeBED = config["rangeBED"]

dict_chrom_regions = {}
for k in config["regionFiles"]:
	regionFile = config["regionFiles"][k]
	chrom = k[3:]
	with open(regionFile) as fin:
		for line in fin:
			line = line.strip().split()
			start = int(line[0]) - rangeBED
			if start < 0:
				start = 0
			end = int(line[0]) + rangeBED
			if chrom not in dict_chrom_regions:
				dict_chrom_regions[chrom] = []
			dict_chrom_regions[chrom].append("%s_%s_%s" % (k, start, end))

wildcard_constraints:
    REGION="chr[0-9]+_\d+_\d+",
    CHROM="[0-9]+",


def all_chroms_and_regions(wildcards):
	rtn = []
	pat = "RelateOutput/chr{CHROM}/{REGION}/{REGION}.relateOut.mut"
	for CHROM in dict_chrom_regions.keys():
		for REGION in dict_chrom_regions[CHROM]:
			rtn.append( pat.format(CHROM=CHROM, REGION=REGION) )
	return(rtn)

	

rule dummy:
	input: "finished.txt"



rule RMintermediateFiles:
	input: "RelateOutput/allTargets.relateOut.txt"
	output: "finished.txt"
	params: sge_opts="-l mfree=4G -l h_rt=48:00:00"
	shell:
		" rm -rf phasedVCF/; touch {output} ;"
		" rm -rf phasedVCFmerged/chr22/*/*.rename.simpleI*; touch {output}"

rule gatherRelateOut:
	input: 
		all_chroms_and_regions
	output: "RelateOutput/allTargets.relateOut.txt"
	params: sge_opts="-l mfree=4G -l h_rt=48:00:00"
	shell:
		" python ./scripts/pullTargetOut_fromMUTfile.py {PATH2RelateOut} > {output} "


rule Relate:
	input:
		hap = "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.relate.mask.haps",
		dist = "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.relate.mask.dist",
		sample= "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.relate.sample"
	output: "RelateOutput/chr{CHROM}/{REGION}/{REGION}.relateOut.anc", 
			"RelateOutput/chr{CHROM}/{REGION}/{REGION}.relateOut.mut", 
	params: sge_opts="-l mfree=20G -l h_rt=48:00:00"
	run:
		chrom, start, end = list({wildcards.REGION})[0].split("_")
		FILE_RelateGeneMap = dict_chrom_RelateGeneMap[chrom[3:]]
		shell(" cd RelateOutput/chr{wildcards.CHROM}/{wildcards.REGION} ; {PATH_RELATE}/Relate --mode All -m {mu} -N {Ne} --map ../../../%s --haps ../../../{input.hap} --sample ../../../{input.sample} --dist ../../../{input.dist}  --memory 20 -o {wildcards.REGION}.relateOut ; cd ../../../ ; touch {output[0]}" % FILE_RelateGeneMap)


rule relate_mask:
	input:
		hap = "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.relate.ancestral.haps",
		sample= "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.relate.sample"
	output:
		hap = "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.relate.mask.haps",
		dist = "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.relate.mask.dist"
	params:
		sge_opts="-l mfree=4G -l h_rt=120:00:00 ", prefix="{REGION}"
	run:
		chrom, start, end = list({wildcards.REGION})[0].split("_")
		FILE_mask = dict_chrom_mask[chrom[3:]]
		shell(" x=`basename {input.hap} .ancestral.haps `; cd phasedVCFmerged/chr{wildcards.CHROM}/{wildcards.REGION} ; {PATH_RELATE}/RelateFileFormats --mode FilterHapsUsingMask --haps ../../../{input.hap} --sample ../../../{input.sample} --mask ../../../%s -o ${{x}}.mask ;" % FILE_mask)


rule relate_ancestralFlip:
	input:
		hap = "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.relate.biallelic.haps",
		sample= "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.relate.sample"
	output:
		hap = "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.relate.ancestral.haps"
	params:
		sge_opts="-l mfree=4G -l h_rt=120:00:00 ", prefix="{REGION}"
	run:
		chrom, start, end = list({wildcards.REGION})[0].split("_")
		FILE_ancestral = dict_chrom_ancestral[chrom[3:]]
		shell(" x=`basename {input.hap} .biallelic.haps `; cd phasedVCFmerged/chr{wildcards.CHROM}/{wildcards.REGION} ; {PATH_RELATE}/RelateFileFormats --mode FlipHapsUsingAncestor --haps ../../../{input.hap} --sample ../../../{input.sample} --ancestor ../../../%s -o ${{x}}.ancestral ;" % FILE_ancestral)


rule relate_rmNonbiallelicSNV:
	input:
		hap= "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.relate.haps",
		sample= "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.relate.sample"
	output:
		hap = "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.relate.biallelic.haps"
	params:
		sge_opts="-l mfree=4G -l h_rt=120:00:00 ", prefix="{REGION}"
	shell:
		" x=`basename {input.hap} .haps `; cd phasedVCFmerged/chr{wildcards.CHROM}/{wildcards.REGION} ; {PATH_RELATE}/RelateFileFormats --mode RemoveNonBiallelicSNPs --haps ../../../{input.hap} -o ${{x}}.biallelic ;"


rule convertVCF2Relate:
	input:
		vcf= "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.phased.rename.simpleINDEL.vcf.gz",
		tbi= "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.phased.rename.simpleINDEL.vcf.gz.tbi",
#		reject= "rejLiftOver/{REGION}.phased.rename.rejLiftOver.vcf"
	output:
		hap= "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.relate.haps",
		sample= "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.relate.sample"
	params:
		sge_opts="-l mfree=4G -l h_rt=120:00:00", prefix="{REGION}"
	shell:
		" x=`basename {input.vcf} .phased.rename.simpleINDEL.vcf.gz `; cd phasedVCFmerged/chr{wildcards.CHROM}/{wildcards.REGION} ; {PATH_RELATE}/RelateFileFormats --mode ConvertFromVcf --haps ${{x}}.relate.haps --sample ${{x}}.relate.sample -i ../../../phasedVCFmerged/chr{wildcards.CHROM}/{wildcards.REGION}/${{x}}.phased.rename.simpleINDEL ;"


rule simplifyINDEL:
	input:
		vcf= "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.phased.rename.subset.vcf.gz",
		tbi= "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.phased.rename.subset.vcf.gz.tbi",
	output:
		vcf= "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.phased.rename.simpleINDEL.vcf.gz",
		tbi= "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.phased.rename.simpleINDEL.vcf.gz.tbi",
	params:
		sge_opts="-l mfree=4G -l h_rt=120:00:00", prefix="{REGION}"
	shell:
		" zcat {input.vcf} | python scripts/simplifyAllelesInVCF.py - | bgzip -c > {output.vcf}; "
		" tabix -p vcf {output.vcf} "


rule mergeANDpullSamples:
	input:
		vcf= expand("phasedVCF/chr{{CHROM}}/{{REGION}}/{NAME}.{{REGION}}.phased.rename.vcf.gz", NAME=DATA),
		tbi= expand("phasedVCF/chr{{CHROM}}/{{REGION}}/{NAME}.{{REGION}}.phased.rename.vcf.gz.tbi", NAME=DATA),
	output:
		vcf= "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.phased.rename.subset.vcf.gz",
		tbi= "phasedVCFmerged/chr{CHROM}/{REGION}/{REGION}.phased.rename.subset.vcf.gz.tbi"
	params:
		sge_opts="-l mfree=4G -l h_rt=120:00:00 -pe serial 3"
	run:
		if len(input.vcf) > 1:
			shell(" bcftools merge --threads 3 -Ov {input.vcf} | bcftools view -Ov --samples-file {FILE_SAMPLES} | bcftools +missing2ref - -- -p | bgzip -c > {output.vcf} ; tabix -p vcf {output.vcf} ")
		else:
			shell(" bcftools view -Ov --samples-file {FILE_SAMPLES} {input.vcf} | bcftools +missing2ref - -- -p | bgzip -c > {output.vcf} ; tabix -p vcf {output.vcf} ")

rule renameCHR:
	input: 
		vcf="phasedVCF/chr{CHROM}/{REGION}/{NAME}.{REGION}.phased.vcf.gz",
		tbi="phasedVCF/chr{CHROM}/{REGION}/{NAME}.{REGION}.phased.vcf.gz.tbi",
	output:
		vcf="phasedVCF/chr{CHROM}/{REGION}/{NAME}.{REGION}.phased.rename.vcf.gz",
		tbi="phasedVCF/chr{CHROM}/{REGION}/{NAME}.{REGION}.phased.rename.vcf.gz.tbi",
	params:
		sge_opts="-l mfree=4G -l h_rt=120:00:00 -pe serial 3"
	shell:
		" bcftools annotate --threads 3 --rename-chrs {SUB_CHR_NAME} {input.vcf}| bgzip -c > {output.vcf}; "
		" tabix -p vcf {output.vcf} ;"

rule phase:
	input:
		vcf = ancient(_get_unphasedVCF)
	output: 
		vcf="phasedVCF/chr{CHROM}/{REGION}/{NAME}.{REGION}.phased.vcf.gz",
		tbi="phasedVCF/chr{CHROM}/{REGION}/{NAME}.{REGION}.phased.vcf.gz.tbi",
	params:
		sge_opts="-l mfree=8G -l h_rt=200:00:00 -pe serial 3"
	run:
		chrom, start, end = wildcards.REGION.split("_")
		region = "%s:%s-%s" % (chrom, start, end)
		ref_VCF = dict_chrom_refVCF[chrom]
		geneMap = dict_chrom_geneMap[chrom]
		shell(""" java -Xmx8G -jar {beagle} gt={input.vcf} map={geneMap} nthreads=3 chrom={region} out=phasedVCF/chr{wildcards.CHROM}/{wildcards.REGION}/{wildcards.NAME}.{wildcards.REGION}.phased """)
		shell(""" tabix -p vcf {output.vcf} """)



