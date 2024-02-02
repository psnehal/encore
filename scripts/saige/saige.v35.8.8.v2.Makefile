# This make file will run SAIGE pipeline

# This script uses the reference fai index file (REFFILEINDEX) to
#   create chcuks of size BINSIZE for STEP 2
# Use THREADS= to set the number of threads to use for STEP 1
# Use "make -j" to specify the number of parallel jobs for STEP 2

# Make variables can be replaced on the command line. For example
#  make -f saige.Makefile -j20 THREADS=20 RESPONSE=BMI


OUTDIR = ./
OUTNAME = results.txt.gz
LOGNAME = saige.log

VCFFIELD = DS
PLINKFILE = plink.bed
SAVFILE = chr1.sav
SAMPLEFILE = samples.txt
REFFILE = /Users/snehalpatil/Documents/GithubProjects/ForkEncore/encore/anno/hs38DH.fa
REFFILEINDEX = $(addsuffix .fai, $(REFFILE))
PHENOFILE = pheno.txt
PHENOFILEIDCOL = IND_ID
RESPONSE = ""
RESPONSETYPE = quantitative
INVNORM = FALSE
COVAR = ""
BINSIZE = 1500000
THREADS = 8
STEP1OPT = --memoryChunk=2
STEP2OPT = --minMAF=0.001 --IsOutputAFinCaseCtrl=FALSE

#RLIBPATH = /net/encore1/saige/lib-v35.8.8
RSCRIPT = R_LIBS=/sw/pkgs/arc/encore/R-modules-SAIGE-0.39.4 /sw/pkgs/arc/stacks/gcc/10.3.0/Ropenblas/4.2.0/bin/Rscript
STEP1SCRIPT = $(APPDIR)step1_fitNULLGLMM_v35.R
STEP2SCRIPT = $(APPDIR)step2_SPATests_savvy_v35.R
APPDIR = /sw/pkgs/arc/encore/encore/scripts/saige/
TABIX = /sw/pkgs/arc/encore/tabix/1.9/bin

ifeq ("$(wildcard $(REFFILEINDEX))","")
$(error Cannot read REFFILEINDEX ($(REFFILEINDEX)) in order to create chunks for step 2)
endif

#MAKE BINS CHR.START.END
CHRS = $(shell seq 1 22) X
$(foreach chr, $(CHRS), $(eval CHR$(chr)BINS = $(shell awk -v bin=$(BINSIZE) '/^(chr)?$(chr)\y/ {s=1;e=bin; while(s<$$2) {if (e>$$2) {e=$$2}; print $$1 "." s "." e; s=s+bin; e=e+bin;}}' $(REFFILEINDEX))))
BINS = $(foreach chr, $(CHRS), $(CHR$(chr)BINS))

all: $(OUTDIR)$(OUTNAME) $(OUTDIR)$(OUTNAME).tbi $(OUTDIR)$(LOGNAME).gz

print-%  : ; @echo $* = $($*)

$(OUTDIR)step1.rda:
	$(RSCRIPT) $(STEP1SCRIPT) \
		--plinkFile=$(basename $(PLINKFILE)) \
		--phenoFile=$(PHENOFILE) \
		--phenoCol=$(RESPONSE) \
		--covarColList=$(COVAR) \
		--sampleIDColinphenoFile=$(PHENOFILEIDCOL) \
		--traitType=$(RESPONSETYPE) \
		--outputPrefix=$(basename $@) \
		--invNormalize=$(INVNORM) \
		--nThreads=$(THREADS) \
		$(STEP1OPT) &> $(OUTDIR)$(LOGNAME)

freezeid := 2c7adba4-b4ac-11ee-9161-000c29de3a12


$(OUTDIR)step2.bin.%.txt: $(OUTDIR)step1.rda
	$(RSCRIPT) $(STEP2SCRIPT) \
		--savFile=$(SAVFILE) \
		--sepchr \
		--chrom=$(subst chr,,$(word 1,$(subst ., , $*)))\
		--start=$(word 2, $(subst ., ,$*)) \
		--end=$(word 3, $(subst ., ,$*)) \
		--sampleFile=$(SAMPLEFILE) \
		--GMMATmodelFile=$^ \
		--varianceRatioFile=$(addsuffix .varianceRatio.txt, $(basename $^)) \
		--SAIGEOutputFile=$@ \
		--numLinesOutput=2 \
		$(STEP2OPT) &> $@.log

$(OUTDIR)step2.bin.%.txt: $(OUTDIR)step1.rda
	$(RSCRIPT) $(STEP2SCRIPT) \
		--savFile=$(SAVFILE) \
		--vcfField=$(VCFFIELD)\
		--sepchr \
		--chrom=$(word 1, $(subst ., ,$*)) \
		--start=$(word 2, $(subst ., ,$*)) \
		--end=$(word 3, $(subst ., ,$*)) \
		--sampleFile=$(SAMPLEFILE) \
		--GMMATmodelFile=$^ \
		--varianceRatioFile=$(addsuffix .varianceRatio.txt, $(basename $^)) \
		--SAIGEOutputFile=$@ \
		--numLinesOutput=2 \
		$(STEP2OPT) &> $@.log



STEP2FILES = $(foreach bin,$(BINS),$(OUTDIR)step2.bin.$(bin).txt)
STEP2LOGS = $(foreach bin,$(STEP2FILES),$(bin).log)
$(OUTDIR)$(OUTNAME): $(STEP2FILES)
	awk 'FNR!=1 || NR==1' $^ | tr " " "\t" | bgzip -c > $@

$(OUTDIR)$(OUTNAME).tbi: $(OUTDIR)$(OUTNAME)
	tabix -s1 -b2 -e2 -S1 $^

$(OUTDIR)$(LOGNAME).gz: $(OUTDIR)$(OUTNAME).tbi
	cat $(OUTDIR)$(LOGNAME) $(STEP2LOGS) | gzip -c > $@

.PHONY: clean

clean:
	rm -f $(STEP2FILES)
	rm -f $(STEP2LOGS)
	rm -f $(OUTDIR)$(LOGNAME)
