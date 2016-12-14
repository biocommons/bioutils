# Download assembly info from NCBI and convert to json
# Use:
# 

.PHONY: FORCE
.SUFFIXES:
.DELETE_ON_ERROR:

SHELL:=/bin/bash -o pipefail
PATH:=../../../sbin:${PATH}

update:
	make pull
	make json

# rsyncs all assembly records into pull/ directory
pull: FORCE
	mkdir -p $@
	rsync -L -v --no-motd ftp.ncbi.nlm.nih.gov::genomes/ASSEMBLY_REPORTS/All/{GCF_000001405,GCF_000306695}.*.assembly.txt $@

pull/%.json: pull/%.txt
	assembly-to-json -p ${@D} $<

json:
	for f in pull/*.txt; do assembly-to-json -p. "$$f"; done
