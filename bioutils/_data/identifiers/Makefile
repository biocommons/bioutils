# Download assembly info from NCBI and convert to json
# Use:
# 

.PHONY: FORCE
.SUFFIXES:
.DELETE_ON_ERROR:

SHELL:=/bin/bash -o pipefail
PATH:=../../../sbin:${PATH}

registry.json.gz: FORCE
	curl -s https://identifiers.org/service/registryxml | ./convert-registry | gzip >$@
