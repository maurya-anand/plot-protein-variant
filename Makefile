.PHONY: install build_docker pull_docker

install:
	mkdir ./bin
	wget https://raw.githubusercontent.com/maurya-anand/gene-to-protein-domains/refs/heads/main/gene-to-protein-domains.py -O ./bin/gene-to-protein-domains.py
	chmod a+x ./bin/*.*

CONTAINER := $(shell sed -n '/withName: "PLOT_VARIANTS"/,/\}/p' conf/base.config | sed -n "s/.*container = '\([^']*\)'.*/\1/p" | head -n1)
SIF_NAME := $(shell echo $(CONTAINER) | sed -e 's|/|-|g' -e 's|:|-|').img
SIF_PATH := ./container/$(SIF_NAME)

build_docker:
	docker build -t $(CONTAINER) -f container/Dockerfile .

pull_docker:
	docker pull $(CONTAINER)

./container:
	mkdir -p ./container

# ensure pull_singularity is a normal target depending on the image file
pull_singularity: $(SIF_PATH)

$(SIF_PATH): ./container
	if [ -f "$@" ]; then \
		echo "$@ already exists, skipping"; \
	else \
		singularity pull --name "$@" docker://$(CONTAINER); \
	fi