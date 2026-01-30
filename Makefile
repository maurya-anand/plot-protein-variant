.PHONY: install

install:
	mkdir ./bin
	wget https://raw.githubusercontent.com/maurya-anand/gene-to-protein-domains/refs/heads/main/gene-to-protein-domains.py -O ./bin/gene-to-protein-domains.py
	chmod a+x ./bin/*.*