#!/usr/bin/env python
# -*- coding=utf-8 -*-
import sys
import gzip

def scanning_COSMIC(gene_list,COSMIC_file,ofile):
    """
    Reading through COSMIC file and  geting gene list related items out
    """
    if COSMIC_file.endswith('gz'):
        COSMIC_data = gzip.open(COSMIC_file,'r')
    else:
        COSMIC_data = open(COSMIC_file,'r')
    with COSMIC_data,open(ofile,'w') as odata:
        title = COSMIC_data.readline()
        odata.write(title)
        for line in COSMIC_data:
            genename = line.strip().split('\t')[0].split('_')[0]
            if genename in gene_list:
                odata.write(line)

def reading_in_gene_list(gene_list_file):
    """
    Reading in gene list and return a gene list
    """
    gene_list = []
    with open(gene_list_file,'r') as indata:
        for line in indata:
            gene_list.append(line.strip())
    return gene_list

def main():
    if len(sys.argv) != 4:
        print('Usage: python cosmic_gene_extracting.py gene_list COSMIC_file ofile')
        exit(0)
    gene_list_file = sys.argv[1]
    COSMIC_file = sys.argv[2]
    ofile = sys.argv[3]
    gene_list = reading_in_gene_list(gene_list_file)
    scanning_COSMIC(gene_list,COSMIC_file,ofile)

if __name__ == "__main__":
    main()