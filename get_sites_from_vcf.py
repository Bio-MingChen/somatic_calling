#!/usr/bin/env python
# -*- coding=utf-8 -*-
import sys
import gzip

def make_sites_dict(sites_file):
    """
    extracing chr pos ref alt from sites_file and generating sites_list
    """
    sites_list = []
    with open(sites_file,'r') as indata:
        for line in indata:
            line_list = line.strip().split('\t')
            sites_tup = tuple(line_list)
            sites_list.append(sites_tup)
    return sites_list

def extracting_sites_from_vcf(sites_list,vcf,ofile):
    """
    """
    if vcf.endswith('gz'):
        indata = gzip.open(vcf,'r')
    else:
        indata = open(vcf,'r')
    with indata,\
        open(ofile,'w') as odata:
        for line in indata:
            if not line.startswith('#'):
                line_list = line.strip().split('\t')
                sites_tup = (line_list[0],line_list[1],line_list[3],line_list[4])
                if sites_tup in sites_list:
                    odata.write(line)
            else:
                odata.write(line)

if __name__ == "__main__":
    vcf = sys.argv[1]
    sites_list = make_sites_dict(sys.argv[2])
    ofile = sys.argv[3]
    extracting_sites_from_vcf(sites_list,vcf,ofile)