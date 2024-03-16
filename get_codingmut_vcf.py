#!/usr/bin/env python
# -*- coding=utf-8 -*-

#从COSMIC提供的Coding区的突变的vcf文件CosmicCodingMuts.vcf
#中提取感兴趣的相关基因的位点
import sys


def get_gene_list(gene_list_file):
    """
    reading genes list
    """
    gene_list = []
    with open(gene_list_file,'r') as indata:
        for line in indata:
            gene_list.append(line.strip())
    return gene_list

def vcf_parser(vcf,gene_list,ofile):
    """
    parsing vcf
    """
    with open(vcf,'r') as indata,\
        open(ofile,'w') as odata:
        for line in indata:
            if not line.startswith('#'):
                line_list = line.strip().split('\t')
                # chrom = line_list[0]
                # pos = line_list[1]
                # ref = line_list[3]
                # alt = line_list[4]
                info_dict = info_parser(line_list[7])
                if info_dict['GENE'] in gene_list:
                    odata.write(line)
            else:
                odata.write(line)
                
def info_parser(info):
    """
    parsing info field and transforming to a dict
    """
    info_dict = {}
    for field in info.strip().split(';'):
        #print(field.split('='))
        if '=' in field:
            key,value = field.split('=',1)
        else:
            key = field
            value = field
        info_dict[key] = value
    return info_dict

if __name__ == "__main__":
    vcf = sys.argv[1]
    gene_list = get_gene_list(sys.argv[2])
    ofile = sys.argv[3]
    vcf_parser(vcf,gene_list,ofile)