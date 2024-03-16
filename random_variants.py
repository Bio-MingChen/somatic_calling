#!/usr/bin/env python
# -*- coding=utf-8 -*-

#根据输入的配置文件来生成随机的snp位点,以及af值,以用于somatic的随机插入
import argparse
from random import randint,random

def generate_random_snv(bed_file,ofile,af_max=1,af_min=0.2,random_number=5):
    """
    reading in bed_file and generating random site with randome or indicated AF 
    """
    with open(bed_file,'r') as indata,\
    open(ofile,'w') as odata:
        for line in indata:
            line_list = line.strip().split("\t")
            chrom,start,end = line_list[:3]
            try:
                af = line_list[3]
            except:
                af = round(random(),1)
            random_sites = [randint(int(start),int(end)) for _ in range(random_number)]
            for site in random_sites:
                odata.write("\t".join([str(i) for i in [chrom,site,site,af]]) +'\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bed','-b',required=True,help='用于生成随机数的区间文件')
    parser.add_argument('--ofile','-o',default='random.txt',help='输出文件')
    parser.add_argument('--max_af','-max',default=1,type=int,help='随机生成AF时的上限')
    parser.add_argument('--min_af','-min',default=0.2,type=int,help='随机生成AF时的下限')
    parser.add_argument('--random_number','-num',default=5,type=int,help='每个给定区间的随机数生成数')
    args = vars(parser.parse_args())
    bed_file = args.get('bed')
    ofile = args.get('ofile')
    max_af = args.get('max_af')
    min_af = args.get('min_af')
    random_number = args.get('random_number')
    generate_random_snv(bed_file,ofile,max_af,min_af,random_number)