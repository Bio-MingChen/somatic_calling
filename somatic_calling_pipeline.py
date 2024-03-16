#!/usr/bin/env python
# -*- coding=utf-8 -*-

#Short Somatic Variants Calls Pipeline
import os
import configparser
import argparse
import textwrap
import logging
from pprint import pprint
from glob import glob

def set_logging(name):
    """
    Set basic logging config and name a logger
    """
    logging.basicConfig(filename="Somatic_Calling.log",level=logging.DEBUG,
        filemode="w",format='[%(asctime)s:%(funcName)s:%(name)s:%(levelname)s] %(message)s'
        )

    logger = logging.getLogger(name)
    return logger

logger = set_logging('Somatic_Calling')

script_path = os.path.dirname(os.path.abspath(__file__))
default_db_path = os.path.join(os.path.dirname(script_path),'db')
config_name = os.path.join(script_path,'somatic_config.ini')
#Step1 Generating Config File

def generating_config(args):
    """
    Making config for generating Short Somatic Variants Calling Script
    """
    db_path = args.get('database_path') or default_db_path
    hg37_ref = '/PUBLIC/database/HUMAN/genome/Human/human_g1k_v37_decoy.fasta'
    if 'NJPROJ' in os.environ['HOME']:
        hg37_ref = '/NJPROJ2/DISEASE/Database/genome/Human/human_g1k_v37_decoy.fasta'
    gnomad_common_variants = os.path.join(db_path,'af-only-gnomad.raw.sites.hg19.rm_chr.vcf.gz')
    gnomad_biallelic_snp_variants = os.path.join(db_path,'af-only-gnomad.raw.sites.hg19.rm_chr.biallele.vcf.gz')
    official_PON = os.path.join(db_path,'somatic-b37_Mutect2-exome-panel.vcf')
    hg37_ref_img = os.path.join(db_path,'hg37_reference.fasta.img')
    intervals_dir = os.path.join(db_path,'interval-files')
    hg37_dict = os.path.join(db_path,'somatic-b37_Homo_sapiens_assembly19.dict')
    extract_py = os.path.join(script_path,'GetAnnoFromMerge_pipe4.6.py')
    config = configparser.ConfigParser()
    config['db'] = {
        'gnomad_common_variants':gnomad_common_variants,
        'gnomad_biallelic_snp_variants':gnomad_biallelic_snp_variants,
        'official_PON':official_PON,
        'hg37_ref_img':hg37_ref_img,
        'intervals_dir':intervals_dir,
        'hg37_dict':hg37_dict,
        'hg37_ref':hg37_ref,
        'extract_py':extract_py,
        }
    with open(config_name,'w') as odata:
        config.write(odata)
    logger.info('Config file has been generated to {path}'.format(
        path = config_name
    ))

#Step2 Generating All Scripts Needed for Calling Short Somatic variants
class TitleParser():
    """
    Reading in a file title and than getting indicated column element 
    Case Ignore
    """
    
    def __init__(self,title):
        self.title_list = [i.lower() for i in title.strip().split('\t')]

    def get_field(self,line_list,colname):
        if len(self.title_list) != len(line_list):
            raise Exception("Title length differs with line!")
        try:
            idx = self.title_list.index(str(colname).lower())
        except:
            print('{colname} not in title!'.format(colname=colname))
            return None

        return line_list[idx]

def parse_sample_info(sample_info):
    """
    Parsing the sample info and get the meta information
    samplename
    sample_info_dict example
    Paired samples:
    {
        Tag: {
            T:sample_tumor,
            N:sample_normal,
            Sex: 'M'
        },
         
    }
    Single tumor samples:
    {
        sample:{
            T:sample_tumor,
            Sex: 'F'
        },
        
    }
    """
    sample_info_dict = {}
    with open(sample_info,'r') as indata:
        for line in indata:
            if line.lower().startswith('#familyid'):
                title = line
                title_parse = TitleParser(title)
            elif line.startswith('#'): # Jumping annotation information
                pass
            elif 'type' not in title.lower():
                line_list = line.strip().split('\t')
                sample_name = title_parse.get_field(line_list,'sampleid')
                gender = title_parse.get_field(line_list,'sex')
                sample_type =  'T' # Single Tumor Sample was supplied
                if sample_name not in sample_info_dict:
                    sample_info_dict[sample_name] = {sample_type:sample_name,'Sex':gender}
            else:
                # print('Paired Tumor parse is under developed!')
                # exit(0)
                line_list = line.strip().split('\t')
                sample_name = title_parse.get_field(line_list,'sampleid')
                gender = title_parse.get_field(line_list,'sex')
                sample_type = title_parse.get_field(line_list,'Type')
                tag = title_parse.get_field(line_list,'Tag')
                if tag not in sample_info_dict:
                    sample_info_dict[tag] = {sample_type:sample_name,'Sex':gender}
                else:
                    sample_info_dict[tag][sample_type] = sample_name
        pprint(sample_info_dict)
        return sample_info_dict

def sort_gatk_bam(key):
    """
    Sort the bam files with chromosome order by their name 
    name like sample.8.recal.bam
    """
    chroms = [str(i+1) for i in range(22)] + ['X','Y']
    basename = os.path.basename(key)
    return chroms.index(basename.split('.')[1])

def generate_bam_file(is_gatk,tumor,main_path,gender,):
    """
    gather all Bams used to call somatic mutation for both tumor and normal samples
    """
    if is_gatk:
        tumor_bam_file = []
        bam_dir = os.path.join(main_path,'Mutation',
            '{tumor}.gatk'.format(tumor=tumor),'bychr','*.recal.bam')
        bams = glob(bam_dir)
        if (gender == 'M' and len(bams) == 24) or (gender == 'F' and len(bams) == 23):
            for bam_path in bams:
                if os.path.exists(bam_path):
                    tumor_bam_file.append(bam_path)
                else:
                    raise IOError('Tumor bam path {bam_path} does not exists!'.format(bam_path=bam_path))
        else:
            raise Exception('bam number is {bam_number} when gender is {gender}'.format(
            bam_number = len(bams), gender = gender,
        ))
        tumor_bam_file.sort(key=sort_gatk_bam)
    else:
        tumor_bam_file = os.path.join(main_path, 'Mapping', '{tumor}.{tumor}'.format(tumor=tumor),
                tumor + '.final.bam')
        if not os.path.exists(tumor_bam_file):
            logger.error('Tumor bam path {tumor_bam_file} is NOT exist'.format(tumor_bam_file=tumor_bam_file))
            raise IOError('Tumor bam path {tumor_bam_file} is NOT exist'.format(tumor_bam_file=tumor_bam_file))
        tumor_bam_file = [tumor_bam_file]
    return tumor_bam_file

def generating_scripts(args):
    """
    Making Directories and Generating Scripts 
    """
    logger.info('Initializing basic config...')
    generating_config(args)
    config = configparser.ConfigParser()
    logger.info('Reading in config file form {path}...'.format(path=config_name))
    config.read(config_name)
    sample_info_dict = parse_sample_info(args.get('sample_info'))
    samples = [tag for tag in sample_info_dict]
    main_path = args.get('main_path')
    somatic_dir = os.path.join(main_path,'Somatic')
    hg37_ref = config.get('db','hg37_ref')
    extract_py = config.get('db','extract_py')
    if not os.path.exists(somatic_dir):
        os.mkdir(somatic_dir)
    config_job = []
    filter_script_ofile_list = []
    for tag in sample_info_dict:
        tag_dir = os.path.join(somatic_dir,tag)
        if not os.path.exists(tag_dir):
            os.mkdir(tag_dir)
        tumor = sample_info_dict[tag].get('T')
        normal = sample_info_dict[tag].get('N')
        gender = sample_info_dict[tag].get('Sex')
        is_gatk = args.get('v4.7gatk_mutation')
        tumor_bam_file = generate_bam_file(is_gatk,tumor,main_path,gender)
        other_info = {
            'directory':tag_dir,
            'tag': tag,
            'tumor':tumor,
            'gender':gender,
            'tumor_bam':",".join(tumor_bam_file),
        }       
        config['sample_info'] = other_info
        if normal:
            # normal_bam_file = os.path.join(main_path, 'Mapping', '%s.%s' % (normal,normal), normal + '.final.bam')
            normal_bam_file = generate_bam_file(is_gatk,normal,main_path,gender)
            if not os.path.exists(normal_bam_file):
                logger.warning('Single Tumor without matched normal sample mode')
            else:
                config['sample_info']['normal'] = normal
                config['sample_info']['normal_bam'] = ",".join(normal_bam_file)
        hg37_ref = config.get('db','hg37_ref')
        ofile = os.path.join(tag_dir,'{tag}.somatic_config.ini'.format(tag=tag))
        with open(ofile,'w') as odata:
            config.write(odata)
            logger.info('Generating config file to {ofile}'.format(ofile=ofile))
        logger.info('Generating Somatic calling scripts...')
        sub_config_job,filter_script_ofile = generate_Mutect2_script(config)
        config_job += sub_config_job
        filter_script_ofile_list.append(filter_script_ofile)
        
    config = merge_and_annovar(samples,main_path,hg37_ref,extract_py,config_job,filter_script_ofile_list)
    logger.info('Generating merge and annovar scripts ...')
    make_job(config,main_path)

def bam_map_to_dict(bams,chroms):
    """
    mapping bams to every chroms
    """
    if len(bams) == 24: # gatk split bams
        return {c:b for b,c in zip(bams,chroms)}
    elif len(bams) == 1: # single final bam
        return {c:bams[0] for c in chroms}
    else:
        logger.error('Wrong bams!{bams}'.format(bams=bams))
        print('Wrong bams!{bams}'.format(bams=bams))
        exit(0)

def generate_Mutect2_script(config):
    """
    Generating Mutect script
    """
    tag = config.get('sample_info','tag')
    gender = config.get('sample_info','gender')
    if gender == 'M':
        chroms = [str(i+1) for i in range(22)] + ['X','Y']
    elif gender == 'F':
        chroms = [str(i+1) for i in range(22)] + ['X']
    elif gender == 'U':
        logger.warning("{tag}'s gender is unknown! Y chromosome will be called")
        chroms = [str(i+1) for i in range(22)] + ['X','Y']

    hg37_ref = config.get('db','hg37_ref')
    gnomad_common_variants = config.get('db','gnomad_common_variants')
    official_PON = config.get('db','official_PON')
    directory = config.get('sample_info','directory')
    vcf_ofile_list = []
    bam_list = []
    f1r2_list = []
    ofile_list = []
    config_job = []
    script = textwrap.dedent("""
        set -eo pipefail
        cd {directory}
        echo start `date +"%F %T"`
        gatk --java-options -Xmx5000m Mutect2  \\
            -R {hg37_ref} \\""".format(hg37_ref=hg37_ref,directory=directory))
    tumor_bam = config.get('sample_info','tumor_bam').split(',')
    tumor_bam_dict = bam_map_to_dict(tumor_bam,chroms)
    tumor = config.get('sample_info','tumor')
    normal = config.get('sample_info','normal',fallback=None)
    normal_bam = config.get('sample_info', 'normal_bam',fallback=None)
    if normal:
        normal_bam = config.get('sample_info', 'normal_bam').split(',')
        normal_bam_dict = bam_map_to_dict(normal_bam,chroms)
    for chrom in chroms:
        mutect2_script = script + '\n    -I {tumor_bam} \\\n    -tumor {tumor} \\'.format(
            tumor_bam=tumor_bam_dict[chrom],tumor=tumor)
        if normal_bam:
                mutect2_script = script + '\n    -I {normal_bam}\\\n    -normal {normal} \\'.format(
                    normal_bam=normal_bam_dict[chrom],normal=normal)
        mutect2_script += textwrap.dedent("""
            --independent-mates \\
            -O {tag}.{chrom}.unfiltered.vcf \\
            --bam-output {tag}.{chrom}.bamout.bam \\
            --f1r2-tar-gz {tag}.{chrom}.f1r2.tar.gz \\
            -L {chrom} \\
            --germline-resource {gnomad_common_variants} \\
            -pon {official_PON} 
        echo end `date +"%F %T"`       
        """.format(**locals()))
        ofile = os.path.join(directory,'{tag}.{chrom}.Mutect2.sh'.format(tag=tag,chrom=chrom))
        ofile_list.append(ofile)
        vcf_ofile_list.append('{tag}.{chrom}.unfiltered.vcf'.format(tag=tag,chrom=chrom))
        bam_list.append('{tag}.{chrom}.bamout.bam'.format(tag=tag,chrom=chrom))
        f1r2_list.append('{tag}.{chrom}.f1r2.tar.gz'.format(tag=tag,chrom=chrom))
        with open(ofile,'w') as odata:
            odata.write(mutect2_script)
            logger.info('Mutect2 script {script} generated!'.format(script=ofile))
            config_job.append([ofile,'5G'])
    #merge vcf
    merge_vcf_script = textwrap.dedent(
        """
        set -eo pipefail
        cd {directory}
        echo start `date +"%F %T"`
        ## merging vcf
        gatk --java-options "-Xmx1g" MergeVcfs  \\""".format(directory=directory)
    )
    for ofile in vcf_ofile_list:
        merge_vcf_script += '\n -I %s \\' % (ofile)
    merge_vcf_script += textwrap.dedent('''
            -O {tag}.unfiltered.vcf
        
        ## merging stats
        gatk --java-options "-Xmx1g" MergeMutectStats \\'''.format(tag=tag))
    for ofile in vcf_ofile_list:
        merge_vcf_script += '\n    --stats %s.stats \\' % (ofile)
    merge_vcf_script += textwrap.dedent('''
            -O {tag}.unfiltered.merged.stats 
        echo end `date +"%F %T"`
        '''.format(tag=tag))
    merge_ofile = os.path.join(directory,'{tag}.vcf_merge.sh'.format(tag=tag))
    with open(merge_ofile,'w') as odata:
        odata.write(merge_vcf_script)
        logger.info('merge script {script} generated!'.format(script=merge_ofile))
        rely_scripts = ",".join([(i) for i in ofile_list])
        config_job.append([merge_ofile,'1G',rely_scripts])

    # Gather bam files
    gather_bam_script = textwrap.dedent(
        """
        set -eo pipefail
        cd {directory}
        echo start `date +"%F %T"`
        ## merge bams
        gatk --java-options "-Xmx5g" GatherBamFiles \\
            -R {hg37_ref} \\""".format(hg37_ref=hg37_ref,directory=directory)
    )
    for bam in bam_list:
        gather_bam_script += '\n    -I {bam} \\'.format(bam=bam)
    gather_bam_script += textwrap.dedent(
        """
            -O {tag}.unsorted.out.bam

        ## sorting bam
        gatk --java-options "-Xmx5g" SortSam -I {tag}.unsorted.out.bam \\
            -O {tag}.sorted.out.bam \\
            --SORT_ORDER coordinate \\
            -VALIDATION_STRINGENCY LENIENT

        ## making index
        gatk --java-options "-Xmx500m" BuildBamIndex  \\
            -I {tag}.sorted.out.bam -VALIDATION_STRINGENCY LENIENT
        echo end `date +"%F %T"`
        """.format(tag=tag)
    )
    gather_bam_script_ofile = os.path.join(directory,'{tag}.merge_and_sort_bam.sh'.format(tag=tag))
    with open(gather_bam_script_ofile,'w') as odata:
        odata.write(gather_bam_script)
        logger.info('Gathering bam and sorting script generated!')
        config_job.append([gather_bam_script_ofile,'5G',rely_scripts])
    # Get Pileup Summaries
    gnomad_biallelic_snp_variants = config.get('db','gnomad_biallelic_snp_variants')
    #intervals_dir = config.get('db','intervals_dir')
    tumor_pileup_ofile = generating_get_pileup_summary_script(gnomad_biallelic_snp_variants,
        chroms,tumor_bam_dict,'tumor',directory,config_job,hg37_ref,tag)
    if normal_bam:
        normal_pileup_ofile = generating_get_pileup_summary_script(gnomad_biallelic_snp_variants,
        chroms,normal_bam_dict,'normal',directory,config_job,hg37_ref,tag)
    
    # LearnReadOrientationModel
    model_script = textwrap.dedent('''
        set -eo pipefail
        cd {directory}
        echo start `date +"%F %T"`
        gatk --java-options "-Xmx1g" LearnReadOrientationModel \\'''.format(directory=directory))
    for f1r2 in f1r2_list:
        model_script += '\n    -I {f1r2} \\'.format(f1r2=f1r2)
    model_script += textwrap.dedent('''
            -O '{tag}.artifact-priors.tar.gz'
        echo end `date +"%F %T"`
        '''.format(tag=tag))
    model_script_ofile = os.path.join(directory,'{tag}.LearnReadOrientationModel.sh'.format(tag=tag))
    with open(model_script_ofile,'w') as odata:
        odata.write(model_script)
        logger.info('Learn Read Orientation Model Script {script} generated!'.format(script=model_script_ofile))
    config_job.append([model_script_ofile,'1G',(merge_ofile)])
    
    ## Gather Pileup Summaries
    hg37_dict = config.get('db','hg37_dict')
    type_lists = ['tumor','normal'] if normal_bam else ['tumor']
    gather_ofile_rely = []
    for sample_type in type_lists:
        if sample_type == 'tumor':
            pileup_rely = tumor_pileup_ofile
        else:
            pileup_rely = normal_pileup_ofile
        gather_summaries_script = textwrap.dedent(
            """
            set -eo pipefail
            cd {directory}
            echo start `date +"%F %T"`
            gatk --java-options "-Xmx1g" GatherPileupSummaries \\
                --sequence-dictionary {hg37_dict} \\""".format(hg37_dict=hg37_dict,directory=directory))
        # for idx in range(len(tumor_pileup_ofile.split(','))):
        for chrom in chroms:
            gather_summaries_script += '\n  -I {tag}.{sample_type}.{chrom}.pileups.table \\'.format(
                tag=tag,sample_type=sample_type,chrom=chrom)
        gather_summaries_script += textwrap.dedent(
            """
              -O {tag}.{sample_type}.pileups.table 
            echo end `date +"%F %T"`
            """.format(tag=tag,sample_type=sample_type))
        gather_ofile = os.path.join(directory,'{tag}.{sample_type}.gather_pileup.sh'.format(tag=tag,sample_type=sample_type))
        with open(gather_ofile,'w') as odata:
            odata.write(gather_summaries_script)
        logger.info('Gather pileup summaries script {} generated'.format(gather_ofile))
        config_job.append([gather_ofile,'1G',pileup_rely])
        gather_ofile_rely.append(gather_ofile)

    ## Calculate Contamination
    calculate_contamination_script = textwrap.dedent(
        """
        set -eo pipefail
        cd {directory}
        echo start `date +"%F %T"`
        gatk --java-options "-Xmx1g" CalculateContamination \\
            -I {tag}.tumor.pileups.table \\""".format(tag=tag,directory=directory))

    if normal_bam:
        calculate_contamination_script += '\n    --matched normal.pileups.table \\'
    calculate_contamination_script += textwrap.dedent(
        """
            -O {tag}.contamination.table
        echo end `date +"%F %T"`
        """.format(tag=tag))
    calculate_contamination_script_ofile = os.path.join(directory,'{tag}.CalculateContamination.sh'.format(tag=tag))
    with open(calculate_contamination_script_ofile,'w') as odata:
        odata.write(calculate_contamination_script)
        logger.info('Calculating Contamination Script {script} generated!'.format(script=calculate_contamination_script_ofile))
        config_job.append([calculate_contamination_script_ofile,'1G',",".join(gather_ofile_rely)])

    ## Filter Shell
    hg37_ref_img = config.get('db','hg37_ref_img')
    filter_script = textwrap.dedent(
        """
        set -eo pipefail
        cd {directory}
        echo start `date +"%F %T"`
        ## Marking variants
        gatk --java-options "-Xmx500m" FilterMutectCalls \\
            -R {hg37_ref} \\
            -V  {tag}.unfiltered.vcf \\
            -O {tag}.mark.vcf \\
            --ob-priors {tag}.artifact-priors.tar.gz \\
            -stats {tag}.unfiltered.merged.stats \\
            --filtering-stats {tag}.filtering.stats \\
            --contamination-table {tag}.contamination.table
        
        ## extracting PASS variants with more than 4 reads covering
        bcftools filter -i 'INFO/DP>=4&&%FILTER=="PASS"' {tag}.mark.vcf -o {tag}.filtered.vcf
        
        ## split snp and indel
        bcftools filter -i '%TYPE==\"snp\"' {tag}.filtered.vcf |bgzip -c > {tag}.snp.filtered.vcf.gz
        bcftools filter -i '%TYPE==\"indel\"' {tag}.filtered.vcf |bgzip -c > {tag}.indel.filtered.vcf.gz
        bcftools index -t {tag}.snp.filtered.vcf.gz
        bcftools index -t {tag}.indel.filtered.vcf.gz

        ## filter alignment artifacts
        #gatk --java-options "-Xmx500m" FilterAlignmentArtifacts \\
        #    -R {hg37_ref} \\
        #    -V {tag}.filtered.vcf \\
        #    -I {tag}.sorted.out.bam \\
        #    --bwa-mem-index-image {hg37_ref_img} \\
        #    -O {tag}.filtered.rm_artifacts.vcf.gz
        ## extracting PASS variants
        #bcftools filter -i '%FILTER=="PASS"' {tag}.filtered.rm_artifacts.vcf.gz -Oz -o {tag}.filtered.final.vcf.gz

        ## split snp and indel variants
        #bcftools filter -i '%TYPE=="snp"' {tag}.filtered.final.vcf.gz -Oz -o {tag}.snp.filtered.final.vcf.gz
        #bcftools filter -i '%TYPE=="indel"' {tag}.filtered.final.vcf.gz -Oz -o {tag}.indel.filtered.final.vcf.gz
        echo end `date +"%F %T"`
        """.format(**locals())
    )
    filter_script_ofile = os.path.join(directory,'{tag}.filter.sh'.format(tag=tag))
    with open(filter_script_ofile,'w') as odata:
        odata.write(filter_script)
        logger.info('Filter script {script} generated!'.format(script=filter_script_ofile))
    filter_rely = ",".join([(calculate_contamination_script_ofile),
        (model_script_ofile)])
    config_job.append([filter_script_ofile,'1G',filter_rely])
    
    return config_job,filter_script_ofile

def make_job(config_job,main_path):
    config_job_ofile = os.path.join(main_path,'Somatic.config_job')
    with open(config_job_ofile,'w') as odata:
        for line in config_job:
            odata.write("\t".join(line) + '\n')
        logger.info('Configeration {config_job_ofile} used to make job generated!'.format(
            config_job_ofile=config_job_ofile
        ))
    os.system('cd {main_path} && makejob {config_job_ofile} --outfile Somatic.SomaticCalling.job \
--queue "disease1.q,all.q"'.format(**locals()))

def sort_interval(key):
    """
    sort interval by first number
    """
    return int(key.split('-')[0])

def generating_get_pileup_summary_script(gnomad_biallelic_snp_variants,chroms,bams_dict,tag,directory,config_job,hg37_ref,sample):
    # if isinstance(intervals_dir,list) and isinstance(bam,dict): 
    # intervals_file = sorted(os.listdir(intervals_dir),key=sort_interval)
    files_list = []
    # for idx,interval in enumerate(intervals_file):
    #     interval_path = os.path.join(intervals_dir,interval)
    for chrom in chroms:
        bam = bams_dict[chrom]
        get_pileup_summary_script = textwrap.dedent(
            """
            set -eo pipefail
            cd {directory}
            echo start `date +"%F %T"`
            gatk --java-options "-Xmx5g" GetPileupSummaries \\
                -R {hg37_ref} \\
                -I {bam} \\
                -V {gnomad_biallelic_snp_variants} \\
                --interval-set-rule INTERSECTION \\
                -L {chrom} \\
                -O {sample}.{tag}.{chrom}.pileups.table
            echo end `date +"%F %T"`
            """.format(**locals())
        )
        get_pileup_summary_ofile = os.path.join(directory,'{sample}.{tag}.{chrom}.GetPileupSummaries.sh'.format(
            chrom=chrom,tag=tag,sample=sample)) 
        with open(get_pileup_summary_ofile,'w') as odata:
            odata.write(get_pileup_summary_script)
            logger.info('GetPileupSummries script {} generated!'.format(get_pileup_summary_ofile))
            files_list.append(get_pileup_summary_ofile)
            config_job.append([get_pileup_summary_ofile,'5G'])
    #print(files_list)
    return ",".join(files_list)

def merge_and_annovar(samples,main_path,hg37_ref,extract_py,config_job,filter_script_ofile_list):
    """
    One sample directly generate annovar file
    More than one: merge vcf file and annotate with 
    annovar and extract the result to each samples
    """
    if len(samples) == 1:
        sample_dir = os.path.join(main_path,'Somatic',samples)
        for var in ['snp','indel']:
            annovar_shell = os.path.join(sample_dir,'work.{samples}.annovar.sh'.format(samples=samples))
            vcf = os.path.join(sample_dir,'{samples}.{var}.filtered.vcf.gz'.format(samples=samples,var=var))
            with open(annovar_shell,'w') as odata:
                annovar_script = textwrap.dedent(
                    '''
                    set -eo pipefail
                    cd {sample_dir}
                    echo start `date +"%F %T"`
                    # annovar with cosmic90
                    pyannovar {vcf} -p 'GeneName,refGene,\
Gencode,cpgIslandExt,cytoBand,wgRna,targetScanS,tfbsConsSites,\
genomicSuperDups,Repeat,avsnp,clinvar_20170905,gwasCatalog,\
1kg_Chinese,1000g_EAS,1000g_ALL,esp6500si_all,\
gnomad_gnome_exome_ALL_AF_AN,gnomad_gnome_exome_EAS_AF_AN,\
NovoDb_WES_2573,NovoDb_WGS_568,dbscsnv,Spidex,\
dbnsfp31a_interpro,dbnsfp293apart,phylopcadd,gerp++gt2,\
mcap,revel,cosmic90'

                    #cut of shared_hom and shared_het
                    mv {var}.annovar.hg19_multianno.xls.gz {var}.annovar.hg19_multianno.xls.bak.gz
                    getcol -e --col shared_hom,shared_het {var}.annovar.hg19_multianno.xls.bak.gz|bgzip -c >\
                    {var}.annovar.hg19_multianno.xls.gz

                    echo end `date +"%F %T"`
                    ''')
                odata.write(annovar_script)
                config_job.append([annovar_shell,'5G',",".join(filter_script_ofile_list)])

    else:
        VCF_path = os.path.join(main_path,'VCF')
        if not os.path.exists(VCF_path):
            os.mkdir(VCF_path)
        #Generate snp and indel merge and annovar shell
        merge_vcf_shell_list = []
        for var in ['snp','indel']:
            merge_vcf_shell = os.path.join(VCF_path,'{var}.merged_annovar.sh'.format(var=var))
            merge_vcf_shell_list.append(merge_vcf_shell)
            merge_vcf =  os.path.join(VCF_path,'{var}.merged.vcf.gz'.format(var=var))
            vcfs_list = [os.path.join(main_path,'Somatic',s,
            '{s}.{var}.filtered.vcf.gz'.format(s=s,var=var)) for s in samples]
            vcfs = " ".join(vcfs_list) 
            with open(merge_vcf_shell,'w') as odata:
                shell = textwrap.dedent('''
                set -eo pipefail
                cd {VCF_path}
                echo start `date +"%F %T"`
                #merge vcf
                bcftools-1.9 merge -Oz -o {merge_vcf} {vcfs}
                
                # annovar with cosmic90
                pyannovar {merge_vcf} -p 'GeneName,refGene,\
Gencode,cpgIslandExt,cytoBand,wgRna,targetScanS,tfbsConsSites,\
genomicSuperDups,Repeat,avsnp,clinvar_20170905,gwasCatalog,\
1kg_Chinese,1000g_EAS,1000g_ALL,esp6500si_all,\
gnomad_gnome_exome_ALL_AF_AN,gnomad_gnome_exome_EAS_AF_AN,\
NovoDb_WES_2573,NovoDb_WGS_568,dbscsnv,Spidex,\
dbnsfp31a_interpro,dbnsfp293apart,phylopcadd,gerp++gt2,\
mcap,revel,cosmic90'
                
                #cut of shared_hom and shared_het
                mv {var}.annovar.hg19_multianno.xls.gz {var}.annovar.hg19_multianno.xls.bak.gz
                getcol -e --col shared_hom,shared_het {var}.annovar.hg19_multianno.xls.bak.gz|bgzip -c >\
                {var}.annovar.hg19_multianno.xls.gz

                echo end `date +"%F %T"`
                '''.format(**locals()))
                odata.write(shell)
                config_job.append([merge_vcf_shell,'5G',",".join(filter_script_ofile_list)])
        # Generate extract shell for each samples
        for s in samples:
            for var in ['snp','indel']:
                directory = os.path.join(main_path,'Somatic',s)
                extract_shell = os.path.join(directory,'work.{s}.extract.{var}.sh'.format(s=s,var=var))
                with open(extract_shell,'w') as odata:
                    extract_script = textwrap.dedent('''
                    set -eo pipefail
                    cd {directory}
                    # norm
                    bcftools-1.9 norm \\
                        -f {hg37_ref}\\
                        -m -both \\
                        -Oz -o {directory}/{s}.gatk.{var}_sn.vcf.gz \\
                        {directory}/{s}.{var}.filtered.vcf.gz

                    python {extract_py} \\
                        -M {main_path}/VCF/{var}.annovar.hg19_multianno.xls.gz \\
                        -nv {directory}/{s}.gatk.{var}_sn.vcf.gz \\
                        -rv {directory}/{s}.{var}.filtered.vcf.gz \\
                        -S {s} -T {var} \\
                        -ref hg19 \\
                        -O {directory}/{s}.gatk.{var}.annovar.hg19_multianno.xls

                    gzip -f {directory}/{s}.gatk.{var}.annovar.hg19_multianno.xls
                    rm -f {directory}/{s}.gatk.{var}_sn.vcf.gz
                    '''.format(**locals()))
                    odata.write(extract_script)
                    config_job.append([extract_shell,'1G',','.join(merge_vcf_shell_list)])
    return config_job

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers()
    subcommand_mkconfig = subparser.add_parser('mkconfig',
        help='Making config for call short somatic mutations')
    subcommand_mkconfig.add_argument('--database_path','-db',
        help='The path to database,default is the upper directory of script')
    subcommand_mkconfig.set_defaults(func=generating_config)

    subcommand_mkscript = subparser.add_parser('mkscript',
        help='Making scripts for call short somatic mutations')
    subcommand_mkscript.add_argument('--sample_info','-info',required=True,
        help='sample info with family meta information,Required')
    subcommand_mkscript.add_argument('--v4.7gatk_mutation','-gatk',action='store_true',
        help='get realigned bams generated by BQSR from Mutation in v4.7 main pipeline')
    subcommand_mkscript.add_argument('--main_path','-path',default = os.environ['PWD'],
        help='The main project directory containing Mapping directory,default current directory')
    subcommand_mkscript.add_argument('--database_path','-db',
        help='The path to database,default is the upper directory of script')
    subcommand_mkscript.set_defaults(func=generating_scripts)
    args = parser.parse_args()
    args.func(vars(args))