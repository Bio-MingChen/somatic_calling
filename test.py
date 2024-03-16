def bam_map_to_dict(bams,chroms):
    """
    mapping bams to every chroms
    """
    if isinstance(bams,list):
        return {c:b for b,c in zip(bams,chroms)}
    return {c:bams for c in chroms}

if __name__ == "__main__":
    bams = ['a','b','c']
    chroms = ['1','2','3']
    print(bam_map_to_dict(bams,chroms))