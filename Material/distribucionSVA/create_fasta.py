import re

def vcf_to_fasta(vcf_file, fasta_file, line_width=80):
    with open(vcf_file, 'r') as vcf, open(fasta_file, 'w') as fasta:
        for line in vcf:
            if line.startswith('#'):
                continue 
            columns = line.strip().split('\t')
            chrom = columns[0]
            pos = columns[1]
            alt = columns[4]


            clean_alt = re.sub(r'[<>[\]*]', '', alt)

            header = f">{chrom}:{pos}"
            fasta.write(f"{header}\n")
            
            
            for start in range(0, len(clean_alt), line_width):
                fasta.write(f"{clean_alt[start:start+line_width]}\n")


vcf_to_fasta('only_ins.vcf', 'seq_ins_pat.fasta')

