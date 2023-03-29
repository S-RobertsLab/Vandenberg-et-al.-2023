import csv
import time
import subprocess
import os
import re

with open('AAAAAAAAAA_discard.txt', 'w', newline='') as o:
    pass

def write_to_discard(list):
    with open('AAAAAAAAAA_discard.txt', 'a', newline='') as o:
        writer = csv.writer(o, delimiter='\t')
        writer.writerow(list)


def uvc_to_universal(input_filename, output_filename):
    with open(input_filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        with open(output_filename, 'w', newline='') as o:
            writer = csv.writer(o, delimiter='\t')
            writer.writerow(['CHROM','START','END','REF','VAR','MUT','VAR_TYPE','SELECTED_ELEMENT','GENOTYPE'])
            header = next(reader)
            with open('AAAAAAAAAA.txt', 'w', newline = '') as test:
                test_writer = csv.writer(test, delimiter='\t')
                for line in reader:
                    selected_element = line[0]
                    genotype = line[1]
                    dose = line[2]
                    photoreact = line[4]
                    chrom = line[6]
                    pos0 = int(line[8])-1
                    pos1 = int(line[8])
                    var_type = line[10]
                    ref = line[12]

                    # skips mitochondrial chromosomes
                    if 'Mito' in chrom:
                        continue

                    # checks for dose 0, if dose 0, it skips
                    if dose == '0': 
                        # print(line)
                        write_to_discard(line)
                        continue

                    # checks for protoreactivated samples, if it is reactivated, it skips
                    if photoreact != '0':
                        # print(line)
                        write_to_discard(line)
                        continue

                    # checks if variant type is SNP and changes it to SNV, if not it skips
                    if var_type == 'SNP':
                        var_type = 'SNV'
                    else: 
                        # print(line)
                        write_to_discard(line)
                        continue

                    # checks if the reference base is alphabetical, if not it skips
                    if not ref.isalpha(): 
                        # print(line)
                        write_to_discard(line)
                        continue

                    # determines which column to use as the variant base by using the column its not in
                    if line[18]== ref: var = line[21]
                    else: var = line[18]
                    
                    # skips the mutation if the variant base is not alphabetical
                    if not var.isalpha(): 
                        # print(line)
                        write_to_discard(line)
                        continue

                    mut = f'{ref}>{var}'
                    test_writer.writerow(line)
                    writer.writerow([chrom,pos0,pos1,ref,var,mut,var_type,selected_element,genotype])

def uvb_to_universal(input_filename, output_filename):
    with open(input_filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        with open(output_filename, 'w', newline='') as o:
            writer = csv.writer(o, delimiter='\t')
            writer.writerow(['CHROM','START','END','REF','VAR','MUT','VAR_TYPE','SELECTED_ELEMENT','GENOTYPE'])
            header = next(reader)
            for line in reader:
                selected_element_geno = line[0]
                genotype = selected_element_geno.split(' ')[1]
                selected_element = selected_element_geno.split(' ')[0]
                chrom = line[2]
                if 'Mito' in chrom: continue
                var_type = line[11]
                if var_type not in ['SNV', 'MNV']: continue
                pos0 = int(line[3])-1
                pos1 = int(line[3])
                if var_type == 'MNV': pos1 += 1
                ref = line[12]
                if len(ref) > 2: continue
                var = line[13]
                mut = f'{ref}>{var}'
                writer.writerow([chrom,pos0,pos1,ref,var,mut,var_type,selected_element,genotype])

def split_var_type(input_filename):
    with open(input_filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        with open(f'SNV_{input_filename}', 'w', newline='') as snv:
            snv_writer = csv.writer(snv, delimiter='\t')
            with open(f'MNV_{input_filename}', 'w', newline='') as mnv:
                mnv_writer = csv.writer(mnv, delimiter='\t')
                header = next(reader)
                snv_writer.writerow(header)
                mnv_writer.writerow(header)
                for line in reader:
                    if line[6] == 'SNV':
                        snv_writer.writerow(line)
                    if line[6] == 'MNV':
                        mnv_writer.writerow(line)

def split_snv_file(input_filename):
    with open(input_filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        with open(f'{os.path.dirname(input_filename)}/SNV_{os.path.basename(input_filename)}', 'w', newline='') as snv:
            snv_writer = csv.writer(snv, delimiter='\t')
            with open(f'{os.path.dirname(input_filename)}/MNV_{os.path.basename(input_filename)}', 'w', newline='') as mnv:
                mnv_writer = csv.writer(mnv, delimiter='\t')
                header = next(reader)
                snv_writer.writerow(header)
                mnv_writer.writerow(header)
                prev_line = None
                # Iterate over the lines in the input file
                for line in reader:
                    # Check if the previous line is not None
                    if prev_line:
                        # Get the value of column 2 from the previous line and the current line
                        prev_start_pos = int(prev_line[1])
                        prev_chrom = prev_line[0]
                        curr_start_pos = int(line[1])
                        curr_chrom = line[0]
                        # Check if the value of column 2 in the previous line is one less than the next line
                        if curr_start_pos == prev_start_pos + 1 and curr_chrom == prev_chrom:
                            # Write both lines to the output file
                            mnv_writer.writerow([curr_chrom, prev_start_pos, line[2], prev_line[3]+line[3], prev_line[4]+line[4], prev_line[3]+line[3]+'>'+prev_line[4]+line[4],'MNV',line[7],line[8]])
                            prev_line = next(reader)
                        else:
                            snv_writer.writerow(prev_line)
                    # Update the previous line
                    prev_line = line
                snv_writer.writerow(line)

def expand_positional_values(bed_filename: str, reference_filename: str, five_prime_increase: int, three_prime_increase: int):
        """
        Takes a .bed file and returns another .bed file that has increased nucleotide context. (positional values remain the same)
        The context can be increased by fivePrimeCount and threePrimeCount.
        Put any input files (mutation and FASTA files) into the 'Input' directory.
        """
        with open(bed_filename, 'r') as bed_file:
            header = bed_file.readline()
            # creates a temporary copy of the .bed mutation file to
            # change positional numbers according to the values in definition.
            with open(bed_filename[:-4] + '_expanded_positions.bed', 'w') as expanded_potition_bed_file:
                print(  f"Adding nucleotide context...")
                print(  f"From {os.path.basename(bed_filename)}, adding {five_prime_increase} "
                        f"to the 5-prime end and adding {three_prime_increase} to the 3-prime end...")
                for line in bed_file:
                    tsv = line.strip().split('\t')
                    # decrease and increase the positional numbers so when
                    # compared to the fasta file, it will return the nucleotides
                    # that are in the newly defined range
                    new_five_prime_position = int(tsv[1]) - five_prime_increase
                    new_three_prime_position = int(tsv[2]) + three_prime_increase
                    # copy the file but use the new positional numbers
                    expanded_potition_bed_file.write(   '\t'.join([tsv[0],
                                                                str(new_five_prime_position),
                                                                str(new_three_prime_position)] + tsv[3:])
                                                        + '\n')
        # terminal portion to compare new .bed file with .fasta input file and return the nuceotides missing
        time.sleep(1) # Bottleneck for terminal progress bar
        print(f"Running bedtools getfasta using {os.path.basename(reference_filename)}...")
        with subprocess.Popen(  args=[f'bedtools getfasta -fi {reference_filename} -bed {bed_filename[:-4]}_expanded_positions.bed -fo {bed_filename[:-4]}_expanded_context.fa -tab'], stdout=subprocess.PIPE,
                                shell=True) as p: # .replace(' ', '\\ ')

            # This for loop should output everything from the commands as the process is running
            for text in p.stdout:
                # I had to make the output a string and specify it's encoding so it doesn't come out all wonky
                print(text)
        # output file from terminal fasta comparison to confirm that they have the same positional data and lines aren't missing.
        with open(bed_filename, 'r') as bed_file:
            header = bed_file.readline()
            with open(bed_filename[:-4] + '_expanded_context.bed', 'w') as bed_file_expanded_context:
                bed_file_expanded_context.write(header.replace('VAR\t', ''))
                with open(bed_filename[:-4] + '_expanded_context.fa', 'r') as context_fasta:
                    print(f"Using intermediate fasta file to add nucleotide context...")
                    it_worked = True
                    for line in bed_file:
                        # takes the expanded nucleotide context and confirms that it processed correctly.
                        tsv = line.strip().split('\t')
                        bed_chrom_num = tsv[0]
                        bed_five_prime = tsv[1]
                        bed_three_prime = tsv[2]
                        bed_nuc = tsv[3]
                        rest_of_line = tsv[5:]
                        rest_of_line='\t'.join(rest_of_line)
                        # sets line with nucleotides as variable
                        fasta_info = context_fasta.readline()
                        check = re.split(':|-|\t', fasta_info, maxsplit=4)
                        fasta_chrom_num = check[0]
                        fasta_five_prime = check[1]
                        fasta_three_prime = check[2]
                        fasta_context = check[3].upper()
                        if bed_chrom_num == fasta_chrom_num and int(bed_five_prime) == (int(fasta_five_prime)+five_prime_increase) and int(bed_three_prime) == (int(fasta_three_prime)-three_prime_increase):
                            if bed_nuc == fasta_context[1:2]:
                                bed_file_expanded_context.write(f'{bed_chrom_num}\t{bed_five_prime}\t{bed_three_prime}\t{fasta_context.strip()}\t{rest_of_line}\n')
                            else:
                                # print(bed_nuc, fasta_context[1:2])
                                print(f"Error in line \'{fasta_info.strip()}\'")
                                print(f'bed file line \"{line.strip()}\"')
                                print(f"Aborting now, ERROR CODE 1.")
                                it_worked = False
                                break
                        else:
                            print(f"Error in line \'{fasta_info.strip()}\'")
                            print(f'bed file line \"{line.strip()}\"')
                            print(f"Aborting now, ERROR CODE 2.")
                            it_worked = False
                            break
        if it_worked:
            print(f"File validation completed! Removing intermediate files...")
            os.remove(bed_filename[:-4] + '_expanded_positions.bed')
            os.remove(bed_filename[:-4] + '_expanded_context.fa')
            print(f"Intermediate files removed! Process complete")
            print(bed_filename[:-4] + '_expanded_context.bed')

def create_sbs_keys(input_file_T, input_file_C):
    with open(input_file_T, 'r') as t:
        with open(input_file_C, 'r') as c:
            dict = {}
            for line in t:
                dict.setdefault(line[0]+'[T>A]'+line[2], 0)
            t.seek(0)
            for line in t:
                dict.setdefault(line[0]+'[T>C]'+line[2], 0)
            t.seek(0)
            for line in t:
                dict.setdefault(line[0]+'[T>G]'+line[2], 0)
            for line in c:
                dict.setdefault(line[0]+'[C>A]'+line[2], 0)
            c.seek(0)
            for line in c:
                dict.setdefault(line[0]+'[C>T]'+line[2], 0)
            c.seek(0)
            for line in c:
                dict.setdefault(line[0]+'[C>G]'+line[2], 0)
            return dict

def create_dnv_keys(input_file):
    with open(input_file, 'r') as f:
        key_dict = {}
        for line in f:
            key_dict.setdefault(line.strip(), 0)
        return key_dict

def reverse_complement(seq: str):
    '''
    Takes a nucleotide string in IUPAC and regular format and returns the reverse complement
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                        "R":"Y", "Y":"R", "S":"S", "W":"W", "K":"M",
                        "M":"K","B":"V", "D":"H", "H":"D", "V":"B",
                        "N":"N"}
    return "".join(complement.get(base, base) for base in reversed(seq))

def create_graph_code(input_filename):
    with open(input_filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        data_dict = {}
        for line in reader:
            genotype = line[7]
            element = line[6]
            data_dict.setdefault(genotype, {}).setdefault(element, create_sbs_keys('/media/cam/Data8/Brittany_Reversion_Paper/random/dumb_file_T.txt','/media/cam/Data8/Brittany_Reversion_Paper/random/dumb_file_C.txt'))
        f.seek(0)
        header = next(reader)
        i = 0
        for line in reader:
            i += 1
            context = line[3]
            mut = line[4]
            genotype = line[7]
            element = line[6]
            counting_key = f'{context[0]}[{mut}]{context[-1]}'
            try:
                data_dict[genotype][element][counting_key] += 1
            except KeyError:
                counting_key = f'{reverse_complement(context)[0]}[{reverse_complement(mut[0])}>{reverse_complement(mut[-1])}]{reverse_complement(context)[-1]}'
        # print(data_dict.keys())
        for dict_geno in data_dict.keys():
            sample_list = list(data_dict[dict_geno].keys())
            # print(sample_list)
            sample_header = ''+'\t'+'\t'.join(sample_list)
            mut_type_list = list(create_sbs_keys('/media/cam/Data8/Brittany_Reversion_Paper/random/dumb_file_T.txt','/media/cam/Data8/Brittany_Reversion_Paper/random/dumb_file_C.txt').keys())
            with open(f'{os.path.dirname(input_filename)}/{dict_geno}_UVC_graph_data.txt', 'w') as o:
                o.write(sample_header+'\n')
                for mut_type in mut_type_list:
                    count_list = []
                    count_list.append(str(mut_type))
                    for sample in sample_list:
                        count_list.append(str(data_dict[dict_geno][sample][mut_type]))
                    line = '\t'.join(count_list)
                    o.write(f'{line}\n')


def create_dnv_graph_code(input_filename):
    with open(input_filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        data_dict = {}
        for line in reader:
            genotype = line[7]
            element = line[6]
            mutation = line[4]
            data_dict.setdefault(genotype, {}).setdefault(element, create_dnv_keys('MNV_Key.txt'))
        f.seek(0)
        header = next(reader)
        for line in reader:
            genotype = line[7]
            element = line[6]
            mutation = line[4]
            try: data_dict[genotype][element][mutation] += 1
            except KeyError:
                try:
                    new_key = reverse_complement(mutation[:2])+'>'+reverse_complement(mutation[-2:])
                    data_dict[genotype][element][new_key] += 1
                except: pass
        for genotype in data_dict.keys():
            sample_list = list(data_dict[genotype].keys())
            sample_header = ''+'\t'+'\t'.join(sample_list)+'\n'
            with open(f'{os.path.dirname(input_filename)}/{genotype}_graph_data_{os.path.basename(input_filename)}', 'w') as o:
                o.write(sample_header)
                for mutation in create_dnv_keys('MNV_Key.txt').keys():
                    count_list = []
                    count_list.append(str(mutation))
                    for sample in data_dict[genotype].keys():
                        count_list.append(str(data_dict[genotype][sample][mutation]))
                    line = '\t'.join(count_list)
                    o.write(f'{line}\n')
