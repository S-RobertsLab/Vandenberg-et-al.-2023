import subprocess

#################################
### CAM'S CODE IMPORTANT INFO ###
### ------------------------- ###
### I use tsv as a variable   ###
### frequently, it means "tab ###
### separated values" and its ###
### a list of the columns in  ###
### each line for files, and  ###
### these get indexed often   ###
### looking like tsv[index].  ###
#################################

def parse_origins(origin_file, output_file):
    '''
    Converts the file from John/Steve and converts it to a BED3 type format
    '''
    with open(origin_file, 'r') as f:
        with open(output_file, 'w') as o:
            f.readline()
            o.write('CHROM\tSTART\tEND\tORIGINS\n')
            # initializes the first line in the file and all of the data
            first_line = f.readline()
            prev_tsv = first_line.strip().split('\t')
            prev_name = prev_tsv[0]
            prev_chrom = int(prev_tsv[2])
            prev_end = int(prev_tsv[4])
            temp_end = 0
            overlap = False
            for line in f:
                # grabs the specific information from the next line in the file
                tsv = line.strip().split('\t')
                name = tsv[0]
                chrom = int(tsv[2])
                start = int(tsv[3])
                end = int(tsv[4])
                # checks if they're on the same chromosome
                if chrom > prev_chrom:
                    overlap = False
                    temp_end = 0
                # checks when no overlap and origin follows
                elif not overlap and chrom == prev_chrom and start > prev_end:
                    # writes from the previous origin to the current one
                    o.write(f'chr{chrom}\t{prev_end}\t{start}\t{prev_name}-{name}\n')
                # if there is no longer overlap, write the farthest position to the start of the next
                elif overlap and chrom == prev_chrom and start > temp_end:
                    o.write(f'chr{chrom}\t{temp_end}\t{start}\t{temp_name}-{name}\n')
                    overlap = False
                    temp_end = 0
                # detects if there is any overlap and holds the farthest endpoint in the overlap
                elif chrom == prev_chrom and prev_end >= start:
                    overlap = True
                    if temp_end == 0:
                        temp_end = prev_end
                        temp_name = prev_name
                    if end > temp_end: temp_end = end
                    continue
                else:
                    print('ERROR!')
                # reassigns the data to the previous line before reading in the new line
                prev_tsv = tsv
                prev_name = name
                prev_chrom = chrom
                prev_end = end

def mut_file_to_bed(mutations_file, output_file):
    '''
    Takes one of Brittany's csv files and converts it to some BED3 type format
    '''
    with open(mutations_file, 'r') as f:
        with open(output_file, 'w') as o:
            # writes the new headers for the file
            o.write('CHROM\tSTART\tEND\tCONTEXT\tMUTATION\n')
            f.readline()
            for line in f:
                tsv = line.strip().split(',')
                # if there is photo-reactivation, it will skip writing that mutation
                if tsv[9] != '0': continue
                # if the mutation was not able to be counted, put an NA for the mutatio position
                if 'NA' in tsv[32]: mutation = 'NA'
                # if there was a mutation, write it from each trinucleotide context
                else: mutation = f'{tsv[29][2]}>{tsv[32][2]}'
                chrom = roman_to_int(tsv[11].strip('"'))
                start = int(tsv[3])+1
                end = int(tsv[4])-1
                context = tsv[29].strip('"')
                o.write(f'chr{chrom}\t{start}\t{end}\t{context}\t{mutation}\n')

def roman_to_int(s):
    """
    converts a string of roman numerals to an integer
    """
    roman = {'I':1,'V':5,'X':10,'L':50,'C':100,'D':500,'M':1000,'IV':4,'IX':9,'XL':40,'XC':90,'CD':400,'CM':900}
    i = 0
    num = 0
    while i < len(s):
        if i+1<len(s) and s[i:i+2] in roman:
            num+=roman[s[i:i+2]]
            i+=2
        else:
            num+=roman[s[i]]
            i+=1
    return num

def asses_location_relative_to_origins(mutation_file, origins_file, overlap_file, output_file):
    '''
    Using an origins file as well as a mutation file in BED3 type format, it intersects them and calculates the percentage it is towards the origin.
    It also uses bedtools intersect and writes a file thats similar to the BED3 mutation file as input
    '''
    with subprocess.Popen(args=[f'bedtools intersect -wa -wb -a {mutation_file} -b {origins_file} > {overlap_file}'],
                            stdout=subprocess.PIPE, shell=True) as p:
        for text in p.stdout:
            print(text)
    with open(overlap_file, 'r') as f:
        with open(output_file, 'w') as o:
            # writes headers
            o.write('CHROM\tSTART\tEND\tCONTEXT\tMUTATION\tPERCENT_BETWEEN_ORIGINS\n')
            for line in f:
                # reads in data from the overlapped sequence
                tsv = line.strip().split('\t')
                chrom = tsv[0]
                start = int(tsv[1])
                end = int(tsv[2])
                context = tsv[3].strip('"')
                mutation = tsv[4]
                # calculates a percentage based value for its relative location between origins, 0 being the start of
                # the upstream origin and 100 being the start of to the downstream origin, with 50 being directly in the middle
                position_percent = (start-int(tsv[6]))/(int(tsv[7])-int(tsv[6]))
                o.write(f'{chrom}\t{start}\t{end}\t{context}\t{mutation}\t{position_percent}\n')

def count_mutations_in_percentages(input_file, context_input, mutation_input, output_file):
    '''
    Takes the previously overlapped file with relative percentages, the mutation context you want, and the mutation type.
    From there it writes pooled mutation counts and percentages for the mutation and the reverse complement
    '''
    with open(input_file, 'r') as f:
        with open(output_file, 'w') as o:
            # a dictionary to translate iupac to nucleotide notation
            iupac_trans = {"R":"AG", "Y":"CT", "S":"GC", "W":"AT", "K":"GT",
                        "M":"AC","B":"CGT", "D":"AGT", "H":"ACT", "V":"ACG",
                        "N":"ACTG", 'A':'A', 'T':'T', 'C':'C', 'G':'G'}
            # finds the reverse complement of the mutation
            rev_comp_mutation = f'{reverse_complement(mutation_input[0])}>{reverse_complement(mutation_input[2])}'
            # this step finds the possible trinucleotide contexts based on iupac notation in a few steps
            # 1. find all possible characters for each position
            context1 = iupac_trans[context_input[0]]
            context2 = iupac_trans[context_input[1]]
            context3 = iupac_trans[context_input[2]]
            # 2. create the lists where each trinuc context will go
            possible_contexts = []
            rev_comp_contexts = []
            # 3. itterate through each possible nuc for each position and add the resulting trinuc to the list
            for letter1 in context1:
                for letter2 in context2:
                    for letter3 in context3:
                        possible_contexts.append(letter1+letter2+letter3)
            # 4. find the reverse complement of each trinuc and add it to the reverse complement list
            for item in possible_contexts:
                rev_comp_contexts.append(reverse_complement(item))
            # 5. create a dictionary with the binned percentage categories for each group
            context_percents = {0:0, 10:0, 20:0, 30:0, 40:0, 50:0, 60:0, 70:0, 80:0, 90:0}
            rev_comp_percents = {0:0, 10:0, 20:0, 30:0, 40:0, 50:0, 60:0, 70:0, 80:0, 90:0}
            f.readline()
            for line in f:
                # read in the data from each file
                tsv = line.strip().split('\t')
                context = tsv[3]
                mutation = tsv[4]
                percent = float(tsv[5])*100
                # determine if the mutation fits the possible context list and is the correct mutation and add one to the percent bin
                if context in possible_contexts and mutation == mutation_input:
                    if percent < 10: context_percents[0] += 1
                    elif percent < 20: context_percents[10] += 1
                    elif percent < 30: context_percents[20] += 1
                    elif percent < 40: context_percents[30] += 1
                    elif percent < 50: context_percents[40] += 1
                    elif percent < 60: context_percents[50] += 1
                    elif percent < 70: context_percents[60] += 1
                    elif percent < 80: context_percents[70] += 1
                    elif percent < 90: context_percents[80] += 1
                    elif percent < 100: context_percents[90] += 1
                # determine if the mutation fits the rev_comp list and is the correct mutation and add one to the percent bin
                elif context in rev_comp_contexts and mutation == rev_comp_mutation:
                    if percent < 10: rev_comp_percents[0] += 1
                    elif percent < 20: rev_comp_percents[10] += 1
                    elif percent < 30: rev_comp_percents[20] += 1
                    elif percent < 40: rev_comp_percents[30] += 1
                    elif percent < 50: rev_comp_percents[40] += 1
                    elif percent < 60: rev_comp_percents[50] += 1
                    elif percent < 70: rev_comp_percents[60] += 1
                    elif percent < 80: rev_comp_percents[70] += 1
                    elif percent < 90: rev_comp_percents[80] += 1
                    elif percent < 100: rev_comp_percents[90] += 1
            # write all the bins to a file and show percentages relative to each bin category
            o.write(f'Region Counts for {mutation_input} mutations in {context_input} contexts\npercent_region_between_origins\tmutation_count\tpercent_of_bin\n')
            for key, value in context_percents.items():
                percentage = (context_percents[key])/(context_percents[key]+rev_comp_percents[key])
                o.write(f'{key}\t{value}\t{percentage}\n')
            o.write('\n')
            o.write(f'Region Counts for {rev_comp_mutation} mutations in {reverse_complement(context_input)} contexts\npercent_region_between_origins\tmutation_count\tpercent_of_bin\n')
            for key, value in rev_comp_percents.items():
                percentage = (rev_comp_percents[key])/(context_percents[key]+rev_comp_percents[key])
                o.write(f'{key}\t{value}\t{percentage}\n')

def reverse_complement(seq: str):
    '''
    Takes a nucleotide string in IUPAC and regular format and returns the reverse complement
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                        "R":"Y", "Y":"R", "S":"S", "W":"W", "K":"M",
                        "M":"K","B":"V", "D":"H", "H":"D", "V":"B",
                        "N":"N"}
    return "".join(complement.get(base, base) for base in reversed(seq))


##################################################################################
#################### EXAMPLE PIPELINE FOR WT DATA SHOWN BELOW ####################
##################################################################################

# # Takes the Saccer3 origins of replication coordinates and convers it to a BED3 format used in the intersection
# parse_origins(  origin_file= 'Saccer3_OriginsOfRep.tsv', 
#                 output_file= 'Saccer3_Origins.bed')

# # Converts Brittany's mutation file to a usable BED3 format used in the intersection
# mut_file_to_bed(mutations_file= 'AllWT_UVC_snv.csv', 
#                 output_file='WT_snp.bed')

# # Takes the mutation output file from the previous step, and the origins file from the first step (not listed)
# # and intersects them creating an overlap file then creates the final percent mutation file
# # asses_location_relative_to_origins( mutation_file= '/media/cam/Data8/Brittany_Reversion_Paper/WT_snp.bed',
#                                     origins_file= '/media/cam/Data8/Brittany_Reversion_Paper/Saccer3_Origins.bed', 
#                                     overlap_file= '/media/cam/Data8/Brittany_Reversion_Paper/WT_Origins_overlap.tsv', 
#                                     output_file= '/media/cam/Data8/Brittany_Reversion_Paper/WT_OriginPercentage.bed')

# # Counts and calculates percents based on the iupac notation context and mutation given.
# count_mutations_in_percentages( input_file= 'WT_OriginPercentage.bed', 
#                                 context_input= 'NTA',
#                                 mutation_input= 'T>A',
#                                 output_file= 'WT_TtoA_NTA.txt')

# count_mutations_in_percentages( input_file= 'WT_OriginPercentage.bed',
#                                 context_input='YCN',
#                                 mutation_input='C>T',
#                                 output_file= 'WT_CtoT_YCN.txt')

# count_mutations_in_percentages( input_file='WT_OriginPercentage.bed',
#                                 context_input='TTN', 
#                                 mutation_input='T>C',
#                                 output_file= 'WT_TtoC_TTN.txt')