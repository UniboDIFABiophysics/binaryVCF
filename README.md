# binaryVCF
python script that converts genetic variants contained in VCF files into a matrix in which
each row represents a variant reported in the provided VCF files at least once (encoded as
chromosome number-position-reference allele- alternate allele) and each column is named after an 
ID assigned to each sample. In each cell of the matrix is reported the number of alternative alleles for
each locus, thus 0 indicates that the variant is not present in the VCF file of the patient whereas 1
indicates its presence in heterozygosity and 2 the presence of the variant in homozygosity.

