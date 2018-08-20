This project uses the Needleman-Wunsch algorithm to get the global alignment between two protein sequences. In this instance those two sequences are the actin binding domains of D. melanogaster Shot and  human ACF7, but this program can be used for other sequences too. The program can either use PAM250 or BLOSUM62 as the scoring matrix for comparing protein sequences. The program does not consider the affire gap problem. The program outputs its results to the console, and you can copy and paste from there.

This repl contains:

main.py: The file containing the code of the program.

PAM250.txt: The PAM250 scoring matrix.

BLOSUM62.txt: The BLOSUM62 scoring matrix.

Shot.txt: The FASTA protein sequence for D. melanogaster Shot

hACF7.txt: The FASTA protein sequence for human ACF7.

misc testing data: A collection of snippets I used to troubleshoot my code that are not relevant to final product.

To use this program on different strings, change the assignments of string1IO and string2IO.