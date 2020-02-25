# BE-tool

while IFS= read -r line; do python find_sgRNAs_BE.py -g $line -b R33A -o ../TF/; done < TF.txt      
