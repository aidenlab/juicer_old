To run the code, you need an input file where each line represent a chromosome and is formatted:

		chr_name    observed_file_name    expected_file_name    norm_file_name    resolution

for chr_name, use only integers or X. observed files need to be binary outputs of juicebox dump. expected and norm are just the normal text files from juicebox dump. 

run the code using the following command:

		peakcallingGPU18_short3.py    input_file_name    output_file_name1    output_file_name2    fdr    p    w

output_file_name1 gets the fdr threshold values printed to it

output_file_name2 gets the enriched pixels

fdr actually should be an int corresponding to 1/max-q-val (i.e. for 1% FDR do 100, for 10%FDR do 10)

p = peak width

w = window


for fdr, p, and w parameters, look back at what we did in the paper and the parameter sweep stuff i sent you.

		peakcallingGPU_listgiven_edges_short_REVISED.py    input_file_name    output_file_name1    output_file_name2    fdr     input_list    p    w

the input_list should be a file with six fields:

		chr1    x1    x2    chr2    y1    y2

i.e. basically a looplist (no header needed); use 23 instead of X (even for mouse) and make sure that the input list is at whatever resolution that you're running HiCCUPS at.