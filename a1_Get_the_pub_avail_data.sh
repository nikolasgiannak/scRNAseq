
# How to retrieve publically available data


wget -nH --cut-dirs=5  --no-clobber --convert-links --random-wait -r -p -E -e robots=off https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154763/suppl/

ls *gz | xargs gunzip

# Convert the log transformed counts (as found in the published available data) to raw counts
# To do that we will use the script from Immunitas that can be found here:

#  https://github.com/immunitastx/recover-counts


## for one file
 ./recover_counts_from_log_normalized_data.py -m 10000 -d CSV GSE154763_ESCA_normalized_expression.csv -o GSE154763_ESCA_raw_counts.txt
 
## for all files
# mamba install parallel 
# regular expression to remove the normalized_expression.csv and rename the output to _raw_counts.txt
ls *expression.csv | parallel --rpl '{%(.+?)} s/$$1$//;' ./recover_counts_from_log_normalized_data.py -m 10000 -d CSV {} -o {%_normalized_expression.csv}_raw_count.txt
 