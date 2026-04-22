# Satellite Insertions
Detection of satellite DNA insertions from 10x Linked-Read whole genome sequencing.

Usage:

python find_sattels.py -i phased_possorted_bam.bam -o my_sample

Required parameters:

| Argument | Description |
| --- | --- | 
| `-i, --input_bam` | Alignment from 10x lariat (by default the file phased_possorted_bam) |
| `-o, --out_prefix` | Prefix for results |

Optional parameters:

| Argument | Description |
| --- | --- | 
| `-t, --min_num_tel` | Minimum number of telomeric repeats TTAGGG within a read to be counted as telomeric (default 4)|
| `-s, --min_num_sat` | Minimum number of telomeric reads within a molecule (default 1)|
| `-m, --max_mm_sat` | Maximum number of allowed mismatches when detecting the 17bp alpha-satellite motif (default 2)|

Python package requirements:

sys\
getopt\
time\
pysam\
collections\
regex
