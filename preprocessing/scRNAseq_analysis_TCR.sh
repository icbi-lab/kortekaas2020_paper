#$ -S /bin/sh

# Options
while getopts i:f:s:n:o: opt;
do
    case $opt in
        i)  id="$OPTARG";;
        f)  fastqs="$OPTARG";;
        s)  sample="$OPTARG";;
	n)  ncores="$OPTARG";;	
	o)  outdir="$OPTARG";;
    esac
done
id="$outdirH"$id"_TCR"
sample=$sample"_TCR"

reference=/data/databases/CellRanger/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0

cellranger vdj --id=$id \
--sample=$sample \
--reference=$reference \
--fastqs=$fastqs \
--localcores=$ncores
