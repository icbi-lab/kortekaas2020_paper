#$ -S /bin/sh

# Options
chemistry=auto
while getopts i:f:s:n:o:c: opt;
do
    case $opt in
        i)  id="$OPTARG";;
        f)  fastqs="$OPTARG";;
        s)  sample="$OPTARG";;
	n)  ncores="$OPTARG";;	
	o)  outdir="$OPTARG";;
	c)  chemistry="$OPTARG";;
    esac
done
id="$outdirH"$id"_GEX"
sample=$sample"_GEX"

transcriptome=/data/databases/CellRanger/refdata-cellranger-GRCh38-3.0.0

cellranger count --id=$id \
--transcriptome=$transcriptome \
--fastqs=$fastqs \
--sample=$sample \
--localcores=$ncores \
--chemistry=$chemistry

