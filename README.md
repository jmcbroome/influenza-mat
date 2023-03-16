# influenza-mat
Code and information for constructing influenzavirus Mutation Annotated Trees (MATs). 

Influenzavirus data is obtained from Entrez using `entrez_download.py`. This produces a series of segment-divided alphainfluenzavirus fasta. These are then aligned with [ViralMSA](https://github.com/niemasd/ViralMSA) wrapping [Minimap2](https://github.com/lh3/minimap2) and converted to VCF with [faToVcf](http://hgdownload.cse.ucsc.edu/admin/exe/). A MAT for each segment is then constructed with [UShER](https://github.com/yatisht/usher). 

A reference genome can be obtained from [NCBI Virus](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/labs/virus/vssi/#/virus?SeqType_s=Genome&VirusLineage_ss=Alphainfluenzavirus,%20taxid:197911). The data in this repository used [this reference](https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/assembly/GCF_001343785.1). The reference fasta is then split with `split_segments.py` to generate segment-level references for alignment.

The following simple bash script implements this pipeline.

```
#insert your email and choice of reference here.
export e=”my_email” ; 
export ref="reference_" ; 
python3 split_segments.py $ref ;
conda create -f env.yml ;
conda activate mat ;
python3 entrez_download.py -e $e ;
for f in segment_*.fasta ; 
  do echo $f ; 
  #align the sequences.
  python3 ViralMSA.py -s $f -r ref_${f} -e $e -o temp ; 
  mv temp/* . ; 
  rm -r temp ; 
  #split the resulting aligned fasta into chunks for processing
  split -l 10000 -d –additional-suffix=.fasta.aln $f.aln ${f%.fasta}_ ;
  #conver to vcf
  for g in ${f%.fasta}_*.aln ;
    do ./faToVcf $g ${g%fasta.aln}vcf ; 
  done ;
  #create the initial tree then consecutively place onto the tree
  usher-sampled -t seed.nwk -v ${f%.fasta}_00.vcf -o ${f%fasta}pb ;
  for v in ${f%.fasta}*vcf ; 
    do usher-sampled -i ${f%fasta}pb -v $v -o ${f%fasta}pb –optimization_radius 0;
  done ;
done ;
gzip *pb ;
```
