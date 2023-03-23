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
conda env create -f env.yml ;
conda activate mat ;
python3 entrez_download.py -e $e ;
for f in segment_*.fasta ; 
  do echo $f ; 
  #align the sequences.
  python3 ViralMSA.py -s $f -r ref_${f} -e $e -o temp ; 
  mv temp/* . ; 
  rm -r temp ; 
  #split the resulting aligned fasta into chunks for processing
  split -l 10000 -d --additional-suffix=.fasta.aln $f.aln ${f%.fasta}_ ;
  #convert to vcf
  for g in ${f%.fasta}_*.aln ;
    do faToVcf $g ${g%fasta.aln}vcf ; 
  done ;
  #create the initial tree then consecutively place onto the tree
  usher-sampled -t seed.nwk -v ${f%.fasta}_00.vcf -o ${f%fasta}pb --optimization_radius 0;
  for v in ${f%.fasta}*vcf ; 
    do usher-sampled -i ${f%fasta}pb -v $v -o ${f%fasta}pb --optimization_radius 0;
  done ;
done ;
gzip *pb ;
```

In some cases divergence of downloaded data is too high to fit to a single reference for each segment, and instead multiple references per segment are required to capture more complete mapping. 

In this case, you can apply minimap2 to a combined reference file containing all your references for each segment.

```
minimap2 -a --MD --cs=long -c -L -x asm20 --sam-hit-only --secondary=no --score-N=0 all_flu_references.fasta all_flu.fasta | samtools view -b | samtools sort | samtools calmd -e -b - all_flu_references.fasta > all_flu.calmd.sorted.bam
```
Once you've created a bam with calmd tagging and masked bases, you can generate diff files for each segment reference.
```
python3 bam_to_diff.py all_flu.calmd.sorted.bam
```
You will need to split the lump reference file at this point into separate reference files for UShER. One simply option is to use [UCSC's faSplit](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/FOOTER).

```
./faSplit byname all_flu_references.fasta split_references/
```
You can then construct individual protobufs for each of the references.
```
for f in *.diff ; do echo $f ; usher-sampled -t seed.nwk --diff $f --ref split_references/${f%diff}fa -o ${f%diff}pb --optimization_radius 0; done
```
It may be useful for visualization or downstream analysis to group multiple trees of the same segment. 

You can obtain a newick tree from your segment fasta using whatever tools you would like. Just ensure that it is a rooted newick
and that the tip labels are the same as the protobuf file names. In my example, I split the references with split_segments.py and then 
uploaded each segment's fasta to the [ClustalW webtool](https://www.genome.jp/tools-bin/clustalw), creating an aligned fasta and then sending it to FastTree (Full) to produce a midpoint-rooted newick.

Once you have a backbone newick and a series of protobufs, you can graft protobufs together onto your backbone with graft_mats.py, which will go looking for the right protobufs to graft to each tip of your backbone tree.
```
python3 graft_mats.py -t segment_1.reference.nwk -f . -d reference_assignments.tsv -o segment_1.grafted.pb
```

This grafted protobuf is compatible with tools for MATs. You can use matUtils extract, summary, and similar tools. Placement with additional samples via UShER, however, may have unexpected results. To convert to taxonium, you can use usher_to_taxonium.

You may also want to merge your assignments tsv into the metadata file, so you can color the view by the reference used for that subtree.

```
import pandas as pd
pd1 = pd.read_csv('metadata.tsv',sep='\t')
pd2 = pd.read_csv('reference_assignments.tsv',sep='\t',names=['accession','reference'])
pd3 = pd1.merge(pd2, on='accession').drop_duplicates("accession")
pd3.to_csv("metadata.assignments.tsv",sep='\t',index=False)
```

```
usher_to_taxonium -i segment_1.complete.pb -o segment_1.jsonl.gz -m metadata.assignments.tsv --key_column accession -c strain,segment,country,collection_date,host,serotype,reference
```
