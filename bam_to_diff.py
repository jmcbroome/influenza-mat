import sys
import pysam
bam = pysam.AlignmentFile(sys.argv[1], "rb")
lengths = bam.lengths
for c, rlen in zip(bam.references, bam.lengths):
    with open(c+".diff",'w+') as dfile:
        print(f"Adding to file {c+'.diff'}...")
        for read in bam.fetch(contig=c):
            print(">"+read.query_name,file=dfile)
            #fill in Ns for all unaligned parts of the sequence
            for i in range(0, read.reference_start):
                print("N",i,sep='\t',file=dfile)
            for qposition, refposition, refbase in read.get_aligned_pairs(with_seq=True):
                if qposition == None:
                    continue
                base = read.query_sequence[qposition]
                #ignore insertions in the read, since we can't use their variational data.
                if refposition == None:
                    continue
                #deletions are filled as N.
                elif base == None:
                    print("N", refposition, sep='\t',file=dfile)
                #mismatch are recorded as they are.
                elif base != refbase:
                    print(base, refposition, sep='\t',file=dfile)
            #fill in Ns again.
            for i in range(read.reference_end, rlen):
                print("N",i,sep='\t',file=dfile)
