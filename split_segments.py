from Bio import SeqIO
import sys
files = {}
segnames = ['ref_segment_' + str(i) for i in range(1,9)]
for f in segnames:
    files[f] = open(f + ".fasta",'w+')
other = open("unsorted.fasta","w+")
for record in SeqIO.parse(sys.argv[1], "fasta"):
    try:
        search = record.description.find("segment")+8
        idx = record.description[search]
        if idx == " ":
            idx = record.description[search+1]
        idx = int(idx)-1
    except ValueError:
        print("Can't figure out segment of record")
        print(record.description)
        print(">"+record.description,file=other)
        print(record.seq,file=other)
        continue
    print(">"+record.description,file=files[segnames[idx]])
    print(record.seq,file=files[segnames[idx]])
for f in files.values():
    f.close()
other.close()