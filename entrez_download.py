from Bio import Entrez
from bs4 import BeautifulSoup as bs
import sys
import argparse
import time

def argparser():
    parser = argparse.ArgumentParser(description="Download influenzavirus sequences through Entrez.")
    parser.add_argument("-e","--email",help="Set your email.",default="jmcbroom@ucsc.edu")
    parser.add_argument("-g","--groupsize",help='Download this many sequences as a time. Default 10000',type=int,default=10000)
    parser.add_argument("-t","--total",default=0,type=int,help='Set the total number of sequences to download. Default is 0, which gets all available.')
    parser.add_argument("-q","--queryterm",help="Set the term for the entrez search. Default 'Alphainfluenzavirus'")
    parser.add_argument('-f','--fields',default=['segment','country','collection_date','host','isolation_source','serotype'],nargs='*',help="Add these fields as columns to the metadata table.")
    parser.add_argument("-s","--store-name",action='store_true',help="Use to additionally store metadata in fasta headers.")
    return parser.parse_args()
args = argparser()
Entrez.email= args.email
groupsize = args.groupsize
metafields = ['strain'] + args.fields
segment_files = [open('segment_'+str(i)+".fasta",'w+') for i in range(1,9)]
other_file = open('other.fasta','w+')
meta_file = open("metadata.tsv","w+")
print("\t".join(metafields),file=meta_file)
total = int(Entrez.read(Entrez.esearch(db='nuccore', retstart=0, retmax=1, term='Alphainfluenzavirus', retmode='xml'))['Count'])
for i in range(0,total,groupsize):
    try:
        start = time.time()
        handle = Entrez.esearch(db='nuccore', retstart=i, retmax=groupsize, term='Alphainfluenzavirus', retmode='xml') 
        uid_record = Entrez.read(handle)
        handle.close()
        print(f"Retrieving {i} to {i+len(uid_record['IdList'])} of {total} total ({uid_record['Count']} available)",file=sys.stderr)
        dhandle=Entrez.efetch(db='nuccore', id=uid_record['IdList'], retmode='xml') 
        record = bs(dhandle.read(),'xml')
        dhandle.close()
        for data in record.find_all("GBSeq"):
            sname = data.find("GBSeq_primary-accession")
            if sname == None:
                print("Skipping sample with no accession...",file=sys.stderr)
                continue
            sname = next(sname.strings)
            #build the description by qualifier.
            #we care about the host, the segment, the serotype, the isolation source, the country, and the date.
            keys = {k:"N/A" for k in metafields}
            keys['strain'] = sname
            for quals in data.find_all("GBQualifier"):
                try:
                    name, value = list(quals.stripped_strings)
                except ValueError:
                    continue
                if name in keys:
                    keys[name] = value
            if keys.get("segment",'N/A') != 'N/A':
                try:
                    target = segment_files[int(keys['segment'])-1]
                except ValueError:
                    print(f"Skipping- malformed (non-integer) segment key for accession {sname}.",file=sys.stderr)
                    continue
            else:
                target = other_file
            seq = data.find("GBSeq_sequence")
            if seq == None:
                print(f"Skipping sample with accession {sname} that has incorrectly tagged or no sequence",file=sys.stderr)
                continue
            seq = next(seq.strings).upper()
            print("\t".join([keys[k] for k in metafields]),file=meta_file)
            if args.store_name:
                description = "; ".join([k + ": " + keys[k] for k in metafields])
                print(">"+sname+" "+description,file=target)
            else:
                print(">"+sname,file=target)
            print(seq,file=target)
        print(f"Completed chunk download in {time.time() - start} seconds.",file=sys.stderr)
    except KeyboardInterrupt:
        print("Halting...")
        break

for s in segment_files:
    s.close()
other_file.close()
meta_file.close()
