import bte
import argparse
import os
def argparser():
    #thanks to ChatGPT for this boilerplate.
    parser = argparse.ArgumentParser(description='Process some files.')
    parser.add_argument('-t','--tree', help='Path to the base Newick tree file. Tips should have names matching the target protobuf files.')
    parser.add_argument('-f','--pb_folder', help='Path to the folder containing the .pb files to graft.')
    parser.add_argument('-g','--gzipped', action='store_true', help='Indicate that the target tree files are gzipped and end with .gz.')
    parser.add_argument('-o','--output', help="Name of the output MAT protobuf. Default is grafted.pb", default='grafted.pb')
    parser.add_argument('-s','--scale',type=float,help="Scale the branch lengths of the base tree by this value. Default 1",default=1)
    parser.add_argument('-d','--dump',help="Set to a file to write a two-column tab-delimited set of sample-subtree identifiers, for adding to metadata.",default=None)
    args = parser.parse_args()
    if not os.path.isfile(args.tree):
        parser.error(f"File '{args.tree}' does not exist.")
    # Check if the pb_folder exists and is a directory
    if not os.path.isdir(args.pb_folder):
        parser.error(f"'{args.pb_folder}' is not a valid directory.")
    # Get a list of all the .pb files in the pb_folder
    if args.gzipped:
        pb_files = [os.path.join(args.pb_folder, f) for f in os.listdir(args.pb_folder) if f.endswith('.pb.gz')]
    else:
        pb_files = [os.path.join(args.pb_folder, f) for f in os.listdir(args.pb_folder) if f.endswith('.pb')]
    # Check if any .pb files were found
    if not pb_files:
        parser.error(f"No .pb files found in '{args.pb_folder}'.")
    return args, pb_files

def main():
    args, pb_files = argparser()
    if args.dump != None:
        dumpf = open(args.dump,'w+')
    basetree = bte.MATree(nwk_file=args.tree)
    #scale the branch lengths.
    for n in basetree.depth_first_expansion():
        n.set_branch_length(n.branch_length * args.scale)
    #proceed to add the data. 
    names_seen = set()
    for l in basetree.get_leaves():
        if args.gzipped:
            ending = 'pb.gz'
        else:
            ending = 'pb'
        target = args.pb_folder + "/" + l.id + '.' + ending
        if target not in pb_files:
            print(f"WARNING: No matching MAT found for tip {l.id}!")
        else:
            subtree = bte.MATree(target)
            print(f"Adding tree {l.id}...")
            for n in subtree.depth_first_expansion():
                if n.parent == None:
                    basetree.create_node((l.id+"_"+n.id), parent_id = l.id, mutations = n.mutations, annotations = [], branch_length = l.branch_length)
                else:
                    if n.is_leaf():
                        name = n.id
                        if name in names_seen:
                            print(f"WARNING: Sample {name} present multiple times on the tree; tagging with *...")
                            name = name+"*"
                        basetree.create_node(name, parent_id = l.id+"_"+n.parent.id, mutations = n.mutations, annotations = [], branch_length = len(n.mutations))
                        names_seen.add(name)
                        if args.dump != None:
                            print(name, l.id, sep='\t', file=dumpf)
                    else:
                        basetree.create_node(l.id+"_"+n.id, parent_id = l.id+"_"+n.parent.id, mutations = n.mutations, annotations = [], branch_length = len(n.mutations))
    basetree.save_pb(args.output)
    if args.dump != None:
        dumpf.close()

if __name__ == "__main__":
    main()
