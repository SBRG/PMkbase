from Bio import Phylo
import json

def newick_to_json(newick_file):
    tree = Phylo.read(newick_file, "newick")
    return tree_to_dict(tree.clade)

def tree_to_dict(clade):
    node = {"name": clade.name if clade.name else ""}
    if clade.clades:
        node["children"] = [tree_to_dict(c) for c in clade.clades]
    return node

# Save this JSON to a file or use it directly in your Flask app
tree_json = newick_to_json("tree.newick")
with open("tree.json", "w") as f:
    json.dump(tree_json, f)
