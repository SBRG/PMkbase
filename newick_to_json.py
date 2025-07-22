from Bio import Phylo
import json

def newick_to_json(newick_file):
    """
    Converts a Newick format tree file to a JSON-compatible dictionary.

    Args:
        newick_file (str): Path to the Newick format file.

    Returns:
        dict: Nested dictionary representing the tree structure.
    """
    tree = Phylo.read(newick_file, "newick")
    return tree_to_dict(tree.clade)

def tree_to_dict(clade):
    """
    Recursively converts a Biopython Clade object to a dictionary.

    Args:
        clade (Bio.Phylo.Newick.Clade): Clade object from Biopython.

    Returns:
        dict: Dictionary with 'name' and optional 'children' keys.
    """
    node = {"name": clade.name if clade.name else ""}
    if clade.clades:
        node["children"] = [tree_to_dict(c) for c in clade.clades]
    return node

# Save this JSON to a file or use it directly in your Flask app
tree_json = newick_to_json("tree.newick")
with open("tree.json", "w") as f:
    json.dump(tree_json, f)
