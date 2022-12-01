import os, sys
import json
import Bio.Phylo
import argparse

def json_to_tree(json_dict, root=True):
    """Returns a Bio.Phylo tree corresponding to the given JSON dictionary exported
    by `tree_to_json`.
    Assigns links back to parent nodes for the root of the tree.
    """
    # Check for v2 JSON which has combined metadata and tree data.
    if root and "meta" in json_dict and "tree" in json_dict:
        json_dict = json_dict["tree"]

    node = Bio.Phylo.Newick.Clade()

    # v1 and v2 JSONs use different keys for strain names.
    if "name" in json_dict:
        node.name = json_dict["name"]
    else:
        node.name = json_dict["strain"]

    if "children" in json_dict:
        # Recursively add children to the current node.
        node.clades = [json_to_tree(child, root=False) for child in json_dict["children"]]

    # Assign all non-children attributes.
    for attr, value in json_dict.items():
        if attr != "children":
            setattr(node, attr, value)

    # Only v1 JSONs support a single `attr` attribute.
    if hasattr(node, "attr"):
        node.numdate = node.attr.get("num_date")
        node.branch_length = node.attr.get("div")

        if "translations" in node.attr:
            node.translations = node.attr["translations"]
    elif hasattr(node, "node_attrs"):
        node.branch_length = node.node_attrs.get("div")

    if root:
        node = annotate_parents_for_tree(node)

    return node

def annotate_parents_for_tree(tree):
    """Annotate each node in the given tree with its parent.
    """
    tree.root.parent = None
    for node in tree.find_clades(order="level"):
        for child in node.clades:
            child.parent = node

    # Return the tree.
    return tree

def collect_args():
    parser = argparse.ArgumentParser(description = "Extract relevant tip attributes")
    parser.add_argument('--json', required=True, type=str, help="input JSON file")
    return parser.parse_args()

if __name__=="__main__":
    params = collect_args()

    json_fh = open(params.json, "r")
    json_dict = json.load(json_fh)
    tree = json_to_tree(json_dict)

    print(
        "seqName" + "\t" + "clade" + "\t" + "Nextclade_pango" + "\t" \
        + "partiallyAliased" + "\t" + "immune_escape" + "\t" + "ace2_binding" + "\t" \
        + "rbd_count" + "\t" + "non_rbd_spike_count" + "\t" + "non_spike_count"
    )

    data = []
    for n in tree.find_clades(order="postorder"):
        node_elements = {}
        name = n.name
        if hasattr(n, 'node_attrs'):
            if 'clade_membership' in n.node_attrs:
                if 'value' in n.node_attrs["clade_membership"]:
                    clade_membership = n.node_attrs["clade_membership"]["value"]
            else:
                clade_membership = "?"
            if 'Nextclade_pango' in n.node_attrs:
                if 'value' in n.node_attrs["Nextclade_pango"]:
                    Nextclade_pango = n.node_attrs["Nextclade_pango"]["value"]
            else:
                Nextclade_pango = "?"
            if 'partiallyAliased' in n.node_attrs:
                if 'value' in n.node_attrs["partiallyAliased"]:
                    partiallyAliased = n.node_attrs["partiallyAliased"]["value"]
            else:
                partiallyAliased = "?"
            if 'immune_escape' in n.node_attrs:
                if 'value' in n.node_attrs["immune_escape"]:
                    immune_escape = n.node_attrs["immune_escape"]["value"]
                    immune_escape = str(round(float(immune_escape), 8))
            else:
                immune_escape = "?"
            if 'ace2_binding' in n.node_attrs:
                if 'value' in n.node_attrs["ace2_binding"]:
                    ace2_binding = n.node_attrs["ace2_binding"]["value"]
                    ace2_binding = str(round(float(ace2_binding), 8))
            else:
                ace2_binding = "?"

        spike_mutations = set()
        non_spike_mutations = set()
        p = n
        while p.parent != None:
            if 'mutations' in p.branch_attrs:
                if 'S' in p.branch_attrs["mutations"]:
                    mut = p.branch_attrs["mutations"]["S"]
                    sites = [x[1:-1] for x in mut]
                    spike_mutations.update(sites)
                genes = list(p.branch_attrs["mutations"].keys())
                if 'nuc' in genes:
                    genes.remove('nuc')
                if 'S' in genes:
                    genes.remove('S')
                for gene in genes:
                    if gene in p.branch_attrs["mutations"]:
                        mut = p.branch_attrs["mutations"][gene]
                        sites = [gene+x[1:-1] for x in mut]
                        non_spike_mutations.update(sites)
            p = p.parent

        # print("spike_mutations", spike_mutations)
        # print("non_spike_mutations", non_spike_mutations)

        # RBD is positions 318 to 510 in spike
        rbd_mutations = [x for x in spike_mutations if int(x) >= 318 and int(x) <= 510]
        # print("rbd_mutations", rbd_mutations)

        # S1 is positions 14 to 685 in spike
        non_rbd_spike_mutations = [x for x in spike_mutations if int(x) < 318 or int(x) > 510]
        # print("non_rbd_spike_mutations", non_rbd_spike_mutations)

        rbd_count = str(len(rbd_mutations))
        non_rbd_spike_count = str(len(non_rbd_spike_mutations))
        non_spike_count = str(len(non_spike_mutations))

        if name == "BA.2":
            immune_escape = "0"
            ace2_binding = "0"
        if name[:4] != "NODE" and name[:8] != "internal" and name[-8:] != "internal" and name != "BA.3" and name != "rec_parent":
            print(
                name + "\t" + clade_membership + "\t" + Nextclade_pango + "\t" \
                + partiallyAliased + "\t" + immune_escape + "\t" + ace2_binding + "\t" \
                + rbd_count + "\t" + non_rbd_spike_count + "\t" + non_spike_count
            )
