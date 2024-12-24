import argparse
import collections
import json
import pandas as pd
import re

def get_pango_mutations(pango_consensus_seqs_json, starting_clades):
    """
    Retrieve mutations for each Pango lineage relative to parent.
    """
    with open(pango_consensus_seqs_json) as f:
        pango_clades = json.load(f)

    def _build_records(c, recs):
        """
        Build records of Pango clade information.
        Recurses over children.
        """
        if c in recs["clade"]:
            return
        recs["clade"].append(c)
        recs["date"].append(pango_clades[c]["designationDate"])
        recs["parent"].append(pango_clades[c]["parent"])
        recs["all_new_muts_from_ref"].append(
            [
                mut
                for field in ["aaSubstitutionsNew", "aaDeletionsNew"]
                for mut in pango_clades[c][field]
                if mut
            ]
        )
        recs["all_new_muts_reverted_from_ref"].append(
            [
                mut
                for field in ["aaSubstitutionsReverted", "aaDeletionsReverted"]
                for mut in pango_clades[c][field]
                if mut
            ]
        )
        recs["spike_muts_from_ref"].append(
            [
                mut.split(":")[1]
                for field in ["aaSubstitutions", "aaDeletions"]
                for mut in pango_clades[c][field]
                if mut and mut.startswith("S:")
            ]
        )
        for c_child in pango_clades[c]["children"]:
            _build_records(c_child, recs)
        return None
            
    records = collections.defaultdict(list)
    for starting_clade in starting_clades:
        _build_records(starting_clade, records)

    pango_df = pd.DataFrame(records).query("clade not in @exclude_clades")
    return pango_df

def consolidate_reversions(r):
    """If there are reversions, combine with new mutations to get actual changes."""
    reverted = r["all_new_muts_reverted_from_ref"]
    new = r["all_new_muts_from_ref"]
    if not reverted:
        return new
    # get as dicts with key: val of (gene, site): (wt, mutant)
    reverted_dict = {
        (m.split(":")[0], m.split(":")[1][1: -1]): (m.split(":")[1][0], m.split(":")[1][-1])
        for m in reverted
    }
    new_dict = {
        (m.split(":")[0], m.split(":")[1][1: -1]): (m.split(":")[1][0], m.split(":")[1][-1])
        for m in new
    }
    muts = []
    for (gene, site), (rev_wt, rev_mutant) in reverted_dict.items():
        if (gene, site) in new_dict:
            new_wt, new_mutant = new_dict[(gene, site)]
            muts.append(f"{gene}:{rev_mutant}{site}{new_mutant}")
            del new_dict[(gene, site)]
        else:
            muts.append(f"{gene}:{rev_mutant}{site}{rev_wt}")
    for (gene, site), (new_wt, new_mutant) in new_dict.items():
        muts.append(f"{gene}:{new_wt}{site}{new_mutant}")
    return muts
        
def get_pango_pair_differences(pango_df, pango_relationships_path, exclude_muts, pair_only_new_muts_are_spike):
    # Alter parent name here with pango-variant relationships file for computing comparisons
    pango_relationships = pd.read_csv(pango_relationships_path, sep="\t")
    pango_relationships = pango_relationships.set_index("variant")
    pango_df["parent"] = pango_df["variant"].map(lambda x: pango_relationships[x].parent)

    pango_pair_df = (
        pango_df
        .query("parent != ''")
        .assign(
            all_new_muts=lambda x: x.apply(consolidate_reversions, axis=1),
            spike_new_muts_list=lambda x: x["all_new_muts"].map(lambda ms: [m.split(":")[1] for m in ms if m[0] == "S"]),
            n_new_spike_muts=lambda x: x["spike_new_muts_list"].map(len),
            n_new_all_muts=lambda x: x["all_new_muts"].map(len),
        )
        .assign(
            spike_new_muts=lambda x: x["spike_new_muts_list"].map(lambda ms: "; ".join(ms)),
            nonspike_new_muts=lambda x: x["all_new_muts"].map(lambda ms: "; ".join(m for m in ms if not m.startswith("S:"))),
            only_new_muts_are_spike=lambda x: x["n_new_spike_muts"] == x["n_new_all_muts"],
        )
        [["clade", "parent", "date", "n_new_spike_muts", "spike_new_muts_list", "spike_new_muts", "nonspike_new_muts", "only_new_muts_are_spike"]]
        .reset_index(drop=True)
    )

    for mut in exclude_muts:
        pango_pair_df = pango_pair_df[~pango_pair_df["spike_new_muts_list"].map(lambda ms: mut in ms)]

    print("Number of clade pairs:")
    print(
        pango_pair_df
        .groupby(["only_new_muts_are_spike", "n_new_spike_muts"])
        .aggregate(n_clade_pairs=pd.NamedAgg("clade", "count"))
    )

    if pair_only_new_muts_are_spike:
        print("Dropping pairs where new mutations are not spike.")
        pango_pair_df = pango_pair_df.query("only_new_muts_are_spike")
    pango_pair_df = pango_pair_df.drop(columns="only_new_muts_are_spike")
    return pango_pair_df

def process_dms_summary(input_path, split_by_rbd, missing_muts, rename_cols, phenotype_cols):
    """
    Read in dms summary to allow computation of dms phenotypes on sequences in future.
    """
    dms_summary = pd.read_csv(input_path).rename(columns=rename_cols)

    # DMS phenotypes of interest
    phenotypes_basic = phenotype_cols
    assert set(phenotypes_basic).issubset(dms_summary.columns), f"{phenotypes_basic=}\n{dms_summary.columns=}"

    if split_by_rbd is False:
        phenotypes = phenotypes_basic
    elif split_by_rbd is True:
        phenotypes = []
        for s in ["RBD", "non-RBD"]:
            for p in phenotypes_basic:
                p_by_rbd = f"{p} ({s})"
                phenotypes.append(p_by_rbd)
    else:
        raise ValueError(f"invalid {split_by_rbd=}")

    # dict that maps site to wildtype in DMS
    dms_wt = dms_summary.set_index("site")["wildtype"].to_dict()

    # dict that maps site to region in DMS
    if "region" in dms_summary.columns:
        site_to_region = dms_summary.set_index("site")["region"].to_dict()
    else:
        site_to_region = {}
        if split_by_rbd:
            raise ValueError(f"{split_by_rbd=} but 'region' not in {dms_summary.columns=}")

    # dicts that maps (site, wt, mutant) to DMS phenotypes for mutations
    dms_data_all_dict = (
        dms_summary
        .set_index(["site", "wildtype", "mutant"])
        [phenotypes_basic]
        .to_dict(orient="index")
    )

    def mut_dms(m, dms_data_dict):
        """Get DMS phenotypes for a mutation."""
        if missing_muts == "drop":
            null_d = {k: pd.NA for k in phenotypes_basic}
        elif missing_muts == "zero":
            null_d = {k: 0.0 for k in phenotypes_basic}
        else:
            raise ValueError(f"invalid {missing_muts=}")
        site = int(m[1: -1])
        if site not in dms_wt:
            d = null_d
        else:
            parent = m[0]
            mut = m[-1]
            wt = dms_wt[site]
            if parent == wt:
                try:
                    d = dms_data_dict[(site, parent, mut)]
                except KeyError:
                    d = null_d
            elif mut == wt:
                try:
                    d = {k: -v for (k, v) in dms_data_dict[(site, mut, parent)].items()}
                except KeyError:
                    d = null_d
            else:
                try:
                    parent_d = dms_data_dict[(site, wt, parent)]
                    mut_d = dms_data_dict[(site, wt, mut)]
                    d = {p: mut_d[p] - parent_d[p] for p in phenotypes_basic}
                except KeyError:
                    d = null_d
        assert list(d) == phenotypes_basic
        return d

    def muts_dms(ms, dms_data_dict):
        """Get DMS phenotypes for a list of mutations."""
        muts_d = {k: 0 for k in phenotypes}
        for m in ms:
            if split_by_rbd:
                try:
                    is_rbd = site_to_region[int(m[1: -1])] == "RBD"
                except KeyError:
                    return {k: pd.NA for k in phenotypes}
            m_d = mut_dms(m, dms_data_dict)
            for k, val in m_d.items():
                if split_by_rbd:
                    if is_rbd:
                        muts_d[f"{k} (RBD)"] += val
                    else:
                        muts_d[f"{k} (non-RBD)"] += val
                else:
                    muts_d[k] += val
        return muts_d
    return phenotypes, dms_data_all_dict, muts_dms


def compute_dms_phenotype_pairs(pango_pair_df, dms_data_all_dict, phenotypes, muts_dms, pango_pair_dms_path):
    """
    Compute DMS phenotype for pango_pairs.
    """

    pango_pair_dms_df = (
        pango_pair_df
        # to add multiple columns: https://stackoverflow.com/a/46814360
        .apply(
            lambda cols: pd.concat(
                [cols, pd.Series(muts_dms(cols["spike_new_muts_list"], dms_data_all_dict))]
            ),
            axis=1,
        )
        .drop(columns="spike_new_muts_list")
        # remove any clade pairs for which we don't have DMS data for all phenotypes
        .query(" and ".join(f"`{p}`.notnull()" for p in phenotypes))
        .reset_index(drop=True)
    )

    pango_pair_dms_df.to_csv(pango_pair_dms_path, float_format="%.5g", index=False, sep="\t")
    return pango_pair_dms_df

def compute_dms_phenotype(
        pango_df: pd.DataFrame, 
        dms_clade: str,
        dms_data_all_dict,
        phenotypes,
        muts_dms,
        exclude_muts,
        clade_dms_path: str | None =None,
        ) -> pd.DataFrame:
    """
    Compute dms phenotype for clades in `pango_df` relative to `dms_clade`.
    """

    def mutations_from(muts, from_muts):
        """Get mutations from another sequence."""
        new_muts = set(muts).symmetric_difference(from_muts)
        assert all(re.fullmatch("[A-Z\-]\d+[A-Z\-]", m) for m in new_muts)
        new_muts_d = collections.defaultdict(list)
        for m in new_muts:
            new_muts_d[int(m[1: -1])].append(m)
        new_muts_list = []
        for _, ms in sorted(new_muts_d.items()):
            if len(ms) == 1:
                m = ms[0]
                if m in muts:
                    new_muts_list.append(m)
                else:
                    assert m in from_muts
                    new_muts_list.append(m[-1] + m[1: -1] + m[0])
            else:
                m, from_m = ms
                if m not in muts:
                    from_m, m = m, from_m
                assert m in muts and from_m in from_muts
                new_muts_list.append(from_m[-1] + m[1: ])
        return new_muts_list 

    dms_clade_spike_muts_from_ref = pango_df.set_index("clade").at[
        dms_clade, "spike_muts_from_ref"
    ]
    clade_df = (
        pango_df
        .assign(
            spike_muts_from_dms_clade_list=lambda x: x["spike_muts_from_ref"].apply(
                mutations_from, args=(dms_clade_spike_muts_from_ref,),
            ),
            spike_muts_from_dms_clade=lambda x: x["spike_muts_from_dms_clade_list"].map(lambda s: "; ".join(s)),
        )
        .query("clade not in @exclude_clades")
        [["clade", "date", "spike_muts_from_dms_clade_list", "spike_muts_from_dms_clade"]]
    )

    for mut in exclude_muts:
        clade_df = clade_df[~clade_df["spike_muts_from_dms_clade_list"].map(lambda ms: mut in ms)]

    clade_dms_df = (
        clade_df
        # to add multiple columns: https://stackoverflow.com/a/46814360
        .apply(
            lambda cols: pd.concat(
                [cols, pd.Series(muts_dms(cols["spike_muts_from_dms_clade_list"], dms_data_all_dict))]
            ),
            axis=1,
        )
        .drop(columns="spike_muts_from_dms_clade_list")
        # remove any clade pairs for which we don't have DMS data for all phenotypes
        .query(" and ".join(f"`{p}`.notnull()" for p in phenotypes))
        .reset_index(drop=True)
    )

    print(f"Saving to {clade_dms_path}")
    clade_dms_df.to_csv(clade_dms_path, float_format="%.5g", index=False, sep="\t")

    print(f"Retained {len(clade_dms_df)} clades with growth and DMS data")
    return clade_dms_df

def main():
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument("--pango-consensus-seqs-json", type=str, required=True, help="Path to JSON containing consensus seqeuences for clades.")
    parser.add_argument("--input-path", type=str, required=True, help="Path to file containing summaries.")
    parser.add_argument("--clade-phenotype-path", type=str, required=True, help="Path to file containing clade phenotypes.")
    parser.add_argument("--pango-relationships-path", type=str, required=True, help="Path to file containing pango-parent relationships.")
    parser.add_argument("--clade-pair-phenotype-path", type=str, required=True, help="Path to file containing clade pair phenotypes.")
    parser.add_argument("--config", type=str, required=True, help="Configuration file containing details for processing phenotypes")
    args = parser.parse_args()

    config = json.loads(args.config)

    # Add args for pango_consensus_json, starting_clades, .etc
    pango_df = get_pango_mutations(pango_consensus_seqs_json=args.pango_consensus_seqs_json, starting_clades=None)

    # Produce mapping to produce phenotypes from summary file
    phenotypes, dms_data_all_dict, muts_dms = process_dms_summary(
        input_path=args.input_path,
        split_by_rbd=config["split_by_rbd"],
        missing_muts=config["missing_muts"],
        rename_cols=config["rename_cols"],
        phenotype_cols=config["phenotype_cols"],
    )

    # Map from reference sequences to DMS phenotype
    compute_dms_phenotype(
        pango_df, 
        dms_clade=config["dms-clade"], 
        dms_data_all_dict=dms_data_all_dict,
        phenotypes=phenotypes,
        muts_dms=muts_dms,
        exclude_muts=config["exclude_muts"],
        clade_dms_path=args.clade_phenotype_path,
    )

    # Create pair dataframe using pango_df and pango_relationships
    pango_pair_df = get_pango_pair_differences(
        pango_df,
        args.pango_relationships_path,
        exclude_muts=config["exclude_muts"],
        pair_only_new_muts_are_spike=False,
    )

    compute_dms_phenotype_pairs(
        pango_pair_df, 
        dms_data_all_dict=None, 
        phenotypes=None, 
        muts_dms=None, 
        pango_pair_dms_path=None)
    return None

if __name__ == "__main__":
    main()