import os, sys
import json
import argparse

def collect_args():
    parser = argparse.ArgumentParser(description = "Extract relevant results from evofr JSON")
    parser.add_argument('--results', required=True, type=str, help="input JSON file")
    return parser.parse_args()

if __name__=="__main__":
    params = collect_args()

    json_fh = open(params.results, "r")
    json_dict = json.load(json_fh)

    locations = json_dict["metadata"]["location"]
    variants = json_dict["metadata"]["variants"]
    json_ga = [d for d in json_dict["data"] if d["site"]=="ga"]

    print("location" + "\t" + "variant" + "\t" + "median_ga" + "\t" + "ga_upper_80" + "\t" + "ga_lower_80")

    for location in locations:
        for variant in variants:

            # default to 1 if absent, this is the pivot variant
            median = next( (d["value"] for d in json_ga if d["location"]==location and d["variant"]==variant and d["ps"]=="median"), 1)
            upper = next( (d["value"] for d in json_ga if d["location"]==location and d["variant"]==variant and d["ps"]=="HDI_80_upper"), 1)
            lower = next( (d["value"] for d in json_ga if d["location"]==location and d["variant"]==variant and d["ps"]=="HDI_80_lower"), 1)
            print(location + "\t" + variant + "\t" + str(median) + "\t" + str(upper) + "\t" + str(lower))
