import os
import requests
import argparse

# Base URL for the files
BASE_URL = "https://data.nextstrain.org/files/workflows/forecasts-ncov"

def download_file(url, dest):
    """Download a file from a URL and save it locally."""
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    try:
        print(f"Downloading {url}...")
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an error for bad status
        with open(dest, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"Saved to {dest}")
    except requests.RequestException as e:
        print(f"Failed to download {url}: {e}")

def provision_file(data_provenance, variant_classification, geo_resolution, output_path):
    """Construct the URL and download the specified file."""
    file_name = f"{geo_resolution}.tsv.gz"
    url = f"{BASE_URL}/{data_provenance}/{variant_classification}/{file_name}"
    print(f"Constructed URL: {url}")
    download_file(url, output_path)

def main():
    parser = argparse.ArgumentParser(description="Provision sequence counts data files for Snakemake workflows.")
    parser.add_argument("--data-provenance", required=True, help="Data provenance (e.g., gisaid or open).")
    parser.add_argument("--variant-classification", required=True, help="Variant-classification (e.g., nextstrain_clades or pango_lineages).")
    parser.add_argument("--geo-resolution", required=True, help="Geographic resolution (e.g., global or usa).")
    parser.add_argument("--output-path", required=True, help="Path to save the downloaded file.")

    args = parser.parse_args()

    # Provision the requested file
    provision_file(
        data_provenance=args.data_provenance,
        variant_classification=args.variant_classification,
        geo_resolution=args.geo_resolution,
        output_path=args.output_path
    )

if __name__ == "__main__":
    main()
