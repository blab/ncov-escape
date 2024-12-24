import argparse
import os

import requests
import yaml

# Default files to download if no `output` key is found
DEFAULT_FILES_TO_DOWNLOAD = {
    "mutation_phenotypes.csv": "results/mutation_phenotypes.csv",
    "mutation_phenotypes_randomized.csv": "results/mutation_phenotypes_randomized.csv",
    "lineage_phenotypes.csv": "results/clade_phenotypes.csv",
    "lineage_phenotypes_randomized.csv": "results/clade_phenotypes_randomized.csv",
}

def download_file(url, local_path):
    """Download a file from a URL and save it locally."""
    os.makedirs(os.path.dirname(local_path), exist_ok=True)
    try:
        print(f"Downloading {url}...")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(local_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"Saved to {local_path}")
    except requests.RequestException as e:
        print(f"Failed to download {url}: {e}")

def parse_config_and_get_files(config_path, base_url):
    """Retrieve file paths from the specified config or use default files."""
    # Download the config file
    config_url = f"{base_url}/{config_path}"
    response = requests.get(config_url)
    response.raise_for_status()
    config = yaml.safe_load(response.text)

    # Extract output files if available, else return default files
    if "output" in config:
        return config["output"]
    else:
        print("No `output` key in config. Falling back to default files.")
        return DEFAULT_FILES_TO_DOWNLOAD

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Download files specified in a config file from a GitHub repository.")
    parser.add_argument("--config-path", required=True, help="Path to the configuration file in the repository.")
    parser.add_argument("--base-url", required=True, help="Base URL of the GitHub repository.")
    parser.add_argument("--output-path", required=True, help="Path to save files.")
    args = parser.parse_args()

    # Get the files to download
    files_to_download = parse_config_and_get_files(args.config_path, args.base_url)

    # Download each file
    for key, relative_path in files_to_download.items():
        file_url = f"{args.base_url}/{relative_path}"
        local_path = os.path.join(args.output_path, key)
        download_file(file_url, local_path)

if __name__ == "__main__":
    main()
