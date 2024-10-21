#!/usr/bin/env python3
# coding: utf-8

import os
import json
import platform
import urllib.request
import tarfile
import zipfile
import sys

def get_os_arch():
    """Determine the operating system and architecture."""
    os_type = platform.system().lower()
    arch = platform.machine().lower()

    if os_type == 'linux':
        return 'linux-x64' if arch in ('x86_64', 'amd64') else 'linux-x86'
    elif os_type == 'darwin':
        return 'osx'
    elif os_type.startswith('win'):
        return 'win64' if arch == 'amd64' else 'win32'
    else:
        print(f"Unsupported OS: {os_type}")
        sys.exit(1)

def ensure_directory_exists(path):
    """Create directory if it does not exist."""
    if not os.path.exists(path):
        os.makedirs(path)

def download_and_extract(url, dest_path, is_zip=False):
    """Download and extract a package from a URL."""
    filename = url.split('/')[-1]
    file_path = os.path.join(dest_path, filename)

    # Download the file
    print(f"Downloading {filename} from {url}...")
    urllib.request.urlretrieve(url, file_path)

    # Extract the file
    print(f"Extracting {filename}...")
    if is_zip:
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(dest_path)
    else:
        with tarfile.open(file_path, 'r:gz') as tar_ref:
            tar_ref.extractall(dest_path)

    # Clean up
    os.remove(file_path)
    print(f"{filename} installed successfully.\n")

def main():
    # Define paths
    current_folder = os.path.abspath(os.path.dirname(__file__))
    base_folder = os.path.abspath(os.path.join(current_folder, os.pardir))
    default_libs_path = os.path.join(base_folder, "addons", "dFBA", "ext")

    # Load package information
    json_packages = os.path.join(base_folder, 'beta', 'fba_packages.json')
    with open(json_packages) as fh:
        packages_dict = json.load(fh)

    # Determine OS and architecture
    arch = get_os_arch()

    # Install packages
    for pkg in ("coin-or", "libsbml"):
        pkg_path = os.path.join(default_libs_path, pkg)

        # Skip download if package folder already exists and is populated
        # to add: create check if coin-or and libsbml folder exist, if they do not not, create them and then download the files
        if os.path.exists(pkg_path) and os.listdir(pkg_path):
            print(f"Package {pkg} is already installed. Skipping...")
            continue

        pkg_dict = packages_dict.get(pkg, {}).get(arch)
        if not pkg_dict:
            print(f"No package information found for {pkg} on {arch}. Skipping...")
            continue

        url = pkg_dict['url']
        is_zip = url.endswith("zip")
        download_and_extract(url, pkg_path, is_zip)

if __name__ == "__main__":
    main()
