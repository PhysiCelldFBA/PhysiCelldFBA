#!/usr/bin/env python3
# coding: utf-8

import sys, os
import json
import argparse

import requests
import hashlib
import tarfile
import zipfile
import platform


def param_parser():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('--pkg', dest="pkg", required=True, help='Available packages', choices=PACKAGES)
    parser.add_argument('--arch', dest="arch", required=True, choices=ARCHS, help='Current arch')
    parser.add_argument('--checksum', dest="checksum", default=False, help='Check file integrity after downloading')
    parser.add_argument('--path', dest='path', default=DEFAULT_LIBS_PATH, help='Default folder destination to install third-party libs if changed, Makefiles you be updated according')
    return parser


PACKAGES = ("coin-or", "libsbml")
ARCHS = ("linux-x64", "linux-x86", "win64", "win32", "osx")
OS_TYPE = platform.system()
SBML_DIR = "libsbml"

current_folder = os.path.abspath(os.path.dirname(__file__))
BASE_FOLDER = os.path.abspath(os.path.join(current_folder, os.pardir))
DEFAULT_LIBS_PATH = os.path.join(BASE_FOLDER, "addons", "dFBA", "ext")


def main():
    parser = param_parser()
    args = parser.parse_args()
   
    json_packages = os.path.join(BASE_FOLDER, 'beta', 'fba_packages.json')
    packages_dict = {}
    with open(json_packages) as fh:
        packages_dict = json.load(fh)
    
    base_path = args.path
    if not os.path.exists(base_path):
        print(f"Creating {base_path} folder", end=" ")
        os.makedirs(base_path)
        print("Ok!")

    pkg_path = os.path.join(base_path, args.pkg)
    if not os.path.exists(pkg_path):
        print(f"Creating {pkg_path} folder", end=" ")
        os.makedirs(pkg_path)
        print("Ok!")

        print("Moving to %s folder" % pkg_path)
        os.chdir(pkg_path)
            
        pkg_dict = packages_dict[args.pkg][args.arch]
        print(f"Fetching package: {args.pkg} ({args.arch})", end="")
        print(f"Downaling from: {pkg_dict['url']}")
        r = requests.get(pkg_dict['url'], allow_redirects=True)

        # if args.checksum:
        #     print("Package cheksum(sha256)", end=" ")
        #     hash_strn = hashlib.sha256(r.content).hexdigest()
        #     assert pkg_dict['sha256'] == hash_strn
        #     print("Ok!")
        
        fname = pkg_dict["version"]
        with open(fname, 'wb') as fh:
            fh.write(r.content)

        print("Extracting package in %s... " % args.pkg, end=" ")
        if fname.endswith("zip"):
            archiver = zipfile.ZipFile(fname, 'r')
            archiver.extractall()
        elif fname.endswith("tar.gz") or fname.endswith("gz"):
            archiver = tarfile.open(fname, "r:gz")
            archiver.extractall()
            archiver.close()
            if args.pkg == 'libsbml':
                old_lib_dir = os.path.commonprefix(archiver.getnames())
                os.rename(old_lib_dir, SBML_DIR)
                if os.path.exists(old_lib_dir):
                    os.rmdir(old_lib_dir)
        os.remove(fname)
            
        print("Ok!")
        print(f"Dependency {args.pkg} retrived and installed correctly :-)\n")
    
    else:
        print(f"Package {args.pkg} already installed. To reinstall them remove the ext folder and run again.")

if __name__ == "__main__":
    main()