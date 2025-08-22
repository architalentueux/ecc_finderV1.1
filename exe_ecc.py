#!/usr/bin/env python3

"""
# August 2025
# If using this pipeline please cite : XXXXXXXXXX
#--------------------------------------------------------------------------+
#
#	ecc_finder is a tool
#       to detect eccDNA using Illumina and ONT sequencing.
#
#--------------------------------------------------------------------------+
#
#	AUTHOR: Update by Archimede PATIPE / ORIGINAL AUTHOR panpan ZHANG
#	CONTACT
#
#	LICENSE:
# 	GNU General Public License, Version 3
#	http://www.gnu.org/licenses/gpl.html
#
#	VERSION: 1.1.0
#
#--------------------------------------------------------------------------+
"""

import sys
import subprocess
#from eccFinder_lib.utilities import get_eccFinder_version
from utilitaires.utilities import log,run_oae,get_eccFinder_version

def main():
    VERSION = get_eccFinder_version()
    CITATION = """
Zhang et al. "Ecc_finder: a robust and accurate tool for detecting extrachromosomal circular DNA (eccDNA) from sequencing data"
    """

    description = """
ecc_Finder: Tool for detecting eccDNA loci using Illumina and ONT sequencing.
Version: %s

usage: ecc_finder.py <command> [options]

    Mapping mode:
      map-sr         Call candidate eccDNA loci from paired-end short reads
      map-ont        Call candidate eccDNA loci from Nanopore long reads

    Assembly mode:
      asm-sr         Assembly from paired-end short reads
      asm-ont        Assembly from Nanopore long reads
example 
   ecc_finder.py asm-ont /media/archimede/ExtremeSSD/Lucas/eccfinder2/PAU29426_S13_sup.fastq directory/S.meliloti2011complete.fai -r direcory/S.meliloti2011complete.fasta
    options:
      
      -v, --version""" % VERSION


    arg_len = len(sys.argv)

    if arg_len == 1:
        print(description)

    if arg_len > 1:
        cmd = sys.argv[1]

        if cmd == "-h" or cmd == "--help":
            print(description)

        elif cmd == "-v" or cmd == "--version":
            print(VERSION)

        elif cmd == "map-sr":
            if sys.argv[2:] == [] or sys.argv[2:] == ['-h'] or sys.argv[2:] == ['-help']:
                subcmd = "python map-sr.py"
                subprocess.call(subcmd, shell=True)
            else:
                subcmd = ["python", "map-sr.py"] + sys.argv[2:]
                subcmd_out = subprocess.run(subcmd, stdout=subprocess.PIPE).stdout.decode()

#        elif cmd == "map-ont":
#            if sys.argv[2:] == [] or sys.argv[2:] == ['-h'] or sys.argv[2:] == ['-help']:
#                subcmd = "python map-ont.py"
#                subprocess.call(subcmd, shell=True)
#            else:
#                subcmd = ["python", "map-ont.py"] + sys.argv[2:]
#                subcmd_out = subprocess.run(subcmd, stdout=subprocess.PIPE).stdout.decode()
        elif cmd == "map-ont":
            if sys.argv[2:] == [] or sys.argv[2:] == ['-h'] or sys.argv[2:] == ['-help']:
                #subcmd = "python test_eccfinder.py"
                subcmd = [sys.executable, "test_eccfinder.py"]
                subprocess.call(subcmd, shell=True)
            else:
                #subcmd = ["python", "test_eccfinder.py"] + sys.argv[2:]
                subcmd = [sys.executable, "test_eccfinder.py"] + sys.argv[2:]
                print("Running:", subcmd)  # Debug
                subcmd_out = subprocess.run(subcmd, stdout=subprocess.PIPE).stdout.decode()

        elif cmd == "asm-sr":
            if sys.argv[2:] == [] or sys.argv[2:] == ['-h'] or sys.argv[2:] == ['-help']:
                subcmd = "python asm-sr.py"
                subprocess.call(subcmd, shell=True)
            else:
                subcmd = ["python", "asm-sr.py"] + sys.argv[2:]
                subcmd_out = subprocess.run(subcmd, stdout=subprocess.PIPE).stdout.decode()

        elif cmd == "asm-ont":
            if sys.argv[2:] == [] or sys.argv[2:] == ['-h'] or sys.argv[2:] == ['-help']:
                subcmd = "python " + "asm-ont.py"
                subprocess.call(subcmd, shell=True)
            else:
                subcmd = ["python ", "asm-ont.py"] + sys.argv[2:]
                subcmd_out = subprocess.run(subcmd, stdout=subprocess.PIPE).stdout.decode()

        else:
            print(description)
            print("\n** unrecognized command: %s **" % cmd)


if __name__ == "__main__":
    main()
