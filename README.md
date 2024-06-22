# annot2fasta - Create FASTA Files from Annotation

## Overview

`annot2fasta` is a versatile Python script designed to simplify the process of converting annotation data into FASTA files, which are commonly used for biological sequence analysis. This script takes annotation files or folders as input and extracts nucleotide sequences based on the provided coordinates in the annotations. By doing so, it streamlines the conversion of annotated sequences into a more usable FASTA format.

## Usage

The `annot2fasta` script provides a range of options for converting annotation data to FASTA format. Here are the command-line parameters and their descriptions:

```shell
annot2fasta -i <path|file> -o <output_folder> [-h] [-v]
```

### Required Parameters

- `-i <path|file>`: Specifies the input annotation file or folder containing multiple annotation files.


### Optional Parameters
- `-o <output_folder>`: Sets the name of the output folder where the resulting FASTA files will be saved (deafault: fasta_files).

- `-h` or `--help`: Displays the help message, which provides information about the tool's usage.
- `-v`: Prints the current version of the `annot2fasta` script.

## Usage Examples

1. To create FASTA files from a single annotation file:
   ```shell
   annot2fasta -i annotation_file.tab -o output_folder
   ```

2. To process an entire folder of annotation files and save the resulting FASTA files in the default "fasta_files" output folder:
   ```shell
   annot2fasta -i annotation_folder
   ```

## Authors

This tool was developed by Giuliana Pola and Arthur Gruber.

## Contact Information

For any questions, bug reports, or feedback, please do not hesitate to reach out to:

- Arthur Gruber: argruber@usp.br
- Giuliana L. Pola: giulianapola@usp.br
