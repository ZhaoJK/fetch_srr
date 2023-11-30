#!/home/cpc/jiakuan.zhao/anaconda3/envs/scanpy/bin/python 

import os
import re
import sys
import glob
import argparse
import hashlib

def parse_args(args=None):
    Description = "Generate nf-core/rnaseq samplesheet from a directory of FastQ files."
    Epilog = "Example usage: python fastq_dir_to_samplesheet.py <FASTQ_DIR> <SAMPLESHEET_FILE>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FASTQ_DIR", help="Folder containing raw FastQ files.")
    parser.add_argument("SAMPLESHEET_FILE", help="Output samplesheet file.")
    parser.add_argument(
        "-st",
        "--strandedness",
        type=str,
        dest="STRANDEDNESS",
        default="unstranded",
        help="Value for 'strandedness' in samplesheet. Must be one of 'unstranded', 'forward', 'reverse'.",
    )
    parser.add_argument(
        "-r1",
        "--read1_extension",
        type=str,
        dest="READ1_EXTENSION",
        default="_R1_001.fastq.gz",
        help="File extension for read 1.",
    )
    parser.add_argument(
        "-r2",
        "--read2_extension",
        type=str,
        dest="READ2_EXTENSION",
        default="_R2_001.fastq.gz",
        help="File extension for read 2.",
    )
    parser.add_argument(
        "-se",
        "--single_end",
        dest="SINGLE_END",
        action="store_true",
        help="Single-end information will be auto-detected but this option forces paired-end FastQ files to be treated as single-end so only read 1 information is included in the samplesheet.",
    )
    parser.add_argument(
        "-sn",
        "--sanitise_name",
        dest="SANITISE_NAME",
        action="store_true",
        help="Whether to further sanitise FastQ file name to get sample id. Used in conjunction with --sanitise_name_delimiter and --sanitise_name_index.",
    )
    parser.add_argument(
        "-sd",
        "--sanitise_name_delimiter",
        type=str,
        dest="SANITISE_NAME_DELIMITER",
        default="_",
        help="Delimiter to use to sanitise sample name.",
    )
    parser.add_argument(
        "-si",
        "--sanitise_name_index",
        type=int,
        dest="SANITISE_NAME_INDEX",
        default=1,
        help="After splitting FastQ file name by --sanitise_name_delimiter all elements before this index (1-based) will be joined to create final sample name.",
    )
    parser.add_argument('--linkfiles', action='store_true', help="If this argument is used, a unique softlink (using the md5sum of the path) for the file will be created in the same directory as the samplesheet. This is useful if samples have been sequenced on multiple flowcells but have the same filenames.")    
    return parser.parse_args(args)


def fastq_dir_to_samplesheet(
    fastq_dir,
    samplesheet_file,
    strandedness="unstranded",
    read1_extension="_R1_001.fastq.gz",
    read2_extension="_R2_001.fastq.gz",
    single_end=False,
    sanitise_name=False,
    sanitise_name_delimiter="_",
    sanitise_name_index=1,
    linkfiles=False,
):
    def sanitize_sample(path, extension):
        """Retrieve sample id from filename"""
        sample = os.path.basename(path).replace(extension, "")
        if sanitise_name:
            sample = sanitise_name_delimiter.join(
                os.path.basename(path).split(sanitise_name_delimiter)[
                    :sanitise_name_index
                ]
            )
        return sample

    def get_fastqs(extension):
        """
        Needs to be sorted to ensure R1 and R2 are in the same order
        when merging technical replicates. Glob is not guaranteed to produce
        sorted results.
        See also https://stackoverflow.com/questions/6773584/how-is-pythons-glob-glob-ordered
        """
        return sorted(
            glob.glob(os.path.join(fastq_dir, f"*{extension}"), recursive=False)
        )

    def linkfile(readfile, samplesheet_file):
        """
        Create softlinks adding md5 suffix
        """
        if len(readfile) == 0:
            return(readfile)
        sdir = os.path.dirname(samplesheet_file)
        # get dirname of the fastq file, that way md5 will be the same for r1 and r2
        fqdir=os.path.dirname(readfile)
        md5string=hashlib.md5(fqdir.encode('utf-8')).hexdigest()
        
        linkbn=os.path.basename(readfile)
        # prepend
        linkfile = os.path.join(sdir, md5string+linkbn)
        print("ln -s "+readfile+" "+linkfile)
        os.symlink(readfile, linkfile)
        return linkfile

    read_dict = {}

    ## Get read 1 files
    for read1_file in get_fastqs(read1_extension):
        sample = sanitize_sample(read1_file, read1_extension)
        if sample not in read_dict:
            read_dict[sample] = {"R1": [], "R2": []}
        read_dict[sample]["R1"].append(read1_file)

    ## Get read 2 files
    if not single_end:
        for read2_file in get_fastqs(read2_extension):
            sample = sanitize_sample(read2_file, read2_extension)
            read_dict[sample]["R2"].append(read2_file)

    ## Write to file
    if len(read_dict) > 0:
        out_dir = os.path.dirname(samplesheet_file)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)

        with open(samplesheet_file, "w") as fout:
            header = ["sample", "fastq_1", "fastq_2", "strandedness"]
            fout.write(",".join(header) + "\n")
            for sample, reads in sorted(read_dict.items()):
                for idx, read_1 in enumerate(reads["R1"]):
                    read_2 = ""
                    if idx < len(reads["R2"]):
                        read_2 = reads["R2"][idx]
                    # Check if softlinks should be created
                    if linkfiles:
                        read_1 = linkfile(read_1,samplesheet_file)
                        read_2 = linkfile(read_2,samplesheet_file)                        
                    sample_title = re.sub("_S\d*_L\d*$", "",sample)
                    sample_info = ",".join([sample_title, read_1, read_2, strandedness])
                    fout.write(f"{sample_info}\n")
    else:
        error_str = (
            "\nWARNING: No FastQ files found so samplesheet has not been created!\n\n"
        )
        error_str += "Please check the values provided for the:\n"
        error_str += "  - Path to the directory containing the FastQ files\n"
        error_str += "  - '--read1_extension' parameter\n"
        error_str += "  - '--read2_extension' parameter\n"
        print(error_str)
        sys.exit(1)


def main(args=None):
    args = parse_args(args)

    strandedness = "unstranded"
    if args.STRANDEDNESS in ["unstranded", "forward", "reverse"]:
        strandedness = args.STRANDEDNESS

    fastq_dir_to_samplesheet(
        fastq_dir=args.FASTQ_DIR,
        samplesheet_file=args.SAMPLESHEET_FILE,
        strandedness=strandedness,
        read1_extension=args.READ1_EXTENSION,
        read2_extension=args.READ2_EXTENSION,
        single_end=args.SINGLE_END,
        sanitise_name=args.SANITISE_NAME,
        sanitise_name_delimiter=args.SANITISE_NAME_DELIMITER,
        sanitise_name_index=args.SANITISE_NAME_INDEX,
        linkfiles=args.linkfiles,
    )


if __name__ == "__main__":
    sys.exit(main())


