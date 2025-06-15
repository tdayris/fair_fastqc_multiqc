# coding: utf-8

"""
This script aggregates fastq_info text in a CSV table
"""

import pandas


def parse_fastq_info(path: str) -> dict[str, str | int]:
    """
    Read and parse a fastq_info file in memory

    Parameters:
    path    (str): Path to fastq_info output result

    Return:
    (dict[str, str|int]): Parsed content of the file
    """
    result = {}
    with open(path, "r") as fqinfo_stream:
        for line in fqinfo_stream:
            print(line)
            if line.startswith("fastq_utils"):
                result["version"] = line[:-1].split(" ")[-1]
            elif line.startswith("Scanning and indexing all reads from"):
                result["library"] = "Single Ended"
            elif line.startswith("Next "):
                result["library"] = "Pair Ended"
            elif line.startswith("Quality encoding range: "):
                quality_encoding = list(map(int, line[:-1].split(" ")[:-2]))
                result["max_quality"] = quality_encoding[1]
                result["min_quality"] = quality_encoding[0]
            elif line.startswith("Quality encoding: "):
                result["phred_encoding"] = line[:-1].split(": ")[-1]
            elif line.startswith("Read length: "):
                read_len = list(sorted(map(int, line[:-1].split(" ")[:-3])))
                result["min_read_length"] = read_len[0]
                result["mean_read_length"] = read_len[1]
                result["max_read_length"] = read_len[2]
        return result.copy()


# Parse all information
info = [parse_fastq_info(path=fqinfo) for fqinfo in snakemake.input]
print(info)

# Format into table
table = pandas.DataFrame.from_records(info)
print(table)
table.to_csv(
    str(snakemake.output),
    sep=",",
    header=True,
    index=False,
)
