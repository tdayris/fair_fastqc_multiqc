# -*- coding: utf-8 -*-

"""Snakemake wrapper for bash copy within the IGR's Flamingo cluster"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


import logging
import os

from random import randint
from snakemake.shell import shell
from tempfile import TemporaryDirectory
from typing import Callable


# Logging behaviour
try:
    logging.basicConfig(filename=snakemake.log[0], filemode="w", level=logging.DEBUG)
except Exception:
    logging.basicConfig(level=logging.DEBUG)
    logging.warning(
        "No logging file was provided in Snakemake rule. "
        "Logging will be written in standard output."
    )


def on_flamingo() -> bool:
    """
    Return True if pipeline is executed at GustaveRoussy's computing cluster: Flamingo
    """
    return os.environ.get("HOSTNAME", "").lower().startswith("flamingo")


def consider_compression(copy_func: Callable) -> None:
    """
    Decorator designed to handle fastq files compression

    Parameters:
    copy_func   (Callable): The function used to copy

    Return: (None)
    """

    def handle_gzip(*args, **kwargs) -> None:
        logging.debug(f"{args=}, {kwargs=}")
        if not kwargs["src"].lower().endswith(".gz"):
            with TemporaryDirectory() as tmpdir:
                logging.info(
                    f"{kwargs['src']} was not compressed. Compressing it after copy/link/recovery."
                )
                dest: str = kwargs["dest"]
                tmp_dest = f"{tmpdir}/unzipped"
                copy_func(
                    src=kwargs["src"],
                    dest=tmp_dest,
                )
                log_cmd: str = snakemake.log_fmt_shell(
                    stdout=False, stderr=True, append=True
                )
                cmd: str = (
                    f"gzip --verbose --force --stdout {tmp_dest} > {dest} {log_cmd}"
                )
                logging.debug(cmd)
                shell(cmd)

        else:
            copy_func(*args, **kwargs)

    return handle_gzip


@consider_compression
def bash_rsync(src: str, dest: str) -> None:
    """
    Perform bash copy using rsync

    Parameters:
    src     (str): Path to source file
    dest    (str): Path to destination file

    Return (None)
    """
    log_cmd: str = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

    # Consider possible temporary directory use with rsync
    with TemporaryDirectory() as bash_copy_or_link_tmpdir:

        # Build and execute command line
        cmd: str = (
            f"rsync --temp-dir={bash_copy_or_link_tmpdir} "
            "--verbose --checksum --recursive "
            "--human-readable --progress --partial"
            f"{src} {dest} {log_cmd}"
        )
        logging.debug(cmd)
        shell(cmd)


@consider_compression
def bash_iget(src: str, dest: str, threads: int = 0) -> None:
    """
    Perform bash recovery from KDI iRODS using iget

    Parameters:
    src     (str): Path to source file
    dest    (str): Path to destination file

    Return (None)
    """
    log_cmd: str = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

    # Build and execute command line
    cmd: str = f"iget -N {threads} -K -r -V {src} {dest} {log_cmd}"
    logging.debug(cmd)
    shell(cmd)


@consider_compression
def bash_ln(src: str, dest: str) -> None:
    """
    Perform bash symbolic linking using ln

    Parameters:
    src     (str): Path to source file
    dest    (str): Path to destination file

    Return (None)
    """
    log_cmd: str = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

    # Build and execute command line
    cmd: str = f"ln --force --relative --symbolic --verbose {src} {dest} {log_cmd}"
    logging.debug(cmd)
    shell(cmd)


def cat_files(dest: str, *src: list[str]) -> None:
    """
    Concatenate multiple files

    Parameters:
    dest    (str)       : Destination file
    src     (list[str]) : List of source files

    Return (None)
    """
    log_cmd: str = snakemake.log_fmt_shell(stdout=False, stderr=True, append=True)
    infiles: str = " ".join(src)

    # Build and execute command line
    cmd: str = f"cat {infiles} > {dest} {log_cmd}"
    logging.debug(cmd)
    shell(cmd)


def make_available(
    src: str, dest: str, cold_storage: tuple[str] | str, irods_prefix: tuple[str] | str
) -> None:
    """
    Decide whether to use rsync, iget, or ln
    to make a file available.

    Parameters:
    dest            (str)       : Destination file
    src             (list[str]) : List of source files
    cold_storage    (tuple[str]): Mounting points not reachable from computing nodes
    irods_prefix    (tuple[str]): Mounting points which require the use of iget

    Return (None)
    """
    if src.lower().startswith(cold_storage):
        bash_rsync(src=src, dest=dest)
    elif src.lower().startswith(irods_prefix):
        bash_iget(src=src, dest=dest)
    else:
        bash_ln(src=src, dest=dest)


def copy_then_concat(
    dest: str,
    cold_storage: tuple[str] | str,
    irods_prefix: tuple[str] | str,
    *src: list[str],
) -> None:
    """
    Individually make each source file(s) available,
    then concatenate them all in one destination file

    Parameters:
    dest            (str)       : Destination file
    src             (list[str]) : List of source files
    cold_storage    (tuple[str]): Mounting points not reachable from computing nodes
    irods_prefix    (tuple[str]): Mounting points which require the use of iget

    Return (None)
    """
    with TemporaryDirectory() as tmpdir:
        outfiles = []
        for path in src:
            tmp_dest = f"{tmpdir}/{os.path.basename(path)}.{randint(0, 100_000_000)}"
            make_available(path, tmp_dest, cold_storage, irods_prefix)
            outfiles.append(tmp_dest)
        cat_files(dest, *outfiles)


def copy_or_concat(
    dest: str, src: str, cold_storage: tuple[str] | str, irods_prefix: tuple[str] | str
) -> None:
    """
    Recognize whether the source is composed of multiple paths, then
    recognize whether the destination is composed of multiple paths, then
    decides whether the files should be concatenated or made available
    one by one.

    Parameters:
    dest            (str)       : Destination file
    src             (list[str]) : List of source files
    cold_storage    (tuple[str]): Mounting points not reachable from computing nodes
    irods_prefix    (tuple[str]): Mounting points which require the use of iget

    Return (None)
    """
    src_sep: str | None = None
    src_len: int = 1
    if "," in src:
        src_sep = ","
    elif ";" in src:
        src_sep = ";"

    if src_sep:
        src: list[str] = src
        src_len = len(src)

    dest_sep: str | None = None
    dest_len: int = 1
    if "," in dest:
        dest_sep = ","
    elif ";" in dest:
        dest_sep = ";"

    if dest_sep:
        dest: list[str] = dest
        dest_len = len(dest)

    if src_len == dest_len == 1:
        dest = dest[0]
        logging.info(f"Making {src} available at {dest}")
        make_available(src, dest, cold_storage, irods_prefix)
    elif src_len == dest_len:
        for source, destination in zip(src, dest):
            logging.info(f"Making {source} available at {destination}")
            make_available(source, destination, cold_storage, irods_prefix)
    elif (src_len > 1) and (dest_len == 1):
        dest = dest[0]
        logging.info(f"Concatenating each {src} in {dest}")
        copy_then_concat(dest, cold_storage, irods_prefix, *src)
    else:
        raise ValueError("Could not determinate how to make files " "available.")


cold_storage: tuple[str] | str = snakemake.params.get(
    "cold_storage",
    (
        "/mnt/isilon",
        "/mnt/archivage",
        "/mnt/install",
        "glustergv0",
        "nfs01",
        "nas01_test",
    ),
)

irods_prefix: tuple[str] | str = snakemake.params.get(
    "irods_prefix", ("/odin/kdi/dataset/")
)


output_directory = os.path.realpath(os.path.dirname(snakemake.output[0]))
log_cmd: str = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
shell(f"mkdir --parent --verbose {output_directory} {log_cmd}")

sources = snakemake.params.get("in_files", snakemake.input)
destinations = snakemake.output

copy_or_concat(
    dest=destinations, src=sources, cold_storage=cold_storage, irods_prefix=irods_prefix
)
