from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory
from urllib.parse import urlsplit
from urllib.request import urlretrieve
import pandas as pd
import pyranges1 as pr
import pysam


def fetch_file(link: str, output_dir: str, filename: str = None) -> Path:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = filename or PurePosixPath(urlsplit(link).path).name
    if not filename:
        raise ValueError("Filename doesn't specified")

    output_path = output_dir / filename

    if output_path.is_file():
        return output_path

    with TemporaryDirectory(dir=output_dir) as temp_dir:
        temp_path = Path(temp_dir) / 'download'
        urlretrieve(link, temp_path)
        temp_path.replace(output_path)

    return output_path


def remove_file(filepath, remove_file=0):
    if remove_file == 1:
        Path(filepath).unlink()


def create_candidate_regions(
    vcf_filepath,
    bed_filepath=None,
    distance=3,
    only_pass_variants=True,
):
    variants = []

    with pysam.VariantFile(vcf_filepath) as vcf:
        if bed_filepath is None:
            variant_records = vcf
        else:
            intervals = pd.read_csv(
                bed_filepath,
                sep="\t",
                header=None,
                usecols=range(3),
                names=["Chromosome", "Start", "End"],
                dtype={"Chromosome": str, "Start": int, "End": int},
            )

            variant_records = (
                variant
                for chrom, start, end in intervals.itertuples(
                    index=False, name=None
                )
                for variant in vcf.fetch(chrom, start, end)
            )

        for variant in variant_records:
            if only_pass_variants and "PASS" not in variant.filter:
                continue

            if variant.stop - variant.start > distance:
                continue

            variants.append((variant.chrom, variant.start, variant.stop))

    df = (
        pd.DataFrame(variants, columns=["Chromosome", "Start", "End"])
        .drop_duplicates(["Chromosome", "Start"])
        .sort_values(["Chromosome", "Start"])
    )

    regions = []

    for chrom, group in df.groupby("Chromosome", sort=False):
        records = list(
            group[["Start", "End"]].itertuples(index=False, name=None)
        )
        last_end = -1

        for i, (start, end) in enumerate(records):
            count = 1
            limit = start + distance

            for j in range(i + 1, len(records)):
                next_start, next_end = records[j]

                if next_start >= limit:
                    break

                if next_end > limit:
                    continue

                end = max(end, next_end)
                count += 1

            # Keep overlapping windows, but skip contained ones.
            if count >= 2 and end > last_end:
                regions.append((chrom, start, end))
                last_end = end

    regions_df = pd.DataFrame(
        regions,
        columns=["Chromosome", "Start", "End"],
    )

    return pr.PyRanges(regions_df)