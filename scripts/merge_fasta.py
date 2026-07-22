#!/usr/bin/env python3

import gzip
import hashlib
import sys
from pathlib import Path
from typing import Iterator, TextIO


def open_text(filename: str) -> TextIO:
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    return open(filename, "rt")


def read_fasta(filename: str) -> Iterator[tuple[str, str, str]]:
    header = None
    sequence_parts: list[str] = []

    with open_text(filename) as handle:
        for line in handle:
            line = line.rstrip()

            if line.startswith(">"):
                if header is not None:
                    name = header.split()[0]
                    yield name, header, "".join(sequence_parts)

                header = line[1:]
                sequence_parts = []
            else:
                sequence_parts.append(line.strip())

    if header is not None:
        name = header.split()[0]
        yield name, header, "".join(sequence_parts)


def main() -> None:
    if len(sys.argv) < 2:
        sys.exit(
            f"Usage: {Path(sys.argv[0]).name} reference1.fa "
            "reference2.fa [...] > merged.fa"
        )

    seen: dict[str, tuple[int, str, str]] = {}

    for filename in sys.argv[1:]:
        for name, header, sequence in read_fasta(filename):
            sequence = sequence.upper()
            digest = hashlib.sha256(sequence.encode()).hexdigest()
            current = (len(sequence), digest, filename)

            if name in seen:
                old_length, old_digest, old_file = seen[name]

                if old_length != len(sequence) or old_digest != digest:
                    sys.exit(
                        f"ERROR: reference name {name!r} has different "
                        f"sequences in {old_file!r} and {filename!r}"
                    )

                print(
                    f"Skipping identical duplicate {name!r} from {filename}",
                    file=sys.stderr,
                )
                continue

            seen[name] = current
            print(f">{header}")

            for start in range(0, len(sequence), 60):
                print(sequence[start : start + 60])


if __name__ == "__main__":
    main()