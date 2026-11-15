"""Recover integer UMI counts from the normalized GSE121861 expression matrix.

GSE121861 (Kumar et al. 2018, six syngeneic mouse tumor models) provides a
cells x genes matrix of normalized values. Within each cell every nonzero
value is an integer multiple of one unit (the per-cell scale factor), so the
counts are value / unit. This script checks that property for every cell and
writes a genes x cells sparse count matrix (Matrix Market) with the cell and
gene tables.

    python data-raw/bayesprism/01_gse121861_counts.py <geo_dir> <out_dir>

<geo_dir> holds the three GSE121861_syngeneic_*.csv.gz files.
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io, sparse

# The values are written with five significant digits, so value / unit can
# differ from the integer count by a relative error of about 1e-5 per number.
TOLERANCE = 2e-4  # largest allowed relative distance of value / unit from an integer


def recover(geo_dir: Path, out_dir: Path) -> None:
    cells = pd.read_csv(geo_dir / "GSE121861_syngeneic_column_data.csv.gz")
    genes = pd.read_csv(geo_dir / "GSE121861_syngeneic_row_data.csv.gz")
    print(f"cell table: {cells.shape[0]} rows; gene table: {genes.shape[0]} rows")

    rows, cols, vals = [], [], []
    units, worst = [], []
    n_cells = 0
    with gzip.open(geo_dir / "GSE121861_syngeneic_expression.csv.gz", "rt") as fh:
        for i, line in enumerate(fh):
            x = np.fromstring(line, sep=",")
            if x.size != genes.shape[0]:
                raise ValueError(f"cell {i}: {x.size} values, expected {genes.shape[0]}")
            nz = np.flatnonzero(x)
            unit = x[nz].min()
            ratio = x[nz] / unit
            counts = np.rint(ratio)
            worst.append((np.abs(ratio - counts) / np.maximum(counts, 1)).max())
            units.append(unit)
            rows.append(nz)
            cols.append(np.full(nz.size, i))
            vals.append(counts.astype(np.int32))
            n_cells += 1
    worst = np.array(worst)
    print(f"cells read: {n_cells}")
    print(f"largest relative distance from an integer per cell: median {np.median(worst):.2e}, max {worst.max():.2e}")
    bad = int((worst > TOLERANCE).sum())
    print(f"cells above the tolerance {TOLERANCE}: {bad}")
    if n_cells != cells.shape[0]:
        raise ValueError("the expression matrix and the cell table have different numbers of cells")
    if bad:
        raise ValueError("some cells are not integer multiples of one unit; counts cannot be recovered")

    m = sparse.csc_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(genes.shape[0], n_cells),
        dtype=np.int32,
    )
    total = np.asarray(m.sum(axis=0)).ravel()
    print(f"UMI per cell: median {np.median(total):.0f}, range {total.min()}-{total.max()}")
    print(f"cells with a minimum nonzero count of 1 (unit = one UMI): all {n_cells}, by construction")

    cells = cells.copy()
    cells["cell"] = [f"c{i:05d}" for i in range(n_cells)]
    cells["unit"] = units
    cells["umi"] = total
    out_dir.mkdir(parents=True, exist_ok=True)
    io.mmwrite(out_dir / "counts.mtx", m, field="integer")
    cells.to_csv(out_dir / "cells.csv", index=False)
    genes.to_csv(out_dir / "genes.csv", index=False)
    print("cells per model:")
    print(cells.groupby("Model").size().to_string())
    print(f"written to {out_dir}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    recover(Path(sys.argv[1]), Path(sys.argv[2]))
