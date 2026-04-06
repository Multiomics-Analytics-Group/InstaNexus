#!/usr/bin/env python3

from pathlib import Path

import pandas as pd

from instanexus.preprocessing import remove_modifications, main as preprocessing_main


def test_remove_modifications():
    assert remove_modifications("A(ox)BC(mod)D") == "ABCD"
    assert remove_modifications("A[UNIMOD:21]BC[UNIMOD:35]D") == "ABCD"
    assert remove_modifications("A(ox)[UNIMOD:21]BC(mod)[UNIMOD:35]D") == "ABCD"
    assert remove_modifications(None) is None
    assert remove_modifications("ACD") == "ACD"
    assert remove_modifications("A(I)BCD") == "ABCD"
    assert remove_modifications("A(ox)B(I)C(mod)D") == "ABCD"
    assert remove_modifications("A(ox)[UNIMOD:21]B(I)C(mod)[UNIMOD:35]D") == "ABCD"
    assert remove_modifications("AI BCD") == "AL BCD"
    assert remove_modifications("A(ox)I B(mod)CD") == "AL BCD"


def test_preprocessing_no_metadata(tmp_path: Path) -> None:
    """Preprocessing must succeed when metadata_json and contaminants_fasta are both None."""
    input_csv = tmp_path / "sample.csv"
    output_csv = tmp_path / "out" / "cleaned.csv"

    df = pd.DataFrame(
        {
            "preds": ["ACDEF", "GHIKL", "MNPQR"],
            "log_probs": [-0.1, -0.2, -0.3],
        }
    )
    df.to_csv(input_csv, index=False)

    preprocessing_main(
        input_csv=str(input_csv),
        metadata_json=None,
        contaminants_fasta=None,
        chain="",
        reference=False,
        conf=None,
        fdr=None,
        output_csv_path=str(output_csv),
    )

    assert output_csv.exists()
    result = pd.read_csv(output_csv)
    assert "cleaned_preds" in result.columns
    # I is normalized to L during remove_modifications
    assert set(result["cleaned_preds"]) == {"ACDEF", "GHLKL", "MNPQR"}
    assert "protease" not in result.columns
    assert "mapped" not in result.columns
