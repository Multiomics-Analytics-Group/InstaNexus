#!/usr/bin/env python3

"""Assembly output must not depend on Python's string hash seed.

Each hash seed runs in its own subprocess (``PYTHONHASHSEED`` is fixed at
interpreter start-up), which runs every assembly target on the same synthetic
input and prints the results as JSON. The outputs must be identical across
seeds, including order.
"""

import json
import os
import subprocess
import sys
import textwrap

import pytest

HASH_SEEDS = ["0", "1", "2", "3", "42", "12345"]

TARGETS = [
    "assemble_contigs_greedy",
    "merge_contigs_greedy",
    "scaffold_iterative_greedy",
    "assemble_contigs_dbg",
    "merge_sequences_dbg",
    "scaffold_iterative_dbg",
    "assemble_contigs_dbgx",
    "refine_using_overlap_graph",
    "scaffold_iterative_dbgx",
    "Assembler[greedy]",
    "Assembler[dbg]",
    "Assembler[dbg_weighted]",
    "Assembler[dbg_weighted+refine]",
    "Assembler[dbgX]",
    "Assembler[fusion]",
    "Assembler[multimodal_dbg]",
    "Assembler[multimodal_dbg+df]",
    "Assembler[hybrid_dbg+df]",
]

_RUNNER = textwrap.dedent(
    """
    import json
    import logging
    import random
    import sys

    import pandas as pd

    from instanexus import assembly as asm

    logging.disable(logging.CRITICAL)

    # Bovine serum albumin (P02769), residues 1-400
    PROTEIN = (
        "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASL"
        "RETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIET"
        "MREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSH"
        "CIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAK"
    )
    AMINO_ACIDS = "ACDEFGHKLMNPQRSTVWY"

    # random.Random with an int seed does not depend on PYTHONHASHSEED
    rng = random.Random(11)
    peptides = []
    for _ in range(150):
        length = rng.randint(7, 15)
        start = rng.randrange(0, len(PROTEIN) - length)
        peptide = PROTEIN[start : start + length]
        if rng.random() < 0.15:
            pos = rng.randrange(length)
            peptide = peptide[:pos] + rng.choice(AMINO_ACIDS) + peptide[pos + 1 :]
        peptides.append(peptide)
    peptides += peptides[:20]

    df = pd.DataFrame(
        {
            "cleaned_preds": peptides,
            "peptide_abundance": [rng.choice([1e5, 1e6, 1e7]) for _ in peptides],
            "ion_match_intensity": [rng.choice([0.2, 0.5, 0.8]) for _ in peptides],
            "instanovo_token_log_probabilities": [
                str([rng.choice([-0.01, -0.1, -0.5]) for _ in range(len(p))]) for p in peptides
            ],
        }
    )

    params = dict(kmer_size=6, min_overlap=4, size_threshold=10, min_weight=2)


    def assembler(mode, **kwargs):
        return asm.Assembler(mode=mode, **{**params, **kwargs})


    greedy_contigs = asm.assemble_contigs_greedy(peptides, 4)
    long_contigs = [c for c in greedy_contigs if len(c) > 10]

    targets = {
        "assemble_contigs_greedy": lambda: greedy_contigs,
        "merge_contigs_greedy": lambda: asm.merge_contigs_greedy(peptides),
        "scaffold_iterative_greedy": lambda: asm.scaffold_iterative_greedy(greedy_contigs, 4, 10, disable_tqdm=True),
        "assemble_contigs_dbg": lambda: asm.assemble_contigs_dbg(
            asm.get_debruijn_edges_from_kmers(asm.get_kmers(peptides, 6))
        ),
        "merge_sequences_dbg": lambda: asm.merge_sequences_dbg(peptides, disable_tqdm=True),
        "scaffold_iterative_dbg": lambda: asm.scaffold_iterative_dbg(long_contigs, 4, 10, disable_tqdm=True),
        "assemble_contigs_dbgx": lambda: [
            c.seq
            for c in asm.assemble_contigs_dbgx(
                asm.filter_low_weight_edges(asm.build_dbg_from_kmers(asm.get_kmers(peptides, 6)), min_weight=2),
                min_length=10,
            )
        ],
        "refine_using_overlap_graph": lambda: asm.refine_using_overlap_graph(long_contigs, 4),
        "scaffold_iterative_dbgx": lambda: asm.scaffold_iterative_dbgx(peptides, kmer_size=6, size_threshold=10),
        "Assembler[greedy]": lambda: assembler("greedy").run(peptides),
        "Assembler[dbg]": lambda: assembler("dbg").run(peptides),
        "Assembler[dbg_weighted]": lambda: assembler("dbg_weighted").run(peptides),
        "Assembler[dbg_weighted+refine]": lambda: assembler("dbg_weighted", refine_rounds=10).run(peptides),
        "Assembler[dbgX]": lambda: assembler("dbgX").run(peptides),
        "Assembler[fusion]": lambda: assembler("fusion").run(peptides),
        "Assembler[multimodal_dbg]": lambda: assembler("multimodal_dbg").run(peptides),
        "Assembler[multimodal_dbg+df]": lambda: assembler("multimodal_dbg").run(peptides, df_full=df),
        "Assembler[hybrid_dbg+df]": lambda: assembler("hybrid_dbg").run(peptides, df_full=df),
    }

    results = {}
    for name in sys.argv[1:]:
        try:
            results[name] = {"ok": list(targets[name]())}
        except Exception as e:
            results[name] = {"error": f"{type(e).__name__}: {e}"}
    json.dump(results, sys.stdout)
    """
)


def _run_with_hash_seed(seed: str) -> dict:
    """Run every target in a fresh interpreter with the given PYTHONHASHSEED.

    Args:
        seed: Value for PYTHONHASHSEED.

    Returns:
        Mapping of target name to its result (``{"ok": [...]}`` or ``{"error": "..."}``).
    """
    env = {**os.environ, "PYTHONHASHSEED": seed, "TQDM_DISABLE": "1"}
    proc = subprocess.run(
        [sys.executable, "-c", _RUNNER, *TARGETS],
        env=env,
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert proc.returncode == 0, f"runner failed with PYTHONHASHSEED={seed}:\n{proc.stderr[-2000:]}"

    return json.loads(proc.stdout)


@pytest.fixture(scope="module")
def results_by_seed() -> dict:
    return {seed: _run_with_hash_seed(seed) for seed in HASH_SEEDS}


@pytest.mark.parametrize("target", TARGETS)
def test_output_independent_of_hash_seed(results_by_seed: dict, target: str) -> None:
    reference_seed = HASH_SEEDS[0]
    reference = results_by_seed[reference_seed][target]
    assert "error" not in reference, f"{target} raised: {reference['error']}"
    assert reference["ok"], f"{target} returned no sequences; the input does not exercise it"

    differing = {}
    for seed in HASH_SEEDS[1:]:
        result = results_by_seed[seed][target]
        if result != reference:
            if "error" in result:
                differing[seed] = result["error"]
            elif sorted(result["ok"]) == sorted(reference["ok"]):
                differing[seed] = "same sequences, different order"
            else:
                differing[seed] = "different sequences"

    assert not differing, f"{target}: output differs from PYTHONHASHSEED={reference_seed}: {differing}"
