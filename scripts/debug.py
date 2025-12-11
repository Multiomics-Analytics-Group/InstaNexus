import sys
from pathlib import Path

import pandas as pd

# Setup path come nello script vero
SCRIPT_DIR = Path("scripts").resolve()
PROJECT_ROOT = SCRIPT_DIR.parent
sys.path.append(str(PROJECT_ROOT))

try:
    from instanexus import helpers, preprocessing, visualization
    from instanexus.assembly import Assembler

    print("✅ Imports OK")
except Exception as e:
    print(f"❌ Import Error: {e}")
    sys.exit(1)

# Carica dati
try:
    df = pd.read_csv("inputs/ma1_cleaned.csv")
    print(f"✅ CSV Loaded. Columns: {df.columns.tolist()}")
except Exception as e:
    print(f"❌ CSV Error: {e}")
    sys.exit(1)

# Prova una run semplice
try:
    print("🚀 Starting Test Run...")
    sequences = df["cleaned_preds"].dropna().tolist()[:100]  # Solo primi 100 per test

    assembler = Assembler(mode="greedy", min_overlap=3, size_threshold=5)
    scaffolds = assembler.run(sequences=sequences, df_full=df)

    print(f"✅ Assembly OK. Scaffolds found: {len(scaffolds)}")

    if len(scaffolds) > 0:
        # Test Mapping (spesso fallisce qui se manca il riferimento)
        meta = helpers.get_sample_metadata("ma1", chain="light", json_path="json/sample_metadata.json")
        protein = preprocessing.normalize_sequence(meta["protein"])

        mapped = visualization.process_protein_contigs_scaffold(scaffolds, protein, max_mismatches=10, min_identity=0.8)
        print(f"✅ Mapping OK. Mapped: {len(mapped)}")

except Exception as e:
    print("\n❌ CRITICAL ERROR DURING EXECUTION:")
    print(e)
    import traceback

    traceback.print_exc()
