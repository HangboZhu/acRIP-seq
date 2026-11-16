#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Collect all acRIP-seq pair statistics into a unified JSON
"""

import json
import argparse
from pathlib import Path
from datetime import datetime


def main():
    parser = argparse.ArgumentParser(description='Collect all acRIP-seq statistics into one JSON')
    parser.add_argument('-i', '--inputs', nargs='+', required=True, help='Input statistics JSON files')
    parser.add_argument('-o', '--output', required=True, help='Output merged JSON file')
    args = parser.parse_args()

    all_samples = {}
    pairs = []

    for json_file in args.inputs:
        json_path = Path(json_file)
        if not json_path.exists():
            continue
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)

        pair_name = json_path.stem.replace('statistics_', '')
        # Merge samples (avoid duplicates)
        for sample_name, sample_data in data.get('samples', {}).items():
            if sample_name not in all_samples:
                all_samples[sample_name] = sample_data

        pairs.append({
            "pair_name": pair_name,
            "ip_sample": data.get("ip_sample", ""),
            "input_sample": data.get("input_sample", ""),
            "peak_count_exomepeak2": data.get("peak_count_exomepeak2", 0),
            "peak_count_bed6": data.get("peak_count_bed6", 0),
            "m6a_motif_peaks": data.get("m6a_motif_peaks", 0),
            "m6a_motif_ratio": data.get("m6a_motif_ratio", 0.0),
        })

    merged = {
        "collection_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "samples": all_samples,
        "pairs": pairs
    }

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(merged, f, indent=2, ensure_ascii=False)

    print(f"Merged statistics saved to: {output_path}")


if __name__ == "__main__":
    main()
