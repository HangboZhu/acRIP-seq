#!/usr/bin/env python3
"""Generate a complete-pairs TSV for one species from the download manifest
plus on-disk size verification, merging multi-run samples (one GSM sequenced
as several SRR runs) per read direction so every sample stays one PE library.

Usage:
  python3 make_pairs_tsv.py --species mouse \
      --manifest /path/download_manifest.full.tsv.bak \
      --out /tmp/mouse_pairs_current.tsv
Output schema matches make_config_from_pairs.py input.
"""
import argparse
import os
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--species", required=True)
ap.add_argument("--manifest", required=True)
ap.add_argument("--dataroot", default="/hpcdisk1/apod/shijiachen/yuliao/4/data/")
ap.add_argument("--out", required=True)
args = ap.parse_args()

rows = [l.rstrip("\r\n").split("\t") for l in open(args.manifest)]
idx = {k: i for i, k in enumerate(rows[0])}
ROOT = args.dataroot

byg = defaultdict(lambda: defaultdict(list))
for c in rows[1:]:
    if c[idx["species"]] == args.species and c[idx["pair_id"]]:
        byg[c[idx["pair_id"]]][c[idx["GSM"]]].append(c)


def complete(cs):
    for c in cs:
        p = os.path.join(ROOT, c[idx["local_path"]])
        try:
            if str(os.stat(p).st_size) != c[idx["file_size"]]:
                return False
        except OSError:
            return False
    return True


def merge_runs(gsm_dir_files):
    """gsm_dir_files: list of (path, direction) for one GSM side; if >2 files
    (multi run), concatenate per direction into merged_runs/.
    Manifest read_direction values are R1/R2 (paired) or SE (single-end)."""
    def norm(d):
        d = str(d).strip().rstrip("\r").upper()
        if d in ("R1", "1", "READ1", "F"):
            return "1"
        if d in ("R2", "2", "READ2", "R"):
            return "2"
        return "SE"

    by_dir = defaultdict(list)
    for p, d in gsm_dir_files:
        by_dir[norm(d)].append(p)

    if set(by_dir) == {"SE"}:  # 单端样本
        return [by_dir["SE"][0]], False
    if set(by_dir) == {"1", "2"} and all(len(v) == 1 for v in by_dir.values()):
        return [by_dir["1"][0], by_dir["2"][0]], False
    if "SE" in by_dir and len(by_dir) > 1:
        raise ValueError(f"样本同时含 SE 与 PE 文件: {[p for p, _ in gsm_dir_files]}")

    d0 = os.path.dirname(gsm_dir_files[0][0])
    md = os.path.join(d0, "merged_runs")
    os.makedirs(md, exist_ok=True)
    srrs = sorted({os.path.basename(p).rsplit("_", 1)[0] for p, _ in gsm_dir_files})
    tag = srrs[0] + "_plus_" + srrs[-1].replace("SRR", "")
    out = []
    for direction in sorted(by_dir):
        dst = os.path.join(md, f"{tag}_{direction}.fastq.gz")
        if not os.path.exists(dst) or os.path.getsize(dst) == 0:
            with open(dst, "wb") as w:
                for src in sorted(by_dir[direction]):
                    with open(src, "rb") as r:
                        while True:
                            chunk = r.read(64 * 1024 * 1024)
                            if not chunk:
                                break
                            w.write(chunk)
        out.append(dst)
    return sorted(out), True


out_rows = ["\t".join(["species", "GSE", "pair_id", "IP_GSM", "Input_GSM", "IP_SRR",
                       "Input_SRR", "IP_fastq", "Input_fastq", "total_size_GB",
                       "pair_semantics", "eligibility"])]
n, merged = 0, []
for pid, gm in sorted(byg.items()):
    if len(gm) != 2 or not all(complete(cs) for cs in gm.values()):
        continue
    ipg = [g for g, cs in gm.items() if cs[0][idx["assay_role"]] == "IP"]
    ing = [g for g, cs in gm.items() if cs[0][idx["assay_role"]] == "Input"]
    if len(ipg) != 1 or len(ing) != 1:
        continue
    n += 1
    fq = {}
    for label, cs in (("IP", gm[ipg[0]]), ("Input", gm[ing[0]])):
        files, did = merge_runs(
            [(os.path.join(ROOT, c[idx["local_path"]]), c[idx["read_direction"]]) for c in cs])
        if did:
            merged.append(f"{pid}.{label}")
        fq[label] = files
    gse = gm[ipg[0]][0][idx["GSE"]]
    gb = sum(int(c[idx["file_size"]]) for cs in gm.values() for c in cs) / 1e9
    out_rows.append("\t".join([
        args.species, gse, pid, ipg[0], ing[0],
        ";".join(c[idx["SRR"]] for c in gm[ipg[0]]),
        ";".join(c[idx["SRR"]] for c in gm[ing[0]]),
        ";".join(fq["IP"]), ";".join(fq["Input"]),
        f"{gb:.3f}", "Input vs IP (anti-ac4C)", "eligible"]))

with open(args.out, "w") as f:
    f.write("\n".join(out_rows) + "\n")
print(f"{args.species}: {n} complete pairs; merged multi-run samples: {merged or 'none'}")
