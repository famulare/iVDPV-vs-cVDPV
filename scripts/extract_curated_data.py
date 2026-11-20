#!/usr/bin/env python3
"""Extract Sabin-anchored VP1 alignments and canonical metadata from enterovirus-genbank-curated.

Reads the promoted release under `<curated-dir>/final/` and writes a compact,
self-contained snapshot into `<output-dir>/`:

  - PV1_vp1.fasta.gz, PV2_vp1.fasta.gz, PV3_vp1.fasta.gz
      Sabin-anchored codon-aligned VP1 nucleotide sequences (906 / 903 / 900 nt),
      with the canonical Sabin reference (AY184219 / AY184220 / AY184221) placed
      explicitly as sequence #1 in each file.
  - poliovirus_vp1_metadata.tsv.gz
      All 26 canonical columns from `sequence_metadata_vouched.tsv.gz` for the
      extracted rows, plus linked GenBank citation metadata (`pubmed_id`,
      `reference_title`) and VP1 alignment QC metrics (coverage, unambiguous nt,
      ambiguities, internal deletion events/codons, frameshift/stop flags).
  - source_manifest.json
      Provenance record pinning the upstream release version, schema version,
      build manifest hashes, alignment hashes, extraction parameters, and output
      SHA-256 digests.

Re-run whenever a new release of `enterovirus-genbank-curated` lands.
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path

SABIN_REFS = {
    "PV1": ("AY184219", 906),
    "PV2": ("AY184220", 903),
    "PV3": ("AY184221", 900),
}

BASES = frozenset("ACGT")
AMBIG_CHARS = frozenset("NRYSWKMBDHV")
GAP_CHARS = frozenset("-.")

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_gzip_deterministic(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        gz.write(content)


def read_fasta_gz(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur_id: str | None = None
    chunks: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(chunks)
                cur_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if cur_id is not None:
            seqs[cur_id] = "".join(chunks)
    return seqs


def read_vp1_coordinates(coords_path: Path) -> dict[str, tuple[str, int, int, int]]:
    out: dict[str, tuple[str, int, int, int]] = {}
    with coords_path.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["region"] == "VP1":
                out[row["serotype"]] = (
                    row["ref_accession"],
                    int(row["start"]),
                    int(row["end"]),
                    int(row["length"]),
                )
    return out


def read_first_references(refs_path: Path, accessions: set[str]) -> dict[str, tuple[str, str]]:
    """Return first non-Direct-Submission (pubmed_id, title) per accession."""
    out: dict[str, tuple[str, str]] = {}
    if not refs_path.exists():
        return out
    with gzip.open(refs_path, "rt", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            acc = row["record_id"].split(".")[0]
            if acc not in accessions or acc in out:
                continue
            title = (row.get("title") or "").strip()
            pmid = (row.get("pubmed_id") or "").strip()
            if title and title.lower() != "direct submission":
                out[acc] = (pmid, title)
    return out


def analyze_vp1_slice(seq: str) -> dict[str, object]:
    """Compute VP1 alignment QC metrics for a Sabin-anchored nucleotide sequence."""
    s = seq.upper().replace("U", "T").replace(".", "-")
    width = len(s)
    if width % 3 != 0:
        raise ValueError(f"VP1 slice length {width} is not divisible by 3")

    invalid = {ch for ch in s if ch not in BASES and ch not in AMBIG_CHARS and ch not in GAP_CHARS}
    if invalid:
        raise ValueError(f"Unsupported VP1 alignment characters: {sorted(invalid)}")

    material_idx = [i for i, ch in enumerate(s) if ch != "-"]
    non_gap_nt = len(material_idx)
    acgt_nt = sum(1 for ch in s if ch in BASES)
    ambig_nt = non_gap_nt - acgt_nt

    if non_gap_nt == 0:
        return {
            "seq": s,
            "vp1_non_gap_nt": 0,
            "vp1_acgt_nt": 0,
            "vp1_ambig_nt": 0,
            "vp1_coverage_frac": 0.0,
            "vp1_span_start_nt": "",
            "vp1_span_end_nt": "",
            "vp1_internal_del_nt": 0,
            "vp1_internal_del_events": 0,
            "vp1_internal_del_codons": 0,
            "vp1_has_frameshift": False,
            "vp1_has_misphased_gap": False,
            "vp1_has_internal_stop": False,
        }

    span_start = material_idx[0]
    span_end = material_idx[-1]

    # Locate contiguous internal deletion runs within [span_start, span_end]
    internal_runs: list[tuple[int, int]] = []
    i = span_start
    while i <= span_end:
        if s[i] == "-":
            j = i
            while j <= span_end and s[j] == "-":
                j += 1
            internal_runs.append((i, j - i))
            i = j
        else:
            i += 1

    internal_del_nt = sum(run_len for _, run_len in internal_runs)
    internal_del_events = len(internal_runs)
    has_frameshift = any(run_len % 3 != 0 for _, run_len in internal_runs)
    has_misphased_gap = any(start % 3 != 0 or run_len % 3 != 0 for start, run_len in internal_runs)

    # Codon-level internal deletion & internal stop checks
    internal_del_codons = 0
    has_internal_stop = False
    for pos in range(0, width, 3):
        codon = s[pos : pos + 3]
        in_span = pos >= span_start and (pos + 2) <= span_end
        if "-" in codon:
            if in_span or (pos <= span_end and (pos + 2) >= span_start and any(
                span_start <= pos + k <= span_end and codon[k] == "-" for k in range(3)
            )):
                internal_del_codons += 1
        elif all(ch in BASES for ch in codon):
            if CODON_TABLE[codon] == "*":
                has_internal_stop = True

    return {
        "seq": s,
        "vp1_non_gap_nt": non_gap_nt,
        "vp1_acgt_nt": acgt_nt,
        "vp1_ambig_nt": ambig_nt,
        "vp1_coverage_frac": round(non_gap_nt / width, 6),
        "vp1_span_start_nt": span_start + 1,
        "vp1_span_end_nt": span_end + 1,
        "vp1_internal_del_nt": internal_del_nt,
        "vp1_internal_del_events": internal_del_events,
        "vp1_internal_del_codons": internal_del_codons,
        "vp1_has_frameshift": has_frameshift,
        "vp1_has_misphased_gap": has_misphased_gap,
        "vp1_has_internal_stop": has_internal_stop,
    }


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--curated-dir",
        type=Path,
        default=repo_root.parent / "enterovirus-genbank-curated",
        help="Path to the enterovirus-genbank-curated checkout (default: sibling repo).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo_root / "data" / "curated",
        help="Destination directory inside iVDPV-vs-cVDPV (default: data/curated).",
    )
    parser.add_argument(
        "--min-vp1-nt",
        type=int,
        default=1,
        help="Minimum non-gap VP1 nucleotides required to retain a sequence (default: 1).",
    )
    args = parser.parse_args()

    final_dir = args.curated_dir / "final"
    if not final_dir.is_dir():
        raise SystemExit(f"Curated final/ directory not found: {final_dir}")

    build_manifest = json.loads((final_dir / "audit" / "build_manifest.json").read_text(encoding="utf-8"))
    coords = read_vp1_coordinates(final_dir / "alignments" / "reference_region_coordinates.tsv")

    meta_path = final_dir / "canonical" / "sequence_metadata_vouched.tsv.gz"
    with gzip.open(meta_path, "rt", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        canonical_cols = list(reader.fieldnames or [])
        vouched_rows = list(reader)

    meta_by_acc = {row["accession"]: row for row in vouched_rows}
    all_accs = set(meta_by_acc)
    refs_by_acc = read_first_references(final_dir / "source" / "normalized_tsv" / "references.tsv.gz", all_accs)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    extracted_meta_rows: list[dict[str, object]] = []
    alignment_provenance: dict[str, dict[str, object]] = {}
    output_hashes: dict[str, str] = {}
    summary_counts: dict[str, dict[str, int]] = {}

    for serotype, (expected_ref, expected_len) in SABIN_REFS.items():
        ref_acc, vp1_start, vp1_end, vp1_len = coords[serotype]
        if ref_acc != expected_ref or vp1_len != expected_len:
            raise ValueError(
                f"{serotype} coordinate mismatch: got ({ref_acc}, {vp1_len}), "
                f"expected ({expected_ref}, {expected_len})"
            )

        aln_path = final_dir / "alignments" / f"{serotype}_unified_aln.fasta.gz"
        prov_path = final_dir / "alignments" / f"{serotype}_unified.provenance.json"
        prov = json.loads(prov_path.read_text(encoding="utf-8"))
        seqs = read_fasta_gz(aln_path)

        if ref_acc not in seqs:
            raise ValueError(f"Sabin reference {ref_acc} missing from {aln_path}")

        ref_full = seqs[ref_acc].upper().replace("U", "T")
        ref_nongap_cols = [idx for idx, ch in enumerate(ref_full) if ch != "-"]
        if len(ref_nongap_cols) != len(ref_full):
            raise ValueError(
                f"{serotype} reference {ref_acc} contains gaps in Sabin-anchored alignment {aln_path}"
            )

        vp1_cols = ref_nongap_cols[vp1_start - 1 : vp1_end]
        if len(vp1_cols) != expected_len:
            raise ValueError(f"{serotype}: sliced {len(vp1_cols)} VP1 columns, expected {expected_len}")
        col_lo, col_hi = vp1_cols[0], vp1_cols[-1] + 1
        assert col_hi - col_lo == expected_len

        # Collect eligible non-engineered vouched records, placing Sabin reference first
        kept_ids: list[str] = []
        skipped_engineered = 0
        skipped_no_vp1 = 0
        skipped_missing_meta = 0
        qc_by_acc: dict[str, dict[str, object]] = {}

        for acc in sorted(seqs):
            row = meta_by_acc.get(acc)
            if row is None:
                skipped_missing_meta += 1
                continue
            if row["engineered_or_construct"] == "TRUE":
                skipped_engineered += 1
                continue
            vp1_raw = seqs[acc][col_lo:col_hi]
            qc = analyze_vp1_slice(vp1_raw)
            if acc != ref_acc and int(qc["vp1_non_gap_nt"]) < args.min_vp1_nt:
                skipped_no_vp1 += 1
                continue
            qc_by_acc[acc] = qc
            if acc != ref_acc:
                kept_ids.append(acc)

        if ref_acc not in qc_by_acc:
            raise ValueError(f"Sabin reference {ref_acc} was excluded from {serotype}!")
        if qc_by_acc[ref_acc]["vp1_non_gap_nt"] != expected_len or qc_by_acc[ref_acc]["vp1_ambig_nt"] != 0:
            raise ValueError(f"Sabin reference {ref_acc} VP1 is not 100% unambiguous ACGT")

        ordered_ids = [ref_acc] + kept_ids

        # Write Sabin-first VP1 FASTA
        fasta_lines: list[str] = []
        for acc in ordered_ids:
            fasta_lines.append(f">{acc}")
            fasta_lines.append(str(qc_by_acc[acc]["seq"]))
        fasta_bytes = ("\n".join(fasta_lines) + "\n").encode("utf-8")
        out_fasta = args.output_dir / f"{serotype}_vp1.fasta.gz"
        write_gzip_deterministic(out_fasta, fasta_bytes)
        output_hashes[out_fasta.name] = sha256_file(out_fasta)

        # Build metadata rows in the exact same order
        cls_counts: collections.Counter[str] = collections.Counter()
        cov95_counts: collections.Counter[str] = collections.Counter()
        indel95_counts: collections.Counter[str] = collections.Counter()
        fs_counts: collections.Counter[str] = collections.Counter()

        for acc in ordered_ids:
            row = dict(meta_by_acc[acc])
            qc = qc_by_acc[acc]
            pmid, title = refs_by_acc.get(acc, ("", ""))
            out_row: dict[str, object] = {col: row[col] for col in canonical_cols}
            out_row.update(
                {
                    "pubmed_id": pmid,
                    "reference_title": title,
                    "vp1_ref_accession": ref_acc,
                    "vp1_ref_nt": expected_len,
                    "is_sabin_reference": "TRUE" if acc == ref_acc else "FALSE",
                    "vp1_non_gap_nt": qc["vp1_non_gap_nt"],
                    "vp1_acgt_nt": qc["vp1_acgt_nt"],
                    "vp1_ambig_nt": qc["vp1_ambig_nt"],
                    "vp1_coverage_frac": qc["vp1_coverage_frac"],
                    "vp1_span_start_nt": qc["vp1_span_start_nt"],
                    "vp1_span_end_nt": qc["vp1_span_end_nt"],
                    "vp1_internal_del_nt": qc["vp1_internal_del_nt"],
                    "vp1_internal_del_events": qc["vp1_internal_del_events"],
                    "vp1_internal_del_codons": qc["vp1_internal_del_codons"],
                    "vp1_has_frameshift": "TRUE" if qc["vp1_has_frameshift"] else "FALSE",
                    "vp1_has_misphased_gap": "TRUE" if qc["vp1_has_misphased_gap"] else "FALSE",
                    "vp1_has_internal_stop": "TRUE" if qc["vp1_has_internal_stop"] else "FALSE",
                }
            )
            extracted_meta_rows.append(out_row)

            cls = row["poliovirus_classification"] or "<blank>"
            cls_counts[cls] += 1
            if float(qc["vp1_coverage_frac"]) >= 0.95:
                cov95_counts[cls] += 1
                if int(qc["vp1_internal_del_nt"]) > 0:
                    indel95_counts[cls] += 1
            if qc["vp1_has_frameshift"]:
                fs_counts[cls] += 1

        alignment_provenance[serotype] = {
            "source_alignment": f"final/alignments/{serotype}_unified_aln.fasta.gz",
            "reference_accession": ref_acc,
            "vp1_start_1based": vp1_start,
            "vp1_end_1based": vp1_end,
            "vp1_length_nt": expected_len,
            "upstream_rows": prov.get("rows"),
            "upstream_population_sha256": prov.get("population_sha256"),
            "upstream_fasta_sha256": prov.get("artifact_sha256", {}).get("fasta"),
            "extracted_rows": len(ordered_ids),
            "skipped_engineered": skipped_engineered,
            "skipped_below_min_vp1_nt": skipped_no_vp1,
            "skipped_missing_meta": skipped_missing_meta,
        }
        summary_counts[serotype] = {
            "extracted_total": len(ordered_ids),
            "cov_ge_95_by_class": dict(sorted(cov95_counts.items())),
            "cov_ge_95_with_internal_del_by_class": dict(sorted(indel95_counts.items())),
            "frameshift_by_class": dict(sorted(fs_counts.items())),
        }

    # Write combined metadata TSV
    extra_cols = [
        "pubmed_id",
        "reference_title",
        "vp1_ref_accession",
        "vp1_ref_nt",
        "is_sabin_reference",
        "vp1_non_gap_nt",
        "vp1_acgt_nt",
        "vp1_ambig_nt",
        "vp1_coverage_frac",
        "vp1_span_start_nt",
        "vp1_span_end_nt",
        "vp1_internal_del_nt",
        "vp1_internal_del_events",
        "vp1_internal_del_codons",
        "vp1_has_frameshift",
        "vp1_has_misphased_gap",
        "vp1_has_internal_stop",
    ]
    all_out_cols = canonical_cols + extra_cols
    meta_lines: list[str] = ["\t".join(all_out_cols)]
    for row in extracted_meta_rows:
        meta_lines.append("\t".join(str(row[c]).replace("\t", " ").replace("\n", " ") for c in all_out_cols))
    meta_bytes = ("\n".join(meta_lines) + "\n").encode("utf-8")
    out_meta = args.output_dir / "poliovirus_vp1_metadata.tsv.gz"
    write_gzip_deterministic(out_meta, meta_bytes)
    output_hashes[out_meta.name] = sha256_file(out_meta)

    manifest = {
        "source_dataset": "enterovirus-genbank-curated",
        "release_version": build_manifest.get("release_version"),
        "schema_version": build_manifest.get("schema_version"),
        "upstream_code_sha256": build_manifest.get("code_sha256"),
        "upstream_inputs": build_manifest.get("inputs"),
        "source_files": {
            "sequence_metadata_vouched.tsv.gz": sha256_file(meta_path),
            "reference_region_coordinates.tsv": sha256_file(
                final_dir / "alignments" / "reference_region_coordinates.tsv"
            ),
            "build_manifest.json": sha256_file(final_dir / "audit" / "build_manifest.json"),
        },
        "extraction_parameters": {
            "filter_virus_group": "poliovirus",
            "filter_curation_status": "vouched",
            "exclude_engineered_or_construct": True,
            "min_vp1_nt": args.min_vp1_nt,
        },
        "alignments": alignment_provenance,
        "outputs_sha256": output_hashes,
        "summary": summary_counts,
    }
    manifest_path = args.output_dir / "source_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"Extracted enterovirus-genbank-curated release {manifest['release_version']} -> {args.output_dir}")
    for st in ("PV1", "PV2", "PV3"):
        ap = alignment_provenance[st]
        sc = summary_counts[st]
        print(
            f"  {st} (ref={ap['reference_accession']}, {ap['vp1_length_nt']} nt): "
            f"{ap['extracted_rows']} extracted "
            f"(skipped: eng={ap['skipped_engineered']}, no_vp1={ap['skipped_below_min_vp1_nt']})"
        )
        print(f"    cov>=95%: {sc['cov_ge_95_by_class']}")
        print(f"    cov>=95% internal del: {sc['cov_ge_95_with_internal_del_by_class']}")
        print(f"    frameshifts (all cov): {sc['frameshift_by_class']}")


if __name__ == "__main__":
    main()
