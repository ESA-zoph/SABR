from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time
import urllib.parse
import urllib.request
import json

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Find NCBI nucleotide candidate accessions for unresolved host strains."
    )
    parser.add_argument(
        "--interactions",
        type=Path,
        default=Path("data/curation/phage_host_interactions.tsv"),
    )
    parser.add_argument(
        "--coverage",
        type=Path,
        default=Path("data/curation/accession_linkage_coverage.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/curation/host_accession_candidates.tsv"),
    )
    parser.add_argument("--max-strains", type=int, default=80)
    parser.add_argument("--retmax", type=int, default=3)
    parser.add_argument("--sleep", type=float, default=0.35)
    args = parser.parse_args()

    interactions = pd.read_csv(args.interactions, sep="\t", dtype=str).fillna("")
    coverage = pd.read_csv(args.coverage, sep="\t", dtype=str).fillna("")
    unresolved_ids = set(
        coverage.loc[~coverage["bacterium_genome_linked"].map(_truthy), "interaction_id"].astype(str)
    )
    unresolved = interactions[interactions["interaction_id"].astype(str).isin(unresolved_ids)]
    strains = (
        unresolved[["source_key", "bacterium", "strain"]]
        .drop_duplicates()
        .head(args.max_strains)
    )
    rows = []
    for _, row in strains.iterrows():
        query = _query_for(row["bacterium"], row["strain"])
        ids = _esearch(query, retmax=args.retmax)
        if not ids:
            rows.append(_candidate_row(row, query=query, rank=0))
            time.sleep(args.sleep)
            continue
        summaries = _esummary(ids)
        for rank, uid in enumerate(ids, start=1):
            summary = summaries.get(uid, {})
            rows.append(
                _candidate_row(
                    row,
                    query=query,
                    rank=rank,
                    uid=uid,
                    accession=summary.get("caption", ""),
                    title=summary.get("title", ""),
                    length=summary.get("slen", ""),
                )
            )
        time.sleep(args.sleep)
    output = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False)
    print(f"queried_strains\t{len(strains)}")
    print(f"candidate_rows\t{len(output)}")
    print(f"with_candidates\t{int((output['candidate_rank'] > 0).sum())}")


def _query_for(bacterium: str, strain: str) -> str:
    clean_strain = str(strain).replace("_", " ")
    return f'"{bacterium}"[Organism] AND "{clean_strain}"[All Fields]'


def _esearch(query: str, retmax: int) -> list[str]:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(
        {"db": "nuccore", "term": query, "retmode": "json", "retmax": str(retmax)}
    )
    data = json.loads(urllib.request.urlopen(url, timeout=30).read().decode())
    return data["esearchresult"].get("idlist", [])


def _esummary(ids: list[str]) -> dict[str, dict[str, object]]:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(
        {"db": "nuccore", "id": ",".join(ids), "retmode": "json"}
    )
    data = json.loads(urllib.request.urlopen(url, timeout=30).read().decode())["result"]
    return {uid: data[uid] for uid in data.get("uids", [])}


def _candidate_row(
    row: pd.Series,
    query: str,
    rank: int,
    uid: str = "",
    accession: str = "",
    title: str = "",
    length: object = "",
) -> dict[str, object]:
    return {
        "source_key": row["source_key"],
        "bacterium": row["bacterium"],
        "strain": row["strain"],
        "query": query,
        "candidate_rank": rank,
        "ncbi_uid": uid,
        "candidate_accession": accession,
        "candidate_title": title,
        "candidate_length": length,
        "review_status": "needs_review" if rank else "no_candidate",
        "notes": "",
    }


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


if __name__ == "__main__":
    main()
