"""新しく出た論文の候補を集めて、関連度の高い順に JSON に書き出す。

集め先：
  - arXiv（API）：physics.chem-ph ほか
  - Crossref（API）：主な雑誌、ChemRxiv（DOI の接頭辞 10.26434）

すでに紹介した論文（docs/papers/ の front matter にある id）は除く。
使い方：python scripts/fetch_candidates.py --out _work/candidates.json
"""

import argparse
import datetime as dt
import html
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
USER_AGENT = "qclab-daily/0.1 (https://github.com/yugo0323/qc-notes)"

JOURNALS = {
    "J. Chem. Theory Comput.": "1549-9618",
    "J. Chem. Phys.": "0021-9606",
    "J. Phys. Chem. A": "1089-5639",
    "J. Phys. Chem. Lett.": "1948-7185",
    "J. Chem. Inf. Model.": "1549-9596",
    "J. Am. Chem. Soc.": "0002-7863",
    "Inorg. Chem.": "0020-1669",
    "ACS Catal.": "2155-5435",
    "Organometallics": "0276-7333",
    "Chem. Sci.": "2041-6520",
    "Phys. Chem. Chem. Phys.": "1463-9076",
    "Dalton Trans.": "1477-9226",
    "npj Comput. Mater.": "2057-3960",
}
CHEMRXIV_PREFIX = "10.26434"

ARXIV_QUERIES = [
    "cat:physics.chem-ph",
    'cat:physics.comp-ph AND (abs:"quantum chemistry" OR abs:"density functional" OR abs:"electronic structure")',
    '(cat:cs.LG OR cat:cond-mat.mtrl-sci) AND (abs:"interatomic potential" OR abs:"quantum chemistry" OR abs:"density functional" OR abs:"molecular property")',
]

# (正規表現, 重み)。題名と要旨を小文字にして調べる。
KEYWORDS = [
    (r"spin[- ]flip|mixed[- ]reference|\bmrsf\b|spin[- ]complete|spin[- ]adapted", 6),
    (r"multi[- ]?reference|multiconfiguration|\bcasscf\b|\bcaspt2\b|\bnevpt2\b|\bdmrg\b|selected configuration interaction", 5),
    (r"spin contamination|spin[- ]state|spin crossover|spin[- ]crossing|two[- ]state reactivity|minimum energy crossing", 5),
    (r"hydrogen atom transfer|metal hydride|\bmhat\b|bond dissociation (free )?energ", 5),
    (r"transition[- ]metal|organometallic|\b(iron|cobalt|manganese|chromium|ruthenium|palladium|nickel)\b", 3),
    (r"\bdlpno\b|coupled[- ]cluster|ccsd\(t\)|local correlation", 4),
    (r"time[- ]dependent density functional|\btddft\b|linear response|conical intersection", 3),
    (r"exchange[- ]correlation functional|density functional approximation|double[- ]hybrid|range[- ]separated", 3),
    (r"machine[- ]learn\w* (interatomic )?potential|neural network potential|interatomic potential|foundation model|universal model|equivariant|\bomol25\b|\buma\b|\bmace\b|\bnequip\b", 4),
    (r"density functional|\bdft\b|ab initio|electronic structure|quantum chemi", 2),
    (r"catalys|catalyst", 1),
    (r"benchmark", 1),
]
# 計算・理論の論文であることの目印。どれも含まない雑誌論文は除く（実験だけの論文を除くため）。
COMPUTATIONAL = re.compile(
    r"density functional|\bdft\b|ab initio|quantum chemi|electronic structure|calculation|computational|"
    r"simulation|machine learn|neural network|theory|theoretical|coupled[- ]cluster|wave ?function"
)
NEGATIVE = [(r"\bmc-?pdft\b|pair[- ]density functional", -6)]


def get(url, retries=3):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except Exception as e:  # ネットワークの一時的な失敗はやり直す
            if i == retries - 1:
                print(f"warn: {url[:100]}... {e}", file=sys.stderr)
                return None
            time.sleep(5 * (i + 1))


def clean(text):
    text = re.sub(r"<[^>]+>", " ", text or "")
    return re.sub(r"\s+", " ", html.unescape(text)).strip()


def score(title, abstract):
    t = f"{title} {abstract}".lower()
    s = sum(w for pat, w in KEYWORDS + NEGATIVE if re.search(pat, t))
    return s


def posted_ids():
    ids = set()
    for p in (ROOT / "docs" / "papers").rglob("*.md"):
        text = p.read_text(encoding="utf-8")
        m = re.match(r"^---\n(.*?)\n---\n", text, re.S)
        if not m:
            continue
        meta = yaml.safe_load(m.group(1)) or {}
        for paper in meta.get("papers", []) or []:
            ids.add(str(paper.get("id", "")).lower())
    return ids


def from_crossref(items, venue_default):
    out = []
    for it in items:
        title = clean(" ".join(it.get("title") or []))
        if not title:
            continue
        authors = [
            " ".join(x for x in (a.get("given"), a.get("family")) if x)
            for a in it.get("author", [])
        ]
        created = it.get("created", {}).get("date-time", "")[:10]
        venue = clean(" ".join(it.get("container-title") or [])) or venue_default
        licenses = [l.get("URL", "") for l in it.get("license", []) or []]
        out.append({
            "id": it["DOI"].lower(),
            "title": title,
            "authors": authors,
            "venue": venue,
            "date": created,
            "abstract": clean(it.get("abstract", "")),
            "url": f"https://doi.org/{it['DOI']}",
            "licenses": licenses,
            "full_text_open": any("creativecommons" in l for l in licenses),
        })
    return out


def crossref_journal(issn, since, name):
    q = urllib.parse.urlencode({
        "filter": f"from-created-date:{since},type:journal-article",
        "rows": 200,
        "select": "DOI,title,author,abstract,created,container-title,license",
    })
    data = get(f"https://api.crossref.org/journals/{issn}/works?{q}")
    if not data:
        return []
    return from_crossref(json.loads(data)["message"]["items"], name)


def crossref_chemrxiv(since):
    q = urllib.parse.urlencode({
        "filter": f"from-created-date:{since},type:posted-content",
        "rows": 300,
        "select": "DOI,title,author,abstract,created,container-title,license",
    })
    data = get(f"https://api.crossref.org/prefixes/{CHEMRXIV_PREFIX}/works?{q}")
    if not data:
        return []
    items = from_crossref(json.loads(data)["message"]["items"], "ChemRxiv")
    for x in items:
        x["venue"] = "ChemRxiv（プレプリント、査読前）"
        x["full_text_open"] = True
    return items


def arxiv(query, since):
    q = urllib.parse.urlencode({
        "search_query": query,
        "sortBy": "submittedDate",
        "sortOrder": "descending",
        "max_results": 200,
    })
    data = get(f"https://export.arxiv.org/api/query?{q}")
    if not data:
        return []
    ns = {"a": "http://www.w3.org/2005/Atom"}
    out = []
    for e in ET.fromstring(data).findall("a:entry", ns):
        published = e.findtext("a:published", "", ns)[:10]
        if published < since:
            continue
        aid = re.sub(r"v\d+$", "", e.findtext("a:id", "", ns).rsplit("/abs/", 1)[-1])
        out.append({
            "id": f"arxiv:{aid}".lower(),
            "title": clean(e.findtext("a:title", "", ns)),
            "authors": [clean(a.findtext("a:name", "", ns)) for a in e.findall("a:author", ns)],
            "venue": "arXiv（プレプリント、査読前）",
            "date": published,
            "abstract": clean(e.findtext("a:summary", "", ns)),
            "url": f"https://arxiv.org/abs/{aid}",
            "licenses": [],
            "full_text_open": True,
        })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="_work/candidates.json")
    ap.add_argument("--days", type=int, default=4)
    ap.add_argument("--top", type=int, default=30)
    args = ap.parse_args()

    since = (dt.date.today() - dt.timedelta(days=args.days)).isoformat()
    found = []
    for name, issn in JOURNALS.items():
        found += crossref_journal(issn, since, name)
        time.sleep(1)
    found += crossref_chemrxiv(since)
    for q in ARXIV_QUERIES:
        found += arxiv(q, since)
        time.sleep(3)  # arXiv API の利用規約に従い、間隔をあける

    done = posted_ids()
    seen, ranked = set(), []
    for x in found:
        key = x["id"]
        if key in seen or key in done:
            continue
        seen.add(key)
        text = f"{x['title']} {x['abstract']}".lower()
        if not x["abstract"]:
            continue  # 要旨のない論文は、内容を確かめられないので候補にしない
        if not x["id"].startswith("arxiv:") and not COMPUTATIONAL.search(text):
            continue
        x["score"] = score(x["title"], x["abstract"])
        if x["score"] >= 4:
            ranked.append(x)
    ranked.sort(key=lambda x: x["date"], reverse=True)
    ranked.sort(key=lambda x: x["score"], reverse=True)
    ranked = ranked[: args.top]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"since": since, "candidates": ranked}, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"collected {len(found)} records, {len(ranked)} candidates -> {out}")
    for x in ranked[:10]:
        print(f"  [{x['score']:2d}] {x['venue'][:24]:24s} {x['title'][:80]}")


if __name__ == "__main__":
    main()
