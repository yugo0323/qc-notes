"""講義ノートの「参考文献」を機械的に確かめる。

参考文献の各項目は、次の形で書く（1 行 1 件）：
  - 著者, "題名", 雑誌名 巻, ページ (年). DOI: [10.xxxx/yyyy](https://doi.org/10.xxxx/yyyy)
  - 著者, "題名" (出版社, 年). ISBN: 978-xxxxxxxxxx
  - 著者, "題名", arXiv:2506.12345 (年).

確かめること：
  1. DOI は Crossref、ISBN は Open Library、arXiv 番号は arXiv に問い合わせて、実在すること。
  2. 「"」で囲んだ題名が、正式な題名と一致すること（大文字小文字・記号は無視。本の副題は省略してよい）。
  3. (年) が正式な年と一致すること（DOI と arXiv のみ。本は版によって年が違うので確かめない）。
  4. 本文に出てくる DOI がすべて参考文献に載っていること。

使い方：python scripts/check_refs.py docs/courses/qc-theory/04-second-quantization.md [...]
"""

import difflib
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

USER_AGENT = "qclab-refs/0.1 (https://github.com/yugo0323/qc-notes)"
_cache = {}


def get(url):
    if url in _cache:
        return _cache[url]
    data = None
    for i in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
                break
        except Exception:
            time.sleep(4 * (i + 1))
    _cache[url] = data
    return data


def norm(t):
    t = re.sub(r"<[^>]+>", " ", t or "")
    t = t.replace("‐", "-").replace("–", "-").replace("—", "-")
    return re.sub(r"[^a-z0-9]+", " ", t.lower()).strip()


def lookup_doi(doi):
    data = get(f"https://api.crossref.org/works/{urllib.parse.quote(doi)}")
    if not data:
        return None
    m = json.loads(data)["message"]
    title = " ".join(m.get("title") or [])
    sub = " ".join(m.get("subtitle") or [])
    year = None
    for key in ("published-print", "published-online", "issued", "published"):
        parts = (m.get(key) or {}).get("date-parts") or [[None]]
        if parts[0][0]:
            year = parts[0][0]
            break
    years = {
        (m.get(k) or {}).get("date-parts", [[None]])[0][0]
        for k in ("published-print", "published-online", "issued", "published")
    } - {None}
    return {"title": title, "subtitle": sub, "year": year, "years": years}


def lookup_isbn(isbn):
    data = get(f"https://openlibrary.org/isbn/{isbn}.json")
    if not data:
        return None
    try:
        d = json.loads(data)
    except json.JSONDecodeError:
        return None
    return {"title": d.get("title", ""), "subtitle": d.get("subtitle", ""), "year": None, "years": set()}


def lookup_arxiv(aid):
    data = get("https://export.arxiv.org/api/query?" + urllib.parse.urlencode({"id_list": aid}))
    if not data:
        return None
    ns = {"a": "http://www.w3.org/2005/Atom"}
    e = ET.fromstring(data).find("a:entry", ns)
    if e is None or e.findtext("a:title", None, ns) is None:
        return None
    year = int(e.findtext("a:published", "0000", ns)[:4])
    return {"title": e.findtext("a:title", "", ns), "subtitle": "", "year": year, "years": {year}}


def title_ok(written, official, subtitle, is_book):
    w, o = norm(written), norm(official)
    if not w:
        return False, 0.0
    candidates = [o, norm(f"{official} {subtitle}"), norm(f"{official}: {subtitle}")]
    ratio = max(difflib.SequenceMatcher(None, w, c).ratio() for c in candidates)
    if is_book and (o.startswith(w) or w.startswith(o)):
        return True, 1.0
    return ratio >= 0.9, ratio


def check(path):
    text = Path(path).read_text(encoding="utf-8")
    errors, notes = [], []
    m = re.search(r"^## 参考文献\s*$(.*?)(?=^## |\Z)", text, re.M | re.S)
    if not m:
        return [f"{path}: 「## 参考文献」の節がありません"], notes
    refs = [l.strip() for l in m.group(1).splitlines() if l.strip().startswith("- ")]
    if not refs:
        errors.append("参考文献が 1 件もありません")
    listed_dois = set()
    for ref in refs:
        title_m = re.search(r'"([^"]+)"', ref)
        written = title_m.group(1) if title_m else ""
        year_m = re.findall(r"\((?:[^()]*?,\s*)?(\d{4})\)", ref)
        written_year = int(year_m[-1]) if year_m else None
        doi_m = re.search(r"DOI:\s*\[(10\.\d{4,9}/[^\]\s]+)\]", ref) or re.search(r"DOI:\s*(10\.\d{4,9}/[^\s\]\)]+)", ref)
        isbn_m = re.search(r"ISBN:\s*([0-9Xx\-]{10,17})", ref)
        arx_m = re.search(r"arXiv:(\d{4}\.\d{4,5})", ref)
        if doi_m:
            key, info, is_book = doi_m.group(1).rstrip("."), lookup_doi(doi_m.group(1).rstrip(".")), False
            listed_dois.add(key.lower())
        elif isbn_m:
            key, info, is_book = isbn_m.group(1), lookup_isbn(isbn_m.group(1).replace("-", "")), True
        elif arx_m:
            key, info, is_book = arx_m.group(1), lookup_arxiv(arx_m.group(1)), False
        else:
            errors.append(f"DOI・ISBN・arXiv 番号のどれもありません：{ref[:80]}")
            continue
        if not info:
            errors.append(f"実在を確認できません：{key}")
            continue
        if not written:
            errors.append(f"題名が「\"」で囲まれていません：{key}")
            continue
        n_before = len(errors)
        ok, ratio = title_ok(written, info["title"], info["subtitle"], is_book)
        if not ok:
            errors.append(f"題名が一致しません（一致度 {ratio:.2f}）：{key}\n    書いた題名：{written}\n    正式な題名：{info['title']} {info['subtitle']}".rstrip())
        if not is_book and info["years"] and written_year not in info["years"]:
            errors.append(f"年が一致しません：{key}（書いた年 {written_year}、正式 {sorted(info['years'])}）")
        if len(errors) == n_before:
            notes.append(f"ok {key}: {info['title'][:70]}")
        time.sleep(0.5)

    body = text[: m.start()]
    for doi in re.findall(r"10\.\d{4,9}/[^\s\)\]>\"'、。，]+", body):
        if doi.rstrip(".,;").lower() not in listed_dois:
            errors.append(f"本文の DOI が参考文献にありません：{doi}")
    return [f"{path}: {e}" for e in errors], notes


if __name__ == "__main__":
    all_errors = []
    for p in sys.argv[1:]:
        errs, notes = check(p)
        all_errors += errs
        for n in notes:
            print(f"{p}: {n}")
    if all_errors:
        print("\n".join(all_errors))
        sys.exit(1)
    print(f"OK: {len(sys.argv) - 1} 件のファイルの参考文献を確認しました")
