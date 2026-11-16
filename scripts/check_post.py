"""論文紹介の記事を、公開の前に機械的に確かめる。

確かめること：
  1. front matter（記事の先頭の --- で囲まれた部分）に、必要な項目がある。
  2. 紹介した論文が実在し、題名が一致する（Crossref または arXiv に問い合わせる）。
  3. 「この紹介について」の注意書きと、必要な見出しがある。
  4. 本文に出てくる DOI と arXiv 番号が、すべて front matter に載っている（確かめていない文献を書かない）。

使い方：python scripts/check_post.py docs/papers/2026/2026-11-18.md
問題があれば内容を表示して、終了コード 1 で終わる。
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

import yaml

USER_AGENT = "qclab-daily/0.1 (https://github.com/yugo0323/qc-notes)"
REQUIRED_HEADINGS = ["## 3 行でいうと", "## 背景", "## 何をしたか", "## 結果", "## 読むのに必要な基礎", "## 批判的に読むと"]
DISCLAIMER = "人による確認を経ずに公開"


def get(url):
    for i in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except Exception:
            time.sleep(5 * (i + 1))
    return None


def norm(t):
    t = re.sub(r"<[^>]+>", " ", t or "")
    return re.sub(r"[^a-z0-9]+", " ", t.lower()).strip()


def official_title(pid):
    if pid.startswith("arxiv:"):
        q = urllib.parse.urlencode({"id_list": pid[6:]})
        data = get(f"https://export.arxiv.org/api/query?{q}")
        if not data:
            return None
        ns = {"a": "http://www.w3.org/2005/Atom"}
        e = ET.fromstring(data).find("a:entry", ns)
        return e.findtext("a:title", None, ns) if e is not None else None
    data = get(f"https://api.crossref.org/works/{urllib.parse.quote(pid)}")
    if not data:
        return None
    return " ".join(json.loads(data)["message"].get("title") or [])


def check(path):
    errors = []
    text = Path(path).read_text(encoding="utf-8")
    m = re.match(r"^---\n(.*?)\n---\n(.*)$", text, re.S)
    if not m:
        return [f"{path}: front matter がありません"]
    meta, body = yaml.safe_load(m.group(1)) or {}, m.group(2)

    for key in ("title", "date", "papers", "tags"):
        if not meta.get(key):
            errors.append(f"front matter に {key} がありません")
    papers = meta.get("papers") or []
    if sum(1 for p in papers if p.get("role") == "main") != 1:
        errors.append("role: main の論文がちょうど 1 本ではありません")
    if len(papers) > 4:
        errors.append("論文が 4 本を超えています")

    ids = set()
    for p in papers:
        pid = str(p.get("id", "")).lower().strip()
        ids.add(pid)
        if not (pid.startswith("10.") or re.match(r"arxiv:\d{4}\.\d{4,5}$", pid)):
            errors.append(f"id の形が正しくありません：{pid}")
            continue
        got = official_title(pid)
        if not got:
            errors.append(f"実在を確認できません：{pid}")
            continue
        ratio = difflib.SequenceMatcher(None, norm(got), norm(p.get("title", ""))).ratio()
        if ratio < 0.97:
            errors.append(f"題名が一致しません（一致度 {ratio:.2f}）：{pid}\n  記事：{p.get('title')}\n  正式：{got}")
        time.sleep(1)

    if DISCLAIMER not in body:
        errors.append("「この紹介について」の注意書きがありません")
    for h in REQUIRED_HEADINGS:
        if h not in body:
            errors.append(f"見出し「{h}」がありません")

    for doi in re.findall(r"10\.\d{4,9}/[^\s\)\]>\"'、。]+", body):
        doi = doi.rstrip(".,;").lower()
        if doi not in ids:
            errors.append(f"front matter にない DOI が本文にあります：{doi}")
    for aid in re.findall(r"arxiv\.org/(?:abs|html|pdf)/(\d{4}\.\d{4,5})", body, re.I):
        if f"arxiv:{aid}" not in ids:
            errors.append(f"front matter にない arXiv 番号が本文にあります：{aid}")
    return [f"{path}: {e}" for e in errors]


if __name__ == "__main__":
    problems = [e for p in sys.argv[1:] for e in check(p)]
    if problems:
        print("\n".join(problems))
        sys.exit(1)
    print(f"OK: {len(sys.argv) - 1} 件の記事を確認しました")
