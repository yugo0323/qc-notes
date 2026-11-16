"""論文紹介の一覧ページと、トップページの PICKUP・NEWS を、記事の front matter から作り直す。

使い方：python scripts/build_index.py
"""

import datetime as dt
import html
import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
PAPERS = ROOT / "docs" / "papers"
HOME = ROOT / "docs" / "index.md"
INDEX = PAPERS / "index.md"


def load_posts():
    posts = []
    for p in PAPERS.rglob("*.md"):
        if p == INDEX:
            continue
        m = re.match(r"^---\n(.*?)\n---\n", p.read_text(encoding="utf-8"), re.S)
        if not m:
            continue
        meta = yaml.safe_load(m.group(1)) or {}
        if not meta.get("date"):
            continue
        meta["date"] = dt.date.fromisoformat(str(meta["date"]))
        meta["rel"] = p.relative_to(PAPERS).as_posix()
        posts.append(meta)
    posts.sort(key=lambda x: (x["date"], x["rel"]), reverse=True)
    return posts


def main_paper(post):
    return next((p for p in post.get("papers", []) if p.get("role") == "main"), {})


def replace_block(text, name, new):
    pat = re.compile(rf"(<!-- {name}:START -->\n).*?(\n<!-- {name}:END -->)", re.S)
    return pat.sub(lambda m: m.group(1) + new + m.group(2), text)


def build_index(posts):
    lines = [
        "# 論文紹介",
        "",
        "電子構造理論、遷移金属錯体の計算、量子化学のための機械学習に関する新しい論文を、毎朝紹介します。",
        "",
        '!!! note "この紹介について"',
        "    紹介文は AI（Claude）が論文の要旨や本文を読んで自動で作成し、人による確認を経ずに公開しています。"
        "紹介した論文の書誌情報（題名・DOI など）は、Crossref または arXiv で実在を機械的に確かめています。"
        "正確な内容は、必ず原論文で確認してください。",
        "",
        "## 一覧",
        "",
    ]
    if not posts:
        lines.append("（まだありません）")
    month = None
    for post in posts:
        ym = post["date"].strftime("%Y 年 %m 月")
        if ym != month:
            lines += ["", f"### {ym}", ""]
            month = ym
        mp = main_paper(post)
        tags = " ".join(f"`{t}`" for t in post.get("tags", []))
        lines.append(f"- {post['date'].strftime('%m/%d')}　[{post['title']}]({post['rel']})　{tags}")
        if mp.get("venue"):
            lines.append(f"    - {mp.get('venue')}")
    INDEX.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_home(posts):
    text = HOME.read_text(encoding="utf-8")
    if posts:
        p = posts[0]
        mp = main_paper(p)
        link = "papers/" + p["rel"][:-3] + "/"
        pickup = (
            '<div class="qc-pickup">\n'
            '<span class="qc-chip qc-chip--paper">論文紹介</span>\n'
            f'<p class="qc-pickup__title"><a href="{link}">{html.escape(p["title"])}</a></p>\n'
            f'<p class="qc-pickup__meta">{p["date"].strftime("%Y.%m.%d")}　{html.escape(str(mp.get("venue", "")))}</p>\n'
            "</div>"
        )
        text = replace_block(text, "PICKUP", pickup)

    # NEWS：論文紹介の新しいもの 6 件と、手で書いたお知らせ（qc-chip--info / course の行）を日付順に混ぜる
    m = re.search(r"<!-- NEWS:START -->\n(.*?)\n<!-- NEWS:END -->", text, re.S)
    manual = [
        line for line in (m.group(1).splitlines() if m else [])
        if line.startswith("<li>") and "qc-chip--paper" not in line
    ]
    rows = []
    for p in posts[:6]:
        link = "papers/" + p["rel"][:-3] + "/"
        rows.append((
            p["date"].strftime("%Y.%m.%d"),
            f'<li><span class="qc-news__date">{p["date"].strftime("%Y.%m.%d")}</span>'
            f'<span class="qc-chip qc-chip--paper">論文紹介</span>'
            f'<a href="{link}">{html.escape(p["title"])}</a></li>',
        ))
    for line in manual:
        d = re.search(r"\d{4}\.\d{2}\.\d{2}", line)
        rows.append((d.group(0) if d else "0000.00.00", line))
    rows.sort(key=lambda r: r[0], reverse=True)
    news = '<ul class="qc-news">\n' + "\n".join(r[1] for r in rows[:8]) + "\n</ul>"
    text = replace_block(text, "NEWS", news)
    HOME.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    posts = load_posts()
    build_index(posts)
    build_home(posts)
    print(f"index: {len(posts)} posts")
