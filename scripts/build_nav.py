"""mkdocs.yml の nav（サイドメニュー）を、講義ノートのファイルから作り直す。

講義ノートは docs/courses/<連載>/NN-<名前>.md に置き、ファイルの最初の「# 」の行を題名に使う。
mkdocs.yml の最後にある nav: 以下を、まるごと書き換える。
使い方：python scripts/build_nav.py
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
COURSES = [
    ("量子化学の理論", "qc-theory"),
    ("遷移金属錯体の計算", "tm-complexes"),
    ("機械学習の基礎", "ml-basics"),
    ("量子化学のための機械学習", "ml-chem"),
    ("計算パッケージの使い方", "packages"),
    ("研究の道具", "tools"),
]


def q(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def main():
    lines = [
        "nav:",
        "  - ホーム: index.md",
        "  - 論文紹介: papers/index.md",
        "  - 講義ノート:",
        "      - courses/index.md",
    ]
    for name, slug in COURSES:
        files = sorted((ROOT / "docs" / "courses" / slug).glob("[0-9][0-9]-*.md"))
        if not files:
            lines.append(f"      - {name}: courses/{slug}/index.md")
            continue
        lines.append(f"      - {name}:")
        lines.append(f"          - courses/{slug}/index.md")
        for f in files:
            m = re.search(r"^# (.+)$", f.read_text(encoding="utf-8"), re.M)
            title = m.group(1).strip() if m else f.stem
            lines.append(f"          - {q(title)}: courses/{slug}/{f.name}")
    lines += ["  - 用語集: glossary.md", "  - このサイトについて: about.md", ""]

    cfg = ROOT / "mkdocs.yml"
    text = cfg.read_text(encoding="utf-8")
    text = text[: text.index("\nnav:") + 1] + "\n".join(lines)
    cfg.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
