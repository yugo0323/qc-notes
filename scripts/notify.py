"""論文紹介や講義ノートが更新されたら、メールで知らせる。

使い方：python scripts/notify.py <前のコミット> <今のコミット>
  - 2 つのコミットの間で追加・変更された docs/papers/<年>/*.md と docs/courses/**/*.md を調べる。
  - 該当がなければ何もしない。
送信には Gmail を使う。宛先とログイン情報は、公開しないように GitHub の secret から環境変数で受け取る：
  MAIL_USERNAME（送信に使う Gmail のアドレス）、MAIL_PASSWORD（Gmail のアプリ パスワード）、MAIL_TO（宛先）
"""

import os
import re
import smtplib
import subprocess
import sys
from email.message import EmailMessage
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
SITE = "https://yugo0323.github.io/qc-notes/"


def changed_files(before, after):
    out = subprocess.run(
        ["git", "diff", "--name-only", "--diff-filter=AM", before, after, "--", "docs/papers", "docs/courses"],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout.split()
    return [f for f in out if f.endswith(".md")]


def url_of(path):
    rel = Path(path).relative_to("docs").with_suffix("")
    rel = rel.parent if rel.name == "index" else rel
    return SITE + (rel.as_posix() + "/" if rel.as_posix() != "." else "")


def split_front_matter(text):
    m = re.match(r"^---\n(.*?)\n---\n(.*)$", text, re.S)
    return (yaml.safe_load(m.group(1)) or {}, m.group(2)) if m else ({}, text)


def describe_post(path):
    meta, body = split_front_matter((ROOT / path).read_text(encoding="utf-8"))
    main = next((p for p in meta.get("papers", []) if p.get("role") == "main"), {})
    three = re.search(r"## 3 行でいうと\n+(.*?)\n## ", body, re.S)
    briefs = [p for p in meta.get("papers", []) if p.get("role") == "brief"]
    lines = [
        f"■ {meta.get('title', path)}",
        f"  原題：{main.get('title', '')}",
        f"  掲載：{main.get('venue', '')}",
        f"  {url_of(path)}",
    ]
    if three:
        lines += ["", "  3 行でいうと："] + [f"  {l}" for l in three.group(1).strip().splitlines()]
    if briefs:
        lines += ["", "  短信："] + [f"  ・{b.get('title', '')}" for b in briefs]
    return meta.get("title", path), "\n".join(lines)


def describe_course(path):
    _, body = split_front_matter((ROOT / path).read_text(encoding="utf-8"))
    h1 = re.search(r"^# (.+)$", body, re.M)
    title = h1.group(1).strip() if h1 else path
    return title, f"■ {title}\n  {url_of(path)}"


def main():
    before, after = sys.argv[1], sys.argv[2]
    files = changed_files(before, after)
    posts = [f for f in files if re.match(r"docs/papers/\d{4}/", f)]
    courses = [f for f in files if f.startswith("docs/courses/")]
    if not posts and not courses:
        print("通知する更新はありません")
        return

    user, password, to = (os.environ.get(k, "") for k in ("MAIL_USERNAME", "MAIL_PASSWORD", "MAIL_TO"))
    if not (user and password and to):
        print("::warning::メールの secret（MAIL_USERNAME、MAIL_PASSWORD、MAIL_TO）が設定されていないので、通知を送りません。")
        return

    parts, subjects = [], []
    if posts:
        described = [describe_post(p) for p in posts]
        subjects.append("今日の論文：" + described[0][0])
        parts.append("【論文紹介】\n\n" + "\n\n".join(d[1] for d in described))
    if courses:
        described = [describe_course(c) for c in courses]
        subjects.append("講義ノート更新：" + "、".join(d[0] for d in described))
        parts.append("【講義ノート】\n\n" + "\n\n".join(d[1] for d in described))

    msg = EmailMessage()
    msg["Subject"] = "[QCLab] " + " ／ ".join(subjects)
    msg["From"] = user
    msg["To"] = to
    msg.set_content(
        "\n\n".join(parts)
        + f"\n\n――\nQCLab：{SITE}\n論文紹介は AI（Claude）が自動で作成し、人による確認を経ずに公開しています。\n"
    )
    with smtplib.SMTP_SSL("smtp.gmail.com", 465, timeout=60) as s:
        s.login(user, password)
        s.send_message(msg)
    print(f"通知を送りました：{msg['Subject']}")


if __name__ == "__main__":
    main()
