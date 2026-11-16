# QCLab（qc-notes）

量子化学の論文紹介と講義ノートのサイト（MkDocs Material で作り、GitHub Pages で公開）。
公開先：https://yugo0323.github.io/qc-notes/

## 手元で見る

```bash
pip install -r requirements.txt
mkdocs serve
```

## 公開

`main` ブランチに push すると、GitHub Actions がサイトを作り直して公開します（`.github/workflows/deploy.yml`）。

## 毎日の論文紹介

`.github/workflows/daily-paper.yml` が毎朝（日本時間 5:17）動きます。

1. `scripts/fetch_candidates.py`：arXiv と Crossref から候補を集める
2. Claude Code（`scripts/daily_prompt.md` の指示）が `docs/papers/<年>/<日付>.md` を書く
3. `scripts/check_post.py`：論文の実在と題名、見出し、注意書きを確かめる
4. `scripts/build_index.py`：一覧とトップページを作り直す
5. コミット、プッシュ、公開

必要な secret：`CLAUDE_CODE_OAUTH_TOKEN`（`claude setup-token` で作る）。
