# qc-notes

量子化学の論文紹介と講義ノートのサイト（MkDocs Material で作り、GitHub Pages で公開）。

## 手元で見る

```bash
pip install -r requirements.txt
mkdocs serve
```

## 公開

`main` ブランチに push すると、GitHub Actions がサイトを作り直して公開します（`.github/workflows/deploy.yml`）。
