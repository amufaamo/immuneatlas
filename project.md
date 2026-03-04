# **ImmunoAtlas — Single-cell Immune Repertoire Database**
## プロジェクト進行ドキュメント（2026-02-23 更新）

---

## **1. プロジェクト概要 (Project Overview)**

**目的:**
NCBI GEO等の公共データおよびローカルデータから取得したシングルセル免疫レパトア（TCR/BCR）データを、Webブラウザ上でインタラクティブに閲覧・探索できるデータベースを構築する。

**最終ゴール:**
全自動パイプラインにより、GEO上の数千サンプルを自動解析・格納し、CDR3配列やV/J遺伝子で逆引き検索できる大規模プラットフォームの構築。

**ターゲット:**
- プロジェクトメンバー（システム構成とUI/UXの検証）
- 共同研究者（デモンストレーション）

---

## **2. ページ・コンポーネント命名規則**

アプリは2つの主要ページで構成される。AIへの指示を明確にするため、以下の名前を統一して使用すること。

| ページ名 | ファイル | 説明 |
|---|---|---|
| **StudyBrowser** | `webapp/index.html` | Study一覧テーブルが表示されるトップページ。Study・Sample・タグ検索が可能。 |
| **Explorer** | `webapp/viewer.html` | 各SampleのUMAP・クローン分布・各種図を表示するインタラクティブビューア。UMAPに限らず将来的にさまざまな図を追加する想定のため "UMAP Viewer" ではなく "Explorer" と呼ぶ。 |

**用語定義:**
- **Study**: 論文・GEOアクセッション単位のデータセット（例: GSE145926）。複数Sampleを含む場合がある。
- **Sample**: 解析単位。GSM番号や個人ID（例: GSE145926_C141）に対応。
- **Merged Study**: 同一StudyのSampleをHarmonyでバッチ補正して統合したもの（例: `?study=gse145926`）。

---

## **3. システムアーキテクチャ**

```
[データ準備 (R / Bash)]          [Webアプリ (Vanilla JS)]
  process_geo_sample.R   →  JSON  →  StudyBrowser (index.html)
  process_singlecell.R   →  JSON  →  Explorer     (viewer.html)
  merge_study_umap.R     →  JSON  →  Explorer (merged mode)
                              ↓
                         server.py (Python HTTP)  ← ローカル配信
```

### 実行環境
- **すべてのR/Bashスクリプト**: Conda環境 `alew` (`conda activate alew`) 上で実行
- **Webサーバ**: `server.py`（Python `http.server` ベース、ポート8001）

### ディレクトリ構成

```
260212_alew_singlecell_tcrbcr_database/
├── project.md                   ← このファイル
├── server.py                    ← ローカルWebサーバ（python3 server.py）
├── data/                        ← 生データ置き場（.h5, .csv 等）
│   └── <SampleID>/
│       ├── filtered_feature_bc_matrix.h5
│       └── filtered_contig_annotations.csv
├── output/                      ← R処理結果（JSON, qs2）を配置
│   ├── umap_data_<SampleID>.json
│   ├── top_clones_<SampleID>.json
│   ├── seurat_<SampleID>.qs2
│   ├── studies.json             ← StudyBrowserが読み込むStudy一覧
│   └── merged/                  ← merge_study_umap.R の出力
│       └── <StudyID>/
│           ├── umap_merged_<StudyID>.json
│           └── top_clones_merged_<StudyID>.json
├── process_singlecell.R         ← ローカルデータ処理スクリプト
├── process_geo_sample.R         ← GEO公共データ処理スクリプト
├── merge_study_umap.R           ← Merged Study生成スクリプト
├── batch_process_gse145926.sh   ← GSE145926の複数Sample一括処理
├── download_gse145926.sh        ← GSE145926データのダウンロード
├── generate_studies_json.sh     ← studies.json 再生成スクリプト
├── run_pipeline.sh              ← 全体パイプライン実行スクリプト
├── logs/                        ← 処理ログ
└── webapp/
    ├── index.html               ← StudyBrowser
    ├── viewer.html              ← Explorer
    ├── app.js                   ← Explorer のメインロジック
    └── style.css                ← 共通スタイル
```

---

## **4. データフォーマット仕様**

### umap_data_\<ID\>.json（セルレベルデータ）

```json
[
  {
    "barcode": "AAACCTGAGCGTTGCC-1",
    "umap1": -3.12,
    "umap2": 5.67,
    "seurat_cluster": "0",
    "cell_type": "CD8 T",
    "tcr_clone": "TRAV12-2_CAVMDSSYKLIF_TRBV20-1_CASSQETQYF",
    "trav": "TRAV12-2",
    "traj": "TRAJ12",
    "trbv": "TRBV20-1",
    "trbj": "TRBJ2-5",
    "cdr3a": "CAVMDSSYKLIF",
    "cdr3b": "CASSQETQYF",
    "sample": "C141",
    "severity": "Moderate"
  },
  ...
]
```

### top_clones_\<ID\>.json（クローン集計データ）

```json
[
  {
    "clone_id": "TRAV12-2_CAVMDSSYKLIF_TRBV20-1_CASSQETQYF",
    "count": 42,
    "proportion": 0.0123,
    "cdr3a": "CAVMDSSYKLIF",
    "cdr3b": "CASSQETQYF",
    "trav": "TRAV12-2",
    "trbv": "TRBV20-1"
  },
  ...
]
```

### studies.json（StudyBrowser読み込み用）

```json
[
  {
    "id": "gse145926_c141",
    "studyGroup": "gse145926",
    "studyGroupTitle": "Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19",
    "title": "Single-cell landscape ... [C141]",
    "description": "論文アブストラクト等の説明文",
    "disease": "COVID-19",
    "severity": "Moderate",
    "species": "Homo sapiens",
    "assay": "10x 5′ GEX + TCR",
    "organ": "BALF",
    "viewerUrl": "/webapp/viewer.html?sample=gse145926_c141",
    "dataUrl": "/output/umap_data_GSE145926_C141.json",
    "qs2Url": "/output/seurat_GSE145926_C141.qs2"
  }
]
```

---

## **5. 技術スタック**

| カテゴリ | 技術/ツール | 備考 |
|---|---|---|
| **実行環境** | Conda (`alew` 環境) | すべてのデータ処理の基盤 |
| **データ取得** | Bash (wget/curl) | 公共データのダウンロード用 |
| **データ処理 (R)** | Seurat v5, scRepertoire, harmony, tidyverse | Conda環境内で実行 |
| **オブジェクト保存** | qs2 (`qs_save`/`qs_read`) | `.qs2` 形式で高速圧縮保存 |
| **データフォーマット** | JSON | Webに渡すための軽量化 |
| **Web UI** | Vanilla HTML/CSS/JavaScript | フレームワーク不使用・シンプル構成 |
| **グラフ描画** | Plotly.js (CDN) | インタラクティブUMAP描画 |
| **Webサーバ** | Python `http.server` (`server.py`) | ローカル開発用（ポート8001） |

---

## **6. フェーズ別ロードマップ**

### ✅ Phase 1: エンドツーエンドのパイプライン開通（完了）

**完了した内容:**
- [x] R処理スクリプト作成（`process_singlecell.R`, `process_geo_sample.R`）
- [x] ローカルデータ（MVP Sample）でSeurat + scRepertoire統合・JSON出力
- [x] GEO公共データ（GSE145926 C141）のダウンロードと処理
- [x] **StudyBrowser** (`index.html`) の実装
  - Studyテーブル（タイトル・Disease・Species・Assay・Organ列）
  - StudyGroup折りたたみ（複数Sampleを1グループに収納）
  - Studies.jsonからの動的ロード（fallback定義あり）
  - テキスト検索フィルタ
  - セル数の非同期カウント表示
- [x] **Explorer** (`viewer.html` + `app.js`) の実装
  - UMAPインタラクティブ表示（Plotly.js）
  - カラーモード: Clusters / TCR / BCR / Sample / Severity
  - TCRの色分け: Clone / TRAV / TRBV / TRAJ / TRBJ / CDR3α / CDR3β / α Len / β Len
  - BCRの色分け: Clone / IGHV / IGLV-IGKV / CDR3H / CDR3L
  - クラスターフィルタパネル（チェックボックス）
  - クローン一覧・検索・クリックでハイライト
  - Sampleフィルタパネル（Merged Study用）
  - Sample切り替えドロップダウン（ヘッダー）
  - セル数・選択数の表示
  - .qs2ダウンロードボタン
  - ホバーツールチップ
- [x] GSE145926 の複数サンプル処理（C141他）
- [x] Merged Study（Harmony バッチ補正 UMAP）: `merge_study_umap.R`
- [x] `server.py` によるローカルサーバ配信（ポート8001）

---

### 🔲 Phase 2: 機能強化と修正

**目標:** Explorerの機能を充実させ、データの信頼性・使いやすさを向上する。ブラウザの軽快な動作は最優先事項として維持する。

#### 優先課題

- [x] **Gene Expression Feature Plot** (Explorer最優先):
  - [x] Seuratの `FeaturePlot()` 相当の機能。指定した遺伝子の各細胞の発現量をUMAP上にグラデーション（低→高: グレー→シアン等）で重ねて表示する。
  - [x] **UI**: 左パネルにView Modeボタン「Gene Expr」を追加。選択するとテキスト入力欄が現れ、遺伝子名（例: `CD8A`, `IL7R`）を入力してEnterで即座に描画更新。
  - **軽量化の設計方針**（動作が重くならないために必ず守ること）:
    1. [x] **事前データ分離**: 遺伝子発現行列はメインの `umap_data_*.json` から切り離し、別ファイル `expr_data_<SampleID>.json` として保存。Gene Feature Plot モードに入ったときに初めて1回だけfetchする（起動時は不要）。
    2. [x] **スパースJSON形式**: 発現量が0の細胞はJSONに含めず、`{ "CD8A": {"BARCODE1": 2.3, "BARCODE3": 1.1, ...} }` のような辞書形式で保存。ファイルサイズを最小化する。
    3. [x] **上位N遺伝子のみ格納**: 全遺伝子ではなく、高発現・高分散な上位2000〜5000遺伝子のみを `expr_data` に含める（Rスクリプト側でフィルタリング）。
    4. [x] **色付けはPlotlyの `marker.color` 配列の差し替えのみ**: UMAPの点の位置は変えず、`Plotly.restyle()` で色配列だけ更新することで再描画コストを最小化する。
    5. [x] **入力デバウンス**: 遺伝子名入力後300ms待ってから描画（高速タイピング中に毎回描画しない）。

- [x] **クローン多様性指標** (Explorer): Shannon指数・Simpson指数をStatisticsパネルに追加表示（JSONから計算可能）
- [ ] **複数クローン同時選択** (Explorer): クローン一覧でチェックボックス複数選択対応・選択クローンを異なる色でハイライト
- [ ] **CSV/PNG エクスポート機能** (Explorer): 表示中のUMAPのPNG保存、選択クローンのCSV出力
- [ ] **サンプル間クローン比較** (Explorer): 同一StudyのSample間でクローンオーバーラップをベン図またはヒートマップで可視化
- [ ] **StudyBrowserのフィルタ強化**: Disease/Organ/Species/Assayでのドロップダウン絞り込みUI
- [ ] **データ品質確認**: 既存JSONの整合性チェック（barcode重複・欠損値等）
- [ ] **Error handling改善**: JSON読み込み失敗時のユーザーフレンドリーなエラー表示
- [ ] **GSE145926 追加サンプル処理**: C142, C143, C144... の順次追加

#### expr_data_\<SampleID\>.json の仕様（Gene Feature Plot用）

```json
{
  "genes": ["CD8A", "CD3D", "IL7R", ...],
  "expr": {
    "CD8A": {
      "AAACCTGAGCGTTGCC-1": 2.31,
      "AACCGCGTCTTGCAAG-1": 0.87
    },
    "CD3D": {
      "AAACCTGAGCGTTGCC-1": 1.52
    }
  }
}
```

- 発現量0のセルはキーごと省略（スパース表現）
- 値は log-normalized された発現量（Seuratnの `data` slot）
- Rスクリプト（`process_geo_sample.R` / `process_singlecell.R`）に `expr_data` 出力ステップを追加して生成する

---

### 🔲 Phase 3: 网羅的なデータ収集とデータベース完成

**目標:** NCBI GEOの公開シングルセル免疫レパトアデータを網羅的に取得し、大規模データベースを構築する。

#### データ収集戦略

1. **ターゲット選定**: GEO検索で "single cell" + "TCR" or "BCR" の全データセットをリストアップ
2. **自動ダウンロード**: GEO APIまたはスクレイピングで対象GSMを特定し自動取得
3. **バッチ処理**: `batch_process_gse145926.sh` を参考に複数Sample一括処理スクリプト作成
4. **品質管理**: 処理失敗サンプルのブラックリスト化・リトライロジック

#### データベーススケールアップ

- [ ] SQLiteまたはPostgreSQLへの移行（JSON → DB）
- [ ] CDR3配列による逆引き検索エンジン
- [ ] VDJdbやIEDBとの照合（抗原特異性の自動付与）
- [ ] 外部ID（PubMed / GEO / SRA）との自動リンク

#### 将来的な高度機能

- [ ] Alluvial plot（サンキーダイアグラム）によるクローン動態追跡
- [ ] サンプル間クローンオーバーラップのヒートマップ・ベン図
- [ ] scArches等を用いた細胞タイプ自動アノテーション
- [ ] Nextflow/Snakemake による完全自動パイプライン化

---

## **7. アプリの新規構築ガイド（AIへの指示用）**

このセクションは、将来AIに一発でアプリを再構築させるための詳細仕様書。

### 7-1. アプリ全体の設計思想

- **フレームワーク不使用**: Vanilla HTML + CSS + JavaScript（Node.js, React, Next.js 等は使わない）
- **CDN依存のみ**: Plotly.js を CDN から読み込む
- **ダークモード固定**: 常時ダークテーマ
- **デザインシステム**: CSS変数（カスタムプロパティ）で色・サイズを管理

### 7-2. CSS変数（style.css で定義）

```css
:root {
  --bg-base:      #0a0f1a;   /* ページ背景 */
  --bg-surface:   #101828;   /* カード・パネル背景 */
  --bg-raised:    #1a2333;   /* 入力フィールド背景 */
  --bg-glass:     rgba(16,24,40,0.8);

  --text-primary:   #e2e8f0;
  --text-secondary: #94a3b8;
  --text-muted:     #64748b;

  --accent-cyan:    #22d3ee;
  --accent-violet:  #a78bfa;
  --accent-emerald: #34d399;
  --accent-amber:   #fbbf24;
  --accent-pink:    #f472b6;

  --border-subtle:  rgba(148,163,184,0.08);
  --border-strong:  rgba(148,163,184,0.2);

  --header-h:     56px;
  --radius-md:    8px;
  --radius-lg:    12px;
  --transition:   0.15s ease;
}
```

### 7-3. StudyBrowser (index.html) 仕様

**概要:**
- ヘッダー: ロゴ（SVGアイコン）+ "ImmunoAtlas MVP v0.2"
- ヒーローセクション: タイトル・説明・統計（Study数・総細胞数）
- 検索バー: テキスト入力でStudy名・Disease・Species等を横断検索
- Studies テーブル:
  - 列構成: `Study（タイトル+説明+セル数）| Disease | Species | Assay | Organ | →`
  - `grid-template-columns: 3fr 1fr 1fr 1fr 1fr 80px`
  - 1 Sampleのみ: 1行で直接 `viewer.html?sample=<id>` へリンク
  - 複数Sample: グループヘッダー行（クリックで展開）＋サブ行
  - ホバー時: 左端にcyan→violetグラデーションのバー表示

**データ読み込み順序:**
1. `/output/studies.json` をfetchして取得
2. 失敗時は `FALLBACK_STUDIES`（ハードコード）を使用
3. 各StudyのJSONから非同期でセル数を取得して表示更新

**studyGroupフィールド:**
- 同じ `studyGroup` 値を持つStudyは1グループとしてまとめて表示
- `studyGroupTitle` でグループの表示タイトルを指定（なければtitleから `[xxx]` を除去）

### 7-4. Explorer (viewer.html + app.js) 仕様

**URLパラメータ:**
- `?sample=<id>` : 単一Sample表示（例: `?sample=gse145926_c141`）
- `?study=<id>` : Merged Study表示（例: `?study=gse145926`）

**画面レイアウト:**
```
[ヘッダー: ロゴ | ← Studies | [Sample切替ドロップダウン] | Cells数 | Selected数 | .qs2 DL]
[左パネル(260px固定)] | [右パネル(残り幅)]
```

**左パネルのセクション（上から）:**
1. **View Mode**: モードボタン群
   - `Clusters` - seurat_clusterで色分け
   - `TCR` - TCRクローン情報で色分け（サブオプション: Clone/TRAV/TRBV/TRAJ/TRBJ/CDR3α/CDR3β/αLen/βLen）
   - `BCR` - BCRクローン情報で色分け（サブオプション: Clone/IGHV/IGLV-IGKV/CDR3H/CDR3L）
   - `Sample` - サンプルで色分け（Merged Studyデータがある場合のみ表示）
   - `Severity` - 重症度で色分け（該当データがある場合のみ表示）
2. **Clusters**: チェックボックスリスト（クラスタの表示/非表示）
3. **Samples**: チェックボックスリスト（Merged Studyのみ表示）
4. **Top Clonotypes**: クローン一覧（検索可能・クリックでハイライト）
5. **Statistics**: Total Cells / TCR+ Cells / BCR+ Cells / Clusters数

**右パネル:**
- Plotly.js でUMAPを描画（WebGLまたはSVG）
- ホバーで各セルの情報をツールチップ表示
- セル選択時にヘッダーの`Selected`カウントを更新

**データ読み込みロジック（app.js）:**
```
1. URLパラメータを解析
2. ?study= の場合:
   - /output/merged/<studyId>/umap_merged_<studyId>.json を fetch
   - /output/merged/<studyId>/top_clones_merged_<studyId>.json を fetch
   - Sample/Severityモードボタンを表示
3. ?sample= の場合:
   - /output/umap_data_<SampleID>.json を fetch（SampleIDは大文字）
   - /output/top_clones_<SampleID>.json を fetch
   - studies.json から同じstudyGroupの他Sampleを取得しドロップダウン生成
4. データをもとにPlotlyでUMAP描画
5. セル数・統計を更新
```

**Sample切り替えドロップダウン（ヘッダー中央）:**
- 現在のSample名をバッジ表示
- 同一studyGroupの他Sampleがある場合はvシェブロンを表示
- クリックでドロップダウン展開→選択で `?sample=<id>` に遷移

### 7-5. R処理スクリプトの仕様

**process_geo_sample.R** (GEO公共データ用):
```
引数: args[1] = SampleID (例: "GSE145926_C141")
      args[2] = data dir (例: "/path/to/data/GSE145926_C141")
      args[3] = output dir (例: "/path/to/output")
処理:
  1. Read10X_h5() または Read10X() でGEXロード
  2. CreateSeuratObject → QC (nFeature > 200, percent.mt < 20)
  3. NormalizeData → FindVariableFeatures → ScaleData → PCA
  4. FindNeighbors → FindClusters → RunUMAP
  5. filtered_contig_annotations.csv を scRepertoire::combineTCR() / combineBCR() でロード
  6. scRepertoire::combineExpression() でSeuratオブジェクトに付与
  7. UMAP座標 + metadata + TCR/BCR情報を JSON で出力: umap_data_<SampleID>.json
  8. 上位50クローンを JSON で出力: top_clones_<SampleID>.json
  9. qs2::qs_save() で Seurat オブジェクト保存: seurat_<SampleID>.qs2
```

**merge_study_umap.R** (Merged Study用):
```
引数: args[1] = StudyID (例: "gse145926")
      args[2] = inputDir (複数の .qs2 があるディレクトリ)
      args[3] = outputDir
処理:
  1. 各 .qs2 をロード → merge()
  2. harmony::RunHarmony() でバッチ補正
  3. FindNeighbors → FindClusters → RunUMAP (dims.use = "harmony")
  4. umap_merged_<studyId>.json を出力（sample・severity列を含む）
  5. top_clones_merged_<studyId>.json を出力
```

### 7-6. サーバ起動方法

```bash
# alew環境は不要（Python標準ライブラリのみ使用）
cd /home/masakazu/Data/260212_alew_singlecell_tcrbcr_database
python3 server.py
# → http://localhost:8001/webapp/index.html でStudyBrowserが開く
# → http://localhost:8001/webapp/viewer.html?sample=gse145926_c141 でExplorerが開く
```

server.py は Python 標準の `http.server` をベースに以下を追加:
- `/output/` パスを `output/` ディレクトリにマッピング
- `/webapp/` パスを `webapp/` ディレクトリにマッピング
- CORS ヘッダーを付与

---

## **8. 既知の問題・注意事項**

- `viewer.html` の `<script src="app.js?v=N">` のバージョン番号は、app.js更新時にブラウザキャッシュを回避するために手動でインクリメントする
- データパスは大文字小文字を区別する: `umap_data_GSE145926_C141.json`（大文字）
- Merged Study の URL は `?study=gse145926`（小文字）、データファイルは `umap_merged_gse145926.json`（小文字）
- 新しいSampleを追加した後は `generate_studies_json.sh` を実行して `studies.json` を再生成する（またはFALLBACK_STUDIESに手動追加）
- R処理では `alew` conda環境が必須: `conda run -n alew Rscript process_geo_sample.R ...`
