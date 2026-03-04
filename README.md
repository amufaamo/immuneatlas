<img width="1699" height="1043" alt="image" src="https://github.com/user-attachments/assets/ad9eee7e-9b6d-41c4-9f62-d6e957914368" />


# ImmunoAtlas

**ImmunoAtlas: Single-cell Immune Repertoire Database**

ImmunoAtlas is a comprehensively integrated database and interactive web platform designed for exploring single-cell immune repertoire (TCR/BCR) data alongside gene expression data. It automatically processes single-cell sequencing datasets (such as from NCBI GEO), performs quality control, integration, and dimension reduction, and provides an intuitive, fast, web-based explorer to visualize cell clusters, clonotypes, and gene expressions.

## 🔥 Key Features

- **Automated Data Processing Pipeline**: R and Bash scripts that automate the downloading and processing of 10x Genomics single-cell immune profiling data using `Seurat` and `scRepertoire`.
- **Merged Study Integration**: Batch correction and integration of multiple samples across a study using `Harmony`.
- **Interactive Web Explorer**: A lightweight, Vanilla JS web application designed for high performance.
  - Interactive UMAP visualization using `Plotly.js`.
  - Color-coding by Clusters, TCR/BCR clonotype information, Sample origins, and Severity.
  - **Gene Expression Feature Plot**: Dynamically colors cells based on specific gene expression levels without performance drops (uses separate sparse JSON data).
  - Top Clonotypes list, statistical summaries, and dynamic data filtering.
- **Study Browser**: A main portal to search, filter, and access processed datasets.

## ⚙️ Technology Stack

- **Data Processing**: R (Seurat v5, scRepertoire, Harmony), Bash
- **Web Interface**: HTML5, CSS3 (Custom Variables/Design System), Vanilla JavaScript
- **Visualization**: Plotly.js
- **Data Transfer**: Highly optimized, sparse JSON format for fast rendering.

## 📂 Repository Structure

```
.
├── webapp/                 # Frontend web application (HTML, CSS, JS)
├── process_singlecell.R    # Seurat/scRepertoire processing for local data
├── process_geo_sample.R    # Automated processing for GEO datasets
├── merge_study_umap.R      # Batch correction and merging of samples (Harmony)
├── batch_process_*.sh      # Automation scripts for dataset downloads and processing
├── generate_expr_data.R    # Generates sparse JSON for gene expression visualization
├── server.py               # Lightweight Python HTTP server for local development
└── project.md              # Detailed project roadmap and architectural documentation
```

*Note: Large raw data files (`.h5`, `.csv`) and processed output files (`.qs2`, `.json`) are excluded from this repository via `.gitignore`.*

## 🚀 How to Run Locally

1. **Clone the repository:**
   ```bash
   git clone https://github.com/amufaamo/immuneatlas.git
   cd immuneatlas
   ```

2. **Start the local server:**
   You do not need to install complex Node.js dependencies. The app runs via a simple Python server.
   ```bash
   python3 server.py
   ```

3. **Open in Browser:**
   - **Study Browser**: `http://localhost:8001/webapp/index.html`
   - **Explorer**: Available via links in the Study Browser (e.g., `http://localhost:8001/webapp/viewer.html?sample=...`)

---
*Developed as a high-performance, accessible immune repertoire analysis platform.*
