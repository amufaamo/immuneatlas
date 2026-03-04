/* ===================================================
   ImmunoAtlas — app.js  v6  (Gene Expression Feature Plot)
   =================================================== */

// ─── Study-level config ──────────────────────────────
const STUDY_CONFIGS = {
    gse145926: {
        dataPath: '/output/umap_data_gse145926_merged.json',
        clonesPath: '/output/top_clones_gse145926_merged.json',
        summaryPath: '/output/sample_summary_gse145926_merged.json',
        qs2Path: '/output/seurat_gse145926_merged.qs2',
        label: 'GSE145926 · Merged (12 samples · 84,114 cells) · Harmony',
        isMerged: true,
    },
};

// ─── Per-Sample Config ───────────────────────────────
const SAMPLES = {
    mvp: {
        dataPath: '/output/umap_data.json',
        clonesPath: '/output/top_clones.json',
        qs2Path: '/output/seurat_MVP_Sample.qs2',
        label: '10x Genomics · 5\u2032 PBMC · Human Diseased',
    },
    gse145926_c141: { dataPath: '/output/umap_data_GSE145926_C141.json', clonesPath: '/output/top_clones_GSE145926_C141.json', qs2Path: '/output/seurat_GSE145926_C141.qs2', label: 'GSE145926 · C141 · COVID-19 BALF (Moderate)' },
    gse145926_c142: { dataPath: '/output/umap_data_GSE145926_C142.json', clonesPath: '/output/top_clones_GSE145926_C142.json', qs2Path: '/output/seurat_GSE145926_C142.qs2', label: 'GSE145926 · C142 · COVID-19 BALF (Moderate)' },
    gse145926_c143: { dataPath: '/output/umap_data_GSE145926_C143.json', clonesPath: '/output/top_clones_GSE145926_C143.json', qs2Path: '/output/seurat_GSE145926_C143.qs2', label: 'GSE145926 · C143 · COVID-19 BALF (Moderate)' },
    gse145926_c144: { dataPath: '/output/umap_data_GSE145926_C144.json', clonesPath: '/output/top_clones_GSE145926_C144.json', qs2Path: '/output/seurat_GSE145926_C144.qs2', label: 'GSE145926 · C144 · COVID-19 BALF (Severe)' },
    gse145926_c145: { dataPath: '/output/umap_data_GSE145926_C145.json', clonesPath: '/output/top_clones_GSE145926_C145.json', qs2Path: '/output/seurat_GSE145926_C145.qs2', label: 'GSE145926 · C145 · COVID-19 BALF (Severe)' },
    gse145926_c146: { dataPath: '/output/umap_data_GSE145926_C146.json', clonesPath: '/output/top_clones_GSE145926_C146.json', qs2Path: '/output/seurat_GSE145926_C146.qs2', label: 'GSE145926 · C146 · COVID-19 BALF (Severe)' },
    gse145926_c51: { dataPath: '/output/umap_data_GSE145926_C51.json', clonesPath: '/output/top_clones_GSE145926_C51.json', qs2Path: '/output/seurat_GSE145926_C51.qs2', label: 'GSE145926 · C51  · COVID-19 BALF (Severe)' },
    gse145926_c52: { dataPath: '/output/umap_data_GSE145926_C52.json', clonesPath: '/output/top_clones_GSE145926_C52.json', qs2Path: '/output/seurat_GSE145926_C52.qs2', label: 'GSE145926 · C52  · COVID-19 BALF (Moderate)' },
    gse145926_c100: { dataPath: '/output/umap_data_GSE145926_C100.json', clonesPath: '/output/top_clones_GSE145926_C100.json', qs2Path: '/output/seurat_GSE145926_C100.qs2', label: 'GSE145926 · C100 · COVID-19 BALF (Moderate)' },
    gse145926_c148: { dataPath: '/output/umap_data_GSE145926_C148.json', clonesPath: '/output/top_clones_GSE145926_C148.json', qs2Path: '/output/seurat_GSE145926_C148.qs2', label: 'GSE145926 · C148 · COVID-19 BALF (Severe)' },
    gse145926_c149: { dataPath: '/output/umap_data_GSE145926_C149.json', clonesPath: '/output/top_clones_GSE145926_C149.json', qs2Path: '/output/seurat_GSE145926_C149.qs2', label: 'GSE145926 · C149 · COVID-19 BALF (Severe)' },
    gse145926_c152: { dataPath: '/output/umap_data_GSE145926_C152.json', clonesPath: '/output/top_clones_GSE145926_C152.json', qs2Path: '/output/seurat_GSE145926_C152.qs2', label: 'GSE145926 · C152 · COVID-19 BALF (Severe)' },
};

// ─── Mode detection ───────────────────────────────────
const _rawSampleKey = window.IMMUNOATLAS_SAMPLE || new URLSearchParams(location.search).get('sample') || 'mvp';
let STUDY_ID = window.IMMUNOATLAS_STUDY || null;   // e.g. 'gse145926'
if (STUDY_ID === 'merged') {
    STUDY_ID = 'gse145926'; // fallback for old URL
}
const IS_STUDY_MODE = !!STUDY_ID && STUDY_ID in STUDY_CONFIGS;
const sampleKey = IS_STUDY_MODE ? ('__study__' + STUDY_ID) : _rawSampleKey;
const SAMPLE_CFG = IS_STUDY_MODE ? STUDY_CONFIGS[STUDY_ID] : (SAMPLES[sampleKey] || SAMPLES.mvp);
const DATA_PATH = SAMPLE_CFG.dataPath;
const CLONES_PATH = SAMPLE_CFG.clonesPath;

// ─── Color Palettes ──────────────────────────────────
const CLUSTER_COLORS = [
    '#22d3ee', '#a78bfa', '#34d399', '#f472b6', '#fbbf24',
    '#60a5fa', '#f87171', '#4ade80', '#c084fc', '#fb923c',
    '#38bdf8', '#e879f9', '#86efac', '#fdba74', '#a5f3fc',
    '#d8b4fe', '#6ee7b7', '#fca5a5', '#93c5fd', '#fde68a',
    '#7dd3fc', '#f9a8d4', '#6ee7f7', '#bbf7d0',
];

// Sample-level colors (12色, 各サンプルに割り当て)
const SAMPLE_COLORS = [
    '#22d3ee', '#a78bfa', '#34d399', '#f472b6', '#fbbf24',
    '#60a5fa', '#f87171', '#4ade80', '#c084fc', '#fb923c',
    '#e879f9', '#86efac',
];

// Severity colors
const SEVERITY_COLORS = {
    Moderate: '#fbbf24',   // amber
    Severe: '#f472b6',   // rose
    Healthy: '#34d399',   // emerald
    Unknown: '#60a5fa',   // blue
};

// ─── Repertoire Field Configs ─────────────────────────
const TCR_FIELDS = [
    { key: 'clone', label: 'Clone', field: 'TCR_CTgene', type: 'cat', title: 'Top Clonotypes', hint: 'TRAV…TRBV, e.g. TRAV1-2' },
    { key: 'trav', label: 'TRAV', field: '_trav', type: 'cat', title: 'Top TRAV Genes', hint: 'e.g. TRAV1-2, TRAV26-1' },
    { key: 'trbv', label: 'TRBV', field: '_trbv', type: 'cat', title: 'Top TRBV Genes', hint: 'e.g. TRBV20-1, TRBV5-5' },
    { key: 'traj', label: 'TRAJ', field: '_traj', type: 'cat', title: 'Top TRAJ Genes', hint: 'e.g. TRAJ6, TRAJ33' },
    { key: 'trbj', label: 'TRBJ', field: '_trbj', type: 'cat', title: 'Top TRBJ Genes', hint: 'e.g. TRBJ2-1, TRBJ1-5' },
    { key: 'cdr3a', label: 'CDR3α', field: '_cdr3a', type: 'cat', title: 'Top CDR3α Seqs', hint: 'motif e.g. CAGQ, CALS' },
    { key: 'cdr3b', label: 'CDR3β', field: '_cdr3b', type: 'cat', title: 'Top CDR3β Seqs', hint: 'motif e.g. CASS, CASF' },
    { key: 'cdr3a_len', label: 'αLen', field: '_cdr3a_len', type: 'num', title: 'CDR3α Length', hint: '' },
    { key: 'cdr3b_len', label: 'βLen', field: '_cdr3b_len', type: 'num', title: 'CDR3β Length', hint: '' },
];
const BCR_FIELDS = [
    { key: 'clone', label: 'Clone', field: 'BCR_CTgene', type: 'cat', title: 'Top Clonotypes', hint: 'IGHV…, e.g. IGHV1-2' },
    { key: 'ighv', label: 'IGHV', field: '_ighv', type: 'cat', title: 'Top IGHV Genes', hint: 'e.g. IGHV1-2, IGHV3-23' },
    { key: 'lv', label: 'IGLV/IGKV', field: '_bcr_lv', type: 'cat', title: 'Top L-chain V', hint: 'e.g. IGKV1-5, IGLV2-14' },
    { key: 'cdr3h', label: 'CDR3H', field: '_cdr3h', type: 'cat', title: 'Top CDR3H Seqs', hint: 'heavy chain CDR3 motif' },
    { key: 'cdr3l', label: 'CDR3L', field: '_cdr3l', type: 'cat', title: 'Top CDR3L Seqs', hint: 'light chain CDR3 motif' },
];

function getActiveFieldCfg() {
    if (state.colorMode === 'tcr') return TCR_FIELDS.find(f => f.key === state.tcrColorBy);
    if (state.colorMode === 'bcr') return BCR_FIELDS.find(f => f.key === state.bcrColorBy);
    return null;
}

// ─── Cell Enrichment (parse V/J genes & CDR3) ────────
function enrichCells(cells) {
    cells.forEach(c => {
        // seurat_clusters を必ず文字列に正規化 (JSONによって数値/文字列が混在するため)
        if (c.seurat_clusters != null) c.seurat_clusters = String(c.seurat_clusters);
        // TCR genes: "TRAV.TRAJ.TRAC[;...]_TRBV.D.TRBJ.TRBC"
        if (c.TCR_CTgene) {
            const [alpha, beta] = c.TCR_CTgene.split('_');
            const aSegs = (alpha || '').split(';')[0].split('.');
            const bSegs = (beta || '').split('.');
            c._trav = aSegs[0] || null;
            c._traj = aSegs[1] || null;
            c._trbv = bSegs[0] || null;
            c._trbj = bSegs[2] || null;   // V . D . J . C
        } else { c._trav = c._traj = c._trbv = c._trbj = null; }

        // TCR CDR3: "cdr3a[;cdr3a2]_cdr3b"
        if (c.TCR_CTaa) {
            const [ap, bp] = c.TCR_CTaa.split('_');
            c._cdr3a = (ap || '').split(';')[0] || null;
            c._cdr3b = bp || null;
            c._cdr3a_len = c._cdr3a ? c._cdr3a.length : null;
            c._cdr3b_len = c._cdr3b ? c._cdr3b.length : null;
        } else { c._cdr3a = c._cdr3b = c._cdr3a_len = c._cdr3b_len = null; }

        // BCR genes: "IGHV.IGHJ.IGHC_IGKV/IGLV.J.C"
        if (c.BCR_CTgene) {
            const [h, l] = c.BCR_CTgene.split('_');
            c._ighv = (h || '').split('.')[0] || null;
            c._bcr_lv = (l || '').split('.')[0] || null;
        } else { c._ighv = c._bcr_lv = null; }

        // BCR CDR3: "cdr3h_cdr3l"
        if (c.BCR_CTaa) {
            const [h, l] = c.BCR_CTaa.split('_');
            c._cdr3h = h || null;
            c._cdr3l = l || null;
        } else { c._cdr3h = c._cdr3l = null; }
    });
}

// ─── App State ───────────────────────────────────────
const state = {
    cells: [], clones: [], sampleSummary: [],
    colorMode: 'cluster',          // 'cluster'|'tcr'|'bcr'|'sample'|'severity'|'gene'
    tcrColorBy: 'clone',
    bcrColorBy: 'clone',
    selectedFeatureValue: null,    // highlighted value
    selectedFeatureField: null,    // cell property field
    hiddenClusters: new Set(),
    hiddenSamples: new Set(),     // Merged mode: hidden sample_ids
    clusterCounts: {}, clusterIds: [],
    sampleIds: [], sampleColorMap: {},   // Merged mode
    cloneSearchQuery: '',
    // Gene Expr mode
    exprData: null,          // 遅延ロード: { genes: [...], expr: { CD8A: { barcode: val } } }
    exprLoading: false,      // fetch中フラグ
    activeGene: '',          // 現在表示中の遺伝子名
    geneColorValues: null,   // 各セルの発現量配列 (cells と同順)
    diversityData: null,     // Diversity metrics { "GSE145926_C141": {shannon: ...} }
    plotMode: 'umap',        // 'umap' or 'diversity'
};

// ─── DOM ─────────────────────────────────────────────
const $ = id => document.getElementById(id);
const loadingOverlay = $('loading-overlay');
const loadingProgress = $('loading-progress');
const totalCellsEl = $('total-cells');
const selectedCellsEl = $('selected-cells');
const statTotalEl = $('stat-total-cells');
const statTCREl = $('stat-tcr-cells');
const statBCREl = $('stat-bcr-cells');
const statClustersEl = $('stat-clusters');
const clusterLegend = $('cluster-legend');
const cloneList = $('clone-list');
const cloneSection = $('clonotype-section');
const clusterSection = $('cluster-section');
const cloneLegendOverlay = $('clone-legend-overlay');
const diversityOverlay = $('diversity-overlay');
const diversityPlotContainer = $('diversity-plot-container');

// ─── Utils ───────────────────────────────────────────
const fmt = n => (n == null ? '—' : n.toLocaleString());
const truncate = (s, n = 50) => (!s ? 'None' : s.length > n ? s.slice(0, n) + '…' : s);
const setProgress = p => { loadingProgress.style.width = p + '%'; };

function highlightMatch(text, q) {
    if (!q) return text;
    const i = text.toLowerCase().indexOf(q.toLowerCase());
    if (i === -1) return text;
    return text.slice(0, i) +
        `<mark class="clone-highlight">${text.slice(i, i + q.length)}</mark>` +
        text.slice(i + q.length);
}

// ─── Data Loading ────────────────────────────────────
async function loadData() {
    setProgress(10);
    const [cr, lr] = await Promise.all([fetch(DATA_PATH), fetch(CLONES_PATH)]);
    setProgress(50);
    const [cells, clones] = await Promise.all([cr.json(), lr.json()]);
    setProgress(90);
    return { cells, clones };
}

// ─── Init ─────────────────────────────────────────────
async function init() {
    try {
        const { cells, clones } = await loadData();
        enrichCells(cells);
        state.cells = cells;
        state.cells.forEach((c, i) => c._index = i);
        state.clones = Array.isArray(clones) ? clones : [];

        // cluster metadata
        cells.forEach(c => {
            const cl = c.seurat_clusters;
            state.clusterCounts[cl] = (state.clusterCounts[cl] || 0) + 1;
        });
        state.clusterIds = Object.keys(state.clusterCounts).sort((a, b) => +a - +b);

        // stats
        const tcrCells = cells.filter(c => c.TCR_CTgene).length;
        const bcrCells = cells.filter(c => c.BCR_CTgene).length;
        totalCellsEl.textContent = fmt(cells.length);
        statTotalEl.textContent = fmt(cells.length);
        statTCREl.textContent = fmt(tcrCells);
        statBCREl.textContent = fmt(bcrCells);
        statClustersEl.textContent = state.clusterIds.length;

        // Fetch studies.json to get Shannon Index etc.
        try {
            const sres = await fetch('/output/studies.json');
            if (sres.ok) {
                const studies = await sres.json();
                const curId = IS_STUDY_MODE ? `${STUDY_ID}_merged` : sampleKey;
                const meta = studies.find(s => s.id === curId);
                if (meta) {
                    if (meta.shannonIndex != null) $('stat-shannon').textContent = meta.shannonIndex.toFixed(3);
                    if (meta.nClonotypes != null) $('stat-clones').textContent = fmt(meta.nClonotypes);
                }
            }
        } catch (e) { console.warn('Metadata fetch failed', e); }

        setProgress(100);

        const badge = $('sample-badge');
        if (badge) badge.textContent = SAMPLE_CFG.label;

        // ── Study-level Merged mode: 追加セットアップ ─────────
        if (IS_STUDY_MODE) {
            // sample_id の一覧とカラーマップを構築
            const sampleSet = new Set(cells.map(c => c.sample_id).filter(Boolean));
            state.sampleIds = [...sampleSet].sort();
            state.sampleIds.forEach((sid, i) => {
                state.sampleColorMap[sid] = SAMPLE_COLORS[i % SAMPLE_COLORS.length];
            });

            // Sample / Severity ボタンを表示・有効化
            const btnSample = $('btn-sample');
            const btnSeverity = $('btn-severity');
            if (btnSample) { btnSample.style.display = ''; btnSample.classList.remove('mode-disabled'); }
            if (btnSeverity) { btnSeverity.style.display = ''; btnSeverity.classList.remove('mode-disabled'); }

            // Sample filter panel を表示
            const sfSection = $('sample-filter-section');
            if (sfSection) sfSection.style.display = '';

            buildSampleFilterLegend();

            // sample_summary を非同期ロード (セル数表示用オプション)
            if (SAMPLE_CFG.summaryPath) {
                fetch(SAMPLE_CFG.summaryPath)
                    .then(r => r.json())
                    .then(data => { state.sampleSummary = data; })
                    .catch(() => { });
            }

            // diversity_summary を非同期ロード (多様性プロット用)
            fetch('/output/diversity_summary.json')
                .then(r => r.json())
                .then(data => { state.diversityData = data; })
                .catch(() => { });
        }

        // TCR/BCR ボタン有効/無効
        const btnTCR = $('btn-tcr');
        const btnBCR = $('btn-bcr');
        const cardTCR = statTCREl ? statTCREl.closest('.stat-card') : null;
        const cardBCR = statBCREl ? statBCREl.closest('.stat-card') : null;

        if (tcrCells === 0 && btnTCR) {
            btnTCR.classList.add('mode-disabled');
            btnTCR.title = 'このサンプルにはTCRデータがありません';
            btnTCR.innerHTML = `<svg width="14" height="14" viewBox="0 0 14 14" fill="none">
              <path d="M2 7L7 2l5 5-5 5z" stroke="currentColor" stroke-width="1.5" fill="none" />
            </svg> TCR`;
            if (cardTCR) cardTCR.classList.add('stat-unavailable');
        }
        if (bcrCells === 0 && btnBCR) {
            btnBCR.classList.add('mode-disabled');
            btnBCR.title = 'このサンプルにはBCRデータがありません';
            btnBCR.innerHTML = `<svg width="14" height="14" viewBox="0 0 14 14" fill="none">
              <rect x="2" y="2" width="10" height="10" rx="2" stroke="currentColor" stroke-width="1.5" fill="none" />
            </svg> BCR`;
            if (cardBCR) cardBCR.classList.add('stat-unavailable');
        }

        buildClusterLegend();
        buildPlot();
        bindControls();

        setTimeout(() => {
            loadingOverlay.classList.add('fade-out');
            setTimeout(() => { loadingOverlay.style.display = 'none'; }, 500);
        }, 300);

    } catch (err) {
        console.error(err);
        document.querySelector('.loading-text').textContent = '⚠ Failed to load data.';
        document.querySelector('.loading-text').style.color = '#f87171';
    }
}

// ─── Cluster Legend ───────────────────────────────────
function buildClusterLegend() {
    clusterLegend.innerHTML = '';
    state.clusterIds.forEach((id, i) => {
        const color = CLUSTER_COLORS[i % CLUSTER_COLORS.length];
        const row = document.createElement('div');
        row.className = 'legend-item-row';
        row.id = `cluster-item-${id}`;
        row.dataset.cluster = id;
        row.innerHTML = `
          <div class="legend-color" style="background:${color}"></div>
          <span class="legend-label">Cluster ${id}</span>
          <span class="legend-count">${fmt(state.clusterCounts[id])}</span>`;
        row.addEventListener('click', () => toggleCluster(id));
        clusterLegend.appendChild(row);
    });
}

function toggleCluster(id) {
    if (state.hiddenClusters.has(id)) {
        state.hiddenClusters.delete(id);
        $(`cluster-item-${id}`).classList.remove('dimmed');
    } else {
        state.hiddenClusters.add(id);
        $(`cluster-item-${id}`).classList.add('dimmed');
    }
    const visibleCells = state.cells.filter(isCellVisible).length;
    totalCellsEl.textContent = fmt(visibleCells);
    updatePlot();
}

// ─── Sample Filter Legend (Merged mode) ──────────────
function buildSampleFilterLegend() {
    const legend = $('sample-filter-legend');
    if (!legend) return;
    legend.innerHTML = '';

    // severity でグループ化 + ソート
    const severityOrder = { Moderate: 0, Severe: 1, Healthy: 2, Unknown: 9 };
    const sorted = [...state.sampleIds].sort((a, b) => {
        const sevA = getSeverityForSample(a);
        const sevB = getSeverityForSample(b);
        return (severityOrder[sevA] ?? 9) - (severityOrder[sevB] ?? 9) || a.localeCompare(b);
    });

    let lastSev = null;
    sorted.forEach(sid => {
        const sev = getSeverityForSample(sid);
        const color = state.sampleColorMap[sid];
        const count = state.cells.filter(c => c.sample_id === sid).length;

        // severity グループヘッダー
        if (sev !== lastSev) {
            const hdr = document.createElement('div');
            hdr.className = 'sample-group-header';
            const sevColor = SEVERITY_COLORS[sev] || '#60a5fa';
            hdr.innerHTML = `<span style="color:${sevColor};font-weight:600;font-size:10px;text-transform:uppercase;letter-spacing:0.6px">${sev}</span>`;
            legend.appendChild(hdr);
            lastSev = sev;
        }

        const row = document.createElement('div');
        row.className = 'legend-item-row sample-filter-row';
        row.id = `sample-item-${sid}`;
        row.dataset.sample = sid;
        row.innerHTML = `
          <div class="legend-color" style="background:${color}"></div>
          <span class="legend-label" style="font-family:'JetBrains Mono',monospace;font-size:11px">${sid}</span>
          <span class="legend-count">${fmt(count)}</span>`;
        row.addEventListener('click', () => toggleSample(sid));
        legend.appendChild(row);
    });
}

function getSeverityForSample(sid) {
    // sampleSummaryからseverityを取得、なければcellsから推定
    const s = state.sampleSummary.find(x => x.sample_id === sid);
    if (s) return s.severity || 'Unknown';
    const cell = state.cells.find(c => c.sample_id === sid && c.severity);
    return cell ? cell.severity : 'Unknown';
}

function toggleSample(sid) {
    if (state.hiddenSamples.has(sid)) {
        state.hiddenSamples.delete(sid);
        $(`sample-item-${sid}`)?.classList.remove('dimmed');
    } else {
        state.hiddenSamples.add(sid);
        $(`sample-item-${sid}`)?.classList.add('dimmed');
    }
    // 表示中のセル数を更新
    const visibleCells = state.cells.filter(isCellVisible).length;
    totalCellsEl.textContent = fmt(visibleCells);
    updatePlot();
}

// ─── Feature List (sidebar) ───────────────────────────
function renderFeatureList() {
    const cfg = getActiveFieldCfg();
    if (!cfg) return;

    // Update section title and search placeholder
    const titleEl = $('clone-section-title');
    if (titleEl) titleEl.textContent = cfg.title;
    const searchEl = $('clone-search');
    if (searchEl && cfg.hint) searchEl.placeholder = `Search… ${cfg.hint}`;

    cloneList.innerHTML = '';

    if (cfg.type === 'num') {
        renderNumericDist(cfg);
        return;
    }

    // Categorical
    const q = state.cloneSearchQuery.trim().toLowerCase();
    const { field } = cfg;

    const counts = {};
    state.cells.forEach(c => { const v = c[field]; if (v != null) counts[v] = (counts[v] || 0) + 1; });
    const sorted = Object.entries(counts).sort((a, b) => b[1] - a[1]);
    const max = sorted.length ? sorted[0][1] : 1;
    const colorMap = {};
    sorted.forEach(([v], i) => { colorMap[v] = CLUSTER_COLORS[i % CLUSTER_COLORS.length]; });

    let entries = q ? sorted.filter(([v]) => v.toLowerCase().includes(q)) : sorted;

    if (entries.length === 0) {
        cloneList.innerHTML = `<div class="clone-empty">No results${q ? ` for <span class="clone-empty-query">"${q}"</span>` : ''}</div>`;
        return;
    }

    entries.forEach(([value, count]) => {
        const pct = (count / max * 100).toFixed(1);
        const rank = sorted.findIndex(([v]) => v === value);
        const color = colorMap[value];
        const card = document.createElement('div');
        card.className = 'clone-card';
        if (state.selectedFeatureValue === value) card.classList.add('active');

        const displayVal = q ? highlightMatch(truncate(value, 50), q) : truncate(value, 50);

        // For clone mode show rank; for gene mode show gene name prominently
        const rankLabel = cfg.key === 'clone'
            ? `#${rank + 1} · ${count.toLocaleString()} cells`
            : `#${rank + 1} · ${count.toLocaleString()} cells`;

        card.innerHTML = `
          <div class="clone-rank">${rankLabel}</div>
          <div class="clone-id ${state.selectedFeatureValue === value ? 'active' : ''}"
               style="color:${color}">${displayVal}</div>
          <div class="clone-bar-wrap">
            <div class="clone-bar-bg">
              <div class="clone-bar" style="width:${pct}%;background:${color}"></div>
            </div>
            <span class="clone-count">${fmt(count)}</span>
          </div>`;
        card.addEventListener('click', () => selectFeature(value, field));
        cloneList.appendChild(card);
    });
}

function renderNumericDist(cfg) {
    const vals = state.cells.map(c => c[cfg.field]).filter(v => v != null);
    if (!vals.length) {
        cloneList.innerHTML = `<div class="clone-empty">No ${cfg.label} data</div>`;
        return;
    }
    const counts = {};
    vals.forEach(v => { counts[v] = (counts[v] || 0) + 1; });
    const entries = Object.entries(counts).sort((a, b) => +a[0] - +b[0]);
    const maxC = Math.max(...entries.map(([, c]) => c));

    cloneList.innerHTML = entries.map(([len, count]) => `
      <div class="clone-card len-card">
        <div class="clone-bar-wrap">
          <span class="len-label">${len} aa</span>
          <div class="clone-bar-bg">
            <div class="clone-bar" style="width:${(count / maxC * 100).toFixed(1)}%;background:var(--accent-cyan)"></div>
          </div>
          <span class="clone-count">${fmt(count)}</span>
        </div>
      </div>`).join('');
}

function selectFeature(value, field) {
    if (state.selectedFeatureValue === value && state.selectedFeatureField === field) {
        state.selectedFeatureValue = null;
        state.selectedFeatureField = null;
        selectedCellsEl.textContent = '0';
        cloneLegendOverlay.classList.add('hidden');
    } else {
        state.selectedFeatureValue = value;
        state.selectedFeatureField = field;
        const n = state.cells.filter(c => c[field] === value).length;
        selectedCellsEl.textContent = fmt(n);
        cloneLegendOverlay.classList.remove('hidden');
    }
    renderFeatureList();
    updatePlot();
}

// ─── Value Counts for a field ──────────────────────────
function valueCounts(field) {
    const m = {};
    state.cells.forEach(c => { const v = c[field]; if (v != null) m[v] = (m[v] || 0) + 1; });
    return Object.entries(m).sort((a, b) => b[1] - a[1]);
}

// ─── Plotly Traces ─────────────────────────────────────
let plotDiv = $('umap-plot');
let plotInitialized = false;

// ─── Visible cells (global filter applied) ───────────
function isCellVisible(c) {
    if (state.hiddenClusters.has(c.seurat_clusters)) return false;
    // IS_STUDY_MODE may be true but some old code might use null sample_ids, but sample_id should match.
    if (state.hiddenSamples.has(c.sample_id)) return false;
    return true;
}

function getVisibleCells() {
    return state.cells;
}

function getTraces() {
    const { colorMode, hiddenClusters, clusterIds } = state;
    const cells = getVisibleCells();

    // ── Gene Expr color mode ──────────────────────────
    if (colorMode === 'gene') {
        return getGeneExprTraces(cells);
    }

    // ── Sample color mode ─────────────────────────────
    if (colorMode === 'sample') {
        return state.sampleIds.map(sid => {
            const color = state.sampleColorMap[sid];
            const sub = cells.filter(c => c.sample_id === sid);
            return {
                type: 'scattergl', mode: 'markers', name: sid,
                x: sub.map(c => c.umap_1), y: sub.map(c => c.umap_2),
                text: sub.map(c => buildHover(c)),
                hovertemplate: '%{text}<extra></extra>',
                marker: {
                    color: sub.map(c => isCellVisible(c) ? color : 'rgba(30,40,60,0.3)'),
                    size: 3, opacity: sub.map(c => isCellVisible(c) ? 0.75 : 0.05), line: { width: 0 }
                },
                showlegend: false,
            };
        });
    }

    // ── Severity color mode ───────────────────────────
    if (colorMode === 'severity') {
        const severities = [...new Set(cells.map(c => c.severity || 'Unknown'))];
        return severities.map(sev => {
            const color = SEVERITY_COLORS[sev] || '#60a5fa';
            const sub = cells.filter(c => (c.severity || 'Unknown') === sev);
            return {
                type: 'scattergl', mode: 'markers', name: sev,
                x: sub.map(c => c.umap_1), y: sub.map(c => c.umap_2),
                text: sub.map(c => buildHover(c)),
                hovertemplate: '%{text}<extra></extra>',
                marker: {
                    color: sub.map(c => isCellVisible(c) ? color : 'rgba(30,40,60,0.3)'),
                    size: 3,
                    opacity: sub.map(c => isCellVisible(c) ? 0.8 : 0.05),
                    line: { width: 0 }
                },
                showlegend: false,
            };
        });
    }

    // ── Cluster color mode ────────────────────────────
    if (colorMode === 'cluster') {
        return clusterIds.map((id, i) => {
            const color = CLUSTER_COLORS[i % CLUSTER_COLORS.length];
            const sub = cells.filter(c => c.seurat_clusters === id);
            return {
                type: 'scattergl', mode: 'markers', name: `Cluster ${id}`,
                x: sub.map(c => c.umap_1), y: sub.map(c => c.umap_2),
                text: sub.map(c => buildHover(c)),
                hovertemplate: '%{text}<extra></extra>',
                marker: {
                    color: sub.map(c => isCellVisible(c) ? color : 'rgba(30,40,60,0.3)'),
                    size: 3, opacity: sub.map(c => isCellVisible(c) ? 0.85 : 0.05), line: { width: 0 }
                },
                showlegend: false,
            };
        });
    }

    const cfg = getActiveFieldCfg();
    if (!cfg) return [];
    if (cfg.type === 'num') return getNumericTraces(cfg);
    return getCategoricalTraces(cfg);
}

function getCategoricalTraces(cfg) {
    const { field } = cfg;
    const sorted = valueCounts(field);
    const colorMap = {};
    sorted.forEach(([v], i) => { colorMap[v] = CLUSTER_COLORS[i % CLUSTER_COLORS.length]; });

    const sel = state.selectedFeatureValue;
    const selF = state.selectedFeatureField;
    const isMatch = c => selF && sel != null && c[selF] === sel;

    if (sel != null) {
        const hi = state.cells.filter(c => isMatch(c));
        const oth = state.cells.filter(c => !isMatch(c));
        return [
            {
                type: 'scattergl', mode: 'markers',
                x: oth.map(c => c.umap_1), y: oth.map(c => c.umap_2),
                text: oth.map(c => buildHover(c, cfg)),
                hovertemplate: '%{text}<extra></extra>',
                marker: { color: 'rgba(30,45,70,0.35)', size: 3, opacity: 0.2, line: { width: 0 } },
                showlegend: false
            },
            {
                type: 'scattergl', mode: 'markers',
                x: hi.map(c => c.umap_1), y: hi.map(c => c.umap_2),
                text: hi.map(c => buildHover(c, cfg)),
                hovertemplate: '%{text}<extra></extra>',
                marker: {
                    color: '#f472b6', size: 5, opacity: 0.95,
                    line: { color: 'rgba(255,255,255,0.3)', width: 0.8 }
                },
                showlegend: false
            },
        ];
    }

    // Group: no-value (dim), top-24 colored, overflow gray
    const topKeys = new Set(sorted.slice(0, 24).map(([v]) => v));
    const noVal = state.cells.filter(c => c[field] == null);
    const inTop = state.cells.filter(c => c[field] != null && topKeys.has(c[field]));
    const others = state.cells.filter(c => c[field] != null && !topKeys.has(c[field]));

    const traces = [];
    if (noVal.length) traces.push({
        type: 'scattergl', mode: 'markers',
        x: noVal.map(c => c.umap_1), y: noVal.map(c => c.umap_2),
        text: noVal.map(c => buildHover(c, cfg)),
        hovertemplate: '%{text}<extra></extra>',
        marker: { color: noVal.map(c => isCellVisible(c) ? 'rgba(20,35,55,0.5)' : 'rgba(30,40,60,0.3)'), size: 2.5, opacity: noVal.map(c => isCellVisible(c) ? 0.15 : 0.05), line: { width: 0 } },
        showlegend: false
    });
    if (others.length) traces.push({
        type: 'scattergl', mode: 'markers',
        x: others.map(c => c.umap_1), y: others.map(c => c.umap_2),
        text: others.map(c => buildHover(c, cfg)),
        hovertemplate: '%{text}<extra></extra>',
        marker: { color: others.map(c => isCellVisible(c) ? 'rgba(100,140,175,0.5)' : 'rgba(30,40,60,0.3)'), size: 3, opacity: others.map(c => isCellVisible(c) ? 0.35 : 0.05), line: { width: 0 } },
        showlegend: false
    });

    // One trace per top value (enables individual hover labels)
    sorted.slice(0, 24).forEach(([val]) => {
        const sub = state.cells.filter(c => c[field] === val);
        traces.push({
            type: 'scattergl', mode: 'markers', name: truncate(val, 30),
            x: sub.map(c => c.umap_1), y: sub.map(c => c.umap_2),
            text: sub.map(c => buildHover(c, cfg)),
            hovertemplate: '%{text}<extra></extra>',
            marker: { color: sub.map(c => isCellVisible(c) ? colorMap[val] : 'rgba(30,40,60,0.3)'), size: 3.5, opacity: sub.map(c => isCellVisible(c) ? 0.9 : 0.05), line: { width: 0 } },
            showlegend: false
        });
    });
    return traces;
}

function getNumericTraces(cfg) {
    const { field } = cfg;
    const withVal = state.cells.filter(c => c[field] != null);
    const noVal = state.cells.filter(c => c[field] == null);
    const traces = [];
    if (noVal.length) traces.push({
        type: 'scattergl', mode: 'markers',
        x: noVal.map(c => c.umap_1), y: noVal.map(c => c.umap_2),
        text: noVal.map(c => buildHover(c, cfg)),
        hovertemplate: '%{text}<extra></extra>',
        marker: { color: 'rgba(20,35,55,0.5)', size: 2.5, opacity: noVal.map(c => isCellVisible(c) ? 0.15 : 0.05) },
        showlegend: false
    });
    if (withVal.length) traces.push({
        type: 'scattergl', mode: 'markers',
        x: withVal.map(c => c.umap_1), y: withVal.map(c => c.umap_2),
        text: withVal.map(c => buildHover(c, cfg)),
        hovertemplate: '%{text}<extra></extra>',
        marker: {
            color: withVal.map(c => c[field]),
            colorscale: 'Viridis', showscale: true,
            colorbar: {
                title: { text: cfg.label, font: { color: '#8ba3c1', size: 11 } },
                tickfont: { color: '#8ba3c1', size: 10 }, thickness: 12, len: 0.5, x: 1.02
            },
            size: 3.5, opacity: withVal.map(c => isCellVisible(c) ? 0.9 : 0.05)
        },
        showlegend: false
    });
    return traces;
}

function buildHover(c, cfg, geneVal) {
    const lines = [`<b>Cluster ${c.seurat_clusters}</b>`];
    // Merged mode では sample / severity を先頭に
    if (c.sample_id) lines.push(`Sample: <b>${c.sample_id}</b>`);
    if (c.severity) lines.push(`Severity: ${c.severity}`);
    lines.push(`UMAP: (${c.umap_1.toFixed(2)}, ${c.umap_2.toFixed(2)})`);
    // Gene Expr mode: 発現量を表示
    if (geneVal !== undefined && state.activeGene) {
        lines.push(`<b>${state.activeGene}</b>: ${geneVal > 0 ? geneVal.toFixed(3) : '0'}`);
    } else if (cfg) {
        const v = c[cfg.field];
        if (v != null) lines.push(`${cfg.label}: ${truncate(String(v), 45)}`);
    }
    if (c.TCR_CTaa) lines.push(`CDR3α: ${c._cdr3a || '—'}  CDR3β: ${c._cdr3b || '—'}`);
    if (c.TCR_CTgene && !cfg && state.colorMode !== 'gene') lines.push(`TCR: ${truncate(c.TCR_CTgene, 45)}`);
    if (c.BCR_CTgene && !cfg && state.colorMode !== 'gene') lines.push(`BCR: ${truncate(c.BCR_CTgene, 45)}`);
    return lines.join('<br>');
}

// ─── Plot ─────────────────────────────────────────────
const plotLayout = {
    paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
    margin: { l: 40, r: 20, t: 30, b: 40 },
    xaxis: {
        title: { text: 'UMAP 1', font: { color: '#4a6080', size: 12 } },
        gridcolor: 'rgba(255,255,255,0.04)', zerolinecolor: 'rgba(255,255,255,0.08)',
        tickfont: { color: '#4a6080', size: 10 }, color: '#4a6080'
    },
    yaxis: {
        title: { text: 'UMAP 2', font: { color: '#4a6080', size: 12 } },
        gridcolor: 'rgba(255,255,255,0.04)', zerolinecolor: 'rgba(255,255,255,0.08)',
        tickfont: { color: '#4a6080', size: 10 }, color: '#4a6080'
    },
    hovermode: 'closest', dragmode: 'pan',
    font: { family: 'Inter, sans-serif', color: '#8ba3c1' },
};
const plotConfig = {
    responsive: true, displayModeBar: true, displaylogo: false,
    modeBarButtonsToRemove: ['select2d', 'lasso2d', 'autoScale2d'],
    toImageButtonOptions: { format: 'png', filename: 'immunoatlas_umap', scale: 2 },
};

function buildPlot() { Plotly.newPlot(plotDiv, getTraces(), plotLayout, plotConfig); plotInitialized = true; }
function updatePlot() { if (!plotInitialized) return; Plotly.react(plotDiv, getTraces(), plotLayout, plotConfig); }

// ─── Controls ──────────────────────────────────────────
function switchToCloneSection() {
    clusterSection.style.display = 'none';
    cloneSection.style.display = '';
    renderFeatureList();
}
function resetSelection() {
    state.selectedFeatureValue = null;
    state.selectedFeatureField = null;
    selectedCellsEl.textContent = '0';
    cloneLegendOverlay.classList.add('hidden');
}

function bindControls() {
    // ── VIEW MODE buttons ──────────────────────────────
    document.querySelectorAll('.mode-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            if (btn.classList.contains('mode-disabled')) return;

            document.querySelectorAll('.mode-btn').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
            state.colorMode = btn.dataset.mode;

            // Revert back to UMAP mode if we switch modes
            if (state.plotMode === 'diversity') {
                togglePlotMode('umap');
            }

            $('tcr-colorby-group').style.display = state.colorMode === 'tcr' ? '' : 'none';
            $('bcr-colorby-group').style.display = state.colorMode === 'bcr' ? '' : 'none';
            $('gene-expr-group').style.display = state.colorMode === 'gene' ? '' : 'none';

            resetSelection();

            // Gene Exprモード
            if (state.colorMode === 'gene') {
                clusterSection.style.display = 'none';
                cloneSection.style.display = 'none';
                // exprDataがまだロードされていなければロード開始
                if (!state.exprData && !state.exprLoading) {
                    loadExprData();
                } else if (state.exprData && state.activeGene) {
                    // 既にデータがあれば即再描画
                    applyGeneColor(state.activeGene);
                } else {
                    updatePlot();
                }
                return;
            }

            // Gene以外のモードに戻したときカラーリセット
            state.activeGene = '';
            state.geneColorValues = null;

            const isCloneMode = state.colorMode === 'tcr' || state.colorMode === 'bcr';

            // パネルの表示 / 非表示切り替え
            clusterSection.style.display = isCloneMode ? 'none' : '';
            const sfSection = $('sample-filter-section');
            if (sfSection && IS_STUDY_MODE) sfSection.style.display = isCloneMode ? 'none' : '';

            if (isCloneMode) {
                switchToCloneSection();
            } else {
                cloneSection.style.display = 'none';
            }
            updatePlot();
        });
    });

    // ── Color-By pills ────────────────────────────────
    document.querySelectorAll('.colorby-pill').forEach(pill => {
        pill.addEventListener('click', () => {
            const mode = pill.dataset.mode;
            const field = pill.dataset.field;
            document.querySelectorAll(`.colorby-pill[data-mode="${mode}"]`)
                .forEach(p => p.classList.remove('active'));
            pill.classList.add('active');
            if (mode === 'tcr') state.tcrColorBy = field;
            else state.bcrColorBy = field;
            resetSelection();
            const cfg = getActiveFieldCfg();
            const si = $('clone-search');
            if (si && cfg) si.placeholder = `Search… ${cfg.hint || cfg.title}`;
            state.cloneSearchQuery = '';
            if (si) si.value = '';
            const scb = $('clone-search-clear');
            if (scb) { scb.style.opacity = '0'; scb.style.pointerEvents = 'none'; }
            renderFeatureList();
            updatePlot();
        });
    });

    // ── Cluster: Select All ───────────────────────────
    $('select-all-clusters').addEventListener('click', () => {
        state.hiddenClusters.clear();
        const clegend = $('cluster-legend');
        if (clegend) clegend.querySelectorAll('.legend-item-row').forEach(el => el.classList.remove('dimmed'));
        totalCellsEl.textContent = fmt(state.cells.filter(isCellVisible).length);
        updatePlot();
    });

    // ── Sample: Select All ────────────────────────────
    const saBtn = $('select-all-samples');
    if (saBtn) {
        saBtn.addEventListener('click', () => {
            state.hiddenSamples.clear();
            document.querySelectorAll('.sample-filter-row').forEach(el => el.classList.remove('dimmed'));
            totalCellsEl.textContent = fmt(state.cells.filter(isCellVisible).length);
            updatePlot();
        });
    }

    // ── Clone search ──────────────────────────────────
    const si = $('clone-search');
    const scb = $('clone-search-clear');
    if (si) {
        si.addEventListener('input', () => {
            state.cloneSearchQuery = si.value;
            if (scb) { scb.style.opacity = si.value ? '1' : '0'; scb.style.pointerEvents = si.value ? 'auto' : 'none'; }
            renderFeatureList();
        });
    }
    if (scb) {
        scb.style.opacity = '0'; scb.style.pointerEvents = 'none';
        scb.addEventListener('click', () => {
            state.cloneSearchQuery = '';
            if (si) si.value = '';
            scb.style.opacity = '0'; scb.style.pointerEvents = 'none';
            renderFeatureList();
        });
    }

    // ── Clear selection ───────────────────────────────
    $('clear-clone').addEventListener('click', () => {
        resetSelection();
        renderFeatureList();
        updatePlot();
    });
}

// ─── Sample Switcher ──────────────────────────────────
function severityStyle(sev) {
    if (sev === 'Severe') return 'background:rgba(244,114,182,0.12);border-color:rgba(244,114,182,0.4);color:#f472b6';
    if (sev === 'Moderate') return 'background:rgba(251,191,36,0.12);border-color:rgba(251,191,36,0.4);color:#fbbf24';
    if (sev === 'Healthy') return 'background:rgba(52,211,153,0.12);border-color:rgba(52,211,153,0.4);color:#34d399';
    return 'background:rgba(99,179,237,0.1);border-color:rgba(99,179,237,0.3);color:#63b3ed';
}

async function initSampleSwitcher() {
    const btn = $('sample-switcher-btn');
    const dropdown = $('sample-dropdown');
    const chevron = $('sample-switcher-chevron');
    if (!btn || !dropdown) return;

    let studies;
    try {
        const res = await fetch('/output/studies.json');
        if (!res.ok) return;
        studies = await res.json();
    } catch (e) { return; }

    // 現在のサンプルと同じ studyGroup の兄弟を探す
    const currentStudy = studies.find(s => s.id === sampleKey);
    if (!currentStudy || !currentStudy.studyGroup) return;

    const siblings = studies.filter(s => s.studyGroup === currentStudy.studyGroup);
    if (siblings.length <= 1) return;

    // チェブロンを表示してボタンをクリッカブルに
    chevron.style.opacity = '1';
    btn.classList.add('has-siblings');

    // サンプルラベル抽出 [C141] → C141
    function sampleLabel(study) {
        const m = study.title.match(/\[([^\]]+)\]\s*$/);
        return m ? m[1] : study.id;
    }

    const groupTitle = currentStudy.studyGroupTitle || 'Samples in this study';

    // ドロップダウン内容を構築
    dropdown.innerHTML =
        `<div class="sample-dropdown-header">${groupTitle}</div>` +
        siblings.map(s => {
            const isCurrent = s.id === sampleKey;
            const label = sampleLabel(s);
            const sevBadge = s.severity
                ? `<span class="sample-dropdown-sev" style="${severityStyle(s.severity)}">${s.severity}</span>`
                : '';
            const checkIcon = isCurrent
                ? `<svg class="sample-dropdown-check" width="14" height="14" viewBox="0 0 14 14" fill="none"><path d="M2.5 7l3.5 3.5 5.5-6" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/></svg>`
                : '';
            return `<a class="sample-dropdown-item${isCurrent ? ' current' : ''}" href="/webapp/viewer.html?sample=${s.id}">
              <div class="sample-dropdown-id">${label}</div>
              <div class="sample-dropdown-meta">${sevBadge}<span class="sample-dropdown-cells" id="switcher-cells-${s.id}">—</span>${checkIcon}</div>
            </a>`;
        }).join('');

    // 現在サンプルのセル数は init() 完了後にセット
    setTimeout(() => {
        const el = $(`switcher-cells-${sampleKey}`);
        if (el && state.cells.length > 0) el.textContent = state.cells.length.toLocaleString() + ' cells';
    }, 1500);

    // トグル
    btn.addEventListener('click', e => {
        e.stopPropagation();
        const isOpen = dropdown.classList.contains('open');
        dropdown.classList.toggle('open', !isOpen);
        chevron.style.transform = isOpen ? '' : 'rotate(180deg)';
    });

    // 外クリックで閉じる
    document.addEventListener('click', () => {
        dropdown.classList.remove('open');
        chevron.style.transform = '';
    });
}

// ─── qs2 Download ─────────────────────────────────────
function initDownloadButton() {
    const btn = $('btn-download-qs2');
    if (!btn) return;
    const qs2Path = SAMPLE_CFG.qs2Path;
    if (!qs2Path) { btn.disabled = true; return; }

    fetch(qs2Path, { method: 'HEAD' })
        .then(res => {
            if (!res.ok) throw new Error();
            btn.disabled = false;
            const sz = res.headers.get('content-length');
            if (sz) btn.title = `Download Seurat object (.qs2)  ${(+sz / 1024 / 1024).toFixed(1)} MB`;
        })
        .catch(() => { btn.disabled = true; btn.title = 'qs2 not found'; btn.classList.add('unavailable'); });

    btn.addEventListener('click', () => {
        if (btn.disabled) return;
        const a = document.createElement('a');
        a.href = qs2Path;
        a.download = qs2Path.split('/').pop();
        document.body.appendChild(a); a.click(); document.body.removeChild(a);
        const orig = btn.innerHTML;
        btn.innerHTML = `<svg width="14" height="14" viewBox="0 0 14 14" fill="none">
          <path d="M7 1v8M4 6l3 3 3-3" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round"/>
          <path d="M2 11h10" stroke="currentColor" stroke-width="1.8" stroke-linecap="round"/>
        </svg> Downloading…`;
        btn.classList.add('downloading');
        setTimeout(() => { btn.innerHTML = orig; btn.classList.remove('downloading'); }, 2500);
    });
}

// ─── Gene Expression Feature Plot ─────────────────────

// exprDataのURL計算 (sample or merged study)
function getExprDataPath() {
    if (IS_STUDY_MODE) {
        // Merged study用: e.g. /output/expr_data_gse145926_merged.json
        return `/output/expr_data_${STUDY_ID.toLowerCase()}_merged.json`;
    }

    // MVPサンプル
    if (sampleKey === 'mvp') {
        return '/output/expr_data_mvp.json';
    }

    // GSE145926 sampleKey例: 'gse145926_c141' → '/output/expr_data_GSE145926_C141.json'
    if (sampleKey.startsWith('gse145926_')) {
        const sid = sampleKey.split('_')[1].toUpperCase(); // C141
        return `/output/expr_data_GSE145926_${sid}.json`;
    }

    return null;
}

// expr_data JSON を1回だけfetch
async function loadExprData() {
    const path = getExprDataPath();
    const statusEl = $('gene-expr-status');
    if (!path) {
        if (statusEl) { statusEl.className = 'gene-expr-status error'; statusEl.textContent = 'Gene expr data not available for merged study.'; }
        return;
    }
    state.exprLoading = true;
    if (statusEl) {
        statusEl.className = 'gene-expr-status loading';
        statusEl.innerHTML = '<span class="spin-sm"></span> Loading gene expression data…';
    }
    try {
        const res = await fetch(path);
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        state.exprData = await res.json();
        state.exprLoading = false;
        if (statusEl) {
            statusEl.className = 'gene-expr-status ok';
            statusEl.textContent = `✓ ${state.exprData.genes.length.toLocaleString()} genes loaded`;
        }
        // 既に遺伝子名が入力されていれば描画
        const geneInput = $('gene-input');
        if (geneInput && geneInput.value.trim()) {
            applyGeneColor(geneInput.value.trim());
        }
    } catch (e) {
        state.exprLoading = false;
        if (statusEl) {
            statusEl.className = 'gene-expr-status error';
            statusEl.textContent = '⚠ expr_data not found. Run process_geo_sample.R first.';
        }
    }
}

// 遺伝子名から発現量配列を構築し、描画を更新する
function applyGeneColor(geneName) {
    const statusEl = $('gene-expr-status');
    const scaleEl = $('gene-expr-scale');
    const scaleMaxEl = $('gene-scale-max');

    if (!state.exprData) {
        if (!state.exprLoading) loadExprData();
        return;
    }

    // 大文字小文字を無視して遺伝子名を検索
    const geneKey = state.exprData.genes.find(
        g => g.toLowerCase() === geneName.toLowerCase()
    );

    if (!geneKey) {
        state.activeGene = '';
        state.geneColorValues = null;
        if (statusEl) {
            statusEl.className = 'gene-expr-status error';
            statusEl.textContent = `"${geneName}" not found in expr data`;
        }
        if (scaleEl) scaleEl.style.display = 'none';

        // Reset to default grey uniformly
        updatePlot();
        return;
    }

    state.activeGene = geneKey;
    const exprDataRaw = state.exprData.expr[geneKey] || {};

    // cells と同順の発現量配列を構築 (O(n))
    state.geneColorValues = new Array(state.cells.length).fill(0);

    // 新形式 (sparse arrays: [[indices], [values]]) または旧形式 (dict) に対応
    if (Array.isArray(exprDataRaw) && exprDataRaw.length === 2 && Array.isArray(exprDataRaw[0])) {
        const indices = exprDataRaw[0];
        const values = exprDataRaw[1];
        for (let i = 0; i < indices.length; i++) {
            if (indices[i] < state.geneColorValues.length) {
                state.geneColorValues[indices[i]] = values[i];
            }
        }
    } else if (Array.isArray(exprDataRaw) && exprDataRaw.length === 2 && typeof exprDataRaw[0] === 'number') {
        // auto_unbox によって要素1つの場合に [index, value] として数値スカラー化されてしまうケース
        const index = exprDataRaw[0];
        const value = exprDataRaw[1];
        if (index < state.geneColorValues.length) {
            state.geneColorValues[index] = value;
        }
    } else {
        // legacy map format { barcode: val }
        state.cells.forEach((c, i) => {
            state.geneColorValues[i] = exprDataRaw[c.barcode] || 0;
        });
    }

    const maxVal = Math.max(...state.geneColorValues, 0.01);

    if (statusEl) {
        statusEl.className = 'gene-expr-status ok';
        const nonzero = state.geneColorValues.filter(v => v > 0).length;
        statusEl.textContent = `${geneKey} · ${nonzero.toLocaleString()} cells expressed`;
    }
    if (scaleEl) {
        scaleEl.style.display = '';
        if (scaleMaxEl) scaleMaxEl.textContent = maxVal.toFixed(2);
    }

    // Plotly.restyle の代わりに確実に WebGL が再描画される Plotly.react (updatePlot) を使用
    updatePlot();
}

// Gene Expr用のPlotlyトレース (シングル描画)
function getGeneExprTraces(cells) {
    const vals = state.geneColorValues || new Array(state.cells.length).fill(0);
    const maxVal = Math.max(...vals, 0.01);

    // 常に1つのトレースを返し、のちの restyle で marker.color を更新できるようにする。
    return [{
        type: 'scattergl', mode: 'markers',
        x: cells.map(c => c.umap_1),
        y: cells.map(c => c.umap_2),
        text: cells.map(c => buildHover(c, null, vals[c._index])),
        hovertemplate: '%{text}<extra></extra>',
        marker: {
            color: cells.map(c => vals[c._index]),
            colorscale: [
                [0, 'rgba(40,60,90,0.3)'],
                [0.000001, 'rgba(40,60,90,0.3)'],
                [0.01, '#1e3a5f'],
                [0.3, '#0e5f5f'],
                [0.6, '#0a9a9a'],
                [0.8, '#17c8e8'],
                [1.0, '#22d3ee'],
            ],
            cmin: 0,
            cmax: maxVal,
            showscale: false,
            size: 3.5,
            opacity: cells.map(c => isCellVisible(c) ? 0.95 : 0.05),
            line: { width: 0 },
        },
        showlegend: false,
    }];
}

// ─── Gene Inputのイベントバインド ────────────────────
function bindGeneInput() {
    const input = $('gene-input');
    const clearBtn = $('gene-input-clear');
    const autocomplete = $('gene-autocomplete');
    if (!input) return;

    let debounceTimer = null;

    input.addEventListener('input', () => {
        const val = input.value.trim();
        // Xクリアボタンの表示/非表示
        if (clearBtn) clearBtn.classList.toggle('visible', val.length > 0);

        // オートコンプリートの表示
        if (autocomplete && state.exprData && val.length > 0) {
            const matches = state.exprData.genes
                .filter(g => g.toLowerCase().includes(val.toLowerCase()))
                .slice(0, 50); // トップ50件制限
            if (matches.length > 0) {
                autocomplete.innerHTML = matches.map(g => `<div class="gene-ac-item">${g}</div>`).join('');
                autocomplete.style.display = 'block';
            } else {
                autocomplete.style.display = 'none';
            }
        } else if (autocomplete) {
            autocomplete.style.display = 'none';
        }

        // デバウンス: 300ms後に描画
        clearTimeout(debounceTimer);
        if (!val) {
            state.activeGene = '';
            state.geneColorValues = null;
            const scaleEl = $('gene-expr-scale');
            const statusEl = $('gene-expr-status');
            if (scaleEl) scaleEl.style.display = 'none';
            if (statusEl) { statusEl.className = 'gene-expr-status'; statusEl.textContent = ''; }
            updatePlot();
            return;
        }
        debounceTimer = setTimeout(() => {
            if (state.colorMode === 'gene') applyGeneColor(val);
        }, 500); // 連続入力を考慮し500msに調整
    });

    if (autocomplete) {
        autocomplete.addEventListener('mousedown', (e) => {
            if (e.target.classList.contains('gene-ac-item')) {
                input.value = e.target.textContent;
                autocomplete.style.display = 'none';
                clearTimeout(debounceTimer);
                if (state.colorMode === 'gene') applyGeneColor(input.value);
            }
        });
    }

    input.addEventListener('blur', () => {
        if (autocomplete) autocomplete.style.display = 'none';
    });

    input.addEventListener('focus', () => {
        if (input.value.trim() && autocomplete && autocomplete.innerHTML) {
            autocomplete.style.display = 'block';
        }
    });

    input.addEventListener('keydown', e => {
        if (e.key === 'Enter') {
            clearTimeout(debounceTimer);
            if (autocomplete) autocomplete.style.display = 'none';
            const val = input.value.trim();
            if (val && state.colorMode === 'gene') applyGeneColor(val);
        }
    });

    if (clearBtn) {
        clearBtn.addEventListener('click', () => {
            input.value = '';
            clearBtn.classList.remove('visible');
            if (autocomplete) autocomplete.style.display = 'none';
            state.activeGene = '';
            state.geneColorValues = null;
            const scaleEl = $('gene-expr-scale');
            const statusEl = $('gene-expr-status');
            if (scaleEl) scaleEl.style.display = 'none';
            if (statusEl) { statusEl.className = 'gene-expr-status'; statusEl.textContent = ''; }
            updatePlot();
            input.focus();
        });
    }
}

// ─── Plot Mode Toggle (UMAP vs Diversity) ─────────────
function bindPlotModeToggle() {
    const btnUmap = $('btn-plot-umap');
    const btnDiv = $('btn-plot-diversity');
    if (!btnUmap || !btnDiv) return;

    btnUmap.addEventListener('click', () => togglePlotMode('umap'));
    btnDiv.addEventListener('click', () => togglePlotMode('diversity'));
}

function togglePlotMode(mode) {
    state.plotMode = mode;
    const btnUmap = $('btn-plot-umap');
    const btnDiv = $('btn-plot-diversity');
    if (btnUmap) btnUmap.classList.toggle('active', mode === 'umap');
    if (btnDiv) btnDiv.classList.toggle('active', mode === 'diversity');

    if (mode === 'umap') {
        $('umap-container').style.display = '';
        $('diversity-overlay').style.display = 'none';
        updatePlot();
    } else if (mode === 'diversity') {
        $('umap-container').style.display = 'none';
        $('diversity-overlay').style.display = 'flex';
        renderDiversityPlot();
    }
}

// Helper to get severity from sample/study configs
function getSeverityForSample(sid) {
    if (!sid) return 'Unknown';
    // Match based on study prefix (e.g. gse145926_c141)
    const lowerSid = sid.toLowerCase();
    const mockSid = `gse145926_${lowerSid}`; // Currently GSE only hardcode fallback mapping

    // Look up in SAMPLES if populated
    if (SAMPLES[mockSid]) {
        return SAMPLES[mockSid].severity || 'Unknown';
    }

    // Search cells arrays as fallback (expensive)
    const sampleCells = state.cells.filter(c => c.sample_id === sid);
    if (sampleCells.length > 0 && sampleCells[0].severity) {
        return sampleCells[0].severity;
    }

    return 'Unknown';
}

function renderDiversityPlot() {
    if (!state.diversityData) {
        diversityPlotContainer.innerHTML = '<div style="color:#94a3b8;padding:2rem;">Diversity data not available yet.</div>';
        return;
    }

    if (!IS_STUDY_MODE) {
        diversityPlotContainer.innerHTML = '<div style="color:#94a3b8;padding:2rem;">Diversity comparison is only available in Merged Study View.</div>';
        return;
    }

    const data = state.diversityData;
    const items = [];

    // Aggregate diversity data for available samples
    for (const sid of state.sampleIds) {
        // e.g., 'C141'
        const divKey = Object.keys(data).find(k => k.includes(sid.toUpperCase()));
        if (divKey && data[divKey].shannon !== null && data[divKey].shannon !== undefined) {
            items.push({
                sample: sid,
                severity: getSeverityForSample(sid),
                shannon: data[divKey].shannon,
                simpson: data[divKey].simpson || 0,
                color: state.sampleColorMap[sid] || '#94a3b8'
            });
        }
    }

    if (items.length === 0) {
        diversityPlotContainer.innerHTML = '<div style="color:#94a3b8;padding:2rem;">No Shannon index data available for these samples.</div>';
        return;
    }

    // Sort by severity then by Shannon index
    const severityOrder = { Severe: 0, Moderate: 1, Healthy: 2, Unknown: 9 };
    items.sort((a, b) => {
        const sevA = severityOrder[a.severity] ?? 9;
        const sevB = severityOrder[b.severity] ?? 9;
        if (sevA !== sevB) return sevA - sevB;
        return b.shannon - a.shannon;
    });

    const traces = [
        {
            x: items.map(d => d.sample),
            y: items.map(d => d.shannon),
            type: 'bar',
            name: 'Shannon Index',
            marker: {
                color: items.map(d => d.color),
                opacity: 0.9,
                line: { color: 'rgba(255,255,255,0.1)', width: 1 }
            },
            text: items.map(d => `${d.severity}<br>Shannon: ${d.shannon.toFixed(2)}`),
            hovertemplate: '<b>%{x}</b><br>%{text}<extra></extra>'
        }
    ];

    const layout = {
        title: { text: 'TCR Clonal Diversity (Shannon Index)', font: { color: '#e2e8f0', size: 14 } },
        paper_bgcolor: 'rgba(0,0,0,0)',
        plot_bgcolor: 'rgba(0,0,0,0)',
        margin: { l: 70, r: 20, t: 70, b: 60 },
        xaxis: {
            tickfont: { color: '#8ba3c1', size: 11 },
            gridcolor: 'rgba(255,255,255,0.04)'
        },
        yaxis: {
            title: { text: 'Shannon Index', font: { color: '#8ba3c1', size: 12 } },
            tickfont: { color: '#8ba3c1', size: 11 },
            gridcolor: 'rgba(255,255,255,0.08)'
        },
        font: { family: 'Inter, sans-serif' },
        showlegend: false
    };

    Plotly.newPlot(diversityPlotContainer, traces, layout, { ...plotConfig, responsive: true });
}

// ─── Entry ────────────────────────────────────────────
document.addEventListener('DOMContentLoaded', () => {
    init();
    initDownloadButton();
    bindGeneInput();
    bindPlotModeToggle();
    // Sample switcher: Study modeでは不要だが個別サンプルモードでは引き続き動作
    if (!IS_STUDY_MODE) initSampleSwitcher();
    // Study modeではバッジをMerged表示に変更
    if (IS_STUDY_MODE) {
        const badge = $('sample-badge');
        if (badge) {
            badge.textContent = 'Merged Study View';
            badge.style.background = 'rgba(167,139,250,0.12)';
            badge.style.borderColor = 'rgba(167,139,250,0.4)';
            badge.style.color = '#a78bfa';
        }

        const modeControls = $('plot-mode-controls');
        if (modeControls) {
            modeControls.style.display = 'flex';
        }
    }
});
