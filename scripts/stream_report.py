#!/usr/bin/env python3
"""Generate the self-contained real-time dashboard for stream_spades.py."""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
import html
import json
import os
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence


HTML_ASSET_DIR = Path(__file__).resolve().parent.parent / "data" / "html"
RESOURCE_SPECS = (
    ("<script src=\"/publicdata/js/vue.global.prod.js\"></script>", "js/vue.global.prod.js", "script"),
    ("<script src=\"/publicdata/js/primevue.min.js\"></script>", "js/primevue.min.js", "script"),
    ("<script src=\"/publicdata/js/aura.js\"></script>", "js/aura.js", "script"),
    ("<link rel=\"stylesheet\" href=\"/publicdata/css/primeicons.css\">", "css/primeicons.css", "primeicons"),
)


def _number(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _integer(value: Any, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _local_timestamp(value: Any) -> str:
    """Convert a stored timestamp to the machine's local timezone."""
    if value in (None, ""):
        return ""
    text = str(value)
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return text
    if parsed.tzinfo is None:
        parsed = parsed.astimezone()
    return parsed.astimezone().isoformat(timespec="seconds")


def _record_timestamp(record: Dict[str, Any], name: str) -> str:
    return _local_timestamp(record.get(name, record.get(f"{name}_utc", "")))


def _resolve_path(value: Any, base_dir: Path) -> Path:
    path = Path(str(value or ""))
    return path if path.is_absolute() else base_dir / path


def _relative_path(value: Any, base_dir: Path) -> str:
    if value in (None, ""):
        return ""
    try:
        return os.path.relpath(_resolve_path(value, base_dir).resolve(), base_dir.resolve())
    except (OSError, ValueError):
        return str(value)


def _is_human_pathogen(row: Dict[str, str]) -> bool:
    value = str(row.get("HUMAN_PATHOGEN", "")).strip().lower()
    return value in {"yes", "true", "1", "y"}


def _taxon_rows(path: Path) -> List[Dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader if str(row.get("NAME", "")).strip()]


def _change_count(history: List[Dict[str, Any]], field: str) -> int:
    return sum(
        history[index][field] != history[index - 1][field]
        for index in range(1, len(history))
    )


def build_payload(state: Dict[str, Any], base_dir: Optional[Path] = None) -> Dict[str, Any]:
    """Build browser data from cumulative stream state.

    ``base_dir`` is both the origin for relative state paths and the directory relative
    to which paths are exposed in the report.
    """
    base_dir = (base_dir or Path.cwd()).resolve()
    records = sorted(state.get("timepoints", []), key=lambda item: int(item["timepoint"]))
    timepoints: List[Dict[str, Any]] = []
    pathogens: Dict[str, Dict[str, Any]] = {}

    for record in records:
        timepoint = int(record["timepoint"])
        result_value = record.get("pathogen_full_tsv", "")
        result_path = _resolve_path(result_value, base_dir)
        completed_at = _record_timestamp(record, "completed_at")
        raw_reads = record.get("raw_reads")
        filtered_reads = record.get("filtered_reads")
        input_path = _relative_path(record.get("input_file", ""), base_dir)
        timepoints.append(
            {
                "timepoint": timepoint,
                "input_file": Path(input_path).name,
                "input_path": input_path,
                "raw_reads": _integer(raw_reads) if raw_reads not in (None, "") else None,
                "filtered_reads": _integer(filtered_reads) if filtered_reads not in (None, "") else None,
                "cumulative_raw_reads": record.get("cumulative_raw_reads"),
                "cumulative_filtered_reads": record.get("cumulative_filtered_reads"),
                "status": str(record.get("status", "unknown")),
                "completed_at": completed_at,
                "pathogen_full_tsv": _relative_path(result_value, base_dir),
                "run_log": _relative_path(record.get("run_log", ""), base_dir),
            }
        )

        for row in _taxon_rows(result_path):
            taxid = str(row.get("TAXID", "")).strip()
            name = str(row.get("NAME", "Unknown species")).strip() or "Unknown species"
            level = str(row.get("LEVEL", "unknown")).strip().lower() or "unknown"
            key = f"{level}:{taxid or name}"
            pathogen = pathogens.setdefault(
                key,
                {
                    "key": key,
                    "taxid": taxid,
                    "name": name,
                    "level": level,
                    "human_pathogen": _is_human_pathogen(row),
                    "history": [],
                },
            )
            pathogen["history"].append(
                {
                    "timepoint": timepoint,
                    "completed_at": completed_at,
                    "read_count": _integer(row.get("READ_COUNT")),
                    "best_sig_cov": _number(row.get("BEST_SIG_COV", row.get("SIG_COV", 0.0))),
                    "sni_score": _number(row.get("SNI_SCORE")),
                }
            )

    latest_timepoint = int(records[-1]["timepoint"]) if records else 0
    pathogen_list: List[Dict[str, Any]] = []
    for pathogen in pathogens.values():
        history = sorted(pathogen["history"], key=lambda item: item["timepoint"])
        pathogen["history"] = history
        pathogen["latest"] = history[-1]
        pathogen["first_seen"] = history[0]["timepoint"]
        pathogen["last_seen"] = history[-1]["timepoint"]
        pathogen["present_latest"] = history[-1]["timepoint"] == latest_timepoint
        pathogen["change_counts"] = {
            "read_count": _change_count(history, "read_count"),
            "best_sig_cov": _change_count(history, "best_sig_cov"),
            "sni_score": _change_count(history, "sni_score"),
        }
        pathogen_list.append(pathogen)
    pathogen_list.sort(
        key=lambda item: (
            not item["present_latest"],
            -item["latest"]["read_count"],
            item["name"].lower(),
        )
    )
    latest_taxa = [item for item in pathogen_list if item["present_latest"]]
    pathogen_species = [
        item
        for item in pathogen_list
        if item["level"] == "species" and item["human_pathogen"]
    ]
    latest_pathogen_species = [item for item in pathogen_species if item["present_latest"]]

    raw_values = [item["raw_reads"] for item in timepoints if item["raw_reads"] is not None]
    filtered_values = [
        item["filtered_reads"] for item in timepoints if item["filtered_reads"] is not None
    ]
    pending = state.get("pending") or {}
    monitor_status = str(state.get("monitor_status", "stopped"))
    if pending:
        current_file = Path(str(pending.get("signature", {}).get("path", ""))).name
        current = {
            "label": f"Processing batch {int(pending.get('timepoint', 0)):06d}",
            "detail": current_file or "Pipeline work in progress",
            "tone": "active" if monitor_status == "running" else "warning",
        }
    elif monitor_status == "running":
        current = {
            "label": "Monitoring for stable read files",
            "detail": "Ready for the next FASTA or FASTQ file",
            "tone": "active",
        }
    else:
        current = {
            "label": "Monitor stopped — results are current",
            "detail": f"Last completed batch: {latest_timepoint:06d}" if records else "No completed batches",
            "tone": "ready",
        }

    configuration = dict(state.get("configuration", {}))
    for field in ("input_dir", "database"):
        configuration[field] = _relative_path(configuration.get(field, ""), base_dir)

    return {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "revision": f"{len(records)}:{latest_timepoint}:{monitor_status}:{pending.get('stage', '')}",
        "configuration": configuration,
        "monitor_status": monitor_status,
        "current": current,
        "summary": {
            "timepoints": len(timepoints),
            "raw_reads": sum(raw_values),
            "filtered_reads": sum(filtered_values),
            "qc_timepoints": len(raw_values),
            "taxa": len(pathogen_list),
            "taxa_latest": len(latest_taxa),
            "pathogen_species": len(pathogen_species),
            "pathogen_species_latest": len(latest_pathogen_species),
        },
        "latest_timepoint": latest_timepoint,
        "timepoints": timepoints,
        "pathogens": pathogen_list,
    }


def _read_resource(path: Path) -> bytes:
    try:
        return path.read_bytes()
    except OSError as error:
        raise RuntimeError(f"Unable to read local report asset {path}: {error}") from error


def _embed_primeicons_font(css: str, stylesheet: Path) -> str:
    font = _read_resource(stylesheet.parent / "fonts" / "primeicons.woff2")
    font_face = (
        "@font-face{font-family:'primeicons';font-display:block;"
        f"src:url('data:font/woff2;base64,{base64.b64encode(font).decode('ascii')}') format('woff2');"
        "font-weight:normal;font-style:normal}"
    )
    css, count = re.subn(r"@font-face\s*\{.*?\}", font_face, css, count=1, flags=re.DOTALL)
    if count != 1:
        raise ValueError("PrimeIcons stylesheet is missing its expected font declaration")
    return css


def _embed_browser_resources(template: str) -> str:
    resources = []
    for marker, relative, kind in RESOURCE_SPECS:
        source = HTML_ASSET_DIR / relative
        content = _read_resource(source).decode("utf-8")
        content = re.sub(r"/\*[#@]\s*sourceMappingURL=.*?\*/", "", content, flags=re.DOTALL)
        content = re.sub(r"^\s*//[#@]\s*sourceMappingURL=.*$", "", content, flags=re.MULTILINE)
        if kind == "primeicons":
            content = _embed_primeicons_font(content, source)
        resources.append({"kind": "script" if kind == "script" else "style", "content": content})
        if marker not in template:
            raise ValueError(f"HTML resource marker not found: {marker}")

    packed = json.dumps({"resources": resources}, separators=(",", ":")).encode("utf-8")
    payload = base64.b64encode(gzip.compress(packed, compresslevel=9, mtime=0)).decode("ascii")
    style_count = sum(item["kind"] == "style" for item in resources)
    loader = f"""<script>
window.browserResourcesReady = (async function() {{
  if (typeof DecompressionStream === 'undefined') throw new Error('This report requires a current browser.');
  const slots = Array.from({{length:{style_count}}}, () => {{ const node=document.createElement('style'); document.head.appendChild(node); return node; }});
  const binary=atob('{payload}'); const bytes=new Uint8Array(binary.length);
  for(let i=0;i<binary.length;i+=1) bytes[i]=binary.charCodeAt(i);
  const stream=new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
  const embedded=JSON.parse(await new Response(stream).text()); let styleIndex=0;
  for(const resource of embedded.resources) {{
    if(resource.kind==='style') slots[styleIndex++].textContent=resource.content;
    else {{ const script=document.createElement('script'); script.textContent=resource.content; document.head.appendChild(script); }}
  }}
}})();
</script>"""
    document = template.replace(RESOURCE_SPECS[0][0], loader, 1)
    for marker, _, _ in RESOURCE_SPECS[1:]:
        document = document.replace(marker, "", 1)
    return document


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>REPORT_TITLE</title>
  <script src="/publicdata/js/vue.global.prod.js"></script>
  <script src="/publicdata/js/primevue.min.js"></script>
  <script src="/publicdata/js/aura.js"></script>
  <link rel="stylesheet" href="/publicdata/css/primeicons.css">
  <style>
    :root { --page:#f5f7fa; --surface:#fff; --ink:#17212b; --muted:#627181; --line:#dce4ea; --teal:#197278; --teal-soft:#e5f3f1; --blue:#3f7cac; --blue-soft:#eaf2f9; --coral:#d95d39; --coral-soft:#fff0eb; --gold:#b98113; --shadow:0 16px 42px rgba(24,33,47,.08); }
    * { box-sizing:border-box; } body { margin:0; min-height:100vh; background:radial-gradient(circle at 5% 0,rgba(25,114,120,.08),transparent 28rem),var(--page); color:var(--ink); font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; }
    [v-cloak] { display:none; } button,input { font:inherit; } .app-shell { width:min(1600px,calc(100% - 32px)); margin:auto; padding:26px 0 42px; }
    .report-header { display:flex; justify-content:space-between; align-items:flex-start; gap:24px; margin-bottom:16px; } .eyebrow { color:var(--teal); font-size:.73rem; font-weight:800; letter-spacing:.09em; text-transform:uppercase; } h1 { margin:5px 0 0; font-size:clamp(1.65rem,3vw,2.35rem); line-height:1.08; } .header-meta,.live-control,.toolbar,.pills { display:flex; align-items:center; flex-wrap:wrap; gap:8px; }
    .header-meta { margin-top:13px; } .meta-pill,.count-pill { display:inline-flex; align-items:center; gap:7px; min-height:30px; border:1px solid var(--line); border-radius:999px; background:rgba(255,255,255,.84); color:var(--muted); font-size:.8rem; font-weight:700; padding:0 11px; } .status-dot { width:8px; height:8px; border-radius:50%; background:#94a3af; } .status-dot.active { background:var(--teal); box-shadow:0 0 0 5px rgba(25,114,120,.12); animation:pulse 1.7s infinite; } .status-dot.warning { background:var(--gold); } @keyframes pulse { 50% { box-shadow:0 0 0 8px rgba(25,114,120,0); } }
    .live-control { justify-content:flex-end; color:var(--muted); font-size:.8rem; font-weight:700; } .live-control label { display:flex; align-items:center; gap:7px; cursor:pointer; } .live-state { color:var(--teal); min-width:70px; text-align:right; }
    .stats-grid { display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:12px; margin-bottom:16px; } .metric-card,.panel { border:1px solid var(--line); border-radius:10px; background:rgba(255,255,255,.96); box-shadow:var(--shadow); } .metric-card { padding:16px; } .metric-label { color:var(--muted); font-size:.72rem; font-weight:800; letter-spacing:.055em; text-transform:uppercase; } .metric-value { margin-top:8px; font-size:1.65rem; font-weight:820; line-height:1; } .metric-note { margin-top:7px; color:var(--muted); font-size:.77rem; }
    .panel { overflow:hidden; } .panel-heading { display:flex; justify-content:space-between; gap:20px; align-items:center; padding:16px 18px; border-bottom:1px solid var(--line); } .panel-title { margin:0; font-size:1.12rem; } .panel-subtitle { margin:4px 0 0; color:var(--muted); font-size:.82rem; } .search-box { position:relative; min-width:270px; } .search-box i { position:absolute; left:12px; top:50%; transform:translateY(-50%); color:#82919f; } .search-input { width:100%; height:38px; border:1px solid var(--line); border-radius:8px; background:#fff; color:var(--ink); padding:0 12px 0 35px; outline:none; } .search-input:focus { border-color:var(--teal); box-shadow:0 0 0 3px rgba(25,114,120,.12); }
    .p-datatable { border:0; } .p-datatable-header { border:0!important; background:#fff!important; padding:12px 16px!important; } .p-datatable-thead>tr>th { background:#f7f9fb!important; color:#596979!important; font-size:.7rem; letter-spacing:.045em; text-transform:uppercase; } .p-datatable-tbody>tr>td { border-color:var(--line)!important; padding:.78rem .9rem!important; } .p-datatable-tbody>tr:hover { background:#f7fbfa!important; } .p-datatable-row-expansion>td { background:#f8fafb!important; padding:0!important; } .name-primary { display:block; font-weight:780; } .name-secondary { display:block; margin-top:3px; color:var(--muted); font-size:.73rem; }
    .value-cell { font-variant-numeric:tabular-nums; font-weight:720; white-space:nowrap; } .pill { display:inline-flex; align-items:center; gap:5px; margin-left:7px; border-radius:999px; padding:3px 7px; font-size:.65rem; font-weight:850; vertical-align:middle; } .pill.read { background:var(--teal-soft); color:#0f5b60; } .pill.coverage { background:var(--blue-soft); color:#315f83; } .pill.sni { background:var(--coral-soft); color:#a8462a; } .pill.zero { opacity:.58; } .pill i { font-size:.61rem; }
    .history { padding:18px; } .history-head { display:flex; justify-content:space-between; gap:18px; align-items:flex-start; margin-bottom:13px; } .history-head strong { font-size:.88rem; } .history-head span { color:var(--muted); font-size:.76rem; } .chart-grid { display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:12px; } .chart-card { min-width:0; border:1px solid var(--line); border-radius:9px; background:#fff; padding:12px; } .chart-head { display:flex; justify-content:space-between; gap:10px; align-items:baseline; } .chart-head strong { font-size:.8rem; } .chart-head span { color:var(--muted); font-size:.71rem; } .chart-svg { display:block; width:100%; height:190px; overflow:visible; } .grid-line { stroke:#e7edf1; stroke-width:1; } .chart-line { fill:none; stroke-width:3; stroke-linecap:round; stroke-linejoin:round; } .chart-point { fill:#fff; stroke-width:3; cursor:crosshair; } .chart-axis { fill:var(--muted); font-size:10px; } .chart-tooltip { min-height:20px; color:var(--muted); font-size:.7rem; text-align:center; }
    .empty { padding:32px; color:var(--muted); text-align:center; } .batch-panel { margin-top:16px; } .batch-table-wrap { overflow:auto; max-height:320px; } .batch-table { width:100%; border-collapse:collapse; font-size:.8rem; } .batch-table th { position:sticky; top:0; background:#f7f9fb; color:var(--muted); font-size:.68rem; letter-spacing:.04em; text-align:left; text-transform:uppercase; padding:10px 13px; } .batch-table td { border-top:1px solid var(--line); padding:10px 13px; } .path-cell { max-width:330px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; } .status-tag { display:inline-flex; border-radius:999px; background:var(--teal-soft); color:#0f5b60; font-size:.69rem; font-weight:800; padding:4px 8px; text-transform:capitalize; } .status-tag.no_alignments { background:#fff6df; color:#825d0d; } footer { display:flex; justify-content:space-between; flex-wrap:wrap; gap:12px; margin-top:15px; color:var(--muted); font-size:.72rem; } footer span { overflow-wrap:anywhere; }
    @media(max-width:900px){.stats-grid{grid-template-columns:repeat(2,1fr)}.chart-grid{grid-template-columns:1fr}.panel-heading,.report-header{align-items:flex-start;flex-direction:column}.live-control{justify-content:flex-start}.search-box{min-width:220px}}
    @media(max-width:560px){.app-shell{width:calc(100% - 18px);padding-top:16px}.stats-grid{grid-template-columns:1fr}.toolbar,.search-box{width:100%}}
  </style>
</head>
<body>
<main id="app" class="app-shell" v-cloak>
  <header class="report-header">
    <div><div class="eyebrow">SPADES-GOTTCHA2</div><h1>REPORT_HEADING</h1><div class="header-meta"><span class="meta-pill"><span class="status-dot" :class="payload.current.tone"></span>{{ payload.current.label }}</span><span class="meta-pill"><i class="pi pi-clock"></i>{{ formatTime(payload.generated_at) }}</span><span class="meta-pill"><i class="pi pi-database"></i>{{ payload.configuration.database || '—' }}</span></div></div>
    <div class="live-control"><label><input type="checkbox" v-model="liveEnabled"> Live updates</label><span class="live-state">{{ liveMessage }}</span><button class="p-button p-component p-button-sm p-button-outlined" type="button" @click="pollNow"><span class="pi pi-refresh"></span><span class="p-button-label">Update</span></button></div>
  </header>
  <section class="stats-grid">
    <article class="metric-card"><div class="metric-label">Latest batch</div><div class="metric-value">{{ payload.latest_timepoint ? `B${payload.latest_timepoint}` : '—' }}</div><div class="metric-note">{{ payload.summary.timepoints }} completed</div></article>
    <article class="metric-card"><div class="metric-label">Taxa in latest result</div><div class="metric-value">{{ formatNumber(payload.summary.taxa_latest) }}</div><div class="metric-note">across reported taxonomic levels</div></article>
    <article class="metric-card"><div class="metric-label">Input reads</div><div class="metric-value">{{ formatNumber(payload.summary.raw_reads) }}</div><div class="metric-note">from {{ payload.summary.qc_timepoints }} QC batches</div></article>
    <article class="metric-card"><div class="metric-label">Post-QC records</div><div class="metric-value">{{ formatNumber(payload.summary.filtered_reads) }}</div><div class="metric-note">available to cumulative profiling</div></article>
  </section>
  <section class="panel">
    <div class="panel-heading"><div><h2 class="panel-title">Latest profiling results</h2><p class="panel-subtitle">Expand a taxon to inspect metric history. Pill counts are the number of value changes.</p></div><span class="count-pill">{{ formatNumber(latestTaxa.length) }} taxa</span></div>
    <p-data-table v-model:expanded-rows="expandedRows" v-model:selection="selectedRows" v-model:filters="filters" :value="latestTaxa" data-key="key" state-storage="local" state-key="spades-stream-results" :global-filter-fields="['name','taxid','level']" filter-display="menu" selection-mode="multiple" :meta-key-selection="false" scrollable scroll-height="540px" paginator :rows="25" :rows-per-page-options="[10,25,50,100]" striped-rows removable-sort>
      <template #header><div class="toolbar"><div class="search-box"><i class="pi pi-search"></i><input class="search-input" v-model="filters.global.value" placeholder="Filter taxon name or TaxID"></div><button class="p-button p-component p-button-sm p-button-text" type="button" @click="clearFilter"><span class="pi pi-filter-slash"></span><span class="p-button-label">Clear filter</span></button></div></template>
      <template #empty><div class="empty">No taxa match the current filter, or no latest result is available yet.</div></template>
      <p-column expander style="width:3.5rem"></p-column>
      <p-column field="name" header="Taxon" sortable style="min-width:250px"><template #body="slot"><span class="name-primary">{{ slot.data.name }}</span><span class="name-secondary">{{ formatLevel(slot.data.level) }} · TaxID {{ slot.data.taxid || 'not available' }}</span></template></p-column>
      <p-column field="latest.read_count" header="READ_COUNT" sortable style="min-width:190px"><template #body="slot"><span class="value-cell">{{ formatNumber(slot.data.latest.read_count) }}</span><change-pill tone="read" label="RC" :count="slot.data.change_counts.read_count"></change-pill></template></p-column>
      <p-column field="latest.best_sig_cov" header="BEST_SIG_COV" sortable style="min-width:210px"><template #body="slot"><span class="value-cell">{{ formatCoverage(slot.data.latest.best_sig_cov) }}</span><change-pill tone="coverage" label="COV" :count="slot.data.change_counts.best_sig_cov"></change-pill></template></p-column>
      <p-column field="latest.sni_score" header="SNI_SCORE" sortable style="min-width:190px"><template #body="slot"><span class="value-cell">{{ formatSni(slot.data.latest.sni_score) }}</span><change-pill tone="sni" label="SNI" :count="slot.data.change_counts.sni_score"></change-pill></template></p-column>
      <p-column field="first_seen" header="Observed" sortable style="min-width:130px"><template #body="slot">B{{ slot.data.first_seen }}–B{{ slot.data.last_seen }}</template></p-column>
      <template #expansion="slot"><div class="history"><div class="history-head"><div><strong>{{ slot.data.name }} metric history</strong><br><span>Charts update in place whenever a new cumulative result arrives.</span></div><span>{{ slot.data.history.length }} observations</span></div><div class="chart-grid"><history-chart title="READ_COUNT" field="read_count" color="#197278" :history="slot.data.history" format="integer"></history-chart><history-chart title="BEST_SIG_COV" field="best_sig_cov" color="#3f7cac" :history="slot.data.history" format="percent"></history-chart><history-chart title="SNI_SCORE" field="sni_score" color="#d95d39" :history="slot.data.history" format="decimal"></history-chart></div></div></template>
    </p-data-table>
  </section>
  <section class="panel batch-panel"><div class="panel-heading"><div><h2 class="panel-title">Analysis batches</h2><p class="panel-subtitle">Paths are relative to this report wherever possible. Times use the local system timezone.</p></div></div><div class="batch-table-wrap"><table class="batch-table"><thead><tr><th>Batch</th><th>Input</th><th>Input reads</th><th>Post-QC</th><th>Status</th><th>Completed</th><th>Pipeline log</th></tr></thead><tbody><tr v-for="item in reversedTimepoints" :key="item.timepoint"><td><strong>B{{ String(item.timepoint).padStart(6,'0') }}</strong></td><td class="path-cell" :title="item.input_path">{{ item.input_file }}</td><td>{{ formatNullable(item.raw_reads) }}</td><td>{{ formatNullable(item.filtered_reads) }}</td><td><span class="status-tag" :class="item.status">{{ item.status.replace('_',' ') }}</span></td><td>{{ formatTime(item.completed_at) }}</td><td class="path-cell" :title="item.run_log">{{ item.run_log || '—' }}</td></tr><tr v-if="!payload.timepoints.length"><td colspan="7" class="empty">No completed batches yet.</td></tr></tbody></table></div></section>
  <footer><span>Input: {{ payload.configuration.input_dir || '—' }}</span><span>SPADES-GOTTCHA2 REALTIME report</span></footer>
</main>
<script id="report-data" type="application/json">__REPORT_DATA__</script>
<script>
window.browserResourcesReady.then(function(){
  const {createApp,ref,computed,watch,nextTick} = Vue;
  const initialPayload=JSON.parse(document.getElementById('report-data').textContent);
  const ChangePill={props:{tone:String,label:String,count:Number},template:`<span class="pill" :class="[tone,{zero:!count}]" :title="count+' value change'+(count===1?'':'s')"><i class="pi pi-history"></i>{{label}} {{count}}</span>`};
  const HistoryChart={props:{title:String,field:String,color:String,history:Array,format:String},setup(props){const hover=ref(null),width=430,height=180,margin={top:14,right:15,bottom:31,left:54};const values=computed(()=>props.history.map(p=>Number(p[props.field]||0)));const domain=computed(()=>{const vals=values.value;if(!vals.length)return [0,1];let min=props.field==='sni_score'?Math.min(...vals):0;let max=Math.max(...vals);if(props.field==='sni_score')min=Math.max(0,min-Math.max((max-min)*.15,.0005));if(max<=min)max=min+(props.format==='integer'?1:.01);return [min,max]});const x=i=>props.history.length<=1?margin.left+(width-margin.left-margin.right)/2:margin.left+i*(width-margin.left-margin.right)/(props.history.length-1);const y=v=>margin.top+(height-margin.top-margin.bottom)*(1-(v-domain.value[0])/(domain.value[1]-domain.value[0]));const points=computed(()=>values.value.map((v,i)=>`${x(i)},${y(v)}`).join(' '));const ticks=computed(()=>Array.from({length:5},(_,i)=>{const v=domain.value[0]+(domain.value[1]-domain.value[0])*i/4;return {v,y:y(v)}}));const fmt=v=>props.format==='integer'?Math.round(v).toLocaleString():props.format==='percent'?`${(v*100).toFixed(2)}%`:Number(v).toFixed(5);const hoverText=computed(()=>hover.value===null?'Hover over a point for details':`Batch ${props.history[hover.value].timepoint} · ${fmt(values.value[hover.value])}`);return{width,height,margin,values,points,ticks,x,y,fmt,hover,hoverText}},template:`<article class="chart-card"><div class="chart-head"><strong>{{title}}</strong><span>{{fmt(values[values.length-1]||0)}}</span></div><svg class="chart-svg" :viewBox="'0 0 '+width+' '+height" role="img" :aria-label="title+' history'"><g v-for="tick in ticks" :key="tick.y"><line class="grid-line" :x1="margin.left" :x2="width-margin.right" :y1="tick.y" :y2="tick.y"></line><text class="chart-axis" :x="margin.left-7" :y="tick.y+3" text-anchor="end">{{fmt(tick.v)}}</text></g><polyline class="chart-line" :stroke="color" :points="points"></polyline><g v-for="(point,index) in history" :key="point.timepoint"><circle class="chart-point" :stroke="color" :cx="x(index)" :cy="y(values[index])" r="4.5" @mouseenter="hover=index" @mouseleave="hover=null"><title>Batch {{point.timepoint}}: {{fmt(values[index])}}</title></circle><text v-if="history.length<=10 || index===0 || index===history.length-1" class="chart-axis" :x="x(index)" :y="height-10" text-anchor="middle">B{{point.timepoint}}</text></g></svg><div class="chart-tooltip">{{hoverText}}</div></article>`};
  const app=createApp({components:{ChangePill,HistoryChart},setup(){const payload=ref(initialPayload);const filters=ref({global:{value:null,matchMode:'contains'}});const selectedRows=ref([]);let savedExpanded={};try{savedExpanded=JSON.parse(localStorage.getItem('spades-stream-expanded')||'{}')}catch(e){}const expandedRows=ref(savedExpanded);const liveEnabled=ref(localStorage.getItem('spades-stream-live')!=='off');const liveMessage=ref('current');const polling=ref(false);const latestTaxa=computed(()=>payload.value.pathogens.filter(item=>item.present_latest));const reversedTimepoints=computed(()=>[...payload.value.timepoints].reverse());watch(expandedRows,value=>localStorage.setItem('spades-stream-expanded',JSON.stringify(value||{})),{deep:true});watch(liveEnabled,value=>localStorage.setItem('spades-stream-live',value?'on':'off'));
    const formatNumber=value=>Number(value||0).toLocaleString();const formatNullable=value=>value===null||value===undefined?'—':Number(value).toLocaleString();const formatCoverage=value=>`${(Number(value||0)*100).toFixed(2)}%`;const formatSni=value=>Number(value||0).toFixed(5);const formatLevel=value=>String(value||'unknown').replace(/_/g,' ').replace(/\b\w/g,char=>char.toUpperCase());const formatTime=value=>{if(!value)return '—';const date=new Date(value);return Number.isNaN(date.valueOf())?String(value):date.toLocaleString(undefined,{dateStyle:'medium',timeStyle:'medium'})};function clearFilter(){filters.value.global.value=null}
    async function applyPayload(next){if(next.generated_at===payload.value.generated_at)return;const selectedKeys=new Set(selectedRows.value.map(item=>item.key));const windowScroll={x:window.scrollX,y:window.scrollY};const scrollers=[...document.querySelectorAll('.p-datatable-table-container,.batch-table-wrap')].map(node=>({node,left:node.scrollLeft,top:node.scrollTop}));payload.value=next;selectedRows.value=next.pathogens.filter(item=>item.present_latest&&selectedKeys.has(item.key));await nextTick();scrollers.forEach(item=>{item.node.scrollLeft=item.left;item.node.scrollTop=item.top});window.scrollTo(windowScroll.x,windowScroll.y);liveMessage.value='updated'}
    function parseReport(text){const parsed=new DOMParser().parseFromString(text,'text/html');const node=parsed.getElementById('report-data');if(!node)throw new Error('Updated report data was not found');return JSON.parse(node.textContent)}
    async function iframePayload(url){return new Promise((resolve,reject)=>{const frame=document.createElement('iframe');frame.hidden=true;const timer=setTimeout(()=>{frame.remove();reject(new Error('Local report update timed out'))},12000);frame.onload=()=>{try{const node=frame.contentDocument&&frame.contentDocument.getElementById('report-data');if(!node)throw new Error('Local report data unavailable');clearTimeout(timer);const data=JSON.parse(node.textContent);frame.remove();resolve(data)}catch(error){clearTimeout(timer);frame.remove();reject(error)}};frame.src=url;document.body.appendChild(frame)})}
    async function pollNow(){if(polling.value)return;polling.value=true;liveMessage.value='checking…';try{const url=new URL(location.href);url.searchParams.set('_stream_poll',Date.now());let next;try{const response=await fetch(url.href,{cache:'no-store'});if(!response.ok)throw new Error(`HTTP ${response.status}`);next=parseReport(await response.text())}catch(error){if(location.protocol!=='file:')throw error;next=await iframePayload(url.href)}await applyPayload(next);if(liveMessage.value!=='updated')liveMessage.value='current'}catch(error){console.warn('Live update failed:',error);liveMessage.value='retrying'}finally{polling.value=false}}
    setInterval(()=>{if(liveEnabled.value)pollNow()},5000);return{payload,filters,selectedRows,expandedRows,liveEnabled,liveMessage,latestTaxa,reversedTimepoints,formatNumber,formatNullable,formatCoverage,formatSni,formatLevel,formatTime,clearFilter,pollNow}}
  });
  app.use(PrimeVue.Config,{theme:{preset:PrimeUIX.Themes.Aura}});app.component('p-data-table',PrimeVue.DataTable);app.component('p-column',PrimeVue.Column);app.mount('#app');
}).catch(function(error){console.error(error);document.getElementById('app').removeAttribute('v-cloak');document.getElementById('app').innerHTML='<div class="empty">Unable to initialize the embedded report interface.</div>'});
</script>
</body></html>
"""


def generate_stream_report(
    state: Dict[str, Any], output: Path, title: Optional[str] = None
) -> Path:
    output = output.resolve()
    prefix = str(state.get("configuration", {}).get("prefix", "SPADES stream"))
    report_title = title or f"{prefix}"
    payload = build_payload(state, output.parent)
    payload_json = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).replace(
        "</", "<\\/"
    )
    document = (
        HTML_TEMPLATE.replace("REPORT_TITLE", html.escape(report_title))
        .replace("REPORT_HEADING", html.escape(report_title))
        .replace("__REPORT_DATA__", payload_json)
    )
    document = _embed_browser_resources(document)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp")
    temporary.write_text(document, encoding="utf-8")
    os.replace(temporary, output)
    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state", type=Path, required=True, help="stream_state.json")
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument("--title")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    with args.state.open(encoding="utf-8") as handle:
        state = json.load(handle)
    prefix = str(state.get("configuration", {}).get("prefix", "stream"))
    output = args.output or args.state.parent / f"{prefix}.stream.html"
    generated = generate_stream_report(state, output, args.title)
    print(f"INFO: Stream report: {os.path.relpath(generated, Path.cwd())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
