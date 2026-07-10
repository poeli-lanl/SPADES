#!/usr/bin/env python3
"""Generate the self-contained real-time dashboard for stream_spades.py."""

from __future__ import annotations

import argparse
import csv
import html
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence


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


def _is_human_pathogen(row: Dict[str, str]) -> bool:
    value = str(row.get("HUMAN_PATHOGEN", "")).strip().lower()
    return value in {"yes", "true", "1", "y"}


def _species_rows(path: Path) -> List[Dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            row
            for row in reader
            if str(row.get("LEVEL", "")).strip().lower() == "species"
            and _is_human_pathogen(row)
        ]


def build_payload(state: Dict[str, Any]) -> Dict[str, Any]:
    records = sorted(state.get("timepoints", []), key=lambda item: int(item["timepoint"]))
    timepoints: List[Dict[str, Any]] = []
    pathogens: Dict[str, Dict[str, Any]] = {}

    for record in records:
        timepoint = int(record["timepoint"])
        result_path = Path(str(record.get("pathogen_full_tsv", "")))
        raw_reads = record.get("raw_reads")
        filtered_reads = record.get("filtered_reads")
        timepoints.append(
            {
                "timepoint": timepoint,
                "input_file": Path(str(record.get("input_file", ""))).name,
                "input_path": str(record.get("input_file", "")),
                "raw_reads": _integer(raw_reads) if raw_reads not in (None, "") else None,
                "filtered_reads": (
                    _integer(filtered_reads) if filtered_reads not in (None, "") else None
                ),
                "cumulative_raw_reads": record.get("cumulative_raw_reads"),
                "cumulative_filtered_reads": record.get("cumulative_filtered_reads"),
                "status": str(record.get("status", "unknown")),
                "completed_at_utc": str(record.get("completed_at_utc", "")),
                "pathogen_full_tsv": str(result_path),
            }
        )

        for row in _species_rows(result_path):
            taxid = str(row.get("TAXID", "")).strip()
            name = str(row.get("NAME", "Unknown species")).strip() or "Unknown species"
            key = taxid or name
            pathogen = pathogens.setdefault(
                key,
                {
                    "taxid": taxid,
                    "name": name,
                    "history": [],
                },
            )
            point = {
                "timepoint": timepoint,
                "completed_at_utc": str(record.get("completed_at_utc", "")),
                "read_count": _integer(row.get("READ_COUNT")),
                "signature_coverage": _number(
                    row.get("SIG_COV", row.get("BEST_SIG_COV", 0.0))
                ),
                "sni_score": _number(row.get("SNI_SCORE")),
            }
            pathogen["history"].append(point)

    pathogen_list = []
    latest_timepoint = int(records[-1]["timepoint"]) if records else 0
    for pathogen in pathogens.values():
        pathogen["history"].sort(key=lambda item: item["timepoint"])
        pathogen["latest"] = pathogen["history"][-1]
        pathogen["first_seen"] = pathogen["history"][0]["timepoint"]
        pathogen["last_seen"] = pathogen["history"][-1]["timepoint"]
        pathogen["present_latest"] = pathogen["last_seen"] == latest_timepoint
        previous = (
            pathogen["history"][-2]
            if len(pathogen["history"]) > 1
            else pathogen["history"][-1]
        )
        pathogen["change"] = {
            "read_count": pathogen["latest"]["read_count"] - previous["read_count"],
            "signature_coverage": (
                pathogen["latest"]["signature_coverage"]
                - previous["signature_coverage"]
            ),
            "sni_score": pathogen["latest"]["sni_score"] - previous["sni_score"],
        }
        pathogen_list.append(pathogen)
    pathogen_list.sort(
        key=lambda item: (
            not item["present_latest"],
            -item["latest"]["read_count"],
            item["name"].lower(),
        )
    )
    latest_pathogens = [item for item in pathogen_list if item["present_latest"]]

    raw_values = [item["raw_reads"] for item in timepoints if item["raw_reads"] is not None]
    filtered_values = [
        item["filtered_reads"] for item in timepoints if item["filtered_reads"] is not None
    ]
    pending = state.get("pending") or {}
    monitor_status = str(state.get("monitor_status", "stopped"))

    if pending:
        current_file = Path(str(pending.get("signature", {}).get("path", ""))).name
        if monitor_status == "running":
            activity_label = f"Processing timepoint {int(pending.get('timepoint', 0)):06d}"
            activity_tone = "active"
            activity_detail = current_file or "Pipeline work in progress"
        else:
            activity_label = f"Timepoint {int(pending.get('timepoint', 0)):06d} interrupted"
            activity_tone = "warning"
            activity_detail = f"{current_file or 'Input'} will be recovered on restart"
    elif monitor_status == "running":
        activity_label = "Monitoring for stable read files"
        activity_tone = "active"
        activity_detail = "Ready to process the next FASTA or FASTQ file"
    else:
        activity_label = "Monitor stopped — results are current"
        activity_tone = "ready"
        activity_detail = (
            f"Last completed timepoint: {records[-1]['timepoint']:06d}"
            if records
            else "No timepoints have completed yet"
        )

    non_profiled = sum(1 for item in timepoints if item["status"] != "profiled")
    overall_tone = "warning" if non_profiled else "ready"
    overall_label = (
        f"{len(timepoints)} timepoint{'s' if len(timepoints) != 1 else ''} complete"
    )
    overall_detail = (
        f"{non_profiled} completed without qualifying alignments"
        if non_profiled
        else "All completed timepoints have cumulative profiles"
    )

    if latest_pathogens:
        finding_label = (
            f"Sequence evidence detected for {len(latest_pathogens)} human-pathogenic species"
            if len(latest_pathogens) != 1
            else "Sequence evidence detected for 1 human-pathogenic species"
        )
        finding_detail = (
            "Detected in the latest cumulative sequence analysis. Review the "
            "supporting evidence and correlate with the clinical and epidemiologic context."
        )
        finding_tone = "detected"
    elif records:
        finding_label = "No human-pathogenic species detected"
        finding_detail = (
            "No annotated human-pathogenic species met the displayed profiling result "
            "criteria in the latest cumulative analysis."
        )
        finding_tone = "clear"
    else:
        finding_label = "No completed analysis yet"
        finding_detail = "The latest finding will appear after the first timepoint completes."
        finding_tone = "pending"

    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "configuration": state.get("configuration", {}),
        "monitor_status": monitor_status,
        "current": {
            "label": activity_label,
            "detail": activity_detail,
            "tone": activity_tone,
        },
        "overall": {
            "label": overall_label,
            "detail": overall_detail,
            "tone": overall_tone,
        },
        "finding": {
            "label": finding_label,
            "detail": finding_detail,
            "tone": finding_tone,
            "latest_timepoint": latest_timepoint,
            "completed_at_utc": (
                str(records[-1].get("completed_at_utc", "")) if records else ""
            ),
        },
        "summary": {
            "timepoints": len(timepoints),
            "raw_reads": sum(raw_values),
            "filtered_reads": sum(filtered_values),
            "qc_timepoints": len(raw_values),
            "pathogen_species": len(pathogen_list),
            "pathogen_species_latest": len(latest_pathogens),
        },
        "timepoints": timepoints,
        "pathogens": pathogen_list,
    }


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>REPORT_TITLE</title>
  <style>
    :root {
      --page-bg: #f4f7f9; --surface: #fff; --surface-muted: #f2f6f6;
      --ink: #17212b; --muted: #5d6c79; --line: #d9e3e8;
      --accent: #156f73; --accent-dark: #0d585c; --accent-soft: #e5f3f1;
      --finding: #a74729; --finding-soft: #fff2ed; --finding-line: #ebbeaF;
      --gold: #a87509; --gold-soft: #fff8e7; --blue: #356f9e; --blue-soft: #edf5fb;
      --shadow: 0 14px 38px rgba(22, 33, 43, .075); --radius: 10px;
    }
    * { box-sizing: border-box; }
    body { min-height: 100vh; margin: 0; background: radial-gradient(circle at 8% 0%, rgba(21,111,115,.08), transparent 27rem), var(--page-bg); color: var(--ink); font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }
    button, input { font: inherit; }
    .app-shell { width: min(1380px, calc(100% - 32px)); margin: 0 auto; padding: 26px 0 42px; }
    .report-header { display: grid; grid-template-columns: minmax(0,1fr) auto; gap: 24px; align-items: start; margin-bottom: 18px; }
    .eyebrow { color: var(--accent-dark); font-size: .74rem; font-weight: 820; letter-spacing: .09em; margin-bottom: 7px; text-transform: uppercase; }
    h1 { margin: 0; font-size: clamp(1.8rem,4vw,2.4rem); font-weight: 790; line-height: 1.08; }
    .subtitle { margin: 9px 0 0; color: var(--muted); font-size: .97rem; line-height: 1.5; }
    .run-id { color: var(--ink); font-weight: 780; }
    .header-meta { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 13px; }
    .meta-pill { display: inline-flex; align-items: center; gap: 8px; min-height: 32px; border: 1px solid var(--line); border-radius: 999px; background: rgba(255,255,255,.85); color: var(--muted); font-size: .82rem; font-weight: 720; padding: 0 12px; }
    .status-dot { width: 8px; height: 8px; border-radius: 50%; background: var(--accent); box-shadow: 0 0 0 4px rgba(21,111,115,.12); }
    .status-dot.active { animation: pulse 1.8s infinite; }
    .status-dot.warning { background: var(--gold); box-shadow: 0 0 0 4px rgba(168,117,9,.15); }
    @keyframes pulse { 50% { box-shadow: 0 0 0 7px rgba(21,111,115,0); } }
    .refresh-control { display:flex; align-items:center; gap:9px; border:1px solid var(--line); border-radius:var(--radius); background:var(--surface); box-shadow:var(--shadow); padding:10px 12px; color:var(--muted); font-size:.82rem; font-weight:700; }
    .refresh-control input { accent-color: var(--accent); }

    .clinical-banner { display:grid; grid-template-columns:auto minmax(0,1fr) minmax(280px,.72fr); gap:17px; align-items:center; border:1px solid var(--finding-line); border-left:5px solid var(--finding); border-radius:var(--radius); background:linear-gradient(110deg,var(--finding-soft),#fff 58%); box-shadow:var(--shadow); padding:19px 20px; margin-bottom:18px; }
    .clinical-banner.clear { border-color:#b9d9d4; border-left-color:var(--accent); background:linear-gradient(110deg,var(--accent-soft),#fff 58%); }
    .clinical-banner.pending { border-color:var(--line); border-left-color:#8a98a4; background:var(--surface); }
    .finding-icon { display:grid; place-items:center; width:48px; height:48px; border-radius:50%; background:#fff; border:1px solid var(--finding-line); color:var(--finding); font-size:1.2rem; font-weight:900; }
    .clear .finding-icon { border-color:#b9d9d4; color:var(--accent-dark); }
    .finding-label { margin:2px 0 0; font-size:1.25rem; font-weight:830; line-height:1.25; }
    .finding-detail { max-width:740px; margin:6px 0 0; color:var(--muted); font-size:.9rem; line-height:1.5; }
    .finding-meta { display:flex; flex-wrap:wrap; gap:8px; margin-top:10px; }
    .finding-chip { border-radius:999px; background:#fff; border:1px solid var(--line); color:var(--muted); font-size:.75rem; font-weight:760; padding:5px 9px; }
    .clinical-caution { border-left:1px solid var(--finding-line); padding-left:17px; color:#69483b; font-size:.82rem; line-height:1.48; }
    .clinical-caution strong { display:block; color:var(--finding); margin-bottom:3px; }

    .section-heading { display:flex; align-items:end; justify-content:space-between; gap:18px; margin:23px 0 11px; }
    .section-heading h2 { margin:0; font-size:1.18rem; }
    .section-heading p { margin:5px 0 0; color:var(--muted); font-size:.89rem; }
    .count-badge { border-radius:999px; background:var(--accent-soft); color:var(--accent-dark); font-size:.78rem; font-weight:820; padding:6px 11px; white-space:nowrap; }
    .pathogen-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(330px,1fr)); gap:14px; }
    .pathogen-card, .panel, .metric-panel { border:1px solid var(--line); border-radius:var(--radius); background:rgba(255,255,255,.96); box-shadow:var(--shadow); }
    .pathogen-card { overflow:hidden; border-top:4px solid var(--finding); transition:border-color .18s ease,box-shadow .18s ease,transform .18s ease; }
    .pathogen-card.historical { border-top-color:#9aa7b1; opacity:.86; }
    .pathogen-card:hover { box-shadow:0 20px 48px rgba(22,33,43,.11); transform:translateY(-1px); }
    .pathogen-card.expanded { grid-column:1/-1; border-color:rgba(21,111,115,.52); opacity:1; }
    .pathogen-toggle { width:100%; border:0; background:transparent; color:inherit; cursor:pointer; padding:18px; text-align:left; }
    .pathogen-top { display:flex; justify-content:space-between; gap:14px; align-items:start; }
    .pathogen-status { display:inline-flex; align-items:center; gap:6px; border-radius:999px; background:var(--finding-soft); color:var(--finding); font-size:.7rem; font-weight:850; letter-spacing:.045em; padding:5px 8px; text-transform:uppercase; }
    .historical .pathogen-status { background:#eef1f3; color:#60707d; }
    .pathogen-name { margin:9px 0 0; font-size:1.18rem; font-weight:830; line-height:1.22; }
    .taxid { margin-top:5px; color:var(--muted); font-size:.77rem; }
    .expand-mark { display:grid; place-items:center; width:32px; height:32px; border:1px solid var(--line); border-radius:8px; color:var(--accent); font-size:1.12rem; transition:transform .2s ease; }
    .expanded .expand-mark { transform:rotate(45deg); }
    .pathogen-metrics { display:grid; grid-template-columns:repeat(3,1fr); gap:8px; margin-top:17px; }
    .mini-metric { border-radius:8px; background:var(--surface-muted); padding:11px; }
    .mini-metric:nth-child(2) { background:var(--blue-soft); }
    .mini-metric:nth-child(3) { background:var(--finding-soft); }
    .mini-label { color:var(--muted); font-size:.67rem; font-weight:820; line-height:1.3; text-transform:uppercase; }
    .mini-value { margin-top:6px; font-size:1.08rem; font-weight:840; }
    .mini-context { margin-top:3px; color:var(--muted); font-size:.7rem; }
    .pathogen-foot { display:flex; flex-wrap:wrap; gap:8px; align-items:center; margin-top:13px; color:var(--muted); font-size:.77rem; }
    .trend-chip { border-radius:999px; background:var(--accent-soft); color:var(--accent-dark); font-weight:780; padding:4px 8px; }
    .pathogen-detail { border-top:1px solid var(--line); background:#fbfcfd; padding:17px; }
    .evidence-note { margin:0 0 13px; color:var(--muted); font-size:.82rem; line-height:1.5; }
    .charts-grid { display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:12px; }
    .chart-box { min-width:0; border:1px solid var(--line); border-radius:9px; background:var(--surface); padding:12px; }
    .chart-title { display:flex; justify-content:space-between; gap:10px; align-items:baseline; margin-bottom:6px; }
    .chart-title strong { font-size:.84rem; }
    .chart-title span { color:var(--muted); font-size:.73rem; }
    .line-chart { display:block; width:100%; height:auto; overflow:visible; }
    .axis-label { fill:var(--muted); font-size:10px; }
    .grid-line { stroke:#e7edf1; stroke-width:1; }
    .axis-line { stroke:#bdc9d2; stroke-width:1; }
    .empty-state { border:1px dashed var(--line); border-radius:var(--radius); background:rgba(255,255,255,.7); color:var(--muted); padding:30px; text-align:center; }

    .operations { margin-top:25px; border-top:1px solid var(--line); padding-top:1px; }
    .metric-label { color:var(--muted); font-size:.72rem; font-weight:820; letter-spacing:.05em; text-transform:uppercase; }
    .stats-grid { display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:12px; margin-bottom:13px; }
    .metric-panel { min-height:104px; padding:16px; box-shadow:0 8px 24px rgba(22,33,43,.05); }
    .metric-value { margin-top:8px; font-size:clamp(1.42rem,2.5vw,1.9rem); font-weight:830; line-height:1; }
    .metric-subtle { margin-top:7px; color:var(--muted); font-size:.78rem; line-height:1.35; }
    .qc-explainer { border:1px solid #d3e2eb; border-radius:8px; background:var(--blue-soft); color:#3f5d72; font-size:.8rem; line-height:1.45; padding:10px 12px; margin-bottom:12px; }
    .panel { overflow:hidden; margin-top:10px; }
    .table-wrap { overflow-x:auto; }
    table { width:100%; border-collapse:collapse; font-size:.84rem; }
    th { background:#f7f9fb; color:var(--muted); font-size:.69rem; font-weight:830; letter-spacing:.045em; padding:11px 14px; text-align:left; text-transform:uppercase; white-space:nowrap; }
    td { border-top:1px solid var(--line); padding:12px 14px; vertical-align:middle; }
    .file-cell { max-width:340px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }
    .status-tag { display:inline-flex; border-radius:999px; background:var(--accent-soft); color:var(--accent-dark); font-size:.72rem; font-weight:820; padding:5px 9px; text-transform:capitalize; }
    .status-tag.no_alignments { background:var(--gold-soft); color:#805b0a; }
    .limitations { display:grid; grid-template-columns:auto 1fr; gap:11px; align-items:start; border:1px solid #ead5cb; border-radius:var(--radius); background:#fffaf7; color:#654c42; margin-top:18px; padding:14px 16px; font-size:.82rem; line-height:1.52; }
    .limitations strong { color:#7b3c28; }
    footer { display:flex; justify-content:space-between; gap:18px; flex-wrap:wrap; margin-top:17px; color:var(--muted); font-size:.75rem; }
    footer > span { min-width:0; max-width:100%; overflow-wrap:anywhere; }
    #database-name { word-break:break-all; }

    @media (max-width:980px) { .clinical-banner { grid-template-columns:auto 1fr; } .clinical-caution { grid-column:1/-1; border-left:0; border-top:1px solid var(--finding-line); padding:12px 0 0; } .stats-grid { grid-template-columns:repeat(2,1fr); } .charts-grid { grid-template-columns:1fr; } }
    @media (max-width:680px) { .app-shell { width:min(100% - 20px,1380px); padding-top:18px; } .report-header { grid-template-columns:1fr; } .refresh-control { justify-self:start; } .clinical-banner { grid-template-columns:1fr; } .finding-icon { width:40px; height:40px; } .stats-grid { grid-template-columns:1fr; } .pathogen-grid { grid-template-columns:1fr; } .pathogen-metrics { grid-template-columns:1fr; } .section-heading { align-items:start; flex-direction:column; } }
  </style>
</head>
<body>
  <main class="app-shell">
    <header class="report-header">
      <div>
        <div class="eyebrow">SPADES-GOTTCHA2</div>
        <h1>Pathogen screening summary</h1>
        <p class="subtitle">Run <span class="run-id" id="run-id">—</span></p>
        <div class="header-meta">
          <span class="meta-pill" id="header-status-pill"><span id="header-status-dot" class="status-dot"></span><span id="header-status">Loading analysis status…</span></span>
          <span class="meta-pill" id="header-overall-pill"><span aria-hidden="true">✓</span><span id="header-overall">Loading cumulative status…</span></span>
          <span class="meta-pill">Latest batch <strong id="latest-batch">—</strong></span>
          <span class="meta-pill">Updated <span id="generated-at">—</span></span>
        </div>
      </div>
      <label class="refresh-control" title="Reload this page to pick up newly generated results"><input id="auto-refresh" type="checkbox" checked> Live updates <span id="refresh-countdown">15s</span></label>
    </header>

    <section id="clinical-banner" class="clinical-banner" aria-live="polite">
      <div class="finding-icon" id="finding-icon">!</div>
      <div><div class="eyebrow">Latest cumulative finding</div><div class="finding-label" id="finding-label">Loading…</div><p class="finding-detail" id="finding-detail"></p><div class="finding-meta"><span class="finding-chip" id="finding-batch">No completed batch</span><span class="finding-chip" id="finding-time">—</span></div></div>
      <aside class="clinical-caution"><strong>Screening result — not a standalone diagnosis</strong>Interpret with the specimen, symptoms, exposure history, epidemiology, and confirmatory laboratory testing according to institutional protocol.</aside>
    </section>

    <div class="section-heading"><div><h2>Organisms requiring review</h2><p>Human-pathogenic species annotated in cumulative results. Findings in the latest analysis appear first.</p></div><span class="count-badge" id="pathogen-count">0 in latest result</span></div>
    <section id="pathogen-grid" class="pathogen-grid" aria-live="polite"></section>

    <section class="operations" aria-label="Analysis and data quality details">
      <div class="section-heading"><div><h2>Analysis and data quality</h2><p>Sequencing input and batch details for laboratory review.</p></div></div>
      <div class="stats-grid">
        <article class="metric-panel"><div class="metric-label">Batches completed</div><div class="metric-value" id="metric-timepoints">0</div><div class="metric-subtle">processed one at a time</div></article>
        <article class="metric-panel"><div class="metric-label">Input reads</div><div class="metric-value" id="metric-raw">0</div><div class="metric-subtle" id="metric-raw-note">from QC reports</div></article>
        <article class="metric-panel"><div class="metric-label">Post-QC read records</div><div class="metric-value" id="metric-filtered">0</div><div class="metric-subtle">after filtering and long-read splitting</div></article>
        <article class="metric-panel"><div class="metric-label">Species in latest result</div><div class="metric-value" id="metric-pathogens">0</div><div class="metric-subtle">with annotated human pathogenicity</div></article>
      </div>
      <div class="qc-explainer"><strong>Why can post-QC records exceed input reads?</strong> Long Nanopore reads can be split into multiple analysis records after filtering. These are processing records, not additional specimens or independent biological observations.</div>

      <div class="section-heading"><div><h2>Analysis batches</h2><p>Input-level QC counts and cumulative profiling status.</p></div></div>
      <section class="panel table-wrap"><table><thead><tr><th>Batch</th><th>Input file</th><th>Input reads</th><th>Post-QC records</th><th>Profile status</th><th>Completed (UTC)</th></tr></thead><tbody id="timepoint-body"></tbody></table></section>
    </section>

    <aside class="limitations"><span aria-hidden="true">ⓘ</span><div><strong>Interpretation limitation.</strong> Detection of sequence aligned to a pathogen signature does not establish organism viability or prove that the organism is the cause of disease. A result should be interpreted with clinical findings, specimen quality, epidemiologic information, and appropriate confirmatory testing.</div></aside>
    <footer><span>Reference database: <span id="database-name">—</span></span><span>Self-contained SPADES screening report</span></footer>
  </main>

  <script>
    const DATA = REPORT_DATA;
    const fmt = new Intl.NumberFormat('en-US');
    const fmtNumber = value => value === null || value === undefined ? '—' : fmt.format(value);
    const fmtCoverage = value => `${(Number(value || 0) * 100).toFixed(2)}%`;
    const fmtSni = value => Number(value || 0).toFixed(5);
    const fmtSniPercent = value => `${(Number(value || 0) * 100).toFixed(3)}% identity`;
    const escapeHtml = value => String(value ?? '').replace(/[&<>'"]/g, char => ({'&':'&amp;','<':'&lt;','>':'&gt;',"'":'&#39;','"':'&quot;'}[char]));
    const shortUtc = value => value ? value.replace('T',' ').replace('+00:00','Z') : '—';
    const signed = (value, digits=0) => `${value > 0 ? '+' : ''}${Number(value).toFixed(digits)}`;
    function setText(id,value) { document.getElementById(id).textContent = value; }

    function renderSummary() {
      setText('run-id', DATA.configuration.prefix || '—');
      setText('header-status', DATA.current.label);
      setText('header-overall', DATA.overall.label);
      document.getElementById('header-status-pill').title = DATA.current.detail;
      document.getElementById('header-overall-pill').title = DATA.overall.detail;
      setText('latest-batch', DATA.finding.latest_timepoint ? `T${DATA.finding.latest_timepoint}` : '—');
      setText('generated-at', shortUtc(DATA.generated_at_utc));
      setText('finding-label', DATA.finding.label); setText('finding-detail', DATA.finding.detail);
      setText('finding-batch', DATA.finding.latest_timepoint ? `Cumulative result · Batch ${DATA.finding.latest_timepoint}` : 'No completed batch');
      setText('finding-time', shortUtc(DATA.finding.completed_at_utc));
      setText('metric-timepoints', fmtNumber(DATA.summary.timepoints));
      setText('metric-raw', fmtNumber(DATA.summary.raw_reads)); setText('metric-filtered', fmtNumber(DATA.summary.filtered_reads));
      setText('metric-pathogens', fmtNumber(DATA.summary.pathogen_species_latest));
      setText('metric-raw-note', `from ${DATA.summary.qc_timepoints} QC batch${DATA.summary.qc_timepoints === 1 ? '' : 'es'}`);
      setText('pathogen-count', `${DATA.summary.pathogen_species_latest} in latest result`);
      setText('database-name', DATA.configuration.database || '—');
      document.getElementById('header-status-dot').className = `status-dot ${DATA.current.tone}`;
      const banner = document.getElementById('clinical-banner'); banner.classList.add(DATA.finding.tone);
      setText('finding-icon', DATA.finding.tone === 'detected' ? '!' : DATA.finding.tone === 'clear' ? '✓' : '…');
    }

    function svgNode(name,attributes={},text='') { const node=document.createElementNS('http://www.w3.org/2000/svg',name); Object.entries(attributes).forEach(([key,value])=>node.setAttribute(key,String(value))); if(text) node.textContent=text; return node; }
    function renderLineChart(container,history,metric,color,formatter) {
      const width=420,height=205,margin={top:18,right:15,bottom:38,left:55}; const plotWidth=width-margin.left-margin.right,plotHeight=height-margin.top-margin.bottom;
      const values=history.map(point=>Number(point[metric]||0)); let min=metric==='sni_score'?Math.max(0,Math.min(...values)-.005):0; let max=Math.max(...values); if(max<=min) max=min+(metric==='read_count'?1:.01);
      const x=index=>history.length===1?margin.left+plotWidth/2:margin.left+index*plotWidth/(history.length-1); const y=value=>margin.top+plotHeight-((value-min)/(max-min))*plotHeight;
      const svg=svgNode('svg',{viewBox:`0 0 ${width} ${height}`,class:'line-chart',role:'img','aria-label':`${metric} by analysis batch`});
      for(let tick=0;tick<=4;tick+=1){const value=min+(max-min)*tick/4,yPos=y(value);svg.appendChild(svgNode('line',{x1:margin.left,x2:width-margin.right,y1:yPos,y2:yPos,class:'grid-line'}));const label=metric==='read_count'?fmt.format(Math.round(value)):metric==='signature_coverage'?`${(value*100).toFixed(1)}%`:value.toFixed(3);svg.appendChild(svgNode('text',{x:margin.left-8,y:yPos+3,'text-anchor':'end',class:'axis-label'},label));}
      svg.appendChild(svgNode('line',{x1:margin.left,x2:width-margin.right,y1:margin.top+plotHeight,y2:margin.top+plotHeight,class:'axis-line'}));
      svg.appendChild(svgNode('polyline',{points:history.map((point,index)=>`${x(index)},${y(values[index])}`).join(' '),fill:'none',stroke:color,'stroke-width':3,'stroke-linecap':'round','stroke-linejoin':'round'}));
      const labelStep=Math.max(1,Math.ceil(history.length/8)); history.forEach((point,index)=>{const circle=svgNode('circle',{cx:x(index),cy:y(values[index]),r:4.5,fill:'#fff',stroke:color,'stroke-width':3});circle.appendChild(svgNode('title',{},`Batch ${point.timepoint}: ${formatter(values[index])}`));svg.appendChild(circle);if(index%labelStep===0||index===history.length-1)svg.appendChild(svgNode('text',{x:x(index),y:height-14,'text-anchor':'middle',class:'axis-label'},`B${point.timepoint}`));}); container.replaceChildren(svg);
    }
    function chartBox(title,latestValue){const box=document.createElement('div');box.className='chart-box';box.innerHTML=`<div class="chart-title"><strong>${escapeHtml(title)}</strong><span>${escapeHtml(latestValue)}</span></div><div class="chart-target"></div>`;return box;}

    function renderPathogens() {
      const grid=document.getElementById('pathogen-grid'); if(!DATA.pathogens.length){grid.innerHTML='<div class="empty-state">No annotated human-pathogenic species have been identified in completed cumulative profiles.</div>';return;}
      const storedOpen=localStorage.getItem('spades-stream-open-taxid');
      DATA.pathogens.forEach(pathogen=>{
        const latest=pathogen.latest,card=document.createElement('article'); card.className=`pathogen-card ${pathogen.present_latest?'':'historical'}`;
        const button=document.createElement('button');button.type='button';button.className='pathogen-toggle';button.setAttribute('aria-expanded','false');
        const status=pathogen.present_latest?'Sequence evidence in latest result':'Historical finding — not in latest result';
        const trend=pathogen.history.length>1?`${signed(pathogen.change.read_count)} supporting alignments since prior batch`:'First observed in this run';
        button.innerHTML=`<div class="pathogen-top"><div><span class="pathogen-status">${escapeHtml(status)}</span><div class="pathogen-name">${escapeHtml(pathogen.name)}</div><div class="taxid">Species · TaxID ${escapeHtml(pathogen.taxid||'not available')}</div></div><span class="expand-mark" aria-hidden="true">+</span></div><div class="pathogen-metrics"><div class="mini-metric"><div class="mini-label">Supporting read alignments</div><div class="mini-value">${fmtNumber(latest.read_count)}</div><div class="mini-context">reads aligned to signature</div></div><div class="mini-metric"><div class="mini-label">Reference signature covered</div><div class="mini-value">${fmtCoverage(latest.signature_coverage)}</div><div class="mini-context">breadth of signature evidence</div></div><div class="mini-metric"><div class="mini-label">SNI score</div><div class="mini-value">${fmtSni(latest.sni_score)}</div><div class="mini-context">${fmtSniPercent(latest.sni_score)}</div></div></div><div class="pathogen-foot"><span class="trend-chip">${escapeHtml(trend)}</span><span>Observed B${pathogen.first_seen}–B${pathogen.last_seen}</span></div>`;
        const detail=document.createElement('div');detail.className='pathogen-detail';detail.hidden=true;detail.innerHTML='<p class="evidence-note"><strong>Evidence trend.</strong> Supporting alignments show the number of reads mapped to the organism signature; signature coverage shows how much of that reference signature is represented; SNI summarizes nucleotide identity. These are analytical evidence measures, not a diagnosis or measure of organism viability.</p>';
        const charts=document.createElement('div');charts.className='charts-grid';const readChart=chartBox('Supporting read alignments',fmtNumber(latest.read_count));const coverageChart=chartBox('Reference signature coverage',fmtCoverage(latest.signature_coverage));const sniChart=chartBox('SNI score',fmtSni(latest.sni_score));charts.append(readChart,coverageChart,sniChart);detail.appendChild(charts);card.append(button,detail);grid.appendChild(card);
        let rendered=false;const setOpen=open=>{detail.hidden=!open;card.classList.toggle('expanded',open);button.setAttribute('aria-expanded',String(open));if(open&&!rendered){renderLineChart(readChart.querySelector('.chart-target'),pathogen.history,'read_count','#156f73',fmtNumber);renderLineChart(coverageChart.querySelector('.chart-target'),pathogen.history,'signature_coverage','#356f9e',fmtCoverage);renderLineChart(sniChart.querySelector('.chart-target'),pathogen.history,'sni_score','#a74729',fmtSni);rendered=true;}if(open)localStorage.setItem('spades-stream-open-taxid',pathogen.taxid||pathogen.name);else if(localStorage.getItem('spades-stream-open-taxid')===(pathogen.taxid||pathogen.name))localStorage.removeItem('spades-stream-open-taxid');};
        button.addEventListener('click',()=>setOpen(detail.hidden));if(storedOpen===(pathogen.taxid||pathogen.name))setOpen(true);
      });
    }
    function renderTimepoints(){const body=document.getElementById('timepoint-body');if(!DATA.timepoints.length){body.innerHTML='<tr><td colspan="6" class="empty-state">No completed analysis batches yet.</td></tr>';return;}[...DATA.timepoints].reverse().forEach(item=>{const row=document.createElement('tr');row.innerHTML=`<td><strong>B${String(item.timepoint).padStart(6,'0')}</strong></td><td class="file-cell" title="${escapeHtml(item.input_path)}">${escapeHtml(item.input_file)}</td><td>${fmtNumber(item.raw_reads)}</td><td>${fmtNumber(item.filtered_reads)}</td><td><span class="status-tag ${escapeHtml(item.status)}">${escapeHtml(item.status.replace('_',' '))}</span></td><td>${escapeHtml(shortUtc(item.completed_at_utc))}</td>`;body.appendChild(row);});}
    function setupRefresh(){const checkbox=document.getElementById('auto-refresh'),countdown=document.getElementById('refresh-countdown');checkbox.checked=localStorage.getItem('spades-stream-auto-refresh')!=='off';let seconds=15;checkbox.addEventListener('change',()=>{localStorage.setItem('spades-stream-auto-refresh',checkbox.checked?'on':'off');seconds=15;});setInterval(()=>{if(!checkbox.checked){countdown.textContent='paused';return;}seconds-=1;countdown.textContent=`${seconds}s`;if(seconds<=0)location.reload();},1000);}
    renderSummary();renderPathogens();renderTimepoints();setupRefresh();
  </script>
</body>
</html>
"""


def generate_stream_report(
    state: Dict[str, Any], output: Path, title: Optional[str] = None
) -> Path:
    prefix = str(state.get("configuration", {}).get("prefix", "SPADES stream"))
    report_title = title or f"{prefix} · Pathogen screening summary"
    payload = build_payload(state)
    payload_json = json.dumps(payload, ensure_ascii=False, separators=(",", ":")).replace(
        "</", "<\\/"
    )
    document = (
        HTML_TEMPLATE.replace("REPORT_TITLE", html.escape(report_title))
        .replace("REPORT_HEADING", html.escape(report_title))
        .replace("REPORT_DATA", payload_json)
    )
    output = output.resolve()
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
    print(f"INFO: Stream report: {generated}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
