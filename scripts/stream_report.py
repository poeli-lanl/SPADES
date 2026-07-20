#!/usr/bin/env python3
"""Generate the self-contained real-time dashboard for stream_spades.py."""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
import html
import json
import math
import os
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence


HTML_ASSET_DIR = Path(__file__).resolve().parent.parent / "data" / "html"
DEFAULT_LEVELS = ("species")
DEFAULT_MIN_READ_COUNT = 0
SNI_THRESHOLDS = {"species": 0.95, "strain": 0.99}
DEFAULT_SNI_THRESHOLD = 0.90
RESOURCE_SPECS = (
    (
        "<script src=\"/publicdata/js/vue.global.prod.js\"></script>",
        "js/vue.global.prod.js",
        "script",
    ),
    (
        "<script src=\"/publicdata/js/primevue.min.js\"></script>",
        "js/primevue.min.js",
        "script",
    ),
    ("<script src=\"/publicdata/js/aura.js\"></script>", "js/aura.js", "script"),
    (
        "<link rel=\"stylesheet\" href=\"/publicdata/css/primeicons.css\">",
        "css/primeicons.css",
        "primeicons",
    ),
)


def _finite_number(value: Any) -> Optional[float]:
    """Parse a finite number while treating common missing tokens as missing."""
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.lower() in {"na", "n/a", "nan", "none", "null", "."}:
        return None
    try:
        parsed = float(text.rstrip("%"))
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def _integer(value: Any) -> Optional[int]:
    parsed = _finite_number(value)
    return int(parsed) if parsed is not None else None


def _normalized_fraction(value: Any) -> Optional[float]:
    """Normalize fraction- and percentage-scaled metrics to the 0..1 range."""
    parsed = _finite_number(value)
    if parsed is None:
        return None
    text = str(value).strip()
    if text.endswith("%") or 1 < abs(parsed) <= 100:
        parsed /= 100.0
    return parsed if 0 <= parsed <= 1 else None


def _stable_taxon_key(level: Any, taxid: Any, name: Any) -> str:
    normalized_level = str(level or "unknown").strip().lower() or "unknown"
    identifier = str(taxid or "").strip() or str(name or "Unknown taxon").strip()
    return f"{normalized_level}:{identifier}"


def _deduplicated_history(points: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Keep one normalized observation per timepoint, preferring the latest record."""
    by_timepoint = {int(point["timepoint"]): point for point in points}
    return [by_timepoint[key] for key in sorted(by_timepoint)]


def _latest_change(history: Sequence[Dict[str, Any]], field: str) -> Optional[float]:
    """Return the unrounded latest-minus-previous delta for one normalized metric."""
    if len(history) < 2:
        return None
    previous = history[-2].get(field)
    current = history[-1].get(field)
    if previous is None or current is None:
        return None
    return current - previous


def _sni_threshold(level: Any) -> float:
    return SNI_THRESHOLDS.get(str(level or "").strip().lower(), DEFAULT_SNI_THRESHOLD)


def _has_pathogenic_info(value: Any) -> bool:
    if value is None:
        return False
    text = str(value).strip()
    return bool(text and text.lower() != "nan")


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
    return value in {"yes", "true", "1", "y"} or _has_pathogenic_info(
        row.get("PATHOGENIC_INFO")
    )


def _taxon_rows(path: Path) -> List[Dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [row for row in reader if str(row.get("NAME", "")).strip()]


def _database_label(value: Any) -> str:
    path = Path(str(value or ""))
    if not path.name:
        return ""
    return f"{path.parent.name}/{path.name}" if path.parent.name else path.name


def build_payload(state: Dict[str, Any], base_dir: Optional[Path] = None) -> Dict[str, Any]:
    """Build browser data from cumulative stream state.

    ``base_dir`` is both the origin for relative state paths and the directory relative
    to which paths are exposed in the report.
    """
    base_dir = (base_dir or Path.cwd()).resolve()
    records_by_timepoint = {
        int(record["timepoint"]): record for record in state.get("timepoints", [])
    }
    records = [records_by_timepoint[key] for key in sorted(records_by_timepoint)]
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
                "raw_reads": _integer(raw_reads),
                "filtered_reads": _integer(filtered_reads),
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
            key = _stable_taxon_key(level, taxid, name)
            human_pathogen = _is_human_pathogen(row)
            pathogen = pathogens.setdefault(
                key,
                {
                    "key": key,
                    "taxid": taxid,
                    "name": name,
                    "level": level,
                    "human_pathogen": human_pathogen,
                    "history": [],
                },
            )
            pathogen["history"].append(
                {
                    "timepoint": timepoint,
                    "completed_at": completed_at,
                    "read_count": _integer(row.get("READ_COUNT")),
                    "best_sig_cov": _normalized_fraction(
                        row.get("BEST_SIG_COV", row.get("SIG_COV"))
                    ),
                    "sni_score": _normalized_fraction(row.get("SNI_SCORE")),
                    "human_pathogen": human_pathogen,
                }
            )

    latest_timepoint = int(records[-1]["timepoint"]) if records else 0
    pathogen_list: List[Dict[str, Any]] = []
    for pathogen in pathogens.values():
        history = _deduplicated_history(pathogen["history"])
        pathogen["history"] = history
        pathogen["latest"] = history[-1]
        pathogen["human_pathogen"] = history[-1]["human_pathogen"]
        pathogen["passes_sni"] = bool(
            history[-1]["sni_score"] is not None
            and history[-1]["sni_score"] >= _sni_threshold(pathogen["level"])
        )
        pathogen["first_seen"] = history[0]["timepoint"]
        pathogen["last_seen"] = history[-1]["timepoint"]
        pathogen["present_latest"] = history[-1]["timepoint"] == latest_timepoint
        pathogen["change"] = {
            "read_count": _latest_change(history, "read_count"),
            "best_sig_cov": _latest_change(history, "best_sig_cov"),
            "sni_score": _latest_change(history, "sni_score"),
        }
        pathogen_list.append(pathogen)
    pathogen_list.sort(
        key=lambda item: (
            not item["present_latest"],
            item["latest"]["read_count"] is None,
            -(item["latest"]["read_count"] or 0),
            item["name"].lower(),
        )
    )
    latest_taxa = [item for item in pathogen_list if item["present_latest"]]
    human_pathogens = [item for item in pathogen_list if item["human_pathogen"]]
    latest_human_pathogens = [item for item in latest_taxa if item["human_pathogen"]]

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
    database = configuration.pop("database", "")
    configuration["input_dir"] = _relative_path(
        configuration.get("input_dir", ""), base_dir
    )
    configuration["database_display"] = _database_label(database)
    return {
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "revision": f"{len(records)}:{latest_timepoint}:{monitor_status}:{pending.get('stage', '')}",
        "configuration": configuration,
        "monitor_status": monitor_status,
        "processing_batch": int(pending.get("timepoint", 0)) if pending else None,
        "current": current,
        "filter_defaults": {
            "levels": list(DEFAULT_LEVELS),
            "min_read_count": DEFAULT_MIN_READ_COUNT,
            "pass_sni": False,
            "human_pathogens_only": True,
        },
        "summary": {
            "timepoints": len(timepoints),
            "raw_reads": sum(raw_values),
            "filtered_reads": sum(filtered_values),
            "qc_timepoints": len(raw_values),
            "taxa": len(pathogen_list),
            "taxa_latest": len(latest_taxa),
            "human_pathogens": len(human_pathogens),
            "human_pathogens_latest": len(latest_human_pathogens),
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
    :root { --page:#f5f7fa; --surface:#fff; --ink:#17212b; --muted:#627181; --line:#dce4ea; --teal:#197278; --teal-dark:#115e63; --teal-soft:#e5f3f1; --blue:#3f7cac; --blue-soft:#eaf2f9; --coral:#d95d39; --coral-soft:#fff0eb; --gold:#b98113; --shadow:0 16px 42px rgba(24,33,47,.08); }
    * { box-sizing:border-box; }
    body { margin:0; min-height:100vh; background:radial-gradient(circle at 5% 0,rgba(25,114,120,.08),transparent 28rem),var(--page); color:var(--ink); font-family:Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; }
    [v-cloak] { display:none; } button,input { font:inherit; } .app-shell { width:min(1600px,calc(100% - 32px)); margin:auto; padding:26px 0 42px; }
    .report-header { display:flex; justify-content:space-between; align-items:flex-start; gap:24px; margin-bottom:16px; }
    .eyebrow { color:var(--teal); font-size:.73rem; font-weight:700; text-transform:uppercase; }
    h1 { margin:5px 0 0; font-size:30px; line-height:1.5; }
    .header-meta,.live-control,.toolbar { display:flex; align-items:center; flex-wrap:wrap; gap:8px; }
    .header-meta { margin-top:13px; }
    .meta-pill,.count-pill,.filter-status { display:inline-flex; align-items:center; gap:7px; min-height:30px; border:1px solid var(--line); border-radius:999px; background:rgba(255,255,255,.84); color:var(--muted); font-size:.8rem; font-weight:700; padding:0 11px; }
    .meta-pill strong { color:var(--ink); }
    .status-dot { width:8px; height:8px; border-radius:50%; background:#94a3af; }
    .status-dot.active { background:var(--teal); box-shadow:0 0 0 5px rgba(25,114,120,.12); animation:pulse 1.7s infinite; }
    .status-dot.warning { background:var(--gold); }
    @keyframes pulse { 50% { box-shadow:0 0 0 8px rgba(25,114,120,0); } }
    .live-control { justify-content:flex-end; color:var(--muted); font-size:.8rem; font-weight:700; }
    .live-control label { display:flex; align-items:center; gap:7px; cursor:pointer; }
    .live-state { color:var(--teal); min-width:70px; text-align:right; }
    .stream-button { display:inline-flex; align-items:center; justify-content:center; gap:7px; min-height:38px; border:1px solid #b9c9cf; border-radius:8px; background:#fff; color:var(--teal-dark); cursor:pointer; font-size:.8rem; font-weight:760; line-height:1; padding:0 13px; transition:background .16s,border-color .16s,box-shadow .16s,color .16s,transform .16s; }
    .stream-button:hover:not(:disabled) { border-color:var(--teal); background:var(--teal-soft); color:var(--teal-dark); transform:translateY(-1px); }
    .stream-button:focus-visible { outline:0; box-shadow:0 0 0 3px rgba(25,114,120,.18); }
    .stream-button:disabled { cursor:not-allowed; opacity:.52; transform:none; }
    .stream-button.primary { border-color:var(--teal); background:var(--teal); color:#fff; }
    .stream-button.primary:hover:not(:disabled) { border-color:var(--teal-dark); background:var(--teal-dark); color:#fff; }
    .stream-button i { font-size:.78rem; }
    .stats-grid { display:grid; grid-template-columns:repeat(4,minmax(0,1fr)); gap:12px; margin-bottom:16px; }
    .metric-card,.panel { border:1px solid var(--line); border-radius:10px; background:rgba(255,255,255,.96); box-shadow:var(--shadow); }
    .metric-card { padding:16px; } .metric-label { color:var(--muted); font-size:.72rem; font-weight:800; letter-spacing:.055em; text-transform:uppercase; }
    .metric-value { margin-top:8px; font-size:1.65rem; font-weight:820; line-height:1; } .metric-note { margin-top:7px; color:var(--muted); font-size:.77rem; }
    .filters-panel { margin-bottom:16px; padding:16px; }
    .filters-heading { display:flex; align-items:flex-start; justify-content:space-between; gap:18px; margin-bottom:13px; }
    .filters-heading h2 { margin:0; font-size:1rem; } .filters-heading p { margin:4px 0 0; color:var(--muted); font-size:.78rem; }
    .filter-status.active { border-color:#badbd6; background:var(--teal-soft); color:var(--teal-dark); }
    .filters-grid { display:grid; grid-template-columns:minmax(230px,1.4fr) minmax(145px,.65fr) minmax(190px,.8fr) minmax(190px,.8fr) minmax(235px,1.2fr) auto; align-items:end; gap:12px; }
    .filter-field { min-width:0; } .filter-label { display:block; margin-bottom:6px; color:#536372; font-size:.72rem; font-weight:780; }
    .metric-input,.search-input { width:100%; height:40px; border:1px solid #cbd7dd; border-radius:8px; background:#fff; color:var(--ink); outline:0; }
    .metric-input { padding:0 11px; } .metric-input:focus,.search-input:focus { border-color:var(--teal); box-shadow:0 0 0 3px rgba(25,114,120,.12); }
    .search-box { position:relative; } .search-box i { position:absolute; left:12px; top:50%; transform:translateY(-50%); color:#82919f; }
    .search-input { padding:0 12px 0 35px; }
    .switch-control { display:flex; align-items:center; min-height:40px; gap:9px; color:#465664; font-size:.78rem; font-weight:700; cursor:pointer; }
    .switch-control input { width:34px; height:18px; accent-color:var(--teal); cursor:pointer; }
    .filter-actions { display:flex; align-items:center; justify-content:flex-end; min-height:40px; }
    .p-multiselect { width:100%; min-height:40px; border-color:#cbd7dd!important; border-radius:8px!important; }
    .panel { overflow:hidden; } .panel-heading { display:flex; justify-content:space-between; gap:20px; align-items:center; padding:16px 18px; border-bottom:1px solid var(--line); }
    .panel-title { margin:0; font-size:1.12rem; } .panel-subtitle { margin:4px 0 0; color:var(--muted); font-size:.82rem; }
    .p-datatable { border:0; } .p-datatable-thead>tr>th { background:#f7f9fb!important; color:#596979!important; font-size:.7rem; letter-spacing:.045em; text-transform:uppercase; }
    .p-datatable-tbody>tr>td { border-color:var(--line)!important; padding:.78rem .9rem!important; } .p-datatable-tbody>tr:hover { background:#f7fbfa!important; }
    .p-datatable-table-container { overflow-x:auto; }
    .p-datatable-row-expansion>td { height:auto!important; overflow:visible!important; background:#f8fafb!important; padding:0!important; vertical-align:top!important; } .name-primary { display:flex; align-items:center; gap:7px; font-weight:780; }
    .name-secondary { display:block; margin-top:3px; color:var(--muted); font-size:.73rem; } .pathogen-dot { width:7px; height:7px; border-radius:50%; background:var(--coral); }
    .value-cell { font-variant-numeric:tabular-nums; font-weight:720; white-space:nowrap; }
    .change-pill { display:inline-flex; align-items:center; gap:4px; margin-left:7px; border-radius:999px; padding:3px 7px; font-size:.66rem; font-weight:850; vertical-align:middle; }
    .change-pill.read { background:var(--teal-soft); color:#0f5b60; } .change-pill.coverage { background:var(--blue-soft); color:#315f83; } .change-pill.sni { background:var(--coral-soft); color:#a8462a; }
    .change-pill i { font-size:.61rem; }
    .history { width:100%; min-width:0; height:auto; overflow:visible; padding:20px; } .history-head { display:flex; justify-content:space-between; gap:18px; align-items:flex-start; margin-bottom:13px; }
    .history-head strong { font-size:.88rem; } .history-head span { color:var(--muted); font-size:.76rem; }
    .chart-grid { display:grid; grid-template-columns:repeat(3,minmax(280px,1fr)); gap:12px; width:100%; overflow-x:auto; padding-bottom:2px; } .chart-card { min-width:280px; border:1px solid var(--line); border-radius:9px; background:#fff; padding:12px; }
    .chart-head { display:flex; justify-content:space-between; gap:10px; align-items:baseline; } .chart-head strong { font-size:.8rem; } .chart-head span { color:var(--muted); font-size:.71rem; }
    .chart-svg { display:block; width:100%; height:220px; overflow:visible; } .grid-line { stroke:#e7edf1; stroke-width:1; } .chart-line { fill:none; stroke-width:3; stroke-linecap:round; stroke-linejoin:round; }
    .chart-point { fill:#fff; stroke-width:3; cursor:crosshair; } .chart-axis { fill:var(--muted); font-size:10px; } .chart-tooltip { min-height:20px; color:var(--muted); font-size:.7rem; text-align:center; }
    .empty { padding:32px; color:var(--muted); text-align:center; } .batch-panel { margin-top:16px; } .batch-table-wrap { overflow:auto; max-height:320px; }
    .batch-table { width:100%; border-collapse:collapse; font-size:.8rem; } .batch-table th { position:sticky; top:0; background:#f7f9fb; color:var(--muted); font-size:.68rem; letter-spacing:.04em; text-align:left; text-transform:uppercase; padding:10px 13px; }
    .batch-table td { border-top:1px solid var(--line); padding:10px 13px; } .path-cell { max-width:330px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }
    .status-tag { display:inline-flex; border-radius:999px; background:var(--teal-soft); color:#0f5b60; font-size:.69rem; font-weight:800; padding:4px 8px; text-transform:capitalize; }
    .status-tag.no_alignments { background:#fff6df; color:#825d0d; } footer { display:flex; justify-content:space-between; flex-wrap:wrap; gap:12px; margin-top:15px; color:var(--muted); font-size:.72rem; } footer span { overflow-wrap:anywhere; }
    @media(max-width:1200px){.filters-grid{grid-template-columns:repeat(3,minmax(0,1fr))}.filter-actions{justify-content:flex-start}}
    @media(max-width:900px){.stats-grid{grid-template-columns:repeat(2,minmax(0,1fr))}.chart-grid{grid-template-columns:1fr;overflow-x:visible}.chart-card{min-width:0}.panel-heading,.report-header{align-items:flex-start;flex-direction:column}.live-control{justify-content:flex-start}.filters-grid{grid-template-columns:repeat(2,minmax(0,1fr))}}
    @media(max-width:560px){.app-shell{width:calc(100% - 18px);padding-top:16px}.stats-grid{grid-template-columns:1fr}.filters-grid{grid-template-columns:1fr}.filters-heading{flex-direction:column}.filter-actions,.stream-button{width:100%}.live-control{width:100%}.live-state{text-align:left}.history{padding:14px}.history-head{flex-direction:column}}
  </style>
</head>
<body>
<main id="app" class="app-shell" v-cloak>
  <header class="report-header">
    <div>
      <div class="eyebrow">SPADES-GOTTCHA2</div><h1>REPORT_HEADING</h1>
      <div class="header-meta">
        <span class="meta-pill" :title="payload.current.detail"><span class="status-dot" :class="payload.current.tone"></span>Processing batch: <strong>{{ payload.processing_batch ? `B${payload.processing_batch}` : 'stopped' }}</strong></span>
        <span class="meta-pill"><i class="pi pi-check-circle"></i>Batches completed: <strong>{{ formatNumber(payload.summary.timepoints) }}</strong></span>
        <span class="meta-pill"><i class="pi pi-clock"></i>{{ formatTime(payload.generated_at) }}</span>
        <span class="meta-pill"><i class="pi pi-database"></i>{{ payload.configuration.database_display || '—' }}</span>
      </div>
    </div>
    <div class="live-control">
      <label><input type="checkbox" v-model="liveEnabled"> Live updates</label><span class="live-state">{{ liveMessage }}</span>
      <button class="stream-button primary" type="button" :disabled="polling" @click="pollNow"><i class="pi pi-refresh" :class="{'pi-spin':polling}"></i><span>Update</span></button>
    </div>
  </header>
  <section class="stats-grid">
    <article class="metric-card"><div class="metric-label">Human pathogens</div><div class="metric-value">{{ formatNumber(payload.summary.human_pathogens_latest) }}</div><div class="metric-note">taxa detected and requiring review</div></article>
    <article class="metric-card"><div class="metric-label">Taxa in latest result</div><div class="metric-value">{{ formatNumber(totalLatestTaxa.length) }}</div><div class="metric-note">across reported taxonomic levels</div></article>
    <article class="metric-card"><div class="metric-label">Input reads</div><div class="metric-value">{{ formatNumber(payload.summary.raw_reads) }}</div><div class="metric-note">from {{ payload.summary.qc_timepoints }} QC batches</div></article>
    <article class="metric-card"><div class="metric-label">Post-QC records</div><div class="metric-value">{{ formatNumber(payload.summary.filtered_reads) }}</div><div class="metric-note">available to cumulative profiling</div></article>
  </section>
  <section class="panel filters-panel" aria-label="Result filters">
    <div class="filters-heading"><div><h2>Result filters</h2><p>Filters update the table in place and remain active as new batches arrive.</p></div><span class="filter-status" :class="{active:activeFilterCount}"><i class="pi pi-filter"></i>{{ activeFilterCount ? `${activeFilterCount} active` : 'All taxa' }}</span></div>
    <div class="filters-grid">
      <div class="filter-field"><label class="filter-label" for="level-filter">Taxonomic level</label><p-multi-select v-model="selectedLevels" :options="levelOptions" option-label="label" option-value="value" input-id="level-filter" display="chip" placeholder="Select levels"></p-multi-select></div>
      <div class="filter-field"><label class="filter-label" for="read-count-filter">Minimum READ_COUNT</label><input id="read-count-filter" class="metric-input" type="number" min="0" step="1" v-model.number="minReadCount"></div>
      <label class="switch-control" title="Others ≥ 0.90, species ≥ 0.95, strain ≥ 0.99"><input type="checkbox" v-model="passSni"><span>Pass SNI threshold</span></label>
      <label class="switch-control"><input type="checkbox" v-model="humanPathogensOnly"><span>Human pathogens only</span></label>
      <div class="filter-field"><label class="filter-label" for="taxon-search">Search</label><div class="search-box"><i class="pi pi-search"></i><input id="taxon-search" class="search-input" v-model="filters.global.value" placeholder="Taxon name or TaxID"></div></div>
      <div class="filter-actions"><button class="stream-button" type="button" @click="resetFilters"><i class="pi pi-filter-slash"></i><span>Clear Filters</span></button></div>
    </div>
  </section>
  <section class="panel">
    <div class="panel-heading"><div><h2 class="panel-title">Latest profiling results</h2><p class="panel-subtitle">Expand a taxon to inspect normalized metric history. Change pills compare the latest observation with the immediately preceding observation.</p></div><span class="count-pill">{{ formatNumber(filteredTaxa.length) }} of {{ formatNumber(totalLatestTaxa.length) }} taxa</span></div>
    <p-data-table v-model:expanded-rows="expandedRows" v-model:selection="selectedRows" :value="filteredTaxa" data-key="key" state-storage="local" state-key="spades-stream-results-v2" selection-mode="multiple" :meta-key-selection="false" paginator :rows="25" :rows-per-page-options="[10,25,50,100]" striped-rows removable-sort @row-expand="onRowExpand" @row-collapse="onRowCollapse">
      <template #empty><div class="empty">No taxa match the active filters.</div></template>
      <p-column expander style="width:3.5rem"></p-column>
      <p-column field="name" header="Taxon" sortable style="min-width:250px"><template #body="slot"><span class="name-primary"><span v-if="slot.data.human_pathogen" class="pathogen-dot" title="Human pathogen"></span>{{ slot.data.name }}</span><span class="name-secondary">{{ formatLevel(slot.data.level) }} · TaxID {{ slot.data.taxid || 'not available' }}</span></template></p-column>
      <p-column field="latest.read_count" header="READ_COUNT" sortable style="min-width:180px"><template #body="slot"><span class="value-cell">{{ formatNullable(slot.data.latest.read_count) }}</span><change-pill v-if="slot.data.change.read_count > 0" tone="read" metric="READ_COUNT" format="integer" :value="slot.data.change.read_count"></change-pill></template></p-column>
      <p-column field="latest.best_sig_cov" header="BEST_SIG_COV" sortable style="min-width:205px"><template #body="slot"><span class="value-cell">{{ formatCoverage(slot.data.latest.best_sig_cov) }}</span><change-pill v-if="slot.data.change.best_sig_cov > 0" tone="coverage" metric="BEST_SIG_COV" format="percent" :value="slot.data.change.best_sig_cov"></change-pill></template></p-column>
      <p-column field="latest.sni_score" header="SNI_SCORE" sortable style="min-width:185px"><template #body="slot"><span class="value-cell">{{ formatSni(slot.data.latest.sni_score) }}</span><change-pill v-if="slot.data.change.sni_score > 0" tone="sni" metric="SNI_SCORE" format="decimal" :value="slot.data.change.sni_score"></change-pill></template></p-column>
      <p-column field="first_seen" header="First Observed" sortable style="min-width:130px"><template #body="slot">B{{ slot.data.first_seen }}</template></p-column>
      <template #expansion="slot"><div class="history"><div class="history-head"><div><strong>{{ slot.data.name }} metric history</strong><br><span>Missing values remain gaps; charts update in place with each cumulative result.</span></div><span>{{ slot.data.history.length }} observations</span></div><div class="chart-grid"><history-chart title="READ_COUNT" field="read_count" color="#197278" :history="slot.data.history" format="integer"></history-chart><history-chart title="BEST_SIG_COV" field="best_sig_cov" color="#3f7cac" :history="slot.data.history" format="percent"></history-chart><history-chart title="SNI_SCORE" field="sni_score" color="#d95d39" :history="slot.data.history" format="decimal"></history-chart></div></div></template>
    </p-data-table>
  </section>
  <section class="panel batch-panel"><div class="panel-heading"><div><h2 class="panel-title">Analysis batches</h2><p class="panel-subtitle">Paths are relative to this report wherever possible. Times use the local system timezone.</p></div></div><div class="batch-table-wrap"><table class="batch-table"><thead><tr><th>Batch</th><th>Input</th><th>Input reads</th><th>Post-QC</th><th>Status</th><th>Completed</th><th>Pipeline log</th></tr></thead><tbody><tr v-for="item in reversedTimepoints" :key="item.timepoint"><td><strong>B{{ String(item.timepoint).padStart(6,'0') }}</strong></td><td class="path-cell" :title="item.input_path">{{ item.input_file }}</td><td>{{ formatNullable(item.raw_reads) }}</td><td>{{ formatNullable(item.filtered_reads) }}</td><td><span class="status-tag" :class="item.status">{{ item.status.replace('_',' ') }}</span></td><td>{{ formatTime(item.completed_at) }}</td><td class="path-cell" :title="item.run_log">{{ item.run_log || '—' }}</td></tr><tr v-if="!payload.timepoints.length"><td colspan="7" class="empty">No completed batches yet.</td></tr></tbody></table></div></section>
  <footer><span>Input: {{ payload.configuration.input_dir || '—' }}</span><span>SPADES-GOTTCHA2 REALTIME REPORT</span></footer>
</main>
<script id="report-data" type="application/json">__REPORT_DATA__</script>
<script>
window.browserResourcesReady.then(function(){
  const {createApp,ref,computed,watch,nextTick} = Vue;
  const initialPayload=JSON.parse(document.getElementById('report-data').textContent);
  const framePollToken=new URLSearchParams(location.hash.slice(1)).get('_spades_stream_poll');
  if(window.parent!==window&&framePollToken){window.parent.postMessage({type:'spades-stream-payload',token:framePollToken,payload:initialPayload},'*');return}
  const formatMetric=(value,format)=>{if(value===null||value===undefined||!Number.isFinite(Number(value)))return '—';const number=Number(value);if(format==='integer')return Math.round(number).toLocaleString();if(format==='percent')return `${(number*100).toFixed(2)}%`;return number.toFixed(5)};
  const ChangePill={props:{tone:String,metric:String,value:Number,format:String},setup(props){const text=computed(()=>formatMetric(props.value,props.format));const description=computed(()=>`${props.metric} increased by ${text.value} from the immediately preceding recorded timepoint`);return{text,description}},template:`<span class="change-pill" :class="tone" :title="description" :aria-label="description"><i class="pi pi-arrow-up" aria-hidden="true"></i>{{text}}</span>`};
  const HistoryChart={props:{title:String,field:String,color:String,history:Array,format:String},setup(props){const hover=ref(null),width=430,height=180,margin={top:14,right:15,bottom:31,left:54};const values=computed(()=>props.history.map(point=>{const value=point[props.field];return value===null||value===undefined||!Number.isFinite(Number(value))?null:Number(value)}));const finiteValues=computed(()=>values.value.filter(value=>value!==null));const domain=computed(()=>{const vals=finiteValues.value;if(!vals.length)return [0,1];let min=props.field==='sni_score'?Math.min(...vals):0;let max=Math.max(...vals);if(props.field==='sni_score')min=Math.max(0,min-Math.max((max-min)*.15,.0005));if(max<=min)max=min+(props.format==='integer'?1:.01);return [min,max]});const x=index=>props.history.length<=1?margin.left+(width-margin.left-margin.right)/2:margin.left+index*(width-margin.left-margin.right)/(props.history.length-1);const y=value=>margin.top+(height-margin.top-margin.bottom)*(1-(value-domain.value[0])/(domain.value[1]-domain.value[0]));const segments=computed(()=>{const result=[];let current=[];values.value.forEach((value,index)=>{if(value===null){if(current.length)result.push(current);current=[];return}current.push(`${x(index)},${y(value)}`)});if(current.length)result.push(current);return result});const ticks=computed(()=>Array.from({length:5},(_,index)=>{const value=domain.value[0]+(domain.value[1]-domain.value[0])*index/4;return{value,y:y(value)}}));const fmt=value=>formatMetric(value,props.format);const latestText=computed(()=>fmt(values.value[values.value.length-1]));const hoverText=computed(()=>hover.value===null?'Hover over a point for details':`Batch ${props.history[hover.value].timepoint} · ${fmt(values.value[hover.value])}`);return{width,height,margin,values,segments,ticks,x,y,fmt,hover,hoverText,latestText}},template:`<article class="chart-card"><div class="chart-head"><strong>{{title}}</strong><span>{{latestText}}</span></div><svg class="chart-svg" :viewBox="'0 0 '+width+' '+height" role="img" :aria-label="title+' history'"><g v-for="tick in ticks" :key="tick.y"><line class="grid-line" :x1="margin.left" :x2="width-margin.right" :y1="tick.y" :y2="tick.y"></line><text class="chart-axis" :x="margin.left-7" :y="tick.y+3" text-anchor="end">{{fmt(tick.value)}}</text></g><polyline v-for="(segment,index) in segments" :key="index" class="chart-line" :stroke="color" :points="segment.join(' ')"></polyline><g v-for="(point,index) in history" :key="point.timepoint"><circle v-if="values[index]!==null" class="chart-point" :stroke="color" :cx="x(index)" :cy="y(values[index])" r="4.5" @mouseenter="hover=index" @mouseleave="hover=null"><title>Batch {{point.timepoint}}: {{fmt(values[index])}}</title></circle><text v-if="history.length<=10||index===0||index===history.length-1" class="chart-axis" :x="x(index)" :y="height-10" text-anchor="middle">B{{point.timepoint}}</text></g></svg><div class="chart-tooltip">{{hoverText}}</div></article>`};
  const app=createApp({components:{ChangePill,HistoryChart},setup(){
    const payload=ref(initialPayload),filters=ref({global:{value:null}}),selectedRows=ref([]),polling=ref(false);
    const expansionStorageKey='spades-stream-expansion-v3';let expansionState={expanded:{},choices:{}};try{const saved=JSON.parse(localStorage.getItem(expansionStorageKey)||'{}');if(saved&&typeof saved==='object')expansionState={expanded:saved.expanded&&typeof saved.expanded==='object'?saved.expanded:{},choices:saved.choices&&typeof saved.choices==='object'?saved.choices:{}}}catch(error){}
    const expandedRows=ref({...expansionState.expanded}),liveEnabled=ref(localStorage.getItem('spades-stream-live')!=='off'),liveMessage=ref('current');
    const defaults=JSON.parse(JSON.stringify(initialPayload.filter_defaults));let stored={};try{stored=JSON.parse(localStorage.getItem('spades-stream-filters-v2')||'{}')}catch(error){}
    const selectedLevels=ref(Array.isArray(stored.levels)?stored.levels:[...defaults.levels]);const minReadCount=ref(Number.isFinite(Number(stored.min_read_count))?Number(stored.min_read_count):defaults.min_read_count);const passSni=ref(typeof stored.pass_sni==='boolean'?stored.pass_sni:defaults.pass_sni);const humanPathogensOnly=ref(typeof stored.human_pathogens_only==='boolean'?stored.human_pathogens_only:defaults.human_pathogens_only);if(typeof stored.search==='string')filters.value.global.value=stored.search;
    const rankOrder=['superkingdom','kingdom','phylum','class','order','family','genus','species','strain'];const totalLatestTaxa=computed(()=>payload.value.pathogens.filter(item=>item.present_latest));const levelOptions=computed(()=>[...new Set([...(payload.value.filter_defaults.levels||[]),...(selectedLevels.value||[]),...payload.value.pathogens.map(item=>item.level)])].sort((a,b)=>{const ai=rankOrder.indexOf(a),bi=rankOrder.indexOf(b);return(ai<0?99:ai)-(bi<0?99:bi)||a.localeCompare(b)}).map(value=>({value,label:formatLevel(value)})));
    const filteredTaxa=computed(()=>{const selected=new Set(selectedLevels.value||[]),minimum=Number.isFinite(Number(minReadCount.value))?Number(minReadCount.value):0,query=String(filters.value.global.value||'').trim().toLowerCase();return totalLatestTaxa.value.filter(item=>selected.has(item.level)&&(item.latest.read_count!==null&&item.latest.read_count>=minimum)&&(!passSni.value||item.passes_sni)&&(!humanPathogensOnly.value||item.human_pathogen)&&(!query||`${item.name} ${item.taxid} ${item.level}`.toLowerCase().includes(query)))});
    const activeFilterCount=computed(()=>{let count=0;const allLevels=new Set(levelOptions.value.map(item=>item.value));if(selectedLevels.value.length!==allLevels.size||selectedLevels.value.some(level=>!allLevels.has(level)))count+=1;if(Number(minReadCount.value||0)>0)count+=1;if(passSni.value)count+=1;if(humanPathogensOnly.value)count+=1;if(String(filters.value.global.value||'').trim())count+=1;return count});const reversedTimepoints=computed(()=>[...payload.value.timepoints].reverse());
    const captureScroll=()=>({windowX:window.scrollX,windowY:window.scrollY,items:[...document.querySelectorAll('.p-datatable-table-container,.batch-table-wrap')].map(node=>({node,left:node.scrollLeft,top:node.scrollTop}))});const restoreScroll=snapshot=>nextTick(()=>{snapshot.items.forEach(item=>{item.node.scrollLeft=item.left;item.node.scrollTop=item.top});window.scrollTo(snapshot.windowX,snapshot.windowY)});
    function persistExpansion(){localStorage.setItem(expansionStorageKey,JSON.stringify(expansionState))}
    function restoreCanonicalExpansion(){expandedRows.value={...expansionState.expanded}}
    function setExpansionChoice(event,isExpanded){const key=event&&event.data&&event.data.key;if(!key)return;expansionState.choices[key]=isExpanded;if(isExpanded)expansionState.expanded[key]=true;else delete expansionState.expanded[key];restoreCanonicalExpansion();persistExpansion()}
    function onRowExpand(event){setExpansionChoice(event,true)}function onRowCollapse(event){setExpansionChoice(event,false)}
    watch(liveEnabled,value=>localStorage.setItem('spades-stream-live',value?'on':'off'));
    watch([selectedLevels,minReadCount,passSni,humanPathogensOnly,()=>filters.value.global.value],()=>{const snapshot=captureScroll();localStorage.setItem('spades-stream-filters-v2',JSON.stringify({levels:selectedLevels.value,min_read_count:minReadCount.value,pass_sni:passSni.value,human_pathogens_only:humanPathogensOnly.value,search:filters.value.global.value||''}));nextTick(restoreCanonicalExpansion);restoreScroll(snapshot)},{deep:true});
    const formatNumber=value=>Number(value||0).toLocaleString();const formatNullable=value=>value===null||value===undefined?'—':Number(value).toLocaleString();const formatCoverage=value=>formatMetric(value,'percent');const formatSni=value=>formatMetric(value,'decimal');function formatLevel(value){return String(value||'unknown').replace(/_/g,' ').replace(/\b\w/g,char=>char.toUpperCase())}const formatTime=value=>{if(!value)return '—';const date=new Date(value);return Number.isNaN(date.valueOf())?String(value):date.toLocaleString(undefined,{dateStyle:'medium',timeStyle:'medium'})};
    function resetFilters(){selectedLevels.value=[...defaults.levels];minReadCount.value=defaults.min_read_count;passSni.value=defaults.pass_sni;humanPathogensOnly.value=defaults.human_pathogens_only;filters.value.global.value=null}
    async function applyPayload(next){if(next.generated_at===payload.value.generated_at)return;const selectedKeys=new Set(selectedRows.value.map(item=>item.key)),snapshot=captureScroll();payload.value=next;selectedRows.value=next.pathogens.filter(item=>item.present_latest&&selectedKeys.has(item.key));await nextTick();restoreCanonicalExpansion();restoreScroll(snapshot);liveMessage.value='updated'}
    function parseReport(text){const parsed=new DOMParser().parseFromString(text,'text/html'),node=parsed.getElementById('report-data');if(!node)throw new Error('Updated report data was not found');return JSON.parse(node.textContent)}
    async function iframePayload(sourceUrl){return new Promise((resolve,reject)=>{const frame=document.createElement('iframe'),token=`${Date.now()}-${Math.random().toString(36).slice(2)}`,url=new URL(sourceUrl);frame.hidden=true;url.hash=`_spades_stream_poll=${encodeURIComponent(token)}`;const cleanup=()=>{clearTimeout(timer);window.removeEventListener('message',onMessage);frame.remove()};const finish=(callback,value)=>{cleanup();callback(value)};const onMessage=event=>{const message=event.data;if(event.source!==frame.contentWindow||!message||message.type!=='spades-stream-payload'||message.token!==token)return;if(!message.payload||typeof message.payload!=='object'){finish(reject,new Error('Local report returned invalid data'));return}finish(resolve,message.payload)};const timer=setTimeout(()=>finish(reject,new Error('Local report update timed out')),12000);window.addEventListener('message',onMessage);frame.src=url.href;document.body.appendChild(frame)})}
    async function pollNow(){if(polling.value)return;polling.value=true;liveMessage.value='checking…';try{const url=new URL(location.href);url.searchParams.set('_stream_poll',Date.now());let next;try{const response=await fetch(url.href,{cache:'no-store'});if(!response.ok)throw new Error(`HTTP ${response.status}`);next=parseReport(await response.text())}catch(error){if(location.protocol!=='file:')throw error;next=await iframePayload(url.href)}await applyPayload(next);if(liveMessage.value!=='updated')liveMessage.value='current'}catch(error){console.warn('Live update failed:',error);liveMessage.value='retrying'}finally{polling.value=false}}
    setInterval(()=>{if(liveEnabled.value)pollNow()},5000);return{payload,filters,selectedRows,expandedRows,polling,liveEnabled,liveMessage,selectedLevels,minReadCount,passSni,humanPathogensOnly,totalLatestTaxa,filteredTaxa,levelOptions,activeFilterCount,reversedTimepoints,formatNumber,formatNullable,formatCoverage,formatSni,formatLevel,formatTime,resetFilters,onRowExpand,onRowCollapse,pollNow}
  }});
  app.use(PrimeVue.Config,{theme:{preset:PrimeUIX.Themes.Aura}});app.component('p-data-table',PrimeVue.DataTable);app.component('p-column',PrimeVue.Column);app.component('p-multi-select',PrimeVue.MultiSelect);app.mount('#app');
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
