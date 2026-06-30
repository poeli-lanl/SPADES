#!/usr/bin/env python3
import pandas as pd
import json
import argparse
import sys
import os
import minify_html

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>GOTTCHA2: TSVFILENAME</title>
  
  <script src="/publicdata/js/vue.global.prod.js"></script>
  <script src="/publicdata/js/primevue.min.js"></script>
  <script src="/publicdata/js/aura.js"></script>
  <script src="/publicdata/js/plotly-3.0.1.min.js" charset="utf-8"></script>

  <link href="/publicdata/css/bootstrap.min.css" rel="stylesheet">
  <link rel="stylesheet" href="/publicdata/css/primeicons.css">

  <style>
    :root {
        --page-bg: #f6f8fb;
        --surface: #ffffff;
        --surface-muted: #f0f5f3;
        --ink: #18212f;
        --muted: #617080;
        --line: #dce5ea;
        --accent: #197278;
        --accent-soft: #e6f3f1;
        --accent-strong: #115e63;
        --coral: #d95d39;
        --gold: #d89f2a;
        --blue: #3f7cac;
        --shadow: 0 18px 45px rgba(24, 33, 47, 0.08);
        --radius: 8px;
    }

    * {
        box-sizing: border-box;
    }

    body {
        min-height: 100vh;
        margin: 0;
        background: var(--page-bg);
        color: var(--ink);
        font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
        letter-spacing: 0;
    }

    button,
    input {
        font: inherit;
    }

    [v-cloak] {
        display: none;
    }

    .app-shell {
        width: min(1680px, calc(100% - 32px));
        margin: 0 auto;
        padding: 28px 0 40px;
    }

    .report-header {
        display: grid;
        grid-template-columns: minmax(0, 1fr) auto;
        gap: 24px;
        align-items: start;
        margin-bottom: 18px;
    }

    .eyebrow {
        color: var(--accent-strong);
        font-size: 0.76rem;
        font-weight: 700;
        margin-bottom: 7px;
        text-transform: uppercase;
    }

    .report-header h1 {
        max-width: 100%;
        color: var(--text);
        margin: 0;
        font-size: 30px;
        font-weight: 760;
        line-height: 1.12;
    }

    .header-meta,
    .header-actions,
    .panel-tools,
    .table-toolbar,
    .table-toolbar-actions {
        display: flex;
        align-items: center;
        gap: 10px;
        flex-wrap: wrap;
    }

    .header-meta {
        margin-top: 14px;
    }

    .meta-pill,
    .table-count,
    .level-pill {
        display: inline-flex;
        align-items: center;
        gap: 7px;
        min-height: 30px;
        border: 1px solid var(--line);
        border-radius: 999px;
        background: rgba(255, 255, 255, 0.78);
        color: var(--muted);
        font-size: 0.83rem;
        font-weight: 700;
        padding: 0 12px;
        white-space: nowrap;
    }

    .header-actions {
        justify-content: flex-end;
    }

    .stats-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(165px, 1fr));
        gap: 12px;
        margin-bottom: 16px;
    }

    .metric-panel {
        min-height: 96px;
        border: 1px solid var(--line);
        border-radius: var(--radius);
        background: var(--surface);
        box-shadow: 0 8px 25px rgba(24, 33, 47, 0.05);
        padding: 16px;
    }

    .metric-label {
        color: var(--muted);
        font-size: 0.77rem;
        font-weight: 700;
        text-transform: uppercase;
    }

    .metric-value {
        margin-top: 7px;
        color: var(--ink);
        font-size: clamp(1.28rem, 2vw, 1.8rem);
        font-weight: 800;
        line-height: 1;
    }

    .metric-subtle {
        margin-top: 6px;
        color: var(--muted);
        font-size: 0.83rem;
    }

    .control-panel,
    .panel {
        border: 1px solid var(--line);
        border-radius: var(--radius);
        background: rgba(255, 255, 255, 0.9);
        box-shadow: var(--shadow);
    }

    .control-panel {
        position: sticky;
        z-index: 10;
        top: 0;
        margin-bottom: 16px;
        padding: 14px;
        backdrop-filter: blur(14px);
    }

    .controls-grid {
        display: grid;
        grid-template-columns: minmax(220px, 1.25fr) repeat(2, minmax(135px, 0.55fr)) minmax(250px, 0.95fr) auto;
        gap: 12px;
        align-items: end;
    }

    .field-group {
        min-width: 0;
    }

    .control-label {
        display: block;
        margin-bottom: 6px;
        color: var(--muted);
        font-size: 0.78rem;
        font-weight: 800;
        text-transform: uppercase;
    }

    .metric-input,
    .search-input {
        width: 100%;
        min-height: 42px;
        border: 1px solid var(--line);
        border-radius: 7px;
        background: var(--surface);
        color: var(--ink);
        padding: 0.48rem 0.72rem;
        transition: border-color 150ms ease, box-shadow 150ms ease;
    }

    .metric-input:focus,
    .search-input:focus {
        border-color: var(--accent);
        box-shadow: 0 0 0 3px rgba(25, 114, 120, 0.16);
        outline: none;
    }

    .switch-stack {
        display: grid;
        gap: 8px;
    }

    .filter-switch {
        min-height: 42px;
        display: flex;
        align-items: center;
        gap: 9px;
        border: 1px solid var(--line);
        border-radius: 7px;
        background: var(--surface);
        padding: 0 12px 0 2.55rem;
    }

    .filter-switch .form-check-input {
        cursor: pointer;
    }

    .filter-switch .form-check-label {
        color: var(--ink);
        cursor: pointer;
        font-size: 0.92rem;
        font-weight: 700;
    }

    .control-actions {
        display: flex;
        justify-content: flex-end;
        gap: 8px;
    }

    .icon-button {
        min-height: 42px;
        display: inline-flex;
        align-items: center;
        justify-content: center;
        gap: 8px;
        border: 1px solid var(--line);
        border-radius: 7px;
        background: var(--surface);
        color: var(--ink);
        font-weight: 800;
        padding: 0 13px;
        transition: transform 150ms ease, border-color 150ms ease, background 150ms ease;
        white-space: nowrap;
    }

    .icon-button:hover {
        border-color: rgba(25, 114, 120, 0.45);
        background: var(--accent-soft);
        transform: translateY(-1px);
    }

    .icon-button--primary {
        border-color: var(--accent);
        background: var(--accent);
        color: #ffffff;
    }

    .icon-button--primary:hover {
        border-color: var(--accent-strong);
        background: var(--accent-strong);
        color: #ffffff;
    }

    .workspace-grid {
        display: grid;
        gap: 16px;
    }

    .panel {
        overflow: hidden;
    }

    .panel-heading {
        display: flex;
        justify-content: space-between;
        gap: 16px;
        align-items: center;
        border-bottom: 1px solid var(--line);
        padding: 15px 16px;
    }

    .panel-title {
        margin: 0;
        font-size: 1.02rem;
        font-weight: 800;
    }

    .panel-subtitle {
        margin: 4px 0 0;
        color: var(--muted);
        font-size: 0.88rem;
    }

    #sankey {
        width: 100%;
        height: clamp(440px, 60vh, 790px);
    }

    .legend-strip {
        display: flex;
        gap: 8px;
        align-items: center;
        flex-wrap: wrap;
        justify-content: flex-end;
    }

    .legend-chip {
        display: inline-flex;
        align-items: center;
        gap: 6px;
        color: var(--muted);
        font-size: 0.76rem;
        font-weight: 700;
    }

    .legend-dot {
        width: 9px;
        height: 9px;
        border-radius: 999px;
    }

    .table-panel .p-datatable {
        border: 0;
        border-radius: 0;
    }

    .table-toolbar {
        justify-content: space-between;
        padding: 2px 0;
    }

    .search-box {
        position: relative;
        width: min(360px, 100%);
    }

    .search-box .pi {
        position: absolute;
        top: 50%;
        left: 12px;
        color: var(--muted);
        transform: translateY(-50%);
    }

    .search-box .search-input {
        padding-left: 38px;
    }

    .p-datatable-table {
        font-size: 0.86rem;
    }

    .p-datatable .p-datatable-header {
        border: 0;
        border-bottom: 1px solid var(--line);
        background: var(--surface);
        padding: 13px 16px;
    }

    .p-datatable-table .p-datatable-thead > tr > th {
        border-color: var(--line);
        background: #edf4f1;
        color: #263440;
        font-size: 0.78rem;
        font-weight: 800;
        white-space: nowrap;
    }

    .p-datatable-table .p-datatable-tbody > tr > td {
        border-color: #edf1f4;
        color: #263440;
        font-size: 0.84rem;
        white-space: nowrap;
    }

    .p-datatable-table .p-datatable-tbody > tr:hover {
        background: #f4faf8;
    }

    .p-datatable .p-sortable-column.p-highlight {
        background: #dfeeea;
        color: var(--accent-strong);
    }

    .p-datatable .p-paginator {
        border: 0;
        border-top: 1px solid var(--line);
        background: var(--surface);
        padding: 10px 16px;
    }

    .p-column-filter {
        margin-top: 0.5rem;
    }

    .p-column-filter .p-inputtext {
        font-size: 0.75rem;
        padding: 0.25rem 0.5rem;
    }

    .p-inputtext,
    .p-dropdown,
    .p-multiselect {
        border: 1px solid var(--line);
        border-radius: 7px;
    }

    .p-multiselect {
        width: 100%;
        min-height: 42px;
    }

    .p-multiselect-label {
        padding: 0.55rem 0.75rem;
    }

    .p-selectbutton .p-button {
        border-color: var(--line);
        color: var(--ink);
        font-size: 0.86rem;
        font-weight: 800;
        padding: 0.48rem 0.8rem;
    }

    .p-selectbutton .p-button.p-highlight {
        border-color: var(--accent);
        background: var(--accent);
        color: #ffffff;
    }

    .level-pill {
        min-height: 26px;
        border-color: transparent;
        background: var(--accent-soft);
        color: var(--accent-strong);
        font-size: 0.76rem;
        padding: 0 10px;
    }

    .name-cell {
        display: inline-flex;
        align-items: center;
        gap: 8px;
        max-width: 520px;
    }

    .cell-text {
        overflow: hidden;
        text-overflow: ellipsis;
    }

    .pathogen-badge {
        display: inline-flex;
        align-items: center;
        min-height: 22px;
        border-radius: 999px;
        font-size: 0.68rem;
        font-weight: 800;
        padding: 0 8px;
        cursor: pointer;
    }

    .pathogen-badge--yes {
        background: #fde8df;
        color: #a33d20;
    }

    .pathogen-badge--watch {
        background: #edf1f4;
        color: #52616f;
    }

    .empty-state {
        color: var(--muted);
        font-weight: 700;
        padding: 24px;
        text-align: center;
    }

    .tooltip-inner {
        max-width: 560px;
        padding: 0;
        border-radius: 7px;
        overflow: hidden;
    }

    .tooltip-inner .tooltip-table-wrap {
        max-height: 500px;
        overflow: auto;
        display: block;
        background: #ffffff;
        color: var(--ink);
    }

    .tooltip-inner table {
        width: 100%;
        border-collapse: collapse;
        font-size: 0.8rem;
    }

    .tooltip-inner th,
    .tooltip-inner td {
        border: 1px solid var(--line);
        padding: 0.35rem 0.5rem;
        text-align: left;
        vertical-align: top;
        white-space: normal;
    }

    .tooltip-inner thead th {
        background: #edf4f1;
        color: var(--ink);
        font-weight: 800;
    }

    @media (max-width: 1100px) {
        .report-header,
        .controls-grid {
            grid-template-columns: 1fr;
        }

        .header-actions,
        .control-actions,
        .legend-strip {
            justify-content: flex-start;
        }

        .control-panel {
            position: static;
        }
    }

    @media (max-width: 700px) {
        .app-shell {
            width: min(100% - 20px, 1680px);
            padding-top: 18px;
        }

        .panel-heading,
        .table-toolbar {
            align-items: flex-start;
            flex-direction: column;
        }

        .icon-button {
            width: 100%;
        }

        .header-actions .p-selectbutton {
            width: 100%;
        }

        #sankey {
            height: 520px;
        }
    }
  </style>
</head>
<body>

<main class="app-shell">
    <div id="app" v-cloak>
        <header class="report-header">
            <div>
                <div class="eyebrow">GOTTCHA2 result viewer</div>
                <h1>TSVFILENAME</h1>
                <div class="header-meta">
                    <span class="meta-pill"><i class="pi pi-database"></i>{{ formatNumber(totalRecords) }} rows</span>
                    <span class="meta-pill"><i class="pi pi-sitemap"></i>{{ formatNumber(levelOptions.length) }} levels</span>
                    <span class="meta-pill"><i class="pi pi-filter"></i>{{ activeFilterLabel }}</span>
                </div>
            </div>
            <div class="header-actions">
                <p-select-button
                    v-model="activeView"
                    :options="viewOptions"
                    option-label="label"
                    option-value="value"
                    aria-label="View mode"
                ></p-select-button>
                <button type="button" class="icon-button" @click="resetFilters" title="Reset filters">
                    <i class="pi pi-refresh"></i>
                    <span>Reset</span>
                </button>
            </div>
        </header>

        <section class="stats-grid" aria-label="Filtered summary">
            <div class="metric-panel">
                <div class="metric-label">Pathogen</div>
                <div class="metric-value">{{ formatNumber(summaryStats.pathogens) }}</div>
                <div class="metric-subtle">{{ formatPercent(summaryStats.pathogenRate) }} of visible taxa</div>
            </div>
            <div class="metric-panel">
                <div class="metric-label">Visible taxa</div>
                <div class="metric-value">{{ formatNumber(summaryStats.rows) }}</div>
                <div class="metric-subtle">of {{ formatNumber(totalRecords) }} records</div>
            </div>
            <div class="metric-panel">
                <div class="metric-label">Read count</div>
                <div class="metric-value">{{ formatNumber(summaryStats.reads) }}</div>
                <div class="metric-subtle">filtered total</div>
            </div>
            <div class="metric-panel">
                <div class="metric-label">Mean SNI</div>
                <div class="metric-value">{{ formatDecimal(summaryStats.meanSni, 4) }}</div>
                <div class="metric-subtle">across visible rows</div>
            </div>
        </section>

        <section id="controls" class="control-panel" aria-label="Filters">
            <div class="controls-grid">
                <div class="field-group">
                    <label for="levelFilter" class="control-label">Levels</label>
                    <p-multi-select
                        v-model="selectedLevels"
                        :options="levelOptions"
                        option-label="label"
                        option-value="value"
                        input-id="levelFilter"
                        display="chip"
                        placeholder="Select levels"
                    ></p-multi-select>
                </div>

                <div class="field-group">
                    <label for="readCountMin" class="control-label">Read count &gt;=</label>
                    <input
                        type="number"
                        id="readCountMin"
                        class="metric-input"
                        v-model.number="rcMin"
                        min="RC_MIN"
                        max="RC_MAX"
                    >
                </div>

                <div class="field-group">
                    <label for="sniAdjMin" class="control-label">SNI score &gt;=</label>
                    <input
                        type="number"
                        id="sniAdjMin"
                        class="metric-input"
                        step="0.01"
                        v-model.number="sniMin"
                        min="SNI_MIN"
                        max="SNI_MAX"
                    >
                </div>

                <div class="switch-stack">
                    <div class="form-check form-switch filter-switch">
                        <input type="checkbox" role="switch" id="showValidOnly" class="form-check-input" v-model="showValidOnly">
                        <label for="showValidOnly" class="form-check-label">Qualified taxa only</label>
                    </div>
                    <div class="form-check form-switch filter-switch">
                        <input type="checkbox" role="switch" id="pathogenicOnly" class="form-check-input" v-model="showPathogenicOnly">
                        <label for="pathogenicOnly" class="form-check-label">Pathogen only</label>
                    </div>
                </div>

                <div class="control-actions">
                    <button type="button" class="icon-button" @click="clearSearch" title="Clear search">
                        <i class="pi pi-times"></i>
                        <span>Clear</span>
                    </button>
                    <button type="button" class="icon-button icon-button--primary" @click="exportCsv" title="Export filtered rows">
                        <i class="pi pi-download"></i>
                        <span>CSV</span>
                    </button>
                </div>
            </div>
        </section>

        <div class="workspace-grid">
            <section class="panel chart-panel" v-show="activeView !== 'table'">
                <div class="panel-heading">
                    <div>
                        <h2 class="panel-title">Taxonomic flow</h2>
                        <p class="panel-subtitle">{{ sankeySummary }}</p>
                    </div>
                    <div class="legend-strip" aria-label="Taxonomic level colors">
                        <span class="legend-chip" v-for="item in levelColorLegend" :key="item.level">
                            <span class="legend-dot" :style="{ background: item.color }"></span>
                            {{ item.label }}
                        </span>
                    </div>
                </div>
                <div id="sankey" v-once></div>
            </section>

            <section class="panel table-panel" v-show="activeView !== 'flow'">
                <p-data-table
                    v-model:filters="filters"
                    :value="filteredTableData"
                    :resizable-columns="true"
                    scrollable
                    scrollHeight="520px"
                    paginator
                    :rows="25"
                    :rows-per-page-options="[10, 25, 50, 100]"
                    :global-filter-fields="globalFilterFields"
                    filterDisplay="row"
                    responsive-layout="scroll"
                    striped-rows
                >
                    <template #header>
                        <div class="table-toolbar">
                            <div class="search-box">
                                <i class="pi pi-search"></i>
                                <input
                                    type="text"
                                    v-model="filters.global.value"
                                    placeholder="Search results"
                                    class="search-input"
                                />
                            </div>
                            <div class="table-toolbar-actions">
                                <span class="table-count">{{ formatNumber(filteredTableData.length) }} rows</span>
                                <button type="button" class="icon-button" @click="clearSearch" title="Clear search">
                                    <i class="pi pi-times"></i>
                                    <span>Clear</span>
                                </button>
                            </div>
                        </div>
                    </template>
                    <template #empty>
                        <div class="empty-state">No records match the current filters.</div>
                    </template>
                    <template #loading>
                        <div class="empty-state">Loading records...</div>
                    </template>

                    <p-column
                        field="LEVEL"
                        header="LEVEL"
                        :sortable="true"
                        frozen
                        alignFrozen="left"
                        style="min-width: 135px;"
                    >
                        <template #body="slotProps">
                            <span class="level-pill">{{ formatLevelLabel(slotProps.data['LEVEL']) }}</span>
                        </template>
                    </p-column>

                    <p-column
                        v-for="column in columns"
                        :key="column.field"
                        :field="column.field"
                        :class="column.class"
                        :header="column.header"
                        :sortable="true"
                        :showFilterMenu="false"
                        :filterField="column.field"
                    >
                        <template #body="slotProps">
                            <template v-if="column.field === 'NAME'">
                                <span class="name-cell">
                                    <span class="cell-text">{{ formatCellValue(slotProps.data[column.field], column.type) }}</span>
                                    <span
                                        v-if="getPathogenBadge(slotProps.data)"
                                        class="pathogen-badge"
                                        :class="getPathogenBadge(slotProps.data).class"
                                        v-bs-tooltip="{ title: formatPathogenicTooltip(slotProps.data['PATHOGENIC_INFO']) }"
                                    >
                                        {{ getPathogenBadge(slotProps.data).label }}
                                    </span>
                                </span>
                            </template>
                            <template v-else>
                                <span>{{ formatCellValue(slotProps.data[column.field], column.type) }}</span>
                            </template>
                        </template>
                    </p-column>
                </p-data-table>
            </section>
        </div>
    </div>

</main>
<script src="/publicdata/js/bootstrap.bundle.min.js"></script>

<script>
// Embed data
const records = RECORDS; // Records from Python
const levels = LEVELS_JSON;
const defaultLevels = DEFAULT_LEVELS_JSON;
const initialRcMin = RC_DEFAULT;
const initialSniMin = SNI_DEFAULT;
const reportFilename = REPORT_FILENAME_JSON;

const levelColorMap = {
  superkingdom: '#197278',
  kingdom: '#3f7cac',
  phylum: '#d89f2a',
  class: '#6a994e',
  order: '#d95d39',
  family: '#5e548e',
  genus: '#2a9d8f',
  species: '#e76f51',
  strain: '#8d6e63'
};
const fallbackLevelColors = ['#197278', '#3f7cac', '#d89f2a', '#6a994e', '#d95d39', '#5e548e', '#2a9d8f', '#e76f51'];
const plotConfig = {
  responsive: true,
  displaylogo: false,
  modeBarButtonsToRemove: ['lasso2d', 'select2d']
};

let vueApp; // Vue app instance
function initVueApp() {
    const { createApp, ref, computed, watch, nextTick } = Vue;

    const app = createApp({
        setup() {
            const tableData = ref(records);
            const filters = ref({
                global: { value: null, matchMode: 'contains' }
            });
            const activeView = ref('all');
            const viewOptions = [
                { label: 'All', value: 'all' },
                { label: 'Flow', value: 'flow' },
                { label: 'Table', value: 'table' }
            ];
            const globalFilterFields = records.length > 0
                ? Object.keys(records[0]).filter(key => key !== 'PATHOGENIC_INFO')
                : [];

            if (records.length > 0) {
                Object.keys(records[0]).forEach(key => {
                    if (key === 'PATHOGENIC_INFO') return; // hidden column for display
                    filters.value[key] = { value: null, matchMode: 'contains' };
                });
            }

            const availableLevels = computed(() => {
                return levels.map(level => ({
                    label: formatLevelLabel(level),
                    value: level
                }));
            });

            const selectedLevels = ref(
                defaultLevels && defaultLevels.length > 0 ? [...defaultLevels] : [...levels]
            );
            const rcMin = ref(initialRcMin);
            const sniMin = ref(initialSniMin);
            const showValidOnly = ref(true);
            const showPathogenicOnly = ref(false);

            const filteredTableData = computed(() => {
                return tableData.value;
            });

            const columns = computed(() => {
                if (records.length === 0) return [];
                return Object.keys(records[0])
                    .filter(key => key !== 'PATHOGENIC_INFO' && key !== 'LEVEL')
                    .map(key => ({
                        field: key,
                        header: key,
                        type: typeof records[0][key] === 'number' ? 'numeric' : 'text',
                        class: typeof records[0][key] === 'number' ? 'text-end' : ''
                    }));
            });

            const totalRecords = computed(() => records.length);
            const summaryStats = computed(() => {
                const rows = tableData.value || [];
                const reads = rows.reduce((total, row) => {
                    const value = Number(row.READ_COUNT);
                    return total + (Number.isFinite(value) ? value : 0);
                }, 0);
                const sniValues = rows
                    .map(row => Number(row.SNI_SCORE))
                    .filter(value => Number.isFinite(value));
                const meanSni = sniValues.length
                    ? sniValues.reduce((total, value) => total + value, 0) / sniValues.length
                    : 0;
                const pathogens = rows.filter(row => {
                    const humanPathogen = String(row.HUMAN_PATHOGEN || '').trim().toLowerCase();
                    return humanPathogen === 'yes' || hasPathogenicInfo(row.PATHOGENIC_INFO);
                }).length;

                return {
                    rows: rows.length,
                    reads,
                    meanSni,
                    pathogens,
                    pathogenRate: rows.length ? pathogens / rows.length : 0
                };
            });
            const activeFilterLabel = computed(() => {
                const selectedCount = selectedLevels.value && selectedLevels.value.length
                    ? selectedLevels.value.length
                    : levels.length;
                const pieces = [
                    `${selectedCount} levels`,
                    `RC >= ${formatNumber(rcMin.value)}`,
                    `SNI >= ${formatDecimal(sniMin.value, 2)}`
                ];
                if (showValidOnly.value) pieces.push('qualified');
                if (showPathogenicOnly.value) pieces.push('pathogen');
                return pieces.join(', ');
            });
            const sankeySummary = computed(() => {
                return `${formatNumber(summaryStats.value.rows)} taxa with ${formatNumber(summaryStats.value.reads)} reads`;
            });
            const levelColorLegend = computed(() => {
                return availableLevels.value.map((item, index) => ({
                    ...item,
                    level: item.value,
                    color: getLevelColor(item.value, 1, index)
                }));
            });

            function formatCellValue(value, type) {
                if (value === null || value === undefined) return '';

                if (type === 'numeric' && typeof value === 'number') {
                    if (value % 1 === 0) {
                        return value.toLocaleString();
                    } else {
                        return value.toLocaleString(undefined, {
                            minimumFractionDigits: 0,
                            maximumFractionDigits: 4
                        });
                    }
                }

                return String(value);
            }

            function formatNumber(value) {
                const number = Number(value);
                return Number.isFinite(number) ? number.toLocaleString() : '0';
            }

            function formatDecimal(value, digits = 2) {
                const number = Number(value);
                return Number.isFinite(number) ? number.toFixed(digits) : '0';
            }

            function formatPercent(value) {
                const number = Number(value);
                return Number.isFinite(number)
                    ? number.toLocaleString(undefined, { style: 'percent', maximumFractionDigits: 1 })
                    : '0%';
            }

            function formatLevelLabel(level) {
                if (level === null || level === undefined) return 'Unknown';
                const value = String(level);
                return value.charAt(0).toUpperCase() + value.slice(1);
            }

            function updateTableData(newData) {
                tableData.value = newData;
            }

            function getPathogenBadge(row) {
                if (!row || typeof row !== 'object') return null;
                const value = row['HUMAN_PATHOGEN'];
                if (value === null || value === undefined) return null;
                const normalized = String(value).trim().toLowerCase();
                if (normalized === 'yes') {
                    return { class: 'pathogen-badge--yes', label: 'Pathogen' };
                }
                if (normalized === 'no') {
                    return { class: 'pathogen-badge--watch', label: 'Info' };
                }
                return null;
            }

            function formatPathogenicTooltip(raw) {
                if (raw === null || raw === undefined) return '';
                const content = String(raw);
                if (!content.trim()) return '';
                return `<div class="tooltip-table-wrap">${content}</div>`;
            }

            function clearSearch() {
                if (filters.value.global) {
                    filters.value.global.value = null;
                }
            }

            function resetFilters() {
                selectedLevels.value = defaultLevels && defaultLevels.length > 0 ? [...defaultLevels] : [...levels];
                rcMin.value = initialRcMin;
                sniMin.value = initialSniMin;
                showValidOnly.value = true;
                showPathogenicOnly.value = false;
                clearSearch();
                nextTick(() => {
                    if (typeof updateAll === 'function') {
                        updateAll();
                    }
                });
            }

            function exportCsv() {
                const rows = filteredTableData.value || [];
                const fields = records.length > 0
                    ? Object.keys(records[0]).filter(key => key !== 'PATHOGENIC_INFO')
                    : [];
                if (!rows.length || !fields.length) return;

                const csvRows = [
                    fields.map(escapeCsvValue).join(','),
                    ...rows.map(row => fields.map(field => escapeCsvValue(row[field])).join(','))
                ];
                const blob = new Blob([csvRows.join('\n')], { type: 'text/csv;charset=utf-8;' });
                const url = URL.createObjectURL(blob);
                const link = document.createElement('a');
                const baseName = reportFilename.replace(/\.[^/.]+$/, '') || 'gottcha2-results';
                link.href = url;
                link.download = `${baseName}-filtered.csv`;
                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);
                URL.revokeObjectURL(url);
            }

            watch([selectedLevels, rcMin, sniMin, showValidOnly, showPathogenicOnly], () => {
                if (typeof updateAll === 'function') {
                    updateAll();
                }
            });
            watch(activeView, () => {
                nextTick(() => {
                    if (typeof updateAll === 'function') {
                        updateAll();
                    }
                });
            });

            return {
                tableData,
                filters,
                activeView,
                viewOptions,
                globalFilterFields,
                selectedLevels,
                rcMin,
                sniMin,
                showValidOnly,
                showPathogenicOnly,
                levelOptions: availableLevels,
                filteredTableData,
                columns,
                totalRecords,
                summaryStats,
                activeFilterLabel,
                sankeySummary,
                levelColorLegend,
                formatCellValue,
                formatNumber,
                formatDecimal,
                formatPercent,
                formatLevelLabel,
                updateTableData,
                getPathogenBadge,
                formatPathogenicTooltip,
                clearSearch,
                resetFilters,
                exportCsv,
            };
        },
    });

    app.use(PrimeVue.Config, {
        theme: {
            preset: PrimeUIX.Themes.Aura
        }
    });

    app.component('p-data-table', PrimeVue.DataTable);
    app.component('p-column', PrimeVue.Column);
    app.component('p-select-button', PrimeVue.SelectButton);
    app.component('p-multi-select', PrimeVue.MultiSelect);
    app.directive('bs-tooltip', {
        mounted(el, binding) {
            const title = binding.value?.title || '';
            if (!title) return;
            el._tooltipInstance = new bootstrap.Tooltip(el, {
                title,
                html: true,
                sanitize: false,
                container: 'body'
            });
        },
        updated(el, binding) {
            if (el._tooltipInstance) {
                el._tooltipInstance.dispose();
                el._tooltipInstance = null;
            }
            const title = binding.value?.title || '';
            if (!title) return;
            el._tooltipInstance = new bootstrap.Tooltip(el, {
                title,
                html: true,
                sanitize: false,
                container: 'body'
            });
        },
        unmounted(el) {
            if (el._tooltipInstance) {
                el._tooltipInstance.dispose();
                el._tooltipInstance = null;
            }
        }
    });

    const mountedApp = app.mount('#app');
    vueApp = mountedApp;
    return mountedApp;
}

function getFilterState() {
  if (vueApp) {
    return {
      selectedLevels: Array.isArray(vueApp.selectedLevels) ? vueApp.selectedLevels : [],
      rcMin: Number(vueApp.rcMin),
      sniMin: Number(vueApp.sniMin),
      showValidOnly: !!vueApp.showValidOnly,
      showPathogenicOnly: !!vueApp.showPathogenicOnly
    };
  }

  const levelSelect = document.getElementById('levelFilter');
  const selectedLevels = levelSelect
    ? Array.from(levelSelect.selectedOptions).map(o => o.value)
    : [];
  const rcInput = document.getElementById('readCountMin');
  const sniInput = document.getElementById('sniAdjMin');
  const validOnlyInput = document.getElementById('showValidOnly');
  const pathogenicOnlyInput = document.getElementById('pathogenicOnly');

  return {
    selectedLevels,
    rcMin: rcInput ? Number(rcInput.value) : 0,
    sniMin: sniInput ? Number(sniInput.value) : 0,
    showValidOnly: validOnlyInput ? validOnlyInput.checked : true,
    showPathogenicOnly: pathogenicOnlyInput ? pathogenicOnlyInput.checked : false
  };
}

function hasPathogenicInfo(info) {
  if (info === null || info === undefined) return false;
  if (typeof info === 'number') return Number.isFinite(info);
  if (typeof info === 'string') {
    const value = info.trim();
    return value !== '' && value.toLowerCase() !== 'nan';
  }
  return Boolean(info);
}

function escapeCsvValue(value) {
  if (value === null || value === undefined) return '';
  const text = String(value).replace(/\r?\n/g, ' ');
  if (/[",\n]/.test(text)) {
    return `"${text.replace(/"/g, '""')}"`;
  }
  return text;
}

function normalizeLevel(level) {
  return String(level || '').trim().toLowerCase();
}

function hexToRgba(hex, alpha) {
  const value = hex.replace('#', '');
  const red = parseInt(value.substring(0, 2), 16);
  const green = parseInt(value.substring(2, 4), 16);
  const blue = parseInt(value.substring(4, 6), 16);
  return `rgba(${red}, ${green}, ${blue}, ${alpha})`;
}

function getLevelColor(level, alpha = 1, fallbackIndex = 0) {
  const normalized = normalizeLevel(level);
  const hex = levelColorMap[normalized] || fallbackLevelColors[fallbackIndex % fallbackLevelColors.length];
  return alpha >= 1 ? hex : hexToRgba(hex, alpha);
}

function isElementVisible(element) {
  return !!(element && element.offsetParent !== null && element.clientWidth > 0 && element.clientHeight > 0);
}

function updateAll() {
  const state = getFilterState();

  // Get selected levels (array)
  let sel = state.selectedLevels;
  if (!Array.isArray(sel) || sel.length === 0) {
    sel = [...levels];
  }

  const rcMin = Number.isFinite(state.rcMin) ? state.rcMin : 0;
  const sniMin = Number.isFinite(state.sniMin) ? state.sniMin : 0;
  const showValidOnly = !!state.showValidOnly;
  const showPathogenicOnly = !!state.showPathogenicOnly;

  // Filter down records
  let currentRecords = records; // Start with all records

  // Apply "NOTE" filter first if checked
  if (showValidOnly) {
    currentRecords = currentRecords.filter(r => {
        // Ensure r.NOTE is a string and check if it includes "filtered" (case-insensitive)
        // If r.NOTE is null, undefined, or not a string, it's considered not to contain "filtered"
        return !(r.NOTE && typeof r.NOTE === 'string' && r.NOTE.toLowerCase().includes("filtered"));
    });
  }

  if (showPathogenicOnly) {
    currentRecords = currentRecords.filter(r => hasPathogenicInfo(r.PATHOGENIC_INFO));
  }

  // Apply other filters for Sankey and Vue table display
  const filteredData = currentRecords.filter(r =>
    sel.includes(r.LEVEL) &&
    r.READ_COUNT >= rcMin &&
    r.SNI_SCORE   >= sniMin
  );

  if (vueApp && vueApp.updateTableData) {
    vueApp.updateTableData(filteredData);
  }

  const sankeyEl = document.getElementById('sankey');
  if (!isElementVisible(sankeyEl)) {
    return;
  }

  const sankeyData = filteredData.filter(r => r.PARENT_NAME != null && r.NAME != null);
  const labels = Array.from(new Set(sankeyData.flatMap(r=>[r.PARENT_NAME,r.NAME])));
  const labelIndex = new Map(labels.map((label, idx) => [label, idx]));
  const nodeLevels = new Map();
  sankeyData.forEach(r => {
    if (r.NAME != null) nodeLevels.set(r.NAME, r.LEVEL);
    if (r.PARENT_NAME != null && !nodeLevels.has(r.PARENT_NAME)) nodeLevels.set(r.PARENT_NAME, r.PARENT_LEVEL || '');
  });

  const source=[], target=[], value=[], linkColor=[], customData=[];
  sankeyData.forEach(r => {
    const sourceIndex = labelIndex.get(r.PARENT_NAME);
    const targetIndex = labelIndex.get(r.NAME);
    const depthValue = Number(r.DEPTH);
    if (sourceIndex !== undefined && targetIndex !== undefined && Number.isFinite(depthValue)) {
        source.push(sourceIndex);
        target.push(targetIndex);
        value.push(depthValue);
        linkColor.push(getLevelColor(r.LEVEL, 0.28, source.length));
        customData.push([r.LEVEL || '', Number(r.READ_COUNT) || 0, Number(r.SNI_SCORE) || 0]);
    }
  });

  const sankeyLayout = {
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
    font: { family: 'Inter, Arial, sans-serif', size: 11, color: '#263440' },
    margin: { l: 22, r: 22, t: 18, b: 18 },
    height: Math.max(sankeyEl.clientHeight, 420),
    width: Math.max(sankeyEl.clientWidth, 320),
    annotations: source.length === 0 ? [{
      text: 'No flow records match the current filters.',
      x: 0.5,
      y: 0.5,
      xref: 'paper',
      yref: 'paper',
      showarrow: false,
      font: { size: 15, color: '#617080' }
    }] : []
  };

  if (source.length === 0) {
    Plotly.react('sankey', [], sankeyLayout, plotConfig);
    return;
  }

  const nodeColor = labels.map((label, index) => getLevelColor(nodeLevels.get(label), 0.9, index));
  Plotly.react('sankey',[{
    type:'sankey',
    orientation:'h',
    arrangement: 'snap',
    node: {
      pad:18,
      thickness:18,
      line:{color:'rgba(24,33,47,0.22)',width:0.7},
      label:labels,
      color: nodeColor,
      textposition:'outside',
      hovertemplate: '%{label}<extra></extra>'
    },
    link: {
      source,target,value,
      color: linkColor,
      customdata: customData,
      hovertemplate: '%{source.label} to %{target.label}<br>Depth: %{value:.4f}<br>Level: %{customdata[0]}<br>Read count: %{customdata[1]:,.0f}<br>SNI: %{customdata[2]:.4f}<extra></extra>'
    }
  }], sankeyLayout, plotConfig);
}

document.addEventListener('DOMContentLoaded', function() {
  initVueApp();

  Plotly.newPlot('sankey', [], {
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
    margin: { l: 22, r: 22, t: 18, b: 18 }
  }, plotConfig);
  updateAll();

  let resizeTimer = null;
  window.addEventListener('resize', function() {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(updateAll, 120);
  });
});
</script>
</body>
</html>
"""

def parse_args():
    p = argparse.ArgumentParser(
        description="Generate an interactive Sankey + Vue DataTable HTML from a TSV."
    )
    p.add_argument(
        '--tsv', '-i',
        required=True,
        help="Input TSV file"
    )
    p.add_argument(
        '--out', '-o',
        help="Output HTML file (default: [input_filename_without_extension].html)"
    )
    p.add_argument(
        '--rc-min', type=int, default=None,
        help="Override initial minimum READ_COUNT filter"
    )
    p.add_argument(
        '--sni-min', type=float, default=None,
        help="Override initial minimum SNI_SCORE filter"
    )
    p.add_argument(
        '-e', '--external', action='store_true',
        help='Use external resources for the HTML visualization'
    )

    return p.parse_args()

def main():
    global HTML_TEMPLATE
    args = parse_args()

    # --- Read TSV ---
    try:
        df = pd.read_csv(args.tsv, sep='\t')
    except Exception as e:
        sys.exit(f"ERROR: Failed to read TSV '{args.tsv}': {e}")

    # --- Determine output HTML filename ---
    if args.out:
        OUTPUT_HTML = args.out
    else:
        base_name = os.path.splitext(os.path.basename(args.tsv))[0]
        OUTPUT_HTML = f"{base_name}.html"

    # --- Ensure required columns exist ---
    required_cols = ['READ_COUNT', 'SNI_SCORE', 'LEVEL', 'PARENT_NAME', 'NAME']
    for col in required_cols:
        if col not in df.columns:
            sys.exit(f"ERROR: Required column '{col}' not found in '{args.tsv}'")

    # --- Normalize numeric columns used in filters/plots ---
    df['READ_COUNT'] = pd.to_numeric(df['READ_COUNT'], errors='coerce')
    df['SNI_SCORE'] = pd.to_numeric(df['SNI_SCORE'], errors='coerce')
    if 'DEPTH' not in df.columns:
        df['DEPTH'] = df['READ_COUNT']
    else:
        df['DEPTH'] = pd.to_numeric(df['DEPTH'], errors='coerce')

    # --- Derive filter bounds ---
    rc_min_val = df['READ_COUNT'].min()
    rc_max_val = df['READ_COUNT'].max()
    sni_min_val = df['SNI_SCORE'].min()
    sni_default = args.sni_min if args.sni_min is not None else 0.9

    rc_min = args.rc_min if args.rc_min is not None else (int(rc_min_val) if pd.notna(rc_min_val) else 0)
    rc_max = int(rc_max_val) if pd.notna(rc_max_val) else 0
    sni_min = args.sni_min if args.sni_min is not None else (float(sni_min_val) if pd.notna(sni_min_val) else 0.0)
    sni_max = 1.0 # SNI is typically between 0 and 1

    tsv_filename = os.path.basename(args.tsv)

    # --- Prepare NOTE column ---
    if 'NOTE' not in df.columns:
        df['NOTE'] = ''
    else:
        df['NOTE'] = df['NOTE'].fillna('').astype(str) # Ensure string type and no NaNs

    # --- Prepare data for JS ---
    records = df.astype(object).where(pd.notna(df), None).to_dict(orient='records')
    levels  = df['LEVEL'].dropna().astype(str).unique().tolist()

    taxonomic_order = [
        'superkingdom', 'kingdom', 'phylum', 'class',
        'order', 'family', 'genus', 'species', 'strain'
    ]
    order_map = {lvl: idx for idx, lvl in enumerate(taxonomic_order)}
    levels.sort(key=lambda lvl: (order_map.get(lvl.lower(), len(order_map)), lvl.lower()))

    # Define default selected levels (fallback to all if none are present)
    default_level_names = {'family', 'genus', 'species'}
    default_levels = [lvl for lvl in levels if lvl.lower() in default_level_names]
    if not default_levels:
        default_levels = levels[:]

    HTML_TEMPLATE = HTML_TEMPLATE.replace('TSVFILENAME', tsv_filename)
    HTML_TEMPLATE = HTML_TEMPLATE.replace('REPORT_FILENAME_JSON', json.dumps(tsv_filename))
    HTML_TEMPLATE = HTML_TEMPLATE.replace('RECORDS', json.dumps(records, allow_nan=False))
    HTML_TEMPLATE = HTML_TEMPLATE.replace('DEFAULT_LEVELS_JSON', json.dumps(default_levels))
    HTML_TEMPLATE = HTML_TEMPLATE.replace('LEVELS_JSON', json.dumps(levels))
    HTML_TEMPLATE = HTML_TEMPLATE.replace('RC_DEFAULT', str(rc_min))
    HTML_TEMPLATE = HTML_TEMPLATE.replace('RC_MIN', str(rc_min))
    HTML_TEMPLATE = HTML_TEMPLATE.replace('RC_MAX', str(rc_max))
    HTML_TEMPLATE = HTML_TEMPLATE.replace('SNI_DEFAULT', str(sni_default))
    HTML_TEMPLATE = HTML_TEMPLATE.replace('SNI_MIN', str(sni_min))
    HTML_TEMPLATE = HTML_TEMPLATE.replace('SNI_MAX', str(sni_max))

    if args.external:
        # Use CDN links for external resources
        html_content = HTML_TEMPLATE.replace(
            '<script src="/publicdata/js/d3.v7.min.js"></script>',
            '<script src="https://d3js.org/d3.v7.min.js"></script>'
        ).replace(
            '<script src="/publicdata/js/vue.global.prod.js"></script>',
            '<script src="https://unpkg.com/vue@3/dist/vue.global.prod.js"></script>'
        ).replace(
            '<script src="/publicdata/js/plotly-3.0.1.min.js" charset="utf-8"></script>',
            '<script src="https://cdn.plot.ly/plotly-3.0.1.min.js" charset="utf-8"></script>'
        ).replace(
            '<script src="/publicdata/js/primevue.min.js"></script>',
            '<script src="https://unpkg.com/primevue@4/umd/primevue.min.js"></script>'
        ).replace(
            '<script src="/publicdata/js/aura.js"></script>',
            '<script src="https://unpkg.com/@primeuix/themes/umd/aura.js"></script>'
        ).replace(
            '<link href="/publicdata/css/bootstrap.min.css" rel="stylesheet">',
            '<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">'
        ).replace(
            '<link rel="stylesheet" href="/publicdata/css/primeicons.css">',
            '<link rel="stylesheet" href="https://unpkg.com/primeicons@7.0.0/primeicons.css">'
        ).replace(
            '<script src="/publicdata/js/bootstrap.bundle.min.js"></script>',
            '<script src="https://cdn.jsdelivr.net/npm/bootstrap@5/dist/js/bootstrap.bundle.min.js"></script>'
        )
    else:
        html_content = HTML_TEMPLATE

    html_content = minify_html.minify(html_content)

    # --- Write output HTML ---
    with open(OUTPUT_HTML, 'w', encoding='utf-8') as f:
        f.write(html_content)
        print(f"INFO: Generated {OUTPUT_HTML}. Open it in your browser.")

if __name__ == '__main__':
    main()
