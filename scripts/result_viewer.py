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
    body {
        font-family: sans-serif;
        margin: 20px;
    }
    h1 {
        text-align: center;
        margin-bottom: 30px;
    }
    #controls {
        margin-bottom: 20px;
    }
    #controls .form-label,
    #controls .form-check-label,
    #controls {
        font-size: 0.9rem;
    }
    #sankey {
        width: 100%;
        min-height: 600px;
        height: 800px;
        margin-bottom: 30px;
    }
    
    /* PrimeVue DataTable Custom Styles */
    .p-datatable-table {
        font-size: 0.875rem;
        border-radius: 0.375rem;
        overflow: hidden;
    }
    .p-datatable-table .p-datatable-thead > tr > th {
        background-color: #f8f9fa;
        border-bottom: 2px solid #dee2e6;
        font-weight: 600;
        font-size: 0.875rem;
        white-space: nowrap;
    }
    .p-datatable-table .p-datatable-tbody > tr > td {
        font-size: 0.875rem;
        white-space: nowrap;
    }
    .p-datatable-table .p-datatable-tbody > tr:hover {
        background-color: #f5f5f5;
    }
    .table-controls {
        display: flex;
        justify-content: space-between;
        align-items: center;
        margin-bottom: 1rem;
        gap: 1rem;
        flex-wrap: wrap;
    }
    .search-input {
        max-width: 300px;
    }
    .highlight {
        background-color: #f0e68c;
        font-weight: bold;
    }
    /* Custom PrimeVue input styling */
    .p-inputtext {
        border: 1px solid #ced4da;
        border-radius: 0.375rem;
        padding: 0.375rem 0.75rem;
    }
    .p-dropdown {
        border: 1px solid #ced4da;
        border-radius: 0.375rem;
    }
    .p-multiselect {
        border: 1px solid #ced4da;
        border-radius: 0.375rem;
    }
    /* SelectButton styling */
    .p-selectbutton .p-button {
        font-size: 0.875rem;
        padding: 0.375rem 0.75rem;
    }
    .p-datatable .p-sortable-column.p-highlight {
        background-color: #020617;
    }
    /* Filter styling */
    .p-column-filter {
        margin-top: 0.5rem;
    }
    .p-column-filter .p-inputtext {
        font-size: 0.75rem;
        padding: 0.25rem 0.5rem;
    }
    .p-column-filter .p-inputnumber {
        font-size: 0.75rem;
    }
    .p-datatable .p-datatable-thead > tr > th .p-column-header-content {
        display: flex;
        flex-direction: column;
        align-items: stretch;
    }
    .pathogen-badge {
        margin-left: 0.5rem;
        cursor: pointer;
    }
    .tooltip-inner {
        max-width: 520px;
        padding: 0;
    }
    .tooltip-inner .tooltip-table-wrap {
        max-height: 500px;
        overflow: auto;
        display: block;
    }
    .tooltip-inner table {
        width: 100%;
        border-collapse: collapse;
        font-size: 0.8rem;
    }
    .tooltip-inner th, .tooltip-inner td {
        border: 1px solid #dee2e6;
        padding: 0.35rem 0.5rem;
        text-align: left;
        vertical-align: top;
        white-space: normal;
    }
    .tooltip-inner thead th {
        background: gray;
        font-weight: 600;
    }
  </style>
</head>
<body>

<div class="container">
    <!-- <h1 class="mb-4">GOTTCHA2: TSVFILENAME</h1> -->

    <div id="app">
        <div id="controls" class="row row-cols-auto g-3 align-items-center justify-content-center mb-4">
            <div class="col">
                <label for="levelFilter" class="form-label mb-0">Levels:</label>
            </div>
            <div class="col">
                <p-multi-select
                    v-model="selectedLevels"
                    :options="levelOptions"
                    option-label="label"
                    option-value="value"
                    input-id="levelFilter"
                    display="chip"
                    placeholder="Select levels"
                    class="w-100"
                ></p-multi-select>
            </div>

            <div class="col">
                <label for="readCountMin" class="form-label mb-0">Read Count ≥</label>
            </div>
            <div class="col">
                <input
                    type="number"
                    id="readCountMin"
                    class="form-control"
                    v-model.number="rcMin"
                    min="RC_MIN"
                    max="RC_MAX"
                    style="width: 90px;"
                >
            </div>

            <div class="col">
                <label for="sniAdjMin" class="form-label mb-0">SNI score ≥</label>
            </div>
            <div class="col">
                <input
                    type="number"
                    id="sniAdjMin"
                    class="form-control"
                    step="0.01"
                    v-model.number="sniMin"
                    min="SNI_MIN"
                    max="SNI_MAX"
                    style="width: 90px;"
                >
            </div>

            <div class="col">
                <div class="form-check">
                    <input type="checkbox" id="showValidOnly" class="form-check-input" v-model="showValidOnly">
                    <label for="showValidOnly" class="form-check-label">Qualified Taxa Only</label>
                </div>
            </div>
            <div class="col">
                <div class="form-check">
                    <input type="checkbox" id="pathogenicOnly" class="form-check-input" v-model="showPathogenicOnly">
                    <label for="pathogenicOnly" class="form-check-label">Pathogen Only</label>
                </div>
            </div>
        </div>

        <div id="sankey" v-once class="mb-4"></div>

        <p-data-table 
            v-model:filters="filters" 
                   :value="filteredTableData"
                   :sortable="true"
                   :resizable-columns="true"
            scrollable
            scrollHeight="400px"
            paginator
            :rows="10"
            :rows-per-page-options="[10, 25, 50, 100]"
            :global-filter-fields="['NAME', 'PARENT_NAME', 'LEVEL', 'TAXID']"
            filterDisplay="row"
            responsive-layout="scroll"
            striped-rows
            class="mt-6"
        >
            <template #header>
                <div class="d-flex justify-content-between align-items-center sm">
                    <div class="d-flex align-items-center">
                        <input 
                            type="text" 
                            v-model="filters.global.value" 
                            placeholder="Search keywords..." 
                            class="form-control form-control-sm"
                            style="width: 200px;"
                        />
                    </div>
                    <div class="d-flex align-items-center">
                        <p-select-button 
                            v-model="selectedLevels" 
                            :options="levelOptions" 
                            option-label="label" 
                            option-value="value"
                            :multiple="true"
                            class="mt-6 sm"
                        />
                    </div>
                </div>
            </template>
            <template #empty> No records found. </template>
            <template #loading> Loading... </template>

            <!-- Frozen First Column -->
            <p-column 
                field="LEVEL" 
                header="LEVEL" 
                :sortable="true"
                :showFilterMenu="true"
                frozen 
                alignFrozen="left"
                style="min-width: 120px;"
            >
                <template #body="slotProps">
                    <span>{{ formatCellValue(slotProps.data['LEVEL'], 'text') }}</span>
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
                        <span>{{ formatCellValue(slotProps.data[column.field], column.type) }}</span>
                        <span
                            v-if="getPathogenBadge(slotProps.data)"
                            class="badge pathogen-badge"
                            :class="getPathogenBadge(slotProps.data).class"
                            v-bs-tooltip="{ title: formatPathogenicTooltip(slotProps.data['PATHOGENIC_INFO']) }"
                        >
                            {{ getPathogenBadge(slotProps.data).label }}
                        </span>
                    </template>
                    <template v-else>
                        <span>{{ formatCellValue(slotProps.data[column.field], column.type) }}</span>
                    </template>
                </template>
                
            </p-column>
        </p-data-table>
    </div>

</div>
<script src="/publicdata/js/bootstrap.bundle.min.js"></script>

<script>
// Embed data
const records = RECORDS; // Records from Python
const levels = LEVELS_JSON;
const defaultLevels = DEFAULT_LEVELS_JSON;
const initialRcMin = RC_DEFAULT;
const initialSniMin = SNI_DEFAULT;

let vueApp; // Vue app instance
function initVueApp() {
    const { createApp, ref, computed, watch } = Vue;

    const app = createApp({
        setup() {
            // Core reactive data
            const tableData = ref(records);
            
            // Initialize filters for PrimeVue DataTable
            const filters = ref({
                global: { value: null, matchMode: 'contains' }
            });
            
            // Initialize column filters dynamically
            if (records.length > 0) {
                Object.keys(records[0]).forEach(key => {
                    if (key === 'PATHOGENIC_INFO') return; // hidden column for display
                    filters.value[key] = { value: null, matchMode: 'contains' };
                });
            }
            
            // Level filter data
            const availableLevels = computed(() => {
                return levels.map(level => ({
                    label: level.charAt(0).toUpperCase() + level.slice(1),
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

            // Filter data by selected levels
            const filteredTableData = computed(() => {
                if (!selectedLevels.value || selectedLevels.value.length === 0) {
                    return tableData.value;
                }
                return tableData.value.filter(row => selectedLevels.value.includes(row.LEVEL));
            });
            
            // Table columns computed from data
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
            
            // Methods
            function formatCellValue(value, type) {
                if (value === null || value === undefined) return '';
                
                if (type === 'numeric' && typeof value === 'number') {
                    // Format numbers with appropriate decimal places
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
            
            function updateTableData(newData) {
                tableData.value = newData;
            }

            function getPathogenBadge(row) {
                if (!row || typeof row !== 'object') return null;
                const value = row['HUMAN_PATHOGEN'];
                if (value === null || value === undefined) return null;
                const normalized = String(value).trim().toLowerCase();
                if (normalized === 'yes') {
                    return { class: 'bg-danger', label: 'Pathogen' };
                }
                if (normalized === 'no') {
                    return { class: 'bg-secondary', label: 'Pathogen' };
                }
                return null;
            }

            function formatPathogenicTooltip(raw) {
                if (raw === null || raw === undefined) return '';
                const content = String(raw);
                if (!content.trim()) return '';
                return `<div class="tooltip-table-wrap">${content}</div>`;
            }

            watch([selectedLevels, rcMin, sniMin, showValidOnly, showPathogenicOnly], () => {
                if (typeof updateAll === 'function') {
                    updateAll();
                }
            });
            
            return {
                tableData,
                filters,
                selectedLevels,
                rcMin,
                sniMin,
                showValidOnly,
                showPathogenicOnly,
                levelOptions: availableLevels,
                filteredTableData,
                columns,
                formatCellValue,
                updateTableData,
                getPathogenBadge,
                formatPathogenicTooltip,
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

  // Rebuild Sankey diagram
  // Filter out records with missing PARENT_NAME or NAME for Sankey to prevent errors
  const sankeyData = filteredData.filter(r => r.PARENT_NAME != null && r.NAME != null);
  const labels = Array.from(new Set(sankeyData.flatMap(r=>[r.PARENT_NAME,r.NAME])));
  const labelIndex = new Map(labels.map((label, idx) => [label, idx]));
  const source=[], target=[], value=[];
  sankeyData.forEach(r => {
    const sourceIndex = labelIndex.get(r.PARENT_NAME);
    const targetIndex = labelIndex.get(r.NAME);
    const depthValue = Number(r.DEPTH);
    // Only push if both source and target labels were found and depth is numeric
    if (sourceIndex !== undefined && targetIndex !== undefined && Number.isFinite(depthValue)) {
        source.push(sourceIndex);
        target.push(targetIndex);
        value.push(depthValue);
    }
  });

  const sankeyLayout = {
    font: { size: 10 },
    margin: { l: 20, r: 20, t: 10, b: 10 },
    height: document.getElementById('sankey').clientHeight,
    width: document.getElementById('sankey').clientWidth
  };

  Plotly.react('sankey',[{
    type:'sankey',
    orientation:'h',
    node: {
      pad:15,
      thickness:20,
      line:{color:'black',width:0.5},
      label:labels,
      textposition:'outside',
    },
    link: { source,target,value }
  }], sankeyLayout);

  // Update Vue table with filtered data
  if (vueApp && vueApp.updateTableData) {
    vueApp.updateTableData(filteredData);
  }
}

document.addEventListener('DOMContentLoaded', function() {
  // Initialize Vue app
  const vueAppInstance = initVueApp();
  
  // Initialize empty Sankey then first draw
  Plotly.newPlot('sankey', []); // Initialize with empty data and layout
  updateAll(); // Call once to populate with initial filters

  // Optional: Resize Plotly chart when window resizes
  window.addEventListener('resize', function() {
    // Check if the Plotly chart div has content (i.e., Plotly initialized)
    if (document.getElementById('sankey').children.length > 0) {
        Plotly.Plots.resize(document.getElementById('sankey'));
    }
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

    # html_content = minify_html.minify(html_content)

    # --- Write output HTML ---
    with open(OUTPUT_HTML, 'w', encoding='utf-8') as f:
        f.write(html_content)
        print(f"INFO: Generated {OUTPUT_HTML}. Open it in your browser.")

if __name__ == '__main__':
    main()
