#!/usr/bin/env python3

import argparse
import os
import json
import re
import pandas as pd
from collections import defaultdict
import minify_html

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Genome Coverage</title>
    <script src="/publicdata/js/d3.v7.min.js"></script>
    <script src="/publicdata/js/vue.global.prod.js"></script>
    <script src="/publicdata/js/primevue.min.js"></script>
    <script src="/publicdata/js/aura.js"></script>
    <link href="/publicdata/css/bootstrap.min.css" rel="stylesheet">
    <link rel="stylesheet" href="/publicdata/css/primeicons.css">
    <style>
        body {
            font-family: sans-serif, 'Segoe UI', Tahoma, Geneva, Verdana;
            margin: 0;
            padding: 0;
            background-color: #fff; /* off-white background */
            color: #4a4a4a;
        }
        
        .container {
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }
        
        header {
            background-color: #a7aca2; /* sage green */
            color: white;
            padding: 15px 20px;
            border-radius: 5px;
            margin-bottom: 20px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }
        
        h1 {
            margin: 0;
            font-size: 24px;
        }
        
        .controls {
            background-color: white;
            padding: 15px;
            border-radius: 5px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            margin-bottom: 20px;
        }
        
        select, button {
            padding: 8px 12px;
            margin-right: 10px;
            border: 1px solid #ddd;
            border-radius: 4px;
            background-color: white;
        }
        
        .btn-primary {
            background-color: #b5a397; /* dusty rose */
            color: white;
            cursor: pointer;
            border: none;
            transition: background-color 0.3s;
        }
        
        .btn-primary:hover {
            background-color: #9c8a7e; /* Darker dusty rose */
        }
        
        .visualization {
            background-color: white;
            padding: 15px;
            border-radius: 5px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            # overflow: hidden;
        }
        
        .tooltip {
            position: absolute;
            background-color: rgba(255, 255, 255, 0.9);
            border: 1px solid #ddd;
            border-radius: 4px;
            padding: 10px;
            pointer-events: none;
            font-size: 12px;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
            z-index: 1000;
        }
        
        .fragment {
            cursor: pointer;
            transition: opacity 0.3s;
        }
        
        .fragment:hover {
            opacity: 0.9;
        }

        .fragment-indicator {
            cursor: pointer;
            transition: opacity 0.3s;
        }
        .fragment-indicator:hover {
            opacity: 0.9;
        }

        .legend {
            margin-top: 20px;
            display: flex;
            justify-content: center;
            align-items: center;
            flex-wrap: wrap;
            font-size: 0.8em;
        }
        
        .legend-item {
            display: flex;
            align-items: center;
            margin: 0 15px 10px 0;
        }
        
        .legend-color {
            width: 15px;
            height: 15px;
            margin-right: 5px;
            border-radius: 2px;
        }
        
        .info-panel {
            margin-top: 20px;
            padding: 15px;
            background-color: #f1f1f1;
            border-radius: 5px;
        }
        
        .axis line, .axis path {
            stroke: #ccc;
        }
        
        .axis text {
            font-size: 12px;
            fill: #666;
        }
        
        .brush .selection {
            stroke: #a7aca2; /* sage green */
            stroke-opacity: 0.6;
            fill: #a7aca2;
            fill-opacity: 0.1;
        }
        
        .annotation-item { 
            margin-bottom: 15px;
            border-bottom: 1px dashed #ddd;
            padding-bottom: 10px;
        }
        
        .annotation-item:last-child {
            border-bottom: none;
        }
        
        #annotation-loading {
            text-align: center;
            padding: 20px;
        }
        
        .card {
            margin-bottom: 20px;
        }
        
        .card-header {
            background-color: #a7aca2; /* sage green */
            color: white;
        }
        
        /* PrimeVue Select custom styles */
        .p-select {
            background-color: white;
            border: 1px solid #ddd;
            border-radius: 4px;
            width: 100%;
        }
        
        .p-select .p-select-label {
            padding: 8px 12px;
            color: #4a4a4a;
            font-weight: 500;
        }
        
        .p-select .p-select-dropdown {
            border: 1px solid #ddd;
            border-radius: 4px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }
        
        .p-select-option {
            padding: 8px 12px;
        }
        
        .p-select-option.p-focus {
            background-color: #a7aca2; /* sage green to match theme */
            color: white;
        }
        
        /* Genome option styles */
        .genome-option {
            display: flex;
            justify-content: space-between;
            align-items: center;
            width: 100%;
        }
        
        .genome-name {
            flex-grow: 1;
        }
        
        .genome-taxid {
            color: #6c757d;
            margin-right: 10px;
            font-size: 0.85em;
        }
        
        .read-count {
            background-color: #a7aca2;
            color: white;
            font-size: 0.8em;
            padding: 2px 6px;
            border-radius: 10px;
            white-space: nowrap;
        }
    </style>
</head>
<body>
    <div id="app" class="container">        
        <div class="visualization controls mb-4">
            <div class="card-body">
                <div class="row">
                    <div class="col-md-8">
                        <p-select
                            v-model="selectedGenome"
                            :options="genomeOptions"
                            option-label="label"
                            option-value="value"
                            placeholder="Search or select a genome..."
                            filter
                            @change="handleGenomeSelect"
                            class="w-100"
                        >
                            <template #option="slotProps">
                                <div class="genome-option">
                                    <span class="genome-name">{{ slotProps.option.name }}</span>
                                    <span class="genome-taxid">({{ slotProps.option.taxid }})</span>
                                    <span class="read-count">{{ slotProps.option.readCount.toLocaleString() }} reads</span>
                                </div>
                            </template>
                        </p-select>
                    </div>
                    <div class="col-md-4">
                        <button @click="resetZoom" class="btn btn-primary">Reset Zoom</button>
                    </div>
                </div>
            </div>
        </div>
        
        <div class="card visualization mb-4">
            <div class="card-body">
                <div id="coverage-plot"></div>
                <div class="legend" id="legend"></div>
            </div>
        </div>
        
        <div class="card" v-if="currentGenome">
            <div class="card-header">
                Genome Signature Information
            </div>
            <div class="card-body">
                <div class="row">
                    <div class="col-md-6">
                        <p><strong>Name:</strong> {{ currentGenome.name }}</p>
                        <p><strong>Taxid:</strong> {{ currentGenome.taxid }}</p>
                        <p><strong>Domain:</strong> {{ currentGenome.superkingdom }}</p>
                        <p><strong>Genome Size:</strong> {{ currentGenome.genomeSize.toLocaleString() }} bp</p>
                    </div>
                    <div class="col-md-6">
                        <p><strong>Signature Level:</strong> {{ currentGenome.db_level }}</p>
                        <p><strong>Total Signature Fragments:</strong> {{ currentGenome.numOfSeq.toLocaleString() }}</p>
                        <p><strong>Total Signature Length:</strong> {{ currentGenome.totalLength.toLocaleString() }} bp</p>
                        
                    </div>
                </div>
            </div>
        </div>
        <div class="card" v-else>
            <div class="card-header">
                Genome Signature Information
            </div>
            <div class="card-body">
                <p class="text-muted">Select a genome to view information</p>
            </div>
        </div>
    </div>

    <!-- NCBI Annotation Query Modal -->
    <div class="modal fade" id="annotationModal" tabindex="-1" aria-labelledby="annotationModalLabel" aria-hidden="true">
        <div class="modal-dialog modal-lg">
            <div class="modal-content">
                <div class="modal-header">
                    <h5 class="modal-title" id="annotationModalLabel">NCBI Nucleotide Annotation Query</h5>
                    <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="Close"></button>
                </div>
                <div class="modal-body">
                    <div id="annotation-loading" class="d-none">
                        <div class="spinner-border" role="status">
                            <span class="visually-hidden">Loading...</span>
                        </div>
                        <p>Fetching annotations...</p>
                    </div>
                    <div id="annotation-results"></div>
                </div>
                <div class="modal-footer">
                    <button type="button" class="btn btn-secondary" data-bs-dismiss="modal">Close</button>
                </div>
            </div>
        </div>
    </div>

    <!-- Bootstrap 5 JS Bundle with Popper -->
    <script src="/publicdata/js/bootstrap.bundle.min.js"></script>
    <script>
        // The data will be injected here
        const genomeData = GENOME_DATA_PLACEHOLDER;
        const coverageData = COVERAGE_DATA_PLACEHOLDER;
        
        // Initialize D3 variables
        let svg, xScale, yScale, xAxis, yAxis, xAxisTop, brush;
        const margin = {top: 50, right: 30, bottom: 50, left: 60};
        const width = 1100 - margin.left - margin.right;
        const height = 400 - margin.top - margin.bottom;
        
        // Initialize Vue app
        const { createApp, ref, computed, onMounted, watch } = Vue;
        
        // Bootstrap modal instance
        let annotationModal;
        
        document.addEventListener('DOMContentLoaded', function() {
            // Initialize Bootstrap modal
            annotationModal = new bootstrap.Modal(document.getElementById('annotationModal'));
        });
        
        // Reactive state variables accessible globally
        let selectedTaxid, currentGenome, currentFragment;
        
        // Register Vue-Select component
        const app = createApp({
            setup() {
                // Reactive state - assign to outer variables for global access
                selectedTaxid = ref('');
                currentGenome = ref(null);
                currentFragment = ref(null);
                const selectedGenome = ref(null);
                
                // Calculate unique taxids and read counts
                const uniqueTaxidsInCoverage = [...new Set(coverageData.map(item => item.genome_taxid))];
                
                // Calculate total read counts for each genome
                const genomeTotalReads = {};
                coverageData.forEach(item => {
                    if (!genomeTotalReads[item.genome_taxid]) {
                        genomeTotalReads[item.genome_taxid] = 0;
                    }
                    genomeTotalReads[item.genome_taxid] += item.numreads;
                });
                
                // Prepare genomes with read counts
                const genomesWithReadCounts = genomeData
                    .filter(genome => uniqueTaxidsInCoverage.includes(genome.taxid))
                    .map(genome => ({
                        ...genome,
                        readCount: genomeTotalReads[genome.taxid] || 0
                    }))
                    .sort((a, b) => b.readCount - a.readCount); // Sort by read count descending
                
                // Format options for vue-select
                const genomeOptions = genomesWithReadCounts.map(genome => ({
                    value: genome.taxid,
                    label: `${genome.name} (${genome.taxid})`,
                    name: genome.name,
                    taxid: genome.taxid,
                    readCount: genome.readCount
                }));
                
                // Methods
                function handleGenomeChange() {
                    console.log("Selected TaxID:", selectedTaxid.value);
                    if (selectedTaxid.value) {
                        currentGenome.value = genomeData.find(g => g.taxid === selectedTaxid.value);
                        updateVisualization();
                    }
                }

                function handleGenomeSelect(eventOrValue) {
                    const value = eventOrValue && typeof eventOrValue === "object" && "value" in eventOrValue
                        ? eventOrValue.value
                        : eventOrValue;

                    if (!value) return;

                    selectedTaxid.value = String(value);
                    selectedGenome.value = selectedTaxid.value;
                    currentGenome.value = genomeData.find(g => String(g.taxid) === selectedTaxid.value);
                    updateVisualization();
                }

                function resetZoom() {
                    if (currentGenome.value) {
                        xScale.domain([0, currentGenome.value.totalLength]);
                        updateChart();
                    }
                }
                
                onMounted(() => {
                    initVisualization();
                    if (genomeOptions.length > 0) {
                        selectedGenome.value = genomeOptions[0].value;
                        selectedTaxid.value = genomeOptions[0].value;
                        currentGenome.value = genomeData.find(g => g.taxid === selectedTaxid.value);
                        updateVisualization();
                    }
                });

                return {
                    selectedTaxid,
                    selectedGenome,
                    genomeOptions,
                    currentGenome,
                    currentFragment,
                    genomesWithReadCounts,
                    handleGenomeChange,
                    handleGenomeSelect,
                    resetZoom
                };
            },
            components: {
                'p-select': PrimeVue.Select,
            }
        });
        
        // Configure PrimeVue with theme
        app.use(PrimeVue.Config, {
            theme: {
                preset: PrimeUIX.Themes.Aura,
                options: {
                    darkModeSelector: '.dark-mode'
                }
            }
        });
        
        // Mount Vue app
        app.mount('#app');
        
        // Create tooltip
        const tooltip = d3.select('body').append('div')
            .attr('class', 'tooltip')
            .style('opacity', 0);
        
        // Initialize the visualization container
        function initVisualization() {
            d3.select("#coverage-plot").html("");
            
            svg = d3.select("#coverage-plot")
                .append("svg")
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom)
                .append("g")
                .attr("transform", `translate(${margin.left},${margin.top})`);
            
            // Create scales
            xScale = d3.scaleLinear().range([0, width]);
            yScale = d3.scaleLinear()
                .domain([0, 0.1]) // Will be adjusted based on data
                .range([height, 0]);
            
            // Create axes
            xAxis = svg.append("g")
                .attr("class", "axis x-axis")
                .attr("transform", `translate(0,${height})`);
            
            xAxisTop = svg.append("g")
                .attr("class", "axis x-axis-top")
                .attr("transform", "translate(0,0)");
            
            yAxis = svg.append("g")
                .attr("class", "axis y-axis");
            
            // Add axis labels
            svg.append("text")
                .attr("transform", `translate(${width/2},${height + 40})`)
                .style("text-anchor", "middle")
                .style("font-size", "14px")
                .style("fill", "#666")
                .text("Total Signature (bp)");
            
            svg.append("text")
                .attr("transform", `translate(${width/2},-30)`)
                .style("text-anchor", "middle")
                .style("font-size", "14px")
                .style("fill", "#666")
                .text("Signature Coverage (%)");
            
            svg.append("text")
                .attr("transform", "rotate(-90)")
                .attr("y", 0 - margin.left)
                .attr("x", 0 - (height / 2) - 5)
                .attr("dy", "1em")
                .style("fill", "#666")
                .style("text-anchor", "middle")
                .style("font-size", "14px")
                .text("Mean Depth (x)");
            
            // Create brush for zooming
            brush = d3.brushX()
                .extent([[0, 0], [width, height]])
                .on("end", brushended);
            
            svg.append("g")
                .attr("class", "brush")
                .call(brush);
        }
        
        // Update the visualization for the selected genome
        function updateVisualization() {
            initVisualization();

            console.log("Updating visualization for genome:", currentGenome.value);
            
            // Access currentGenome through the global reactive reference
            const currentGenomeValue = currentGenome.value;
            if (!currentGenomeValue) return;
            
            const taxid = currentGenomeValue.taxid;
            const genomeFragments = coverageData.filter(d => d.genome_taxid === taxid);
            
            if (genomeFragments.length === 0) {
                d3.select("#coverage-plot").html("<p>No coverage data available for this genome.</p>");
                return;
            }
            
            // Sort fragments by sequence name and position
            genomeFragments.sort((a, b) => {
                if (a.seq_name === b.seq_name) {
                    return a.start_position - b.start_position;
                }
                return a.seq_name.localeCompare(b.seq_name);
            });
            
            // Set x domain based on total fragments
            xScale.domain([0, currentGenomeValue.totalLength]);
            
            // Set y domain based on maximum mean depth with some padding
            const maxMeanDepth = d3.max(genomeFragments, d => d.meandepth);
            yScale.domain([0, maxMeanDepth * 1.2]); // Add 20% padding
            
            // Create percentage scale for top axis (0-100%)
            const xScalePercent = d3.scaleLinear()
                .domain([0, 100])
                .range([0, width]);
            
            // Update axes
            xAxis.call(d3.axisBottom(xScale));
            xAxisTop.call(d3.axisTop(xScalePercent).tickFormat(d => d + "%"));
            yAxis.call(d3.axisLeft(yScale));
            
            // Create color scale for different sequences
            const seqNames = [...new Set(genomeFragments.map(d => d.seq_name))];
            const colorScale = d3.scaleOrdinal(d3.schemeCategory10)
                .domain(seqNames);
            
            // Create legend
            const legend = d3.select("#legend").html("");
            seqNames.forEach(seq => {
                const item = legend.append("div")
                    .attr("class", "legend-item");
                
                item.append("div")
                    .attr("class", "legend-color")
                    .style("background-color", colorScale(seq));
                
                item.append("span")
                    .text(seq);
            });
            
            // Calculate positions - concatenate fragments without gaps
            let position = 0;
            genomeFragments.forEach(fragment => {
                fragment.x_start = position;
                const fragmentLength = fragment.end_position - fragment.start_position + 1;
                fragment.x_end = position + fragmentLength;
                position = fragment.x_end; // Concatenate signature fragments without gaps
            });
            
            // Draw fragments - uncovered portion in gray, covered portion in color
            // First draw the entire fragment in gray
            svg.selectAll(".fragment-background")
                .data(genomeFragments)
                .enter()
                .append("rect")
                .attr("class", "fragment-background")
                .attr("x", d => xScale(d.x_start))
                .attr("y", d => yScale(d.meandepth))
                .attr("width", d => Math.max(2, xScale(d.x_end) - xScale(d.x_start)))
                .attr("height", d => height - yScale(d.meandepth))
                .attr("fill", "#eee")
                .attr("rx", 7)  // Add rounded corners with 7px radius
                .attr("ry", 7); // Add rounded corners with 7px radius
            
            // Add signature fragment indicators on top of histogram
            svg.selectAll(".fragment-indicator")
                .data(genomeFragments)
                .enter()
                .append("rect")
                .attr("class", "fragment-indicator")
                .attr("x", d => xScale(d.x_start))
                .attr("y", 0) // Position at the top below the axis
                .attr("width", d => Math.max(2, xScale(d.x_end) - xScale(d.x_start)))
                .attr("height", 15)
                .attr("fill", d => colorScale(d.seq_name))
                .attr("opacity", d => d.coverage / 100) // Opacity based on coverage percentage
                .attr("rx", 4)  // Add rounded corners
                .attr("ry", 4)
                .on("mouseover", function(event, d) {
                    tooltip.transition()
                        .duration(100)
                        .style("opacity", .95);
                    tooltip.html(`
                        <strong>Sequence:</strong> ${d.seq_name}<br>
                        <strong>Position:</strong> ${d.start_position}-${d.end_position}<br>
                        <strong>Coverage:</strong> ${d.coverage.toFixed(2)}%<br>
                        <strong>Read count:</strong> ${d.numreads.toFixed(0)}<br>
                        <strong>Mean Depth:</strong> ${d.meandepth.toFixed(2)}x<br>
                        <small class="text-muted">Click to view NCBI annotations</small>
                    `)
                        .style("left", (event.pageX + 10) + "px")
                        .style("top", (event.pageY - 28) + "px");
                })
                .on("mouseout", function() {
                    tooltip.transition()
                        .duration(500)
                        .style("opacity", 0);
                })
                .on("click", function(event, d) {
                    // Store current fragment data for NCBI query in Vue's reactive state
                    currentFragment.value = d;
                    
                    // Show the modal with loading indicator
                    document.getElementById('annotation-results').innerHTML = '';
                    document.getElementById('annotation-loading').classList.remove('d-none');
                    annotationModal.show();
                    
                    // Set modal title with fragment info
                    document.getElementById('annotationModalLabel').textContent = 
                        `NCBI Annotations: ${d.seq_name} (${d.start_position}-${d.end_position})`;
                    
                    // Query NCBI for annotations
                    queryNucleotide(d.seq_name, d.start_position, d.end_position);
                });
            
            // Then draw the covered portion in color
            svg.selectAll(".fragment")
                .data(genomeFragments)
                .enter()
                .append("rect")
                .attr("class", "fragment")
                .attr("x", d => xScale(d.x_start))
                .attr("y", d => yScale(d.meandepth))
                .attr("width", d => Math.max(2, xScale(d.x_end) - xScale(d.x_start)) * (d.coverage / 100)) // Width based on coverage percentage
                .attr("height", d => height - yScale(d.meandepth))
                .attr("fill", d => colorScale(d.seq_name))
                .attr("rx", 3)  // Add rounded corners with 3px radius
                .attr("ry", 3)
                .on("mouseover", function(event, d) {
                    tooltip.transition()
                        .duration(100)
                        .style("opacity", .95);
                    tooltip.html(`
                        <strong>Sequence:</strong> ${d.seq_name}<br>
                        <strong>Position:</strong> ${d.start_position}-${d.end_position}<br>
                        <strong>Coverage:</strong> ${d.coverage.toFixed(2)}%<br>
                        <strong>Read count:</strong> ${d.numreads.toFixed(0)}<br>
                        <strong>Mean Depth:</strong> ${d.meandepth.toFixed(2)}x<br>
                        <small class="text-muted">Click to view NCBI annotations</small>
                    `)
                        .style("left", (event.pageX + 10) + "px")
                        .style("top", (event.pageY - 28) + "px");
                })
                .on("mouseout", function() {
                    tooltip.transition()
                        .duration(500)
                        .style("opacity", 0);
                })
                .on("click", function(event, d) {
                    // Store current fragment data for NCBI query in Vue's reactive state
                    currentFragment.value = d;
                    
                    // Show the modal with loading indicator
                    document.getElementById('annotation-results').innerHTML = '';
                    document.getElementById('annotation-loading').classList.remove('d-none');
                    annotationModal.show();
                    
                    // Set modal title with fragment info
                    document.getElementById('annotationModalLabel').textContent = 
                        `NCBI Annotations: ${d.seq_name} (${d.start_position}-${d.end_position})`;
                    
                    // Query NCBI for annotations
                    queryNucleotide(d.seq_name, d.start_position, d.end_position);
                });

            
            // Draw mean depth line
            const overallMeanDepth = d3.mean(genomeFragments, d => d.meandepth);
            svg.append("line")
                .attr("x1", 0)
                .attr("y1", yScale(overallMeanDepth))
                .attr("x2", width)
                .attr("y2", yScale(overallMeanDepth))
                .attr("stroke", "red")
                .attr("stroke-width", 3)
                .attr("stroke-opacity", 0.5)
                .attr("stroke-dasharray", "5,5")
                .style("cursor", "pointer")
                .on("mouseover", function(event) {
                    tooltip.transition()
                        .duration(200)
                        .style("opacity", .9);
                    tooltip.html(`<strong>Overall Mean Depth:</strong> ${overallMeanDepth.toFixed(2)}x`)
                        .style("left", (event.pageX + 10) + "px")
                        .style("top", (event.pageY - 28) + "px");
                })
                .on("mouseout", function() {
                    tooltip.transition()
                        .duration(500)
                        .style("opacity", 0);
                });
        }
        
        // NCBI Nucleotide Annotation Query function
        async function queryNucleotide(accession, startPosition, stopPosition) {
            const resultsDiv = document.getElementById('annotation-results');
            const loadingDiv = document.getElementById('annotation-loading');
            
            try {
                // Check if accession is a generic sequence descriptor, not a real accession
                const genericTerms = ["chromosome", "scaffold", "contig", "node", "chr"];
                const accessionLower = accession.toLowerCase();
                
                for (const term of genericTerms) {
                    if (accessionLower.startsWith(term)) {
                        resultsDiv.innerHTML = `<div class="alert alert-info">
                            <i class="bi bi-info-circle me-2"></i>
                            Sequence "${accession}" appears to be a generic sequence identifier rather than an NCBI accession.
                            NCBI annotations are not available for this sequence.
                        </div>`;
                        loadingDiv.classList.add('d-none');
                        return;
                    }
                }
                
                // Set up NCBI API URLs
                const baseUrl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
                const encodedAccession = encodeURIComponent(accession);

                // Step 1: Search for the accession to get the UID in the 'nuccore' database
                const esearchUrl = `${baseUrl}esearch.fcgi?db=nuccore&term=${encodedAccession}&retmode=json`;
                console.log('ESearch URL:', esearchUrl);
                const esearchResponse = await fetch(esearchUrl);
                const esearchData = await esearchResponse.json();

                if (esearchData.esearchresult.idlist && esearchData.esearchresult.idlist.length > 0) {
                    const uid = esearchData.esearchresult.idlist[0];

                    // Step 2: Fetch the nucleotide record in GenBank format
                    const efetchUrl = `${baseUrl}efetch.fcgi?db=nuccore&id=${encodeURIComponent(uid)}&seq_start=${encodeURIComponent(startPosition)}&seq_stop=${encodeURIComponent(stopPosition)}&rettype=gb&retmode=text`;
                    console.log('EFetch URL:', efetchUrl);
                    const efetchResponse = await fetch(efetchUrl);
                    const gbText = await efetchResponse.text();
                    
                    // Step 3: Parse the GenBank format output
                    const lines = gbText.split('\n');
                    let sequenceName = '';
                    let organismName = '';
                    let definition = '';
                    let plasmid = 'No';
                    const annotations = [];
                    let currentFeature = null;
                    let inFeatures = false;
                    
                    // Estimate sequence length from the query range
                    const sequenceLength = stopPosition - startPosition + 1;
                    
                    // Parse the GenBank header information and features
                    for (let i = 0; i < lines.length; i++) {
                        const line = lines[i].trim();
                        
                        // Parse header information
                        if (line.startsWith('LOCUS')) {
                            // Parse LOCUS line for sequence name
                            const parts = line.split(/\s+/);
                            if (parts.length > 1) {
                                sequenceName = parts[1];
                            }
                        }
                        else if (line.startsWith('DEFINITION')) {
                            // Get sequence definition
                            definition = line.substring(10).trim();
                            // Continue to collect multi-line definitions
                            let j = i + 1;
                            while (j < lines.length && !lines[j].trim().startsWith('ACCESSION') && lines[j].trim().length > 0) {
                                definition += ' ' + lines[j].trim();
                                j++;
                            }

                            plasmid = definition.toLowerCase().includes('plasmid') ? 'Yes' : 'No';
                        }
                        else if (line.startsWith('ORGANISM')) {
                            // Get organism name
                            organismName = line.substring(10).trim();
                        } 
                        // Start of features section
                        else if (line === 'FEATURES             Location/Qualifiers') {
                            inFeatures = true;
                            continue;
                        } 
                        // End of features section
                        else if (line.startsWith('ORIGIN')) {
                            inFeatures = false;
                            break;
                        }
                        // Parse features
                        else if (inFeatures) {
                            // Feature lines start with 5 spaces then the feature type (like "gene", "CDS", etc.)
                            if (line.match(/^[a-zA-Z]/) && !line.includes('Location/Qualifiers')) {
                                // Save the previous feature if it exists
                                if (currentFeature) {
                                    annotations.push(currentFeature);
                                }
                                
                                // Get feature type (first word)
                                const featureParts = line.split(/\s+/);
                                
                                // Create new feature
                                currentFeature = {
                                    type: featureParts[0],
                                    location: featureParts[1],
                                    qualifiers: {}
                                };
                            }
                            // Qualifier lines start with 21 spaces then /qualifier=
                            else if (line.match(/^\//) && currentFeature) {
                                // Extract qualifier name and value
                                const qualParts = line.split('=');
                                
                                if (qualParts.length > 1) {
                                    // Handle qualifier with value like /gene="lacZ"
                                    const qualifierName = qualParts[0].replace('/', ''); // Remove leading '/'
                                    let qualifierValue = qualParts[1];

                                    // Handle quoted values
                                    if (qualifierValue.startsWith('"')) {
                                        qualifierValue = qualifierValue.substring(1);
                                        
                                        // Check if the value ends with quote or continues on next lines
                                        if (qualifierValue.endsWith('"')) {
                                            qualifierValue = qualifierValue.substring(0, qualifierValue.length - 1);
                                        } else {
                                            // Multi-line qualifier value
                                            let j = i + 1;
                                            while (j < lines.length) {
                                                const nextLine = lines[j].trim();
                                                if (nextLine.endsWith('"')) {
                                                    qualifierValue += ' ' + nextLine.substring(0, nextLine.length - 1);
                                                    i = j; // Skip these lines in the main loop
                                                    break;
                                                } else {
                                                    qualifierValue += ' ' + nextLine;
                                                }
                                                j++;
                                            }
                                        }
                                    }
                                    
                                    // Add to feature's qualifiers. Skip if the qualifier name is translation
                                    if (qualifierName !== 'translation') {
                                        currentFeature.qualifiers[qualifierName] = qualifierValue;
                                    }
                                }
                            }
                        }
                    }
                    
                    // Add the last feature if it exists
                    if (currentFeature) {
                        annotations.push(currentFeature);
                    }

                    // Generate result HTML - sequence details in clean format
                    let output = `<div class="card mb-3 shadow-sm">
                        <div class="card-header text-white" style="background-color: #a7aca2;">
                            Sequence Details
                        </div>
                        <div class="card-body">
                            <div class="row">
                                <div class="col-md-6">
                                    <p><strong>Accession:</strong> ${accession}</p>
                                    <p><strong>Name:</strong> ${sequenceName || 'N/A'}</p>
                                    <p><strong>Organism:</strong> ${organismName || 'N/A'}</p>
                                </div>
                                <div class="col-md-6">
                                    <p><strong>Length:</strong> ${sequenceLength.toLocaleString()} bp</p>
                                    <p><strong>Location:</strong> ${startPosition.toLocaleString()} - ${stopPosition.toLocaleString()}</p>
                                    <p><strong>Plasmid:</strong> ${plasmid || 'N/A'}</p>
                                </div>
                            </div>
                            <div class="row mt-2">
                                <div class="col-12">
                                    <p><strong>Definition:</strong> ${definition || 'N/A'}</p>
                                </div>
                            </div>
                        </div>
                    </div>`;

                    if (annotations.length > 0) {
                        // Group annotations by type for better organization
                        const groupedAnnotations = {};
                        annotations.forEach(anno => {
                            if (!groupedAnnotations[anno.type]) {
                                groupedAnnotations[anno.type] = [];
                            }
                            groupedAnnotations[anno.type].push(anno);
                        });
                        
                        output += `<div class="d-flex justify-content-between align-items-center mb-3">
                            <h5 class="mb-0">Annotations (${annotations.length})</h5>
                            <span class="badge" style="background-color: #b5a397; color: white">${startPosition.toLocaleString()}-${stopPosition.toLocaleString()} bp</span>
                        </div>`;
                        
                        // Create tabs for feature types
                        const featureTypes = Object.keys(groupedAnnotations);
                        
                        if (featureTypes.length > 1) {
                            // Add tab navigation
                            output += `<ul class="nav nav-tabs mb-3" id="annotationTabs" role="tablist">`;
                            featureTypes.forEach((type, index) => {
                                const isActive = index === 0 ? 'active' : '';
                                output += `
                                <li class="nav-item" role="presentation">
                                    <button class="nav-link ${isActive}" id="tab-${type}" data-bs-toggle="tab" 
                                    data-bs-target="#content-${type}" type="button" role="tab" aria-selected="${index === 0}">
                                        ${type} (${groupedAnnotations[type].length})
                                    </button>
                                </li>`;
                            });
                            output += `</ul>`;
                            
                            // Add tab content
                            output += `<div class="tab-content" id="annotationTabContent">`;
                            featureTypes.forEach((type, index) => {
                                const isActive = index === 0 ? 'show active' : '';
                                output += `
                                <div class="tab-pane fade ${isActive}" id="content-${type}" role="tabpanel" aria-labelledby="tab-${type}">`;
                                
                                // Feature cards in each tab
                                groupedAnnotations[type].forEach(anno => {
                                    output += createAnnotationCard(anno);
                                });
                                
                                output += `</div>`;
                            });
                            output += `</div>`;
                        } else {
                            // If only one feature type, no need for tabs
                            const type = featureTypes[0];
                            const annos = groupedAnnotations[type];
                            
                            annos.forEach(anno => {
                                output += createAnnotationCard(anno);
                            });
                        }
                    } else {
                        output += `<div class="alert alert-info mt-3">
                            <i class="bi bi-info-circle me-2"></i>
                            No specific annotations found overlapping the range ${startPosition.toLocaleString()}-${stopPosition.toLocaleString()}.
                        </div>`;
                    }
                    
                    // Helper function to create an annotation card
                    function createAnnotationCard(anno) {
                        let result = `
                        <div class="card mb-3 shadow-sm annotation-card">
                            <div class="card-body">
                                <div class="d-flex justify-content-between align-items-start">
                                    <h5 class="text-dark">${anno.type}</h5>
                                    <span class="badge bg-light text-dark">${anno.location || anno.start + '..' + anno.end}</span>
                                </div>`;
                        
                        // Display important qualifiers at the top
                        const priorityQualifiers = ['gene', 'product', 'note', 'function', 'protein_id'];
                        
                        if (Object.keys(anno.qualifiers).length > 0) {
                            // First show priority qualifiers
                            result += `<div class="annotation-qualifiers mt-2">`;
                            
                            // Show priority qualifiers first
                            priorityQualifiers.forEach(key => {
                                if (anno.qualifiers[key]) {
                                    const value = anno.qualifiers[key];
                                    if (value === true) {
                                        result += `<div class="mb-1"><span class="badge bg-secondary">${key}</span></div>`;
                                    } else {
                                        result += `<div class="mb-1"><strong>${key}:</strong> ${value}</div>`;
                                    }
                                }
                            });
                            
                            // Then show all other qualifiers
                            for (const [key, value] of Object.entries(anno.qualifiers)) {
                                // Skip already shown priority qualifiers
                                if (priorityQualifiers.includes(key)) continue;
                                
                                if (value === true) {
                                    result += `<div class="mb-1"><span class="badge bg-secondary">${key}</span></div>`;
                                } else {
                                    result += `<div class="mb-1"><strong>${key}:</strong> ${value}</div>`;
                                }
                            }
                            
                            result += `</div>`;
                        }
                        
                        result += `</div></div>`;
                        return result;
                    }

                    // Update results
                    resultsDiv.innerHTML = output;
                } else {
                    resultsDiv.innerHTML = `<div class="alert alert-danger">
                        Nucleotide sequence with accession "${accession}" not found.
                    </div>`;
                }
            } catch (error) {
                console.error('Error fetching data:', error);
                resultsDiv.innerHTML = `<div class="alert alert-danger">
                    An error occurred: ${error.message}. Please check the accession number and your internet connection.
                </div>`;
            } finally {
                loadingDiv.classList.add('d-none');
            }
        }
        
        // Update genome information panel
        function updateGenomeInfo() {
            if (!currentGenome) return;
            
            const infoDiv = document.getElementById('genome-info');
            infoDiv.innerHTML = `
                <div class="row">
                    <div class="col-md-6">
                        <p><strong>Name:</strong> ${currentGenome.name}</p>
                        <p><strong>Taxid:</strong> ${currentGenome.taxid}</p>
                        <p><strong>Domain:</strong> ${currentGenome.superkingdom}</p>
                    </div>
                    <div class="col-md-6">
                        <p><strong>Number of Sequences:</strong> ${currentGenome.numOfSeq}</p>
                        <p><strong>Total Length:</strong> ${currentGenome.totalLength.toLocaleString()}</p>
                        <p><strong>Genome Size:</strong> ${currentGenome.genomeSize.toLocaleString()}</p>
                    </div>
                </div>
            `;
        }
        
        // Handle brush zoom
        function brushended(event) {
            if (!event.selection) return;
            const [x0, x1] = event.selection.map(xScale.invert);
            xScale.domain([x0, x1]);
            
            // Remove brush
            svg.select(".brush").call(brush.move, null);
            
            // Update chart with new domain
            updateChart();
        }
        
        // Reset zoom
        function resetZoom() {
            const currentGenomeValue = currentGenome.value;
            if (!currentGenomeValue) return;
            xScale.domain([0, currentGenomeValue.totalLength]);
            updateChart();
        }
        
        // Update chart after zoom
        function updateChart() {
            // Update x-axis
            svg.select(".x-axis")
                .transition()
                .duration(750)
                .call(d3.axisBottom(xScale));
            
            // Update top x-axis (percentage scale needs to be recalculated based on current zoom)
            const currentDomain = xScale.domain();
            const totalLength = currentGenome.value.totalLength;
            const startPercent = (currentDomain[0] / totalLength) * 100;
            const endPercent = (currentDomain[1] / totalLength) * 100;
            
            const xScalePercentZoomed = d3.scaleLinear()
                .domain([startPercent, endPercent])
                .range([0, width]);
            
            svg.select(".x-axis-top")
                .transition()
                .duration(750)
                .call(d3.axisTop(xScalePercentZoomed).tickFormat(d => d.toFixed(1) + "%"));
            
            // Update fragments - all components including indicators
            svg.selectAll(".fragment-background")
                .transition()
                .duration(750)
                .attr("x", d => xScale(d.x_start))
                .attr("width", d => Math.max(2, xScale(d.x_end) - xScale(d.x_start)));
            
            svg.selectAll(".fragment-indicator")
                .transition()
                .duration(750)
                .attr("x", d => xScale(d.x_start))
                .attr("width", d => Math.max(2, xScale(d.x_end) - xScale(d.x_start)));
                
            svg.selectAll(".fragment")
                .transition()
                .duration(750)
                .attr("x", d => xScale(d.x_start))
                .attr("width", d => Math.max(2, xScale(d.x_end) - xScale(d.x_start)) * (d.coverage / 100));
        }
    </script>
</body>
</html>
"""

def parse_coverage_file(coverage_file):
    """Parse the coverage file and extract relevant information using pandas."""
    # Read the coverage file with pandas
    # Define column names based on the file format
    col_names = ['rname', 'startpos', 'endpos', 'numreads', 'covbases', 'coverage', 'meandepth', 'meanbaseq', 'meanmapq']
    df = pd.read_csv(coverage_file, sep='\t', comment='#', names=col_names, header=None)
    
    # Parse the rname field which has format: [seq_name]|[start_position]|[stop_position]|[genom_taxid]
    df[['seq_name', 'start_position', 'end_position', 'genome_taxid']] = df['rname'].str.rstrip('|').str.split('|', expand=True)
    
    # Convert data types for numerical columns
    numeric_cols = ['startpos', 'endpos', 'numreads', 'covbases', 'coverage', 'meandepth', 'meanbaseq']
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric)
    
    # Handle meanmapq - set to 0 if it's not a valid number
    df['meanmapq'] = pd.to_numeric(df['meanmapq'], errors='coerce').fillna(0)
    
    # Remove rows with missing data in essential columns
    df = df.dropna(subset=['seq_name', 'start_position', 'end_position', 'genome_taxid'])
    
    # Convert DataFrame to list of dictionaries for compatibility with the existing code
    coverage_data = df.to_dict('records')
    
    return coverage_data, df['genome_taxid'].unique().tolist()

def _stats_has_header(stats_file: str) -> bool:
    with open(stats_file, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            first = line.rstrip("\n").split("\t", 1)[0].strip().lower()
            return first in {"db_level", "db level", "level"}
    return False

def _read_stats(stats_file: str) -> pd.DataFrame:
    has_header = _stats_has_header(stats_file)
    STATS_COLUMNS = ['DB_level', 'Name', 'Taxid', 'Domain', 'NumOfSeq', 'Max', 'Min', 'TotalLength', 'GenomeSize', 'Note']
    df = pd.read_csv(
        stats_file,
        sep="\t",
        comment="#",
        names=STATS_COLUMNS,
        header=0 if has_header else None,
        low_memory=False,
        dtype={"Taxid": "string"},
    )

    if "Note" not in df:
        df["Note"] = ""

    return df

def parse_stats_file(stats_file, uniq_taxid_list):
    df = _read_stats(stats_file)
    uniq_taxids = {str(t) for t in uniq_taxid_list}
    df["Taxid"] = df["Taxid"].astype(str)
    df = df[df["Taxid"].isin(uniq_taxids)]
    
    # Convert to dictionary format for the JavaScript
    genome_data = []
    
    for _, row in df.iterrows():
        genome_data.append({
            'db_level': row['DB_level'],
            'name': row['Name'],  # This will be the full name as read by pandas
            'taxid': row['Taxid'],
            'superkingdom': row['Domain'],
            'numOfSeq': int(row['NumOfSeq']),
            'max': int(row['Max']),
            'min': int(row['Min']),
            'totalLength': int(row['TotalLength']),
            'genomeSize': int(row['GenomeSize']) if 'GenomeSize' in df.columns and not pd.isna(row['GenomeSize']) else 0
        })
    
    return genome_data

def generate_html(coverage_data, genome_data, output_file, use_external):
    """Generate the HTML file with the visualization."""
    # Convert data to JSON for JavaScript
    coverage_json = json.dumps(coverage_data)
    genome_json = json.dumps(genome_data)

    if use_external:
        # Use CDN links for external resources
        html_content = HTML_TEMPLATE.replace(
            '<script src="/publicdata/js/d3.v7.min.js"></script>',
            '<script src="https://d3js.org/d3.v7.min.js"></script>'
        ).replace(
            '<script src="/publicdata/js/vue.global.prod.js"></script>',
            '<script src="https://unpkg.com/vue@3/dist/vue.global.prod.js"></script>'
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

    # Replace placeholders in the HTML template
    html_content = html_content.replace('GENOME_DATA_PLACEHOLDER', genome_json)
    html_content = html_content.replace('COVERAGE_DATA_PLACEHOLDER', coverage_json)
    
    minified_html = minify_html.minify(html_content, keep_comments=False, minify_css=True, minify_js=True, remove_processing_instructions=True)

    # Write the HTML file
    with open(output_file, 'w') as f:
        f.write(minified_html)
    
def main():
    """Main function to process files and generate visualization."""
    parser = argparse.ArgumentParser(description='Generate an HTML-based genome coverage visualization')
    parser.add_argument('-c', '--coverage', required=True, help='Path to the coverage file')
    parser.add_argument('-s', '--stats', required=True, help='Path to the stats file')
    parser.add_argument('-o', '--output', default='coverage_visualization.html', 
                       help='Path to the output HTML file (default: coverage_visualization.html)')
    parser.add_argument('-e', '--external', action='store_true',
                       help='Use external resources for the HTML visualization')
    
    args = parser.parse_args()
    
    # Process files
    coverage_data, uniq_taxid_list = parse_coverage_file(args.coverage)
    genome_data = parse_stats_file(args.stats, uniq_taxid_list)
    
    # Generate HTML
    generate_html(coverage_data, genome_data, args.output, args.external)

    print(f"Visualization generated: {os.path.abspath(args.output)}")

if __name__ == "__main__":
    main()
