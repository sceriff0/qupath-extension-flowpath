# QuPath Extension: FlowPath

[![QuPath](https://img.shields.io/badge/QuPath-%E2%89%A50.7.0-blue.svg)](https://qupath.github.io/)
[![Java](https://img.shields.io/badge/Java-25-orange.svg)](https://jdk.java.net/25/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Interactive tree-based cell phenotyping for [QuPath](https://qupath.github.io/). Build hierarchical marker gates with live visualization, drag thresholds, draw 2D regions, and see cells update in real time on multiplexed imaging data (CODEX, MIBI, mIF).

Designed to work with the [mirage](https://github.com/sceriff0/mirage) pipeline for end-to-end multiplexed image analysis — from raw images to cell phenotypes.

<!-- Add screenshots here -->

## Features

- **Hierarchical gating** — Multi-level gate tree (e.g., CD45+ > CD3+ > CD8+ = "T cytotoxic")
- **Multiple gate types** — Threshold (1D), quadrant (2D dual-threshold), polygon, rectangle, and ellipse gates
- **Live histogram & scatter plots** — Per-marker histogram with draggable threshold; 2D scatter plots with interactive gate drawing and crosshair visualization for quadrant gates
- **Real-time preview** — Cells update color and phenotype within ~100ms
- **Raw / Z-score toggle** — Switch value modes per gate
- **Quality filters** — Pre-gating QC with min+max for area, eccentricity, solidity, total intensity, and perimeter
- **Outlier exclusion** — Per-gate percentile clipping with scatter plot axis zoom to clipped range
- **Undo / Redo** — Snapshot-based undo stack (Ctrl+Z / Ctrl+Shift+Z)
- **Drag-and-drop reordering** — Rearrange gates in the tree by dragging between branches
- **Save / Load / Export** — Gate trees as JSON (backward-compatible), phenotypes as CSV

## Installation

### From QuPath Extension Manager

1. Open QuPath → `Extensions` → `Manage extensions` → `Manage extension catalogs` → `Add`
2. Paste: `https://raw.githubusercontent.com/sceriff0/flowpath-catalog/main/catalog.json`
3. Go back to `Manage extensions` → find **FlowPath** → click `+` to install
4. Restart QuPath

### From Release JAR

Download the latest JAR from [Releases](../../releases) and drag it onto QuPath.

### Build from Source

```bash
git clone https://github.com/sceriff0/qupath-extension-flowpath.git
cd qupath-extension-flowpath
./gradlew build
# JAR at build/libs/FlowPath-1.0.0.jar → drag onto QuPath
```

## Quick Start

1. Open your pyramidal OME-TIFF in QuPath
2. Import cell detections (GeoJSON with marker intensities, e.g., from [mirage](https://github.com/sceriff0/mirage))
3. `Extensions` > `FlowPath` (or `Ctrl+G`)
4. Set quality filters to remove segmentation artifacts (area, eccentricity, solidity, perimeter, total intensity)
5. Add root gate → pick gate type (threshold, quadrant, polygon, rectangle, ellipse)
6. For threshold gates: select channel → drag threshold on histogram
7. For 2D gates: select X/Y channels → draw region on scatter plot
8. Build hierarchy: add child gates to branches for sub-gating
9. Name leaf nodes (e.g., "T cytotoxic", "Stroma")
10. Export CSV

## Output Formats

**FlowPath JSON** (`flowpath.json`) — Saves the full gate hierarchy, thresholds, colors, and quality filter settings. Load to reproduce gating.

**Phenotype CSV** (`gate_pheno.csv`) — One row per cell with phenotype name and per-marker +/- status:

```csv
cell_id,phenotype,CD45,CD3,CD8,PANCK
0,T cytotoxic,+,+,+,-
1,Tumor,-,-,-,+
```

Excluded cells (QC-filtered or outlier-excluded) are omitted when "Exclude from CSV" is enabled.

## Citation

If you use this tool in your research, please cite:

> FlowPath: Interactive tree-based cell phenotyping for QuPath. (2026). https://github.com/sceriff0/qupath-extension-flowpath

## License

MIT License. See [LICENSE](LICENSE).
