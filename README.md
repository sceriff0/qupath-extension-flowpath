# QuPath Extension: Gate Tree

[![QuPath](https://img.shields.io/badge/QuPath-%E2%89%A50.7.0-blue.svg)](https://qupath.github.io/)
[![Java](https://img.shields.io/badge/Java-25-orange.svg)](https://jdk.java.net/25/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Interactive tree-based cell phenotyping gating for [QuPath](https://qupath.github.io/). Build hierarchical marker gates with live histogram visualization, drag thresholds, and see cells update in real time on multiplexed imaging data (CODEX, MIBI, mIF).

Designed to work with the [mirage](https://github.com/sceriff0/mirage) pipeline for end-to-end multiplexed image analysis — from raw images to cell phenotypes.

![Screenshot](screenshot.png)

## Features

- **Hierarchical gating** — Multi-level gate tree (e.g., CD45+ > CD3+ > CD8+ = "T cytotoxic")
- **Live histogram** — Per-marker histogram with draggable threshold and percentile clipping
- **Real-time preview** — Cells update color and phenotype within ~100ms
- **Raw / Z-score toggle** — Switch value modes per gate
- **Quality filters** — Pre-gating QC: area, eccentricity, solidity, total intensity
- **Outlier exclusion** — Per-gate percentile clipping with cell exclusion
- **Save / Load / Export** — Gate trees as JSON, phenotypes as CSV

## Installation

### From QuPath Extension Manager

1. Open QuPath → `Extensions` → `Manage extensions` → `Manage extension catalogs` → `Add`
2. Paste: `https://github.com/sceriff0/qupath-extension-gatetree`
3. Go back to `Manage extensions` → find **Gate Tree** → click `+` to install
4. Restart QuPath

### From Release JAR

Download the latest JAR from [Releases](../../releases) and drag it onto QuPath.

### Build from Source

```bash
git clone https://github.com/sceriff0/qupath-extension-gatetree.git
cd qupath-extension-gatetree
./gradlew build
# JAR at build/libs/qupath-extension-gatetree-0.2.0.jar → drag onto QuPath
```

## Quick Start

1. Open your pyramidal OME-TIFF in QuPath
2. Import cell detections (GeoJSON with marker intensities, e.g., from [mirage](https://github.com/sceriff0/mirage))
3. `Extensions` > `Gate Tree` (or `Ctrl+G`)
4. Set quality filters to remove segmentation artifacts
5. Add root gate → select channel → drag threshold on histogram
6. Build hierarchy: "Add Gate to +" for sub-gating positive populations
7. Name leaf nodes (e.g., "T cytotoxic", "Stroma")
8. Export CSV

## Output Formats

**Gate Tree JSON** (`gate_tree.json`) — Saves the full gate hierarchy, thresholds, colors, and quality filter settings. Load to reproduce gating.

**Phenotype CSV** (`gate_pheno.csv`) — One row per cell with phenotype name and per-marker +/- status:

```csv
cell_id,phenotype,CD45,CD3,CD8,PANCK
0,T cytotoxic,+,+,+,-
1,Tumor,-,-,-,+
```

Excluded cells (QC-filtered or outlier-excluded) are omitted when "Exclude from CSV" is enabled.

## Citation

If you use this tool in your research, please cite:

> Gate Tree: Interactive tree-based cell phenotyping gating for QuPath. (2025). https://github.com/sceriff0/qupath-extension-gatetree

## License

MIT License. See [LICENSE](LICENSE).
