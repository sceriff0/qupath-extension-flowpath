package qupath.ext.quflow.io;

import qupath.ext.quflow.engine.GatingEngine;
import qupath.ext.quflow.model.CellIndex;
import qupath.ext.quflow.model.GateNode;
import qupath.ext.quflow.model.GateTree;
import qupath.ext.quflow.model.MarkerStats;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Exports cell phenotype assignments to a CSV file.
 * <p>
 * Each row represents a single cell. Columns include identity, spatial coordinates,
 * geometry measurements, and per-marker triplets (raw intensity, z-score, gating sign).
 */
public class PhenotypeCsvExporter {

    private PhenotypeCsvExporter() {
        // static utility class
    }

    /**
     * Export phenotype assignments to CSV with raw intensities, z-scores, and geometry.
     *
     * @param file   destination CSV file
     * @param index  the cell index containing objects and marker data
     * @param result the gating assignment result
     * @param tree   the gate tree used for the assignment
     * @param stats  per-marker statistics for z-score computation (may be null)
     * @throws IOException if writing fails
     */
    public static void export(File file, CellIndex index, GatingEngine.AssignmentResult result,
                              GateTree tree, MarkerStats stats) throws IOException {

        // 1. Collect all unique marker channels in the gate tree (preserving order)
        List<String> markerColumns = collectMarkerChannels(tree);

        // 2. Build a lookup from phenotype name -> marker sign map
        Map<String, Map<String, String>> phenotypeMarkerSigns = new LinkedHashMap<>();
        for (GateNode root : tree.getRoots()) {
            traceMarkerSigns(root, new LinkedHashMap<>(), phenotypeMarkerSigns);
        }

        boolean excludeFiltered = tree.getQualityFilter().isExcludeFromCsv();
        String[] phenotypes = result.getPhenotypes();
        boolean[] excluded = result.getExcluded();

        try (BufferedWriter writer = new BufferedWriter(new FileWriter(file))) {
            // Write header
            writer.write("cell_id,phenotype,centroid_x,centroid_y,area,perimeter,eccentricity,solidity");
            for (String marker : markerColumns) {
                String safe = escapeCsv(marker);
                writer.write("," + safe + "_raw");
                writer.write("," + safe + "_zscore");
                writer.write("," + safe + "_sign");
            }
            writer.newLine();

            // Write one row per non-excluded cell
            int n = index.getSize();
            for (int i = 0; i < n; i++) {
                if (excludeFiltered && excluded[i]) {
                    continue;
                }

                String phenotype = phenotypes[i] != null ? phenotypes[i] : "";

                // Identity
                writer.write(String.valueOf(i));
                writer.write(',');
                writer.write(escapeCsv(phenotype));

                // Spatial + geometry (empty if NaN)
                writer.write(',' + fmt(index.getCentroidX(i)));
                writer.write(',' + fmt(index.getCentroidY(i)));
                writer.write(',' + fmt(index.getArea(i)));
                writer.write(',' + fmt(index.getPerimeter(i)));
                writer.write(',' + fmt(index.getEccentricity(i)));
                writer.write(',' + fmt(index.getSolidity(i)));

                // Per marker: raw, zscore, sign
                Map<String, String> signs = phenotypeMarkerSigns.get(phenotype);

                for (String marker : markerColumns) {
                    int mIdx = index.getMarkerIndex(marker);
                    double raw = (mIdx >= 0) ? index.getMarkerValues(mIdx)[i] : Double.NaN;

                    double zscore;
                    if (Double.isNaN(raw) || stats == null || stats.getStd(marker) <= 1e-10) {
                        zscore = Double.NaN;
                    } else {
                        zscore = stats.toZScore(marker, raw);
                    }

                    String sign = (signs != null && signs.containsKey(marker)) ? signs.get(marker) : "";

                    writer.write(',' + fmt(raw));
                    writer.write(',' + fmt(zscore));
                    writer.write(',' + sign);
                }
                writer.newLine();
            }
        }
    }

    /** Format a double for CSV; NaN → empty string. */
    private static String fmt(double val) {
        return Double.isNaN(val) ? "" : String.format("%.4f", val);
    }

    /**
     * Collect all unique marker channel names used across the gate tree,
     * preserving the order of first encounter (depth-first).
     */
    private static List<String> collectMarkerChannels(GateTree tree) {
        Set<String> seen = new LinkedHashSet<>();
        for (GateNode root : tree.getRoots()) {
            collectChannelsRecursive(root, seen);
        }
        return new ArrayList<>(seen);
    }

    private static void collectChannelsRecursive(GateNode node, Set<String> seen) {
        if (node.getChannel() != null) {
            seen.add(node.getChannel());
        }
        for (GateNode child : node.getPositiveChildren()) {
            collectChannelsRecursive(child, seen);
        }
        for (GateNode child : node.getNegativeChildren()) {
            collectChannelsRecursive(child, seen);
        }
    }

    /**
     * Recursively trace all paths through the gate tree. At each leaf (a node
     * with no children on a given branch), record the accumulated marker signs
     * keyed by the leaf phenotype name.
     */
    private static void traceMarkerSigns(GateNode node,
                                         Map<String, String> currentSigns,
                                         Map<String, Map<String, String>> result) {
        // Positive branch
        Map<String, String> positivePath = new LinkedHashMap<>(currentSigns);
        positivePath.put(node.getChannel(), "+");

        if (node.getPositiveChildren().isEmpty()) {
            result.put(node.getPositiveName(), new LinkedHashMap<>(positivePath));
        } else {
            for (GateNode child : node.getPositiveChildren()) {
                traceMarkerSigns(child, positivePath, result);
            }
        }

        // Negative branch
        Map<String, String> negativePath = new LinkedHashMap<>(currentSigns);
        negativePath.put(node.getChannel(), "-");

        if (node.getNegativeChildren().isEmpty()) {
            result.put(node.getNegativeName(), new LinkedHashMap<>(negativePath));
        } else {
            for (GateNode child : node.getNegativeChildren()) {
                traceMarkerSigns(child, negativePath, result);
            }
        }
    }

    /**
     * Escape a value for CSV output. Wraps the value in double quotes if it
     * contains a comma, double quote, or newline.
     */
    private static String escapeCsv(String value) {
        if (value == null)
            return "";
        if (value.contains(",") || value.contains("\"") || value.contains("\n")) {
            return "\"" + value.replace("\"", "\"\"") + "\"";
        }
        return value;
    }
}
