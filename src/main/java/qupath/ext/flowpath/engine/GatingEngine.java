package qupath.ext.flowpath.engine;

import qupath.ext.flowpath.model.CellIndex;
import qupath.ext.flowpath.model.GateNode;
import qupath.ext.flowpath.model.GateTree;
import qupath.ext.flowpath.model.MarkerStats;
import qupath.ext.flowpath.model.QualityFilter;
import qupath.lib.roi.interfaces.ROI;

import java.util.Arrays;
import java.util.List;

/**
 * Core gating logic that walks a {@link GateTree} and assigns phenotype labels
 * and colors to every cell in a {@link CellIndex}.
 */
public final class GatingEngine {

    private GatingEngine() {
        // static utility class
    }

    /**
     * Result of running the gating engine over all cells.
     */
    public static final class AssignmentResult {
        private final String[] phenotypes;
        private final boolean[] excluded;
        private final int[] colors;

        AssignmentResult(String[] phenotypes, boolean[] excluded, int[] colors) {
            this.phenotypes = phenotypes;
            this.excluded = excluded;
            this.colors = colors;
        }

        /** Phenotype label per cell, {@code null} for excluded cells. */
        public String[] getPhenotypes() {
            return phenotypes;
        }

        /** {@code true} for cells removed by the quality filter or outlier exclusion. */
        public boolean[] getExcluded() {
            return excluded;
        }

        /** Packed RGB color per cell, 0 for excluded cells. */
        public int[] getColors() {
            return colors;
        }
    }

    /**
     * Assign phenotypes to every cell by walking the gate tree.
     * Delegates to {@link #assignAll(GateTree, CellIndex, MarkerStats, boolean, boolean[])}
     * with no ROI mask.
     */
    public static AssignmentResult assignAll(GateTree tree, CellIndex index, MarkerStats stats, boolean useZScore) {
        return assignAll(tree, index, stats, useZScore, null);
    }

    /**
     * Assign phenotypes to every cell by walking the gate tree.
     *
     * @param tree      the gate tree (roots + quality filter)
     * @param index     columnar cell data
     * @param stats     per-marker statistics (mean, std, percentiles)
     * @param useZScore if {@code true}, thresholds are compared against z-scored values
     * @param roiMask   optional boolean mask where {@code true} means the cell is inside the ROI;
     *                  {@code null} means no ROI filtering
     * @return assignment result with phenotypes, exclusion flags, and colors
     */
    public static AssignmentResult assignAll(GateTree tree, CellIndex index, MarkerStats stats,
                                              boolean useZScore, boolean[] roiMask) {
        int n = index.size();
        String[] phenotypes = new String[n];
        boolean[] excluded = new boolean[n];
        int[] colors = new int[n];

        // 1. Initialize all as Unclassified
        for (int i = 0; i < n; i++) {
            phenotypes[i] = "Unclassified";
        }

        // 2. Apply quality filter
        QualityFilter qf = tree.getQualityFilter();
        if (qf != null) {
            for (int i = 0; i < n; i++) {
                if (!qf.passes(index.getArea(i), index.getEccentricity(i),
                        index.getSolidity(i), index.getTotalIntensity(i))) {
                    excluded[i] = true;
                    phenotypes[i] = null;
                }
            }
        }

        // 2b. Apply ROI mask
        if (roiMask != null) {
            for (int i = 0; i < n; i++) {
                if (!roiMask[i]) {
                    excluded[i] = true;
                    phenotypes[i] = null;
                }
            }
        }

        // 3. Walk gate tree for non-excluded cells
        // Reset transient counts on all nodes before walking
        List<GateNode> roots = tree.getRoots();
        resetCounts(roots);

        for (int i = 0; i < n; i++) {
            if (excluded[i]) {
                continue;
            }
            walkRoots(roots, i, index, stats, useZScore, phenotypes, excluded, colors);
        }

        // Null out phenotypes for cells that got excluded during outlier checks
        for (int i = 0; i < n; i++) {
            if (excluded[i]) {
                phenotypes[i] = null;
                colors[i] = 0;
            }
        }

        return new AssignmentResult(phenotypes, excluded, colors);
    }

    /**
     * Compute a boolean mask indicating which cells pass the quality filter.
     *
     * @param index  columnar cell data
     * @param filter quality filter criteria
     * @return boolean array where {@code true} means the cell passes
     */
    public static boolean[] computeQualityMask(CellIndex index, QualityFilter filter) {
        int n = index.size();
        boolean[] mask = new boolean[n];
        for (int i = 0; i < n; i++) {
            mask[i] = filter.passes(
                    index.getArea(i),
                    index.getEccentricity(i),
                    index.getSolidity(i),
                    index.getTotalIntensity(i));
        }
        return mask;
    }

    /**
     * Compute a boolean mask indicating which cells fall inside the given ROI.
     * If {@code roi} is {@code null}, all cells pass.
     */
    public static boolean[] computeRoiMask(CellIndex index, ROI roi) {
        int n = index.size();
        boolean[] mask = new boolean[n];
        if (roi == null) {
            Arrays.fill(mask, true);
            return mask;
        }
        for (int i = 0; i < n; i++) {
            mask[i] = roi.contains(index.getCentroidX(i), index.getCentroidY(i));
        }
        return mask;
    }

    /**
     * Combine two boolean masks with logical AND. Both arrays must have the same length.
     */
    public static boolean[] combineMasks(boolean[] a, boolean[] b) {
        boolean[] result = new boolean[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = a[i] && b[i];
        }
        return result;
    }

    // ---- private helpers ----

    private static void walkRoots(List<GateNode> roots, int cellIdx,
                                  CellIndex index, MarkerStats stats, boolean useZScore,
                                  String[] phenotypes, boolean[] excluded, int[] colors) {
        for (GateNode root : roots) {
            if (excluded[cellIdx]) {
                return;
            }
            walkNode(root, cellIdx, index, stats, useZScore, phenotypes, excluded, colors);
        }
    }

    private static void walkNode(GateNode node, int cellIdx,
                                 CellIndex index, MarkerStats stats, boolean useZScore,
                                 String[] phenotypes, boolean[] excluded, int[] colors) {
        String channel = node.getChannel();
        int markerIdx = index.getMarkerIndex(channel);
        if (markerIdx < 0) {
            return;
        }

        double rawValue = index.getMarkerValues(markerIdx)[cellIdx];

        // Outlier exclusion based on percentile clip bounds
        if (node.isExcludeOutliers()) {
            double lo = stats.getPercentileValue(channel, node.getClipPercentileLow());
            double hi = stats.getPercentileValue(channel, node.getClipPercentileHigh());
            if (rawValue < lo || rawValue > hi) {
                excluded[cellIdx] = true;
                return;
            }
        }

        double compareValue = useZScore ? stats.toZScore(channel, rawValue) : rawValue;
        double threshold = node.getThreshold();

        if (compareValue >= threshold) {
            node.setPosCount(node.getPosCount() + 1);
            assignBranch(node.getPositiveName(), node.getPositiveColor(),
                    node.getPositiveChildren(), cellIdx, index, stats,
                    useZScore, phenotypes, excluded, colors);
        } else {
            node.setNegCount(node.getNegCount() + 1);
            assignBranch(node.getNegativeName(), node.getNegativeColor(),
                    node.getNegativeChildren(), cellIdx, index, stats,
                    useZScore, phenotypes, excluded, colors);
        }
    }

    private static void assignBranch(String name, int color, List<GateNode> children,
                                      int cellIdx, CellIndex index, MarkerStats stats,
                                      boolean useZScore, String[] phenotypes,
                                      boolean[] excluded, int[] colors) {
        phenotypes[cellIdx] = name;
        colors[cellIdx] = color;
        if (children != null && !children.isEmpty()) {
            for (GateNode child : children) {
                if (excluded[cellIdx]) return;
                walkNode(child, cellIdx, index, stats, useZScore, phenotypes, excluded, colors);
            }
        }
    }

    private static void resetCounts(List<GateNode> nodes) {
        if (nodes == null) {
            return;
        }
        for (GateNode node : nodes) {
            node.setPosCount(0);
            node.setNegCount(0);
            resetCounts(node.getPositiveChildren());
            resetCounts(node.getNegativeChildren());
        }
    }
}
