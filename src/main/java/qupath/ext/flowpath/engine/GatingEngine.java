package qupath.ext.flowpath.engine;

import qupath.ext.flowpath.model.CellIndex;
import qupath.ext.flowpath.model.GateNode;
import qupath.ext.flowpath.model.GateTree;
import qupath.ext.flowpath.model.MarkerStats;
import qupath.ext.flowpath.model.QualityFilter;

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
     *
     * @param tree      the gate tree (roots + quality filter)
     * @param index     columnar cell data
     * @param stats     per-marker statistics (mean, std, percentiles)
     * @param useZScore if {@code true}, thresholds are compared against z-scored values
     * @return assignment result with phenotypes, exclusion flags, and colors
     */
    public static AssignmentResult assignAll(GateTree tree, CellIndex index, MarkerStats stats, boolean useZScore) {
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
        if (node.isHideOutliers()) {
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
