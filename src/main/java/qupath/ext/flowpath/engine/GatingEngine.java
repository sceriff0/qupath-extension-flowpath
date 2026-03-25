package qupath.ext.flowpath.engine;

import qupath.ext.flowpath.model.Branch;
import qupath.ext.flowpath.model.CellIndex;
import qupath.ext.flowpath.model.GateNode;
import qupath.ext.flowpath.model.GateTree;
import qupath.ext.flowpath.model.MarkerStats;
import qupath.ext.flowpath.model.QualityFilter;
import qupath.ext.flowpath.model.EllipseGate;
import qupath.ext.flowpath.model.PolygonGate;
import qupath.ext.flowpath.model.QuadrantGate;
import qupath.ext.flowpath.model.RectangleGate;
import qupath.lib.objects.PathObject;
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
                        index.getSolidity(i), index.getTotalIntensity(i), 0.0)) {
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
                    index.getTotalIntensity(i),
                    0.0);
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
            PathObject obj = index.getObject(i);
            ROI cellRoi = (obj != null) ? obj.getROI() : null;
            mask[i] = (cellRoi != null) && roi.contains(cellRoi.getCentroidX(), cellRoi.getCentroidY());
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
        if (!node.isEnabled()) return;
        if (node instanceof QuadrantGate qg) {
            walkQuadrantNode(qg, cellIdx, index, stats, useZScore, phenotypes, excluded, colors);
        } else if (node instanceof PolygonGate || node instanceof RectangleGate || node instanceof EllipseGate) {
            walk2DNode(node, cellIdx, index, stats, useZScore, phenotypes, excluded, colors);
        } else {
            walkThresholdNode(node, cellIdx, index, stats, useZScore, phenotypes, excluded, colors);
        }
    }

    private static void walkThresholdNode(GateNode node, int cellIdx,
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

        Branch branch;
        if (compareValue >= threshold) {
            branch = node.getBranches().get(0); // positive
        } else {
            branch = node.getBranches().get(1); // negative
        }
        branch.setCount(branch.getCount() + 1);
        assignBranch(branch, cellIdx, index, stats, useZScore, phenotypes, excluded, colors);
    }

    private static void walkQuadrantNode(QuadrantGate gate, int cellIdx,
                                          CellIndex index, MarkerStats stats, boolean useZScore,
                                          String[] phenotypes, boolean[] excluded, int[] colors) {
        int markerIdxX = index.getMarkerIndex(gate.getChannelX());
        int markerIdxY = index.getMarkerIndex(gate.getChannelY());
        if (markerIdxX < 0 || markerIdxY < 0) {
            return;
        }

        double rawX = index.getMarkerValues(markerIdxX)[cellIdx];
        double rawY = index.getMarkerValues(markerIdxY)[cellIdx];

        // Outlier exclusion on X channel
        if (gate.isExcludeOutliers()) {
            double loX = stats.getPercentileValue(gate.getChannelX(), gate.getClipPercentileLow());
            double hiX = stats.getPercentileValue(gate.getChannelX(), gate.getClipPercentileHigh());
            if (rawX < loX || rawX > hiX) {
                excluded[cellIdx] = true;
                return;
            }
            double loY = stats.getPercentileValue(gate.getChannelY(), gate.getClipPercentileLow());
            double hiY = stats.getPercentileValue(gate.getChannelY(), gate.getClipPercentileHigh());
            if (rawY < loY || rawY > hiY) {
                excluded[cellIdx] = true;
                return;
            }
        }

        double compareX = useZScore ? stats.toZScore(gate.getChannelX(), rawX) : rawX;
        double compareY = useZScore ? stats.toZScore(gate.getChannelY(), rawY) : rawY;

        int quadrant = gate.evaluateQuadrant(compareX, compareY);
        Branch branch = gate.getBranches().get(quadrant);
        branch.setCount(branch.getCount() + 1);
        assignBranch(branch, cellIdx, index, stats, useZScore, phenotypes, excluded, colors);
    }

    private static void walk2DNode(GateNode node, int cellIdx,
                                      CellIndex index, MarkerStats stats, boolean useZScore,
                                      String[] phenotypes, boolean[] excluded, int[] colors) {
        List<String> channels = node.getChannels();
        if (channels.size() < 2) return;
        String chX = channels.get(0);
        String chY = channels.get(1);
        int mxIdx = index.getMarkerIndex(chX);
        int myIdx = index.getMarkerIndex(chY);
        if (mxIdx < 0 || myIdx < 0) return;

        double rawX = index.getMarkerValues(mxIdx)[cellIdx];
        double rawY = index.getMarkerValues(myIdx)[cellIdx];

        // Outlier exclusion (same semantics as threshold/quadrant gates)
        if (node.isExcludeOutliers()) {
            double loX = stats.getPercentileValue(chX, node.getClipPercentileLow());
            double hiX = stats.getPercentileValue(chX, node.getClipPercentileHigh());
            if (rawX < loX || rawX > hiX) { excluded[cellIdx] = true; return; }
            double loY = stats.getPercentileValue(chY, node.getClipPercentileLow());
            double hiY = stats.getPercentileValue(chY, node.getClipPercentileHigh());
            if (rawY < loY || rawY > hiY) { excluded[cellIdx] = true; return; }
        }

        double vx = useZScore ? stats.toZScore(chX, rawX) : rawX;
        double vy = useZScore ? stats.toZScore(chY, rawY) : rawY;

        boolean inside;
        if (node instanceof PolygonGate pg) {
            inside = pg.contains(vx, vy);
        } else if (node instanceof RectangleGate rg) {
            inside = rg.contains(vx, vy);
        } else if (node instanceof EllipseGate eg) {
            inside = eg.contains(vx, vy);
        } else {
            inside = false;
        }

        Branch branch = inside ? node.getBranches().get(0) : node.getBranches().get(1);
        branch.setCount(branch.getCount() + 1);
        assignBranch(branch, cellIdx, index, stats, useZScore, phenotypes, excluded, colors);
    }

    private static void assignBranch(Branch branch, int cellIdx,
                                      CellIndex index, MarkerStats stats,
                                      boolean useZScore, String[] phenotypes,
                                      boolean[] excluded, int[] colors) {
        phenotypes[cellIdx] = branch.getName();
        colors[cellIdx] = branch.getColor();
        if (!branch.getChildren().isEmpty()) {
            for (GateNode child : branch.getChildren()) {
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
            for (Branch branch : node.getBranches()) {
                branch.setCount(0);
                resetCounts(branch.getChildren());
            }
        }
    }
}
