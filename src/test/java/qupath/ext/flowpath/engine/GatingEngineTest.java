package qupath.ext.flowpath.engine;

import org.junit.jupiter.api.Test;
import qupath.ext.flowpath.engine.GatingEngine.AssignmentResult;
import qupath.ext.flowpath.model.CellIndex;
import qupath.ext.flowpath.model.GateNode;
import qupath.ext.flowpath.model.GateTree;
import qupath.ext.flowpath.model.MarkerStats;
import qupath.ext.flowpath.model.QualityFilter;
import qupath.lib.objects.PathObject;
import qupath.lib.objects.PathObjects;
import qupath.lib.regions.ImagePlane;
import qupath.lib.roi.ROIs;
import qupath.lib.roi.interfaces.ROI;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class GatingEngineTest {

    // ---- helper ----

    /**
     * Build a synthetic CellIndex with the given markers, marker values, and area values.
     * Each cell is placed at coordinates (i*10, i*10).
     *
     * @param markers      list of marker names
     * @param markerValues [markerIndex][cellIndex] intensity values
     * @param areas        per-cell area values, or null to default to 100.0
     */
    private static CellIndex buildIndex(List<String> markers, double[][] markerValues, double[] areas) {
        int nCells = markerValues[0].length;
        List<PathObject> cells = new ArrayList<>();
        for (int i = 0; i < nCells; i++) {
            ROI roi = ROIs.createPointsROI(i * 10.0, i * 10.0, ImagePlane.getDefaultPlane());
            PathObject obj = PathObjects.createDetectionObject(roi);
            for (int m = 0; m < markers.size(); m++) {
                obj.getMeasurements().put(markers.get(m), markerValues[m][i]);
            }
            obj.getMeasurements().put("area", areas != null ? areas[i] : 100.0);
            cells.add(obj);
        }
        return CellIndex.build(cells, markers);
    }

    /**
     * Build a CellIndex where cells are placed at specific (x, y) coordinates.
     */
    private static CellIndex buildIndexWithCoords(List<String> markers, double[][] markerValues,
                                                   double[] xs, double[] ys) {
        int nCells = markerValues[0].length;
        List<PathObject> cells = new ArrayList<>();
        for (int i = 0; i < nCells; i++) {
            ROI roi = ROIs.createPointsROI(xs[i], ys[i], ImagePlane.getDefaultPlane());
            PathObject obj = PathObjects.createDetectionObject(roi);
            for (int m = 0; m < markers.size(); m++) {
                obj.getMeasurements().put(markers.get(m), markerValues[m][i]);
            }
            obj.getMeasurements().put("area", 100.0);
            cells.add(obj);
        }
        return CellIndex.build(cells, markers);
    }

    /**
     * Create an all-true quality mask of the given size.
     */
    private static boolean[] allTrueMask(int n) {
        boolean[] mask = new boolean[n];
        Arrays.fill(mask, true);
        return mask;
    }

    // ---- tests ----

    @Test
    void assignAllBasicPositiveNegativeSplit() {
        // 10 cells, CD45 values 1..10, threshold 5.5 (raw), no z-score
        List<String> markers = List.of("CD45");
        double[][] values = { {1, 2, 3, 4, 5, 6, 7, 8, 9, 10} };
        CellIndex index = buildIndex(markers, values, null);
        boolean[] mask = allTrueMask(10);
        MarkerStats stats = MarkerStats.compute(index, mask);

        GateNode gate = new GateNode("CD45", 5.5);
        gate.setThresholdIsZScore(false);

        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate);

        AssignmentResult result = GatingEngine.assignAll(tree, index, stats, false);

        int posCount = 0;
        int negCount = 0;
        for (String p : result.getPhenotypes()) {
            assertNotNull(p);
            if (p.equals("CD45+")) posCount++;
            else if (p.equals("CD45-")) negCount++;
        }
        // Values >= 5.5: 6,7,8,9,10 = 5 positive; 1,2,3,4,5 = 5 negative
        assertEquals(5, posCount, "Expected 5 cells positive (>= 5.5)");
        assertEquals(5, negCount, "Expected 5 cells negative (< 5.5)");
    }

    @Test
    void assignAllWithZScore() {
        // 10 cells with values 1..10, mean=5.5, z-score threshold=0 splits at the mean
        List<String> markers = List.of("CD45");
        double[][] values = { {1, 2, 3, 4, 5, 6, 7, 8, 9, 10} };
        CellIndex index = buildIndex(markers, values, null);
        boolean[] mask = allTrueMask(10);
        MarkerStats stats = MarkerStats.compute(index, mask);

        GateNode gate = new GateNode("CD45", 0.0);
        gate.setThresholdIsZScore(true);

        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate);

        AssignmentResult result = GatingEngine.assignAll(tree, index, stats, true);

        // mean = 5.5, so values >= 5.5 have z >= 0 -> positive
        // Values 6,7,8,9,10 -> positive (5); values 1,2,3,4,5 -> negative (5)
        int posCount = 0;
        int negCount = 0;
        for (String p : result.getPhenotypes()) {
            assertNotNull(p);
            if (p.equals("CD45+")) posCount++;
            else if (p.equals("CD45-")) negCount++;
        }
        assertEquals(5, posCount, "Expected 5 positive cells with z-score >= 0");
        assertEquals(5, negCount, "Expected 5 negative cells with z-score < 0");
    }

    @Test
    void qualityFilterExcludesByArea() {
        // 10 cells, areas 10,20,...,100. QF minArea=50 should exclude first 4.
        List<String> markers = List.of("CD45");
        double[][] values = { {1, 2, 3, 4, 5, 6, 7, 8, 9, 10} };
        double[] areas = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
        CellIndex index = buildIndex(markers, values, areas);

        QualityFilter qf = new QualityFilter();
        qf.setMinArea(50);

        boolean[] mask = GatingEngine.computeQualityMask(index, qf);

        int passing = 0;
        int excluded = 0;
        for (boolean b : mask) {
            if (b) passing++;
            else excluded++;
        }
        // Areas < 50: 10,20,30,40 = 4 excluded; areas >= 50: 50..100 = 6 passing
        assertEquals(6, passing, "Expected 6 cells passing with area >= 50");
        assertEquals(4, excluded, "Expected 4 cells excluded with area < 50");
    }

    @Test
    void roiMaskExcludesCellsOutsideRegion() {
        // Place 6 cells at known coordinates; rectangle ROI covers (0,0)-(55,55)
        List<String> markers = List.of("CD45");
        double[][] values = { {1, 2, 3, 4, 5, 6} };
        double[] xs = {5, 15, 25, 50, 80, 100};
        double[] ys = {5, 15, 25, 50, 80, 100};
        CellIndex index = buildIndexWithCoords(markers, values, xs, ys);

        // Create a rectangular ROI from (0,0) to (55,55)
        ROI rectRoi = ROIs.createRectangleROI(0, 0, 55, 55, ImagePlane.getDefaultPlane());

        boolean[] mask = GatingEngine.computeRoiMask(index, rectRoi);

        // Cells inside (55x55): (5,5), (15,15), (25,25), (50,50) = 4 inside
        // Cells outside: (80,80), (100,100) = 2 outside
        assertTrue(mask[0], "Cell at (5,5) should be inside ROI");
        assertTrue(mask[1], "Cell at (15,15) should be inside ROI");
        assertTrue(mask[2], "Cell at (25,25) should be inside ROI");
        assertTrue(mask[3], "Cell at (50,50) should be inside ROI");
        assertFalse(mask[4], "Cell at (80,80) should be outside ROI");
        assertFalse(mask[5], "Cell at (100,100) should be outside ROI");
    }

    @Test
    void outlierExclusionMarksOutliersExcluded() {
        // 100 cells: 98 with values 1..98, plus 2 outliers at -100 and +200
        List<String> markers = List.of("CD45");
        double[] vals = new double[100];
        for (int i = 0; i < 98; i++) {
            vals[i] = i + 1; // 1..98
        }
        vals[98] = -100.0; // outlier low
        vals[99] = 200.0;  // outlier high
        double[][] values = { vals };

        CellIndex index = buildIndex(markers, values, null);
        boolean[] mask = allTrueMask(100);
        MarkerStats stats = MarkerStats.compute(index, mask);

        GateNode gate = new GateNode("CD45", 50.0);
        gate.setThresholdIsZScore(false);
        gate.setExcludeOutliers(true);
        gate.setClipPercentileLow(1.0);
        gate.setClipPercentileHigh(99.0);

        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate);

        AssignmentResult result = GatingEngine.assignAll(tree, index, stats, false);

        // The two outlier cells should be excluded
        boolean[] excluded = result.getExcluded();
        assertTrue(excluded[98], "Low outlier (-100) should be excluded");
        assertTrue(excluded[99], "High outlier (+200) should be excluded");
        assertNull(result.getPhenotypes()[98], "Excluded cell phenotype should be null");
        assertNull(result.getPhenotypes()[99], "Excluded cell phenotype should be null");
    }

    @Test
    void nestedGatingOnlySeesParentPositiveCells() {
        // 10 cells: CD45 values 1..10, CD3 values 10..1 (reversed)
        // Root gate: CD45 threshold=5.5 (raw) -> cells 6-10 are CD45+
        // Child gate on positive branch: CD3 threshold=3.5 (raw)
        // Of cells 6-10, CD3 values are 5,4,3,2,1 -> CD3+ for cells 6,7 (values 5,4 >= 3.5)
        List<String> markers = List.of("CD45", "CD3");
        double[][] values = {
            {1, 2, 3, 4, 5, 6, 7, 8, 9, 10},  // CD45
            {10, 9, 8, 7, 6, 5, 4, 3, 2, 1}    // CD3
        };
        CellIndex index = buildIndex(markers, values, null);
        boolean[] mask = allTrueMask(10);
        MarkerStats stats = MarkerStats.compute(index, mask);

        GateNode root = new GateNode("CD45", 5.5);
        root.setThresholdIsZScore(false);

        GateNode childGate = new GateNode("CD3", 3.5);
        childGate.setThresholdIsZScore(false);
        root.getPositiveChildren().add(childGate);

        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(root);

        AssignmentResult result = GatingEngine.assignAll(tree, index, stats, false);
        String[] phenotypes = result.getPhenotypes();

        // Cells 0-4 (CD45 values 1-5): CD45- (no children on negative branch)
        for (int i = 0; i < 5; i++) {
            assertEquals("CD45-", phenotypes[i], "Cell " + i + " should be CD45-");
        }

        // Cells 5-6 (CD45 values 6,7 -> CD45+, CD3 values 5,4 >= 3.5 -> CD3+)
        assertEquals("CD3+", phenotypes[5], "Cell 5 should be CD3+ (CD3 value=5)");
        assertEquals("CD3+", phenotypes[6], "Cell 6 should be CD3+ (CD3 value=4)");

        // Cells 7-9 (CD45 values 8,9,10 -> CD45+, CD3 values 3,2,1 < 3.5 -> CD3-)
        assertEquals("CD3-", phenotypes[7], "Cell 7 should be CD3- (CD3 value=3)");
        assertEquals("CD3-", phenotypes[8], "Cell 8 should be CD3- (CD3 value=2)");
        assertEquals("CD3-", phenotypes[9], "Cell 9 should be CD3- (CD3 value=1)");
    }

    @Test
    void missingChannelSkipsGate() {
        // Gate on "NONEXISTENT" channel. All cells should remain Unclassified.
        List<String> markers = List.of("CD45");
        double[][] values = { {1, 2, 3, 4, 5} };
        CellIndex index = buildIndex(markers, values, null);
        boolean[] mask = allTrueMask(5);
        MarkerStats stats = MarkerStats.compute(index, mask);

        GateNode gate = new GateNode("NONEXISTENT", 0.0);
        gate.setThresholdIsZScore(false);

        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate);

        AssignmentResult result = GatingEngine.assignAll(tree, index, stats, false);

        for (int i = 0; i < 5; i++) {
            assertEquals("Unclassified", result.getPhenotypes()[i],
                    "Cell " + i + " should be Unclassified when channel is missing");
        }
    }

    @Test
    void emptyTreeLeavesAllUnclassified() {
        List<String> markers = List.of("CD45");
        double[][] values = { {1, 2, 3, 4, 5} };
        CellIndex index = buildIndex(markers, values, null);
        boolean[] mask = allTrueMask(5);
        MarkerStats stats = MarkerStats.compute(index, mask);

        GateTree tree = new GateTree();
        tree.setQualityFilter(null);

        AssignmentResult result = GatingEngine.assignAll(tree, index, stats, false);

        for (int i = 0; i < 5; i++) {
            assertEquals("Unclassified", result.getPhenotypes()[i],
                    "Cell " + i + " should be Unclassified with empty tree");
        }
        for (boolean ex : result.getExcluded()) {
            assertFalse(ex, "No cells should be excluded with empty tree");
        }
    }

    @Test
    void allCellsExcludedByQualityFilter() {
        // QF with impossible minArea should exclude everything
        List<String> markers = List.of("CD45");
        double[][] values = { {1, 2, 3, 4, 5} };
        CellIndex index = buildIndex(markers, values, null); // default area=100

        QualityFilter qf = new QualityFilter();
        qf.setMinArea(Double.MAX_VALUE);

        GateNode gate = new GateNode("CD45", 3.0);
        gate.setThresholdIsZScore(false);

        GateTree tree = new GateTree();
        tree.setQualityFilter(qf);
        tree.addRoot(gate);

        boolean[] mask = allTrueMask(5);
        MarkerStats stats = MarkerStats.compute(index, mask);

        AssignmentResult result = GatingEngine.assignAll(tree, index, stats, false);

        for (int i = 0; i < 5; i++) {
            assertTrue(result.getExcluded()[i], "Cell " + i + " should be excluded");
            assertNull(result.getPhenotypes()[i], "Excluded cell phenotype should be null");
        }
    }

    @Test
    void combineMasksAndsCorrectly() {
        boolean[] a = {true, true, false, false};
        boolean[] b = {true, false, true, false};
        boolean[] result = GatingEngine.combineMasks(a, b);

        boolean[] expected = {true, false, false, false};
        assertArrayEquals(expected, result, "combineMasks should AND element-wise");
    }

    @Test
    void computeQualityMaskHandlesNaN() {
        // Cells with NaN area should pass since QualityFilter.passes() skips NaN comparisons
        List<String> markers = List.of("CD45");
        double[][] values = { {1, 2, 3} };
        double[] areas = {50.0, Double.NaN, 150.0};
        CellIndex index = buildIndex(markers, values, areas);

        QualityFilter qf = new QualityFilter();
        qf.setMinArea(100);

        boolean[] mask = GatingEngine.computeQualityMask(index, qf);

        // area=50 < 100 -> fail; area=NaN -> pass (NaN check skips); area=150 >= 100 -> pass
        assertFalse(mask[0], "Cell with area=50 should fail minArea=100");
        assertTrue(mask[1], "Cell with NaN area should pass (NaN comparison bypassed)");
        assertTrue(mask[2], "Cell with area=150 should pass minArea=100");
    }
}
