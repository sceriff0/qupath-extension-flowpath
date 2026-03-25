package qupath.ext.flowpath.io;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import qupath.ext.flowpath.engine.GatingEngine;
import qupath.ext.flowpath.engine.GatingEngine.AssignmentResult;
import qupath.ext.flowpath.model.*;
import qupath.lib.objects.PathObject;
import qupath.lib.objects.PathObjects;
import qupath.lib.regions.ImagePlane;
import qupath.lib.roi.ROIs;
import qupath.lib.roi.interfaces.ROI;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests that close coverage gaps in the CSV pipeline. Every test verifies
 * actual numeric values in the CSV output, not just phenotype names.
 */
class CsvCorrectnessTest {

    @TempDir
    Path tempDir;

    // ========== Helpers ==========

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

    private static boolean[] allTrueMask(int n) {
        boolean[] mask = new boolean[n];
        Arrays.fill(mask, true);
        return mask;
    }

    private static List<String> parseCsvLine(String line) {
        List<String> fields = new ArrayList<>();
        StringBuilder current = new StringBuilder();
        boolean inQuotes = false;
        for (int i = 0; i < line.length(); i++) {
            char c = line.charAt(i);
            if (inQuotes) {
                if (c == '"') {
                    if (i + 1 < line.length() && line.charAt(i + 1) == '"') {
                        current.append('"'); i++;
                    } else { inQuotes = false; }
                } else { current.append(c); }
            } else {
                if (c == '"') { inQuotes = true; }
                else if (c == ',') { fields.add(current.toString()); current.setLength(0); }
                else { current.append(c); }
            }
        }
        fields.add(current.toString());
        return fields;
    }

    record CsvResult(List<String> header, List<List<String>> rows, AssignmentResult result) {
        String val(int row, String col) {
            int idx = header.indexOf(col);
            assertTrue(idx >= 0, "Column '" + col + "' not in header: " + header);
            return rows.get(row).get(idx);
        }
    }

    private CsvResult run(GateTree tree, CellIndex index, MarkerStats stats,
                          boolean useZScore, String name) throws IOException {
        AssignmentResult result = GatingEngine.assignAll(tree, index, stats, useZScore);
        File f = tempDir.resolve(name).toFile();
        PhenotypeCsvExporter.export(f, index, result, tree, stats);
        List<String> lines = Files.readAllLines(f.toPath());
        List<String> hdr = parseCsvLine(lines.get(0));
        List<List<String>> rows = new ArrayList<>();
        for (int i = 1; i < lines.size(); i++) rows.add(parseCsvLine(lines.get(i)));
        return new CsvResult(hdr, rows, result);
    }

    // ========== Gap 10: Numeric raw values verified ==========

    @Test
    void rawValuesInCsvMatchInput() throws IOException {
        // Use both markers in gates so both appear in CSV columns
        List<String> markers = List.of("CD45", "CD3");
        double[][] values = { {2.5, 7.3}, {1.1, 9.9} };
        CellIndex index = buildIndex(markers, values, null);
        MarkerStats stats = MarkerStats.compute(index, allTrueMask(2));

        GateNode gate1 = new GateNode("CD45", 5.0);
        gate1.setThresholdIsZScore(false);
        GateNode gate2 = new GateNode("CD3", 5.0);
        gate2.setThresholdIsZScore(false);
        gate1.getPositiveChildren().add(gate2);

        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate1);

        CsvResult csv = run(tree, index, stats, false, "raw.csv");

        // Cell 0: CD45=2.5, CD3=1.1
        assertEquals("2.5000", csv.val(0, "CD45_raw"));
        assertEquals("1.1000", csv.val(0, "CD3_raw"));
        // Cell 1: CD45=7.3, CD3=9.9
        assertEquals("7.3000", csv.val(1, "CD45_raw"));
        assertEquals("9.9000", csv.val(1, "CD3_raw"));
    }

    @Test
    void ungatedMarkersAppearInCsv() throws IOException {
        // 3 markers but only CD45 is gated — CD3 and CD8 should still appear as raw+zscore columns
        List<String> markers = List.of("CD45", "CD3", "CD8");
        double[][] values = { {2, 8}, {3, 7}, {4, 6} };
        CellIndex index = buildIndex(markers, values, null);
        MarkerStats stats = MarkerStats.compute(index, allTrueMask(2));

        GateNode gate = new GateNode("CD45", 5.0);
        gate.setThresholdIsZScore(false);
        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate);

        CsvResult csv = run(tree, index, stats, false, "ungated.csv");

        // All 3 markers should have raw columns
        assertEquals("2.0000", csv.val(0, "CD45_raw"));
        assertEquals("3.0000", csv.val(0, "CD3_raw"));
        assertEquals("4.0000", csv.val(0, "CD8_raw"));
        // Sign columns: only CD45 has signs (gated), CD3 and CD8 should be empty
        assertEquals("-", csv.val(0, "CD45_sign"));
        assertEquals("", csv.val(0, "CD3_sign"));
        assertEquals("", csv.val(0, "CD8_sign"));
    }

    // ========== Gap 10 cont: Z-score values verified ==========

    @Test
    void zscoreValuesInCsvAreCorrect() throws IOException {
        // 4 cells, CD45 values [2, 4, 6, 8]. mean=5, std=sqrt(5)≈2.2361
        List<String> markers = List.of("CD45");
        double[][] values = { {2, 4, 6, 8} };
        CellIndex index = buildIndex(markers, values, null);
        boolean[] mask = allTrueMask(4);
        MarkerStats stats = MarkerStats.compute(index, mask);

        GateNode gate = new GateNode("CD45", 0.0);
        gate.setThresholdIsZScore(true);
        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate);

        CsvResult csv = run(tree, index, stats, true, "zscore.csv");

        double mean = stats.getMean("CD45");
        double std = stats.getStd("CD45");

        for (int i = 0; i < 4; i++) {
            double raw = values[0][i];
            double expectedZ = (raw - mean) / std;
            String actual = csv.val(i, "CD45_zscore");
            assertEquals(expectedZ, Double.parseDouble(actual), 0.001,
                "Cell " + i + " z-score should be " + expectedZ);
        }
    }

    // ========== Gap 11: Locale.US decimal format ==========

    @Test
    void csvUsesDotsNotCommasForDecimals() throws IOException {
        List<String> markers = List.of("CD45");
        double[][] values = { {3.14159} };
        CellIndex index = buildIndex(markers, values, null);
        MarkerStats stats = MarkerStats.compute(index, allTrueMask(1));

        GateNode gate = new GateNode("CD45", 0.0);
        gate.setThresholdIsZScore(false);
        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate);

        CsvResult csv = run(tree, index, stats, false, "locale.csv");

        String raw = csv.val(0, "CD45_raw");
        assertTrue(raw.contains("."), "Decimal should use dot, got: " + raw);
        assertFalse(raw.contains(","), "Should not contain comma as decimal separator");
        assertEquals("3.1416", raw); // 4 decimal places, rounded
    }

    // ========== Gap 6: Outlier exclusion for 2D gates ==========

    @Test
    void outlierExclusionWorksForRectangleGate() throws IOException {
        // 100 cells, CD45 and CD3 values 1..100. Two outliers at indices 98,99.
        List<String> markers = List.of("CD45", "CD3");
        double[] cd45 = new double[100];
        double[] cd3 = new double[100];
        for (int i = 0; i < 98; i++) { cd45[i] = i + 1; cd3[i] = i + 1; }
        cd45[98] = -500; cd3[98] = -500; // outlier
        cd45[99] = 500; cd3[99] = 500; // outlier
        double[][] values = { cd45, cd3 };
        CellIndex index = buildIndex(markers, values, null);
        MarkerStats stats = MarkerStats.compute(index, allTrueMask(100));

        RectangleGate gate = new RectangleGate("CD45", "CD3", 0, 50, 0, 50);
        gate.setExcludeOutliers(true);
        gate.setClipPercentileLow(1.0);
        gate.setClipPercentileHigh(99.0);

        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate);

        AssignmentResult result = GatingEngine.assignAll(tree, index, stats, false);
        assertTrue(result.getExcluded()[98], "Cell 98 (outlier -500) should be excluded");
        assertTrue(result.getExcluded()[99], "Cell 99 (outlier 500) should be excluded");
    }

    @Test
    void outlierExclusionWorksForQuadrantGate() throws IOException {
        double[] cd45 = new double[100];
        double[] cd3 = new double[100];
        for (int i = 0; i < 98; i++) { cd45[i] = i + 1; cd3[i] = i + 1; }
        cd45[98] = -500; cd3[98] = 50;
        cd45[99] = 50; cd3[99] = -500;
        double[][] values = { cd45, cd3 };
        CellIndex index = buildIndex(List.of("CD45", "CD3"), values, null);
        MarkerStats stats = MarkerStats.compute(index, allTrueMask(100));

        QuadrantGate gate = new QuadrantGate("CD45", "CD3");
        gate.setThresholdX(50);
        gate.setThresholdY(50);
        gate.setThresholdIsZScore(false);
        gate.setExcludeOutliers(true);
        gate.setClipPercentileLow(1.0);
        gate.setClipPercentileHigh(99.0);

        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate);

        AssignmentResult result = GatingEngine.assignAll(tree, index, stats, false);
        assertTrue(result.getExcluded()[98], "Cell 98 (CD45 outlier) should be excluded");
        assertTrue(result.getExcluded()[99], "Cell 99 (CD3 outlier) should be excluded");
    }

    // ========== Gap 7: Z-score mode for QuadrantGate ==========

    @Test
    void quadrantGateWithZScoreMode() throws IOException {
        // 4 cells. With z-score mode, threshold 0 splits at the mean.
        // CD45 values [1,9,1,9], mean=5, so z>0 for 9, z<0 for 1
        // CD3 values  [1,1,9,9], mean=5, so z>0 for 9, z<0 for 1
        List<String> markers = List.of("CD45", "CD3");
        double[][] values = { {1, 9, 1, 9}, {1, 1, 9, 9} };
        CellIndex index = buildIndex(markers, values, null);
        MarkerStats stats = MarkerStats.compute(index, allTrueMask(4));

        QuadrantGate gate = new QuadrantGate("CD45", "CD3");
        gate.setThresholdX(0.0);
        gate.setThresholdY(0.0);
        gate.setThresholdIsZScore(true);

        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate);

        AssignmentResult result = GatingEngine.assignAll(tree, index, stats, true);

        // Cell 0: CD45=1 (z<0), CD3=1 (z<0) -> NN
        assertEquals("CD45-/CD3-", result.getPhenotypes()[0]);
        // Cell 1: CD45=9 (z>0), CD3=1 (z<0) -> PN
        assertEquals("CD45+/CD3-", result.getPhenotypes()[1]);
        // Cell 2: CD45=1 (z<0), CD3=9 (z>0) -> NP
        assertEquals("CD45-/CD3+", result.getPhenotypes()[2]);
        // Cell 3: CD45=9 (z>0), CD3=9 (z>0) -> PP
        assertEquals("CD45+/CD3+", result.getPhenotypes()[3]);
    }

    // ========== Gap 8: Z-score mode for 2D gates ==========

    @Test
    void rectangleGateWithZScoreMode() throws IOException {
        // Rectangle gate in z-score space: bounds from -0.5 to 0.5 on both axes
        // CD45 values [1,5,9], mean=5, std≈3.27
        // z-scores: -1.22, 0.0, 1.22
        // Only cell 1 (z=0) is inside [-0.5, 0.5]
        List<String> markers = List.of("CD45", "CD3");
        double[][] values = { {1, 5, 9}, {1, 5, 9} };
        CellIndex index = buildIndex(markers, values, null);
        MarkerStats stats = MarkerStats.compute(index, allTrueMask(3));

        RectangleGate gate = new RectangleGate("CD45", "CD3", -0.5, 0.5, -0.5, 0.5);

        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate);

        AssignmentResult result = GatingEngine.assignAll(tree, index, stats, true);

        String insideName = gate.getBranches().get(0).getName();
        String outsideName = gate.getBranches().get(1).getName();

        assertEquals(outsideName, result.getPhenotypes()[0], "Cell 0 z=-1.22 should be outside");
        assertEquals(insideName, result.getPhenotypes()[1], "Cell 1 z=0 should be inside");
        assertEquals(outsideName, result.getPhenotypes()[2], "Cell 2 z=1.22 should be outside");
    }

    // ========== Gap 12: Sign columns for 2D region gates ==========

    @Test
    void regionGateSignsIncludeBothChannels() throws IOException {
        // Rectangle gate on CD45 vs CD3. Inside → both signs "+", outside → both "-"
        List<String> markers = List.of("CD45", "CD3");
        double[][] values = { {5, 20}, {5, 20} };
        CellIndex index = buildIndex(markers, values, null);
        MarkerStats stats = MarkerStats.compute(index, allTrueMask(2));

        RectangleGate gate = new RectangleGate("CD45", "CD3", 2, 8, 2, 8);

        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate);

        CsvResult csv = run(tree, index, stats, false, "region_signs.csv");

        // Cell 0: (5,5) inside → both signs "+"
        assertEquals("+", csv.val(0, "CD45_sign"));
        assertEquals("+", csv.val(0, "CD3_sign"));
        // Cell 1: (20,20) outside → both signs "-"
        assertEquals("-", csv.val(1, "CD45_sign"));
        assertEquals("-", csv.val(1, "CD3_sign"));
    }

    // ========== Gap 12 cont: Quadrant gate signs for both channels ==========

    @Test
    void quadrantGateSignsForBothChannelsInCsv() throws IOException {
        List<String> markers = List.of("CD45", "CD3");
        double[][] values = { {8, 2, 8, 2}, {8, 8, 2, 2} };
        CellIndex index = buildIndex(markers, values, null);
        MarkerStats stats = MarkerStats.compute(index, allTrueMask(4));

        QuadrantGate gate = new QuadrantGate("CD45", "CD3");
        gate.setThresholdX(5.0);
        gate.setThresholdY(5.0);
        gate.setThresholdIsZScore(false);

        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate);

        CsvResult csv = run(tree, index, stats, false, "quad_signs.csv");

        // Cell 0: PP → CD45+, CD3+
        assertEquals("+", csv.val(0, "CD45_sign"));
        assertEquals("+", csv.val(0, "CD3_sign"));
        // Cell 1: NP → CD45-, CD3+
        assertEquals("-", csv.val(1, "CD45_sign"));
        assertEquals("+", csv.val(1, "CD3_sign"));
        // Cell 2: PN → CD45+, CD3-
        assertEquals("+", csv.val(2, "CD45_sign"));
        assertEquals("-", csv.val(2, "CD3_sign"));
        // Cell 3: NN → CD45-, CD3-
        assertEquals("-", csv.val(3, "CD45_sign"));
        assertEquals("-", csv.val(3, "CD3_sign"));
    }

    // ========== Gap 14: Large dataset CSV ==========

    @Test
    void largeDatasetCsvHasCorrectRowCount() throws IOException {
        int n = 500;
        List<String> markers = List.of("CD45");
        double[] vals = new double[n];
        for (int i = 0; i < n; i++) vals[i] = i;
        double[][] values = { vals };
        CellIndex index = buildIndex(markers, values, null);
        MarkerStats stats = MarkerStats.compute(index, allTrueMask(n));

        GateNode gate = new GateNode("CD45", 250.0);
        gate.setThresholdIsZScore(false);

        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate);

        CsvResult csv = run(tree, index, stats, false, "large.csv");
        assertEquals(n, csv.rows.size(), "All 500 cells should be in CSV");

        // Verify first and last rows have correct raw values
        assertEquals("0.0000", csv.val(0, "CD45_raw"));
        assertEquals("499.0000", csv.val(n - 1, "CD45_raw"));

        // Verify phenotype split: 0..249 negative, 250..499 positive
        int pos = 0, neg = 0;
        for (int i = 0; i < n; i++) {
            String p = csv.val(i, "phenotype");
            if (p.equals("CD45+")) pos++;
            else if (p.equals("CD45-")) neg++;
        }
        assertEquals(250, pos, "250 cells should be positive (>= 250)");
        assertEquals(250, neg, "250 cells should be negative (< 250)");
    }

    // ========== Gap 15: Serialization round-trip for 2D gates ==========

    @Test
    void polygonGateSerializationRoundTrip() throws IOException {
        PolygonGate gate = new PolygonGate("CD45", "CD3");
        gate.setVertices(List.of(new double[]{0, 0}, new double[]{10, 0}, new double[]{5, 10}));
        gate.getBranches().get(0).setName("Custom Inside");
        gate.setClipPercentileLow(2.0);

        GateTree tree = new GateTree();
        tree.addRoot(gate);
        File f = tempDir.resolve("polygon.json").toFile();
        FlowPathSerializer.save(tree, f);
        GateTree loaded = FlowPathSerializer.load(f);

        assertEquals(1, loaded.getRoots().size());
        assertInstanceOf(PolygonGate.class, loaded.getRoots().get(0));
        PolygonGate lg = (PolygonGate) loaded.getRoots().get(0);
        assertEquals(3, lg.getVertices().size());
        assertEquals(0, lg.getVertices().get(0)[0]);
        assertEquals(10, lg.getVertices().get(1)[0]);
        assertEquals(5, lg.getVertices().get(2)[0]);
        assertEquals("CD45", lg.getChannelX());
        assertEquals("CD3", lg.getChannelY());
        assertEquals("Custom Inside", lg.getBranches().get(0).getName());
        assertEquals(2.0, lg.getClipPercentileLow());
    }

    @Test
    void rectangleGateSerializationRoundTrip() throws IOException {
        RectangleGate gate = new RectangleGate("CD45", "CD3", 1.5, 8.5, 2.0, 9.0);
        gate.getBranches().get(1).setName("Custom Outside");

        GateTree tree = new GateTree();
        tree.addRoot(gate);
        File f = tempDir.resolve("rect.json").toFile();
        FlowPathSerializer.save(tree, f);
        GateTree loaded = FlowPathSerializer.load(f);

        assertInstanceOf(RectangleGate.class, loaded.getRoots().get(0));
        RectangleGate lg = (RectangleGate) loaded.getRoots().get(0);
        assertEquals(1.5, lg.getMinX());
        assertEquals(8.5, lg.getMaxX());
        assertEquals(2.0, lg.getMinY());
        assertEquals(9.0, lg.getMaxY());
        assertEquals("Custom Outside", lg.getBranches().get(1).getName());
    }

    @Test
    void ellipseGateSerializationRoundTrip() throws IOException {
        EllipseGate gate = new EllipseGate("CD45", "CD3", 5.0, 5.0, 3.0, 2.0);

        GateTree tree = new GateTree();
        tree.addRoot(gate);
        File f = tempDir.resolve("ellipse.json").toFile();
        FlowPathSerializer.save(tree, f);
        GateTree loaded = FlowPathSerializer.load(f);

        assertInstanceOf(EllipseGate.class, loaded.getRoots().get(0));
        EllipseGate lg = (EllipseGate) loaded.getRoots().get(0);
        assertEquals(5.0, lg.getCenterX());
        assertEquals(5.0, lg.getCenterY());
        assertEquals(3.0, lg.getRadiusX());
        assertEquals(2.0, lg.getRadiusY());
    }

    // ========== Gap 17: Nested children in 2D gate branches survive round-trip ==========

    @Test
    void nestedChildInPolygonBranchSurvivesRoundTrip() throws IOException {
        PolygonGate parent = new PolygonGate("CD45", "CD3");
        parent.setVertices(List.of(new double[]{0, 0}, new double[]{10, 0}, new double[]{5, 10}));

        GateNode child = new GateNode("CD8", 3.0);
        child.setThresholdIsZScore(false);
        parent.getBranches().get(0).getChildren().add(child);

        GateTree tree = new GateTree();
        tree.addRoot(parent);
        File f = tempDir.resolve("nested_poly.json").toFile();
        FlowPathSerializer.save(tree, f);
        GateTree loaded = FlowPathSerializer.load(f);

        PolygonGate lp = (PolygonGate) loaded.getRoots().get(0);
        assertEquals(1, lp.getBranches().get(0).getChildren().size());
        GateNode lc = lp.getBranches().get(0).getChildren().get(0);
        assertEquals("CD8", lc.getChannel());
        assertEquals(3.0, lc.getThreshold());
    }

    // ========== Gap 4: Z-score NaN when std=0 ==========

    @Test
    void zscoreIsEmptyInCsvWhenAllValuesIdentical() throws IOException {
        // All cells have identical CD45=5.0 → std=0 → zscore should be empty in CSV
        List<String> markers = List.of("CD45");
        double[][] values = { {5, 5, 5} };
        CellIndex index = buildIndex(markers, values, null);
        MarkerStats stats = MarkerStats.compute(index, allTrueMask(3));

        GateNode gate = new GateNode("CD45", 0.0);
        gate.setThresholdIsZScore(false);
        GateTree tree = new GateTree();
        tree.setQualityFilter(null);
        tree.addRoot(gate);

        CsvResult csv = run(tree, index, stats, false, "zero_std.csv");

        // Z-score should be empty (NaN → "") when std is zero
        assertEquals("", csv.val(0, "CD45_zscore"));
        assertEquals("", csv.val(1, "CD45_zscore"));
        // Raw should still be present
        assertEquals("5.0000", csv.val(0, "CD45_raw"));
    }

    // ========== Enabled flag survives serialization ==========

    @Test
    void enabledFlagSerializationRoundTrip() throws IOException {
        GateNode gate = new GateNode("CD45", 3.0);
        gate.setEnabled(false);

        GateTree tree = new GateTree();
        tree.addRoot(gate);
        File f = tempDir.resolve("enabled.json").toFile();
        FlowPathSerializer.save(tree, f);
        GateTree loaded = FlowPathSerializer.load(f);

        assertFalse(loaded.getRoots().get(0).isEnabled());
    }
}
