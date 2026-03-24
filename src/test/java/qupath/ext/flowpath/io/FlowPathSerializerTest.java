package qupath.ext.flowpath.io;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import qupath.ext.flowpath.model.GateNode;
import qupath.ext.flowpath.model.GateTree;
import qupath.ext.flowpath.model.QualityFilter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;

import static org.junit.jupiter.api.Assertions.*;

class FlowPathSerializerTest {

    @TempDir
    Path tempDir;

    @Test
    void saveAndLoadRoundTrip() throws IOException {
        var tree = new GateTree();
        var root = new GateNode("CD45", 1.5);
        root.setPositiveName("Immune+");
        root.setNegativeName("Immune-");
        root.setPositiveColor((0 << 16) | (255 << 8) | 0);
        root.setNegativeColor((128 << 16) | (128 << 8) | 128);
        root.setClipPercentileLow(2.0);
        root.setClipPercentileHigh(98.0);
        root.setExcludeOutliers(true);
        root.setThresholdIsZScore(false);
        tree.addRoot(root);

        File file = tempDir.resolve("test.json").toFile();
        FlowPathSerializer.save(tree, file);
        GateTree loaded = FlowPathSerializer.load(file);

        assertEquals(1, loaded.getRoots().size());
        GateNode loadedRoot = loaded.getRoots().get(0);
        assertEquals("CD45", loadedRoot.getChannel());
        assertEquals(1.5, loadedRoot.getThreshold());
        assertFalse(loadedRoot.isThresholdIsZScore());
        assertEquals("Immune+", loadedRoot.getPositiveName());
        assertEquals("Immune-", loadedRoot.getNegativeName());
        assertEquals((0 << 16) | (255 << 8) | 0, loadedRoot.getPositiveColor());
        assertEquals((128 << 16) | (128 << 8) | 128, loadedRoot.getNegativeColor());
        assertEquals(2.0, loadedRoot.getClipPercentileLow());
        assertEquals(98.0, loadedRoot.getClipPercentileHigh());
        assertTrue(loadedRoot.isExcludeOutliers());
    }

    @Test
    void nestedChildrenRoundTrip() throws IOException {
        var tree = new GateTree();
        var root = new GateNode("CD45", 0.0);
        var child = new GateNode("CD3", 1.0);
        var grandchild = new GateNode("CD8", 2.0);
        child.getPositiveChildren().add(grandchild);
        root.getPositiveChildren().add(child);
        root.getNegativeChildren().add(new GateNode("PANCK", -0.5));
        tree.addRoot(root);

        File file = tempDir.resolve("nested.json").toFile();
        FlowPathSerializer.save(tree, file);
        GateTree loaded = FlowPathSerializer.load(file);

        GateNode lr = loaded.getRoots().get(0);
        assertEquals(1, lr.getPositiveChildren().size());
        assertEquals(1, lr.getNegativeChildren().size());
        assertEquals("CD3", lr.getPositiveChildren().get(0).getChannel());
        assertEquals("PANCK", lr.getNegativeChildren().get(0).getChannel());

        GateNode loadedGrandchild = lr.getPositiveChildren().get(0).getPositiveChildren().get(0);
        assertEquals("CD8", loadedGrandchild.getChannel());
        assertEquals(2.0, loadedGrandchild.getThreshold());
    }

    @Test
    void qualityFilterRoundTrip() throws IOException {
        var tree = new GateTree();
        var qf = new QualityFilter();
        qf.setMinArea(25);
        qf.setMaxArea(500);
        qf.setMinTotalIntensity(100);
        qf.setMaxEccentricity(0.9);
        qf.setMinSolidity(0.5);
        tree.setQualityFilter(qf);
        tree.addRoot(new GateNode("CD45"));

        File file = tempDir.resolve("qf.json").toFile();
        FlowPathSerializer.save(tree, file);
        GateTree loaded = FlowPathSerializer.load(file);

        QualityFilter lqf = loaded.getQualityFilter();
        assertEquals(25, lqf.getMinArea());
        assertEquals(500, lqf.getMaxArea());
        assertEquals(100, lqf.getMinTotalIntensity());
        assertEquals(0.9, lqf.getMaxEccentricity());
        assertEquals(0.5, lqf.getMinSolidity());
    }

    @Test
    void legacyHideOutliersKeyLoadsAsExcludeOutliers() throws IOException {
        // Simulate a v1 JSON with the old "hideOutliers" key
        String json = """
                {
                  "version": 1,
                  "qualityFilter": {},
                  "gates": [
                    {
                      "channel": "CD45",
                      "threshold": 0.0,
                      "thresholdIsZScore": true,
                      "positiveName": "CD45+",
                      "negativeName": "CD45-",
                      "positiveColor": [0, 200, 0],
                      "negativeColor": [128, 128, 128],
                      "hideOutliers": true,
                      "positiveChildren": [],
                      "negativeChildren": []
                    }
                  ]
                }
                """;
        File file = tempDir.resolve("legacy.json").toFile();
        try (var writer = new BufferedWriter(new FileWriter(file))) {
            writer.write(json);
        }

        GateTree loaded = FlowPathSerializer.load(file);
        assertTrue(loaded.getRoots().get(0).isExcludeOutliers());
    }

    @Test
    void futureVersionThrowsIOException() throws IOException {
        String json = """
                {
                  "version": 999,
                  "gates": []
                }
                """;
        File file = tempDir.resolve("future.json").toFile();
        try (var writer = new BufferedWriter(new FileWriter(file))) {
            writer.write(json);
        }

        assertThrows(IOException.class, () -> FlowPathSerializer.load(file));
    }
}
