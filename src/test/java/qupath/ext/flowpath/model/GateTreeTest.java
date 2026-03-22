package qupath.ext.flowpath.model;

import org.junit.jupiter.api.Test;
import java.util.List;
import static org.junit.jupiter.api.Assertions.*;

class GateTreeTest {

    @Test
    void addAndRemoveRoot() {
        var tree = new GateTree();
        var node = new GateNode("CD45");
        tree.addRoot(node);
        assertEquals(1, tree.getRoots().size());
        assertSame(node, tree.getRoots().get(0));

        tree.removeRoot(node);
        assertTrue(tree.getRoots().isEmpty());
    }

    @Test
    void deepCopyCreatesIndependentTree() {
        var tree = new GateTree();
        var root = new GateNode("CD45", 1.0);
        root.getPositiveChildren().add(new GateNode("CD3", 0.5));
        tree.addRoot(root);

        var copy = tree.deepCopy();

        // Modify original
        root.setThreshold(99.0);
        tree.addRoot(new GateNode("PANCK"));

        // Copy should be unaffected
        assertEquals(1, copy.getRoots().size());
        assertEquals(1.0, copy.getRoots().get(0).getThreshold());
    }

    @Test
    void transferCountsWalksInParallel() {
        // Build two identical trees
        var orig = new GateNode("CD45");
        orig.getPositiveChildren().add(new GateNode("CD3"));

        var copyNode = new GateNode("CD45");
        var copyChild = new GateNode("CD3");
        copyNode.getPositiveChildren().add(copyChild);

        // Set counts on copy
        copyNode.setPosCount(100);
        copyNode.setNegCount(200);
        copyChild.setPosCount(30);
        copyChild.setNegCount(70);

        GateTree.transferCounts(List.of(orig), List.of(copyNode));

        assertEquals(100, orig.getPosCount());
        assertEquals(200, orig.getNegCount());
        assertEquals(30, orig.getPositiveChildren().get(0).getPosCount());
        assertEquals(70, orig.getPositiveChildren().get(0).getNegCount());
    }

    @Test
    void transferCountsHandlesMismatchedSizes() {
        var orig = List.of(new GateNode("CD45"));
        var copies = List.of(new GateNode("CD45"), new GateNode("CD3"));

        // Should not throw — just processes the shorter list
        assertDoesNotThrow(() -> GateTree.transferCounts(orig, copies));
    }

    @Test
    void collectLeafNamesAcrossMultipleRoots() {
        var tree = new GateTree();
        tree.addRoot(new GateNode("CD45"));
        tree.addRoot(new GateNode("PANCK"));

        var names = tree.collectLeafNames();
        assertEquals(List.of("CD45+", "CD45-", "PANCK+", "PANCK-"), names);
    }

    @Test
    void defaultQualityFilterIsPresent() {
        var tree = new GateTree();
        assertNotNull(tree.getQualityFilter());
    }
}
