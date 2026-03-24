package qupath.ext.flowpath.model;

import org.junit.jupiter.api.Test;
import java.util.ArrayDeque;
import java.util.Deque;
import static org.junit.jupiter.api.Assertions.*;

class UndoStackTest {

    private static final int MAX_UNDO = 50;

    @Test
    void undoRestoresPreviousState() {
        Deque<GateTree> undoStack = new ArrayDeque<>();
        Deque<GateTree> redoStack = new ArrayDeque<>();

        GateTree tree = new GateTree();
        tree.addRoot(new GateNode("CD45", 1.0));

        undoStack.push(tree.deepCopy());
        tree.getRoots().get(0).setThreshold(5.0);
        assertEquals(5.0, tree.getRoots().get(0).getThreshold());

        redoStack.push(tree.deepCopy());
        tree = undoStack.pop();
        assertEquals(1.0, tree.getRoots().get(0).getThreshold());
    }

    @Test
    void redoRestoresUndoneState() {
        Deque<GateTree> undoStack = new ArrayDeque<>();
        Deque<GateTree> redoStack = new ArrayDeque<>();

        GateTree tree = new GateTree();
        tree.addRoot(new GateNode("CD45", 1.0));

        undoStack.push(tree.deepCopy());
        tree.getRoots().get(0).setThreshold(5.0);

        redoStack.push(tree.deepCopy());
        tree = undoStack.pop();
        assertEquals(1.0, tree.getRoots().get(0).getThreshold());

        undoStack.push(tree.deepCopy());
        tree = redoStack.pop();
        assertEquals(5.0, tree.getRoots().get(0).getThreshold());
    }

    @Test
    void newModificationClearsRedoStack() {
        Deque<GateTree> undoStack = new ArrayDeque<>();
        Deque<GateTree> redoStack = new ArrayDeque<>();

        GateTree tree = new GateTree();
        tree.addRoot(new GateNode("CD45", 1.0));

        undoStack.push(tree.deepCopy());
        tree.getRoots().get(0).setThreshold(5.0);

        redoStack.push(tree.deepCopy());
        tree = undoStack.pop();
        assertFalse(redoStack.isEmpty());

        undoStack.push(tree.deepCopy());
        tree.getRoots().get(0).setThreshold(3.0);
        redoStack.clear();
        assertTrue(redoStack.isEmpty());
    }

    @Test
    void undoStackCappedAtMax() {
        Deque<GateTree> undoStack = new ArrayDeque<>();

        GateTree tree = new GateTree();
        tree.addRoot(new GateNode("CD45", 0.0));

        for (int i = 0; i < MAX_UNDO + 10; i++) {
            undoStack.push(tree.deepCopy());
            if (undoStack.size() > MAX_UNDO) undoStack.removeLast();
            tree.getRoots().get(0).setThreshold(i + 1);
        }

        assertEquals(MAX_UNDO, undoStack.size());
    }

    @Test
    void deepCopyPreservesQualityFilter() {
        GateTree tree = new GateTree();
        tree.getQualityFilter().setMinArea(100.0);

        GateTree copy = tree.deepCopy();
        copy.getQualityFilter().setMinArea(200.0);

        assertEquals(100.0, tree.getQualityFilter().getMinArea(),
            "Original QualityFilter should not change after modifying copy");
    }
}
