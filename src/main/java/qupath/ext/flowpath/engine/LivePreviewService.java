package qupath.ext.flowpath.engine;

import javafx.animation.PauseTransition;
import javafx.application.Platform;
import javafx.util.Duration;
import qupath.ext.flowpath.model.CellIndex;
import qupath.ext.flowpath.model.GateTree;
import qupath.ext.flowpath.model.MarkerStats;
import qupath.ext.flowpath.model.QualityFilter;
import qupath.ext.flowpath.model.ColorUtils;
import qupath.lib.common.ColorTools;
import qupath.lib.images.ImageData;
import qupath.lib.objects.PathObject;
import qupath.lib.objects.classes.PathClass;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Debounced live-preview service that re-runs gating on a background thread
 * whenever gate parameters change and pushes the results back onto the
 * JavaFX Application Thread so the QuPath viewer updates.
 */
public class LivePreviewService {

    private static final long DEBOUNCE_MS = 80;

    private final PauseTransition debounce;
    private final ExecutorService executor;

    private volatile GateTree gateTree;
    private volatile CellIndex cellIndex;
    private volatile MarkerStats markerStats;
    private volatile boolean useZScore;
    private volatile ImageData<?> imageData;
    private volatile boolean[] roiMask;

    /** Optional callback fired after MarkerStats is recomputed (e.g., to refresh UI sliders). */
    private Runnable onStatsRecomputed;

    /** Optional callback fired when gating computation starts (e.g., to show a spinner). */
    private Runnable onUpdateStarted;

    /** Optional callback fired after cell classifications are applied (e.g., to refresh tree counts). */
    private Runnable onUpdateComplete;

    /** Count of excluded cells from the most recent gating run (updated on FX thread). */
    private int lastExcludedCount;

    public LivePreviewService() {
        this.executor = Executors.newSingleThreadExecutor(r -> {
            Thread t = new Thread(r, "flowpath-preview");
            t.setDaemon(true);
            return t;
        });
        this.debounce = new PauseTransition(Duration.millis(DEBOUNCE_MS));
        this.debounce.setOnFinished(e -> submitGatingWork());
    }

    // ---- setters ----

    public void setGateTree(GateTree tree) {
        this.gateTree = tree;
    }

    public void setCellIndex(CellIndex index) {
        this.cellIndex = index;
    }

    public void setMarkerStats(MarkerStats stats) {
        this.markerStats = stats;
    }

    public MarkerStats getMarkerStats() {
        return this.markerStats;
    }

    public void setUseZScore(boolean useZScore) {
        this.useZScore = useZScore;
    }

    public void setImageData(ImageData<?> imageData) {
        this.imageData = imageData;
    }

    public void setRoiMask(boolean[] roiMask) {
        this.roiMask = roiMask;
    }

    public void setOnStatsRecomputed(Runnable onStatsRecomputed) {
        this.onStatsRecomputed = onStatsRecomputed;
    }

    public void setOnUpdateStarted(Runnable onUpdateStarted) {
        this.onUpdateStarted = onUpdateStarted;
    }

    public void setOnUpdateComplete(Runnable onUpdateComplete) {
        this.onUpdateComplete = onUpdateComplete;
    }

    public int getLastExcludedCount() {
        return lastExcludedCount;
    }

    // ---- public API ----

    /**
     * Request a gating update. The actual computation is debounced by 80 ms so
     * rapid successive calls (e.g., slider drags) are coalesced.
     */
    public void requestUpdate() {
        if (Platform.isFxApplicationThread()) {
            debounce.playFromStart();
        } else {
            Platform.runLater(debounce::playFromStart);
        }
    }

    /**
     * Recompute {@link MarkerStats} from the current {@link CellIndex} using
     * the quality mask derived from the gate tree's quality filter, then trigger
     * a gating update.
     */
    public void recomputeStats() {
        if (cellIndex == null || gateTree == null || executor.isShutdown()) {
            return;
        }
        final boolean[] roi = this.roiMask;
        final CellIndex idx = this.cellIndex;
        final QualityFilter qf = gateTree.getQualityFilter().deepCopy();
        executor.submit(() -> {
            boolean[] qualityMask = GatingEngine.computeQualityMask(idx, qf);
            boolean[] mask = roi != null ? GatingEngine.combineMasks(qualityMask, roi) : qualityMask;
            MarkerStats recomputed = MarkerStats.compute(cellIndex, mask);
            this.markerStats = recomputed;
            if (onStatsRecomputed != null) {
                Platform.runLater(onStatsRecomputed);
            }
            requestUpdate();
        });
    }

    /**
     * Shut down the background executor. Call this when the extension window
     * is closed or the image changes.
     */
    public void shutdown() {
        debounce.stop();
        executor.shutdown();
    }

    // ---- internal ----

    private void submitGatingWork() {
        if (executor.isShutdown()) return;
        // Capture references to avoid races
        final GateTree originalTree = this.gateTree;
        final CellIndex index = this.cellIndex;
        final MarkerStats stats = this.markerStats;
        final boolean zScore = this.useZScore;
        final ImageData<?> data = this.imageData;
        final boolean[] roi = this.roiMask;

        if (originalTree == null || index == null || stats == null || data == null) {
            return;
        }

        // Deep-copy the tree so the background thread works on an immutable snapshot
        final GateTree tree = originalTree.deepCopy();

        if (onUpdateStarted != null) {
            Platform.runLater(onUpdateStarted);
        }

        executor.submit(() -> {
            GatingEngine.AssignmentResult result = GatingEngine.assignAll(tree, index, stats, zScore, roi);

            Platform.runLater(() -> {
                // Transfer counts from the snapshot back to the live tree for UI display
                GateTree.transferCounts(originalTree.getRoots(), tree.getRoots());
                applyResult(result, index, data);
            });
        });
    }

    private void applyResult(GatingEngine.AssignmentResult result, CellIndex index, ImageData<?> data) {
        String[] phenotypes = result.getPhenotypes();
        boolean[] excluded = result.getExcluded();
        int[] colors = result.getColors();
        int n = phenotypes.length;

        // Count total excluded cells for status display
        int excCount = 0;
        for (boolean ex : excluded) if (ex) excCount++;
        this.lastExcludedCount = excCount;

        Map<String, PathClass> classCache = new HashMap<>();

        // Build cache and force-update colors (PathClass.fromString caches globally by name,
        // so we must explicitly set the color each time to reflect user changes)
        Map<String, Integer> colorByName = new HashMap<>();
        for (int i = 0; i < n; i++) {
            if (!excluded[i] && phenotypes[i] != null) {
                colorByName.put(phenotypes[i], colors[i]);
            }
        }
        for (var entry : colorByName.entrySet()) {
            int packed = entry.getValue();
            int qupathColor = ColorUtils.toQuPathColor(packed);
            PathClass pc = PathClass.fromString(entry.getKey(), qupathColor);
            pc.setColor(qupathColor);  // Force-update cached PathClass color
            classCache.put(entry.getKey(), pc);
        }

        // Near-invisible PathClass for excluded cells (avoids red "Unclassified" default)
        int excludedColor = ColorTools.packRGB(20, 20, 20);
        PathClass excludedClass = PathClass.fromString("Excluded", excludedColor);
        excludedClass.setColor(excludedColor);

        for (int i = 0; i < n; i++) {
            PathObject obj = index.getObject(i);
            if (obj == null) {
                continue;
            }
            if (excluded[i]) {
                obj.setPathClass(excludedClass);
            } else {
                obj.setPathClass(classCache.get(phenotypes[i]));
            }
        }

        data.getHierarchy().fireHierarchyChangedEvent(this);

        if (onUpdateComplete != null) {
            onUpdateComplete.run();
        }
    }
}
