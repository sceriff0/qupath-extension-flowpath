package qupath.ext.gatetree.engine;

import javafx.animation.PauseTransition;
import javafx.application.Platform;
import javafx.util.Duration;
import qupath.ext.gatetree.model.CellIndex;
import qupath.ext.gatetree.model.GateTree;
import qupath.ext.gatetree.model.MarkerStats;
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

    /** Optional callback fired after MarkerStats is recomputed (e.g., to refresh UI sliders). */
    private Runnable onStatsRecomputed;

    /** Optional callback fired after cell classifications are applied (e.g., to refresh tree counts). */
    private Runnable onUpdateComplete;

    public LivePreviewService() {
        this.executor = Executors.newSingleThreadExecutor(r -> {
            Thread t = new Thread(r, "gatetree-preview");
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

    public void setOnStatsRecomputed(Runnable onStatsRecomputed) {
        this.onStatsRecomputed = onStatsRecomputed;
    }

    public void setOnUpdateComplete(Runnable onUpdateComplete) {
        this.onUpdateComplete = onUpdateComplete;
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
        if (cellIndex == null || gateTree == null) {
            return;
        }
        executor.submit(() -> {
            boolean[] mask = GatingEngine.computeQualityMask(cellIndex, gateTree.getQualityFilter());
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
        executor.shutdownNow();
    }

    // ---- internal ----

    private void submitGatingWork() {
        // Capture references to avoid races
        final GateTree tree = this.gateTree;
        final CellIndex index = this.cellIndex;
        final MarkerStats stats = this.markerStats;
        final boolean zScore = this.useZScore;
        final ImageData<?> data = this.imageData;

        if (tree == null || index == null || stats == null || data == null) {
            return;
        }

        executor.submit(() -> {
            GatingEngine.AssignmentResult result = GatingEngine.assignAll(tree, index, stats, zScore);

            Platform.runLater(() -> applyResult(result, index, data));
        });
    }

    private void applyResult(GatingEngine.AssignmentResult result, CellIndex index, ImageData<?> data) {
        String[] phenotypes = result.getPhenotypes();
        boolean[] excluded = result.getExcluded();
        int[] colors = result.getColors();
        int n = phenotypes.length;

        Map<String, PathClass> classCache = new HashMap<>();

        for (int i = 0; i < n; i++) {
            PathObject obj = index.getObject(i);
            if (obj == null) {
                continue;
            }
            if (excluded[i]) {
                obj.setPathClass(null);
            } else {
                String name = phenotypes[i];
                int packedColor = colors[i];
                PathClass pc = classCache.computeIfAbsent(name, k -> {
                    int r = ColorTools.red(packedColor);
                    int g = ColorTools.green(packedColor);
                    int b = ColorTools.blue(packedColor);
                    return PathClass.fromString(k, ColorTools.packRGB(r, g, b));
                });
                obj.setPathClass(pc);
            }
        }

        data.getHierarchy().fireHierarchyChangedEvent(this);

        if (onUpdateComplete != null) {
            onUpdateComplete.run();
        }
    }
}
