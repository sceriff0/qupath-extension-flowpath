package qupath.ext.flowpath.ui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.scene.control.*;
import javafx.scene.input.KeyCode;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import qupath.ext.flowpath.engine.GatingEngine;
import qupath.ext.flowpath.engine.LivePreviewService;
import qupath.ext.flowpath.io.FlowPathSerializer;
import qupath.ext.flowpath.io.PhenotypeCsvExporter;
import qupath.ext.flowpath.model.*;
import qupath.lib.gui.QuPathGUI;
import qupath.lib.gui.dialogs.Dialogs;
import qupath.lib.images.ImageData;
import qupath.lib.objects.PathObject;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Main panel for the FlowPath extension.
 * SplitPane: TreeView + QualityFilterPane (left), GateEditorPane (right).
 * Toolbar at bottom for save/load/export.
 */
public class FlowPathPane extends BorderPane {

    private final QuPathGUI qupath;
    private final TreeView<Object> treeView;
    private final GateEditorPane editorPane;
    private final QualityFilterPane qualityFilterPane;
    private final LivePreviewService previewService;
    private final Label statusBar;

    private GateTree gateTree;
    private CellIndex cellIndex;
    private MarkerStats markerStats;
    private List<String> markerNames;
    private boolean[] cachedQualityMask;

    public FlowPathPane(QuPathGUI qupath) {
        this.qupath = qupath;
        this.gateTree = new GateTree();
        this.previewService = new LivePreviewService();

        // --- Left side: TreeView + Quality Filter ---
        treeView = new TreeView<>();
        treeView.setCellFactory(tv -> new FlowPathCell());
        treeView.setShowRoot(false);
        treeView.setRoot(new TreeItem<>("Root"));
        treeView.getSelectionModel().selectedItemProperty().addListener((obs, old, sel) -> onTreeSelectionChanged(sel));
        treeView.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.DELETE || e.getCode() == KeyCode.BACK_SPACE) {
                removeSelectedGate();
                e.consume();
            }
        });

        // Add root gate button
        Button addRootBtn = new Button("+ Add Root Gate");
        addRootBtn.setMaxWidth(Double.MAX_VALUE);
        addRootBtn.setOnAction(e -> addRootGate());

        qualityFilterPane = new QualityFilterPane(gateTree.getQualityFilter());
        qualityFilterPane.setOnFilterChanged(filter -> onQualityFilterChanged());

        VBox leftPane = new VBox(4, treeView, addRootBtn, qualityFilterPane);
        VBox.setVgrow(treeView, Priority.ALWAYS);
        leftPane.setPadding(new Insets(4));
        leftPane.setPrefWidth(280);

        // --- Right side: Gate Editor ---
        editorPane = new GateEditorPane();
        editorPane.setOnNodeChanged(node -> onGateNodeChanged());
        editorPane.setOnAddToPositive(this::addGateToPositive);
        editorPane.setOnAddToNegative(this::addGateToNegative);
        editorPane.setOnRemoveGate(this::removeSelectedGate);

        ScrollPane editorScroll = new ScrollPane(editorPane);
        editorScroll.setFitToWidth(true);
        editorScroll.setPrefWidth(420);

        // --- SplitPane ---
        SplitPane splitPane = new SplitPane(leftPane, editorScroll);
        splitPane.setDividerPositions(0.4);
        setCenter(splitPane);

        // --- Status bar ---
        statusBar = new Label("Total: 0 cells | Excluded: 0 | Gates: 0");
        statusBar.setStyle("-fx-font-size: 11; -fx-text-fill: #aaaaaa; -fx-padding: 2 6 2 6;");

        // --- Bottom toolbar ---
        Button saveBtn = new Button("Save JSON");
        saveBtn.setOnAction(e -> saveTree());
        Button loadBtn = new Button("Load JSON");
        loadBtn.setOnAction(e -> loadTree());
        Button exportBtn = new Button("Export CSV");
        exportBtn.setOnAction(e -> exportCsv());

        HBox toolbar = new HBox(8, saveBtn, loadBtn, new Separator(Orientation.VERTICAL), exportBtn);
        toolbar.setPadding(new Insets(6));

        VBox bottomBox = new VBox(statusBar, toolbar);
        setBottom(bottomBox);

        // Style
        setStyle("-fx-background-color: #1e1e1e;");

        // Initialize from current image
        Platform.runLater(this::initializeFromImage);

        // Listen for image changes
        qupath.imageDataProperty().addListener((obs, oldImg, newImg) -> {
            Platform.runLater(this::initializeFromImage);
        });
    }

    /**
     * Build CellIndex and MarkerStats from the currently loaded image's detections.
     */
    private void initializeFromImage() {
        ImageData<?> imageData = qupath.getImageData();
        if (imageData == null) {
            cellIndex = null;
            markerStats = null;
            markerNames = Collections.emptyList();
            editorPane.setChannelNames(markerNames);
            return;
        }

        Collection<PathObject> detections = imageData.getHierarchy().getDetectionObjects();
        if (detections.isEmpty()) {
            Dialogs.showWarningNotification("FlowPath", "No detections found. Import GeoJSON cells first.");
            return;
        }

        // Discover marker names from image channels (OME-TIFF metadata)
        markerNames = discoverMarkerNames(imageData, detections);

        cellIndex = CellIndex.build(detections, markerNames);

        // Compute quality mask and stats
        recomputeQualityMask();
        markerStats = MarkerStats.compute(cellIndex, cachedQualityMask);

        // Update UI
        editorPane.setChannelNames(markerNames);
        editorPane.setCellIndex(cellIndex);
        editorPane.setMarkerStats(markerStats);

        // Update quality filter ranges + check which QC metrics have data
        double maxArea = 0, maxTotalInt = 0;
        boolean hasEccentricity = false, hasSolidity = false;
        for (int i = 0; i < cellIndex.size(); i++) {
            maxArea = Math.max(maxArea, cellIndex.getArea(i));
            maxTotalInt = Math.max(maxTotalInt, cellIndex.getTotalIntensity(i));
            if (cellIndex.getEccentricity(i) > 0) hasEccentricity = true;
            if (cellIndex.getSolidity(i) > 0 && cellIndex.getSolidity(i) < 1.0) hasSolidity = true;
        }
        qualityFilterPane.updateRanges(maxArea, maxTotalInt);
        qualityFilterPane.setAvailableMetrics(hasEccentricity, hasSolidity);
        updateFilteredCount();

        // Setup preview service
        previewService.setCellIndex(cellIndex);
        previewService.setMarkerStats(markerStats);
        previewService.setGateTree(gateTree);
        previewService.setImageData(imageData);
        previewService.setUseZScore(editorPane.isUseZScore());
        previewService.setOnUpdateComplete(this::onPreviewUpdated);
        previewService.setOnStatsRecomputed(() -> {
            editorPane.setMarkerStats(previewService.getMarkerStats());
            recomputeQualityMask();
            updateFilteredCount();
        });

        updateStatusBar();
    }

    /**
     * Discover marker names from the image's channel metadata.
     * Falls back to detection measurements if no image channels are found.
     */
    private List<String> discoverMarkerNames(ImageData<?> imageData, Collection<PathObject> detections) {
        PathObject sample = detections.iterator().next();
        var measurements = sample.getMeasurements();
        Set<String> measurementKeys = (measurements != null) ? measurements.keySet() : Set.of();

        // Primary: get channel names from image metadata, validate against measurements
        var server = imageData.getServer();
        var channels = server.getMetadata().getChannels();
        if (channels != null && !channels.isEmpty()) {
            List<String> validated = new ArrayList<>();
            for (var ch : channels) {
                String name = ch.getName();
                if (name == null || name.isEmpty()) continue;
                if (hasMeasurement(measurementKeys, name)) {
                    validated.add(name);
                }
            }
            if (!validated.isEmpty()) return validated;
        }

        // Fallback: extract from detection measurements (minus morphology fields)
        Set<String> exclude = Set.of(
            "Centroid X µm", "Centroid Y µm",
            "area µm²", "eccentricity", "perimeter", "convex_area",
            "axis_major_length", "axis_minor_length",
            "x", "y", "label", "fov", "cell_size"
        );

        if (measurementKeys.isEmpty()) return Collections.emptyList();

        return measurementKeys.stream()
            .filter(name -> !exclude.contains(name))
            .filter(name -> !name.startsWith("["))
            .filter(name -> !name.startsWith("_"))
            .sorted()
            .collect(Collectors.toList());
    }

    private boolean hasMeasurement(Set<String> keys, String channel) {
        if (keys.contains(channel)) return true;
        String suffix = "] " + channel;
        for (String key : keys) {
            if (key.endsWith(suffix)) return true;
        }
        return false;
    }

    // --- Tree building ---

    private void rebuildTreeView() {
        TreeItem<Object> root = new TreeItem<>("Root");
        for (GateNode gate : gateTree.getRoots()) {
            root.getChildren().add(buildTreeItem(gate));
        }
        treeView.setRoot(root);
        root.setExpanded(true);
        expandAll(root);
    }

    private TreeItem<Object> buildTreeItem(GateNode gate) {
        TreeItem<Object> gateItem = new TreeItem<>(gate);
        gateItem.setExpanded(true);

        // Positive branch
        TreeItem<Object> posItem = new TreeItem<>(new FlowPathCell.BranchItem(gate, true));
        posItem.setExpanded(true);
        for (GateNode child : gate.getPositiveChildren()) {
            posItem.getChildren().add(buildTreeItem(child));
        }
        gateItem.getChildren().add(posItem);

        // Negative branch
        TreeItem<Object> negItem = new TreeItem<>(new FlowPathCell.BranchItem(gate, false));
        negItem.setExpanded(true);
        for (GateNode child : gate.getNegativeChildren()) {
            negItem.getChildren().add(buildTreeItem(child));
        }
        gateItem.getChildren().add(negItem);

        return gateItem;
    }

    private void expandAll(TreeItem<?> item) {
        item.setExpanded(true);
        for (TreeItem<?> child : item.getChildren()) {
            expandAll(child);
        }
    }

    // --- Gate operations ---

    private void addRootGate() {
        if (markerNames == null || markerNames.isEmpty()) {
            Dialogs.showWarningNotification("FlowPath", "No markers available. Load an image with detections first.");
            return;
        }

        String channel = markerNames.get(0);
        GateNode node = new GateNode(channel);
        gateTree.addRoot(node);
        rebuildTreeView();
        requestPreviewUpdate();
    }

    private void addGateToPositive() {
        GateNode selected = getSelectedGateNode();
        if (selected == null || markerNames == null || markerNames.isEmpty()) return;

        String channel = markerNames.get(0);
        GateNode child = new GateNode(channel);
        selected.getPositiveChildren().add(child);
        rebuildTreeView();
        requestPreviewUpdate();
    }

    private void addGateToNegative() {
        GateNode selected = getSelectedGateNode();
        if (selected == null || markerNames == null || markerNames.isEmpty()) return;

        String channel = markerNames.get(0);
        GateNode child = new GateNode(channel);
        selected.getNegativeChildren().add(child);
        rebuildTreeView();
        requestPreviewUpdate();
    }

    private void removeSelectedGate() {
        GateNode selected = getSelectedGateNode();
        if (selected == null) return;

        boolean hasChildren = !selected.getPositiveChildren().isEmpty() || !selected.getNegativeChildren().isEmpty();
        if (hasChildren) {
            boolean confirm = Dialogs.showConfirmDialog("Remove Gate",
                "This gate has child gates. Remove entire subtree?");
            if (!confirm) return;
        }

        // Remove from parent
        if (!gateTree.getRoots().remove(selected)) {
            removeFromTree(gateTree.getRoots(), selected);
        }

        editorPane.setGateNode(null);
        rebuildTreeView();
        requestPreviewUpdate();
    }

    private boolean removeFromTree(List<GateNode> nodes, GateNode target) {
        for (GateNode node : nodes) {
            if (node.getPositiveChildren().remove(target)) return true;
            if (node.getNegativeChildren().remove(target)) return true;
            if (removeFromTree(node.getPositiveChildren(), target)) return true;
            if (removeFromTree(node.getNegativeChildren(), target)) return true;
        }
        return false;
    }

    // --- Selection handling ---

    private void onTreeSelectionChanged(TreeItem<Object> selected) {
        if (selected == null) {
            editorPane.setGateNode(null);
            return;
        }

        Object item = selected.getValue();
        if (item instanceof GateNode node) {
            editorPane.setGateNode(node);
        } else if (item instanceof FlowPathCell.BranchItem branch) {
            editorPane.setGateNode(branch.parentGate);
        } else {
            editorPane.setGateNode(null);
        }
    }

    private GateNode getSelectedGateNode() {
        TreeItem<Object> sel = treeView.getSelectionModel().getSelectedItem();
        if (sel == null) return null;
        Object item = sel.getValue();
        if (item instanceof GateNode node) return node;
        if (item instanceof FlowPathCell.BranchItem branch) return branch.parentGate;
        return null;
    }

    // --- Updates ---

    private void onGateNodeChanged() {
        previewService.setUseZScore(editorPane.isUseZScore());
        requestPreviewUpdate();
    }

    private void onQualityFilterChanged() {
        if (cellIndex == null) return;
        recomputeQualityMask();
        // Recompute stats on background thread, then trigger preview update
        previewService.recomputeStats();
        updateFilteredCount();
    }

    private void recomputeQualityMask() {
        if (cellIndex == null) {
            cachedQualityMask = null;
            return;
        }
        cachedQualityMask = GatingEngine.computeQualityMask(cellIndex, gateTree.getQualityFilter());
    }

    private void updateFilteredCount() {
        if (cellIndex == null || cachedQualityMask == null) return;
        int passing = 0;
        for (boolean b : cachedQualityMask) if (b) passing++;
        int filtered = cellIndex.size() - passing;
        qualityFilterPane.setFilteredCount(filtered, cellIndex.size());
    }

    private void requestPreviewUpdate() {
        previewService.setGateTree(gateTree);
        previewService.requestUpdate();
    }

    private void onPreviewUpdated() {
        // Refresh tree cell counts and status bar
        Platform.runLater(() -> {
            treeView.refresh();
            updateStatusBar();
        });
    }

    private void updateStatusBar() {
        int total = cellIndex != null ? cellIndex.size() : 0;
        int excluded = 0;
        if (cachedQualityMask != null) {
            for (boolean b : cachedQualityMask) if (!b) excluded++;
        }
        int gateCount = countGates(gateTree.getRoots());
        statusBar.setText(String.format("Total: %,d cells | Excluded: %,d | Gates: %d", total, excluded, gateCount));
    }

    private int countGates(List<GateNode> nodes) {
        int count = 0;
        for (GateNode node : nodes) {
            count++;
            count += countGates(node.getPositiveChildren());
            count += countGates(node.getNegativeChildren());
        }
        return count;
    }

    // --- IO ---

    private void saveTree() {
        File file = Dialogs.promptToSaveFile("Save FlowPath", null, "flowpath.json", "JSON", ".json");
        if (file == null) return;
        try {
            FlowPathSerializer.save(gateTree, file);
            Dialogs.showInfoNotification("FlowPath", "Saved to " + file.getName());
        } catch (Exception ex) {
            Dialogs.showErrorMessage("Save Error", ex.getMessage());
        }
    }

    private void loadTree() {
        File file = Dialogs.promptForFile("Load FlowPath", null, "JSON", ".json");
        if (file == null) return;
        try {
            gateTree = FlowPathSerializer.load(file);
            // Sync the quality filter pane to the new filter object
            qualityFilterPane.setFilter(gateTree.getQualityFilter());
            qualityFilterPane.setOnFilterChanged(filter -> onQualityFilterChanged());

            rebuildTreeView();
            onQualityFilterChanged();
            Dialogs.showInfoNotification("FlowPath", "Loaded from " + file.getName());
        } catch (Exception ex) {
            Dialogs.showErrorMessage("Load Error", ex.getMessage());
        }
    }

    private void exportCsv() {
        if (cellIndex == null || markerStats == null || gateTree.getRoots().isEmpty()) {
            Dialogs.showWarningNotification("FlowPath", "No gates defined or no cells loaded.");
            return;
        }

        File file = Dialogs.promptToSaveFile("Export Phenotypes", null, "gate_pheno.csv", "CSV", ".csv");
        if (file == null) return;

        try {
            GatingEngine.AssignmentResult result = GatingEngine.assignAll(
                gateTree, cellIndex, markerStats, editorPane.isUseZScore());
            PhenotypeCsvExporter.export(file, cellIndex, result, gateTree, markerStats);
            Dialogs.showInfoNotification("FlowPath", "Exported " + file.getName());
        } catch (Exception ex) {
            Dialogs.showErrorMessage("Export Error", ex.getMessage());
        }
    }

    /**
     * Clean up resources when the window is closed.
     */
    public void shutdown() {
        previewService.shutdown();
    }
}
