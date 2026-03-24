package qupath.ext.flowpath.ui;

import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.geometry.Orientation;
import javafx.scene.control.*;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyCodeCombination;
import javafx.scene.input.KeyCombination;
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
import qupath.lib.objects.hierarchy.events.PathObjectHierarchyEvent;
import qupath.lib.objects.hierarchy.events.PathObjectHierarchyListener;
import qupath.lib.roi.interfaces.ROI;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Main panel for the FlowPath extension.
 * SplitPane: TreeView + QualityFilterPane (left), GateEditorPane (right).
 * Toolbar at bottom for save/load/export.
 */
public class FlowPathPane extends BorderPane {

    /** Wrapper for ROI filter combo box items — stores UUID for matching, displays name. */
    private record AnnotationItem(String id, String displayName) {
        @Override public String toString() { return displayName; }
    }

    private static final AnnotationItem ALL_CELLS = new AnnotationItem(null, "All cells");

    private final QuPathGUI qupath;
    private final TreeView<Object> treeView;
    private final GateEditorPane editorPane;
    private final QualityFilterPane qualityFilterPane;
    private final ComboBox<AnnotationItem> roiComboBox;
    private final LivePreviewService previewService;
    private final Label statusBar;

    private GateTree gateTree;
    private CellIndex cellIndex;
    private MarkerStats markerStats;
    private List<String> markerNames;
    private boolean[] cachedQualityMask;
    private boolean[] cachedRoiMask;
    private PathObjectHierarchyListener hierarchyListener;
    private ImageData<?> listenerImageData;

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
        addRootBtn.setTooltip(new Tooltip("Add a new top-level gate to the gating hierarchy"));

        // ROI filter
        roiComboBox = new ComboBox<>();
        roiComboBox.getItems().add(ALL_CELLS);
        roiComboBox.getSelectionModel().select(0);
        roiComboBox.setMaxWidth(Double.MAX_VALUE);
        roiComboBox.setOnAction(e -> onRoiSelectionChanged());
        roiComboBox.setStyle("-fx-text-fill: black;");
        roiComboBox.setCellFactory(lv -> new ListCell<>() {
            @Override protected void updateItem(AnnotationItem item, boolean empty) {
                super.updateItem(item, empty);
                setText(empty || item == null ? null : item.displayName());
                setStyle("-fx-text-fill: black;");
            }
        });
        roiComboBox.setButtonCell(new ListCell<>() {
            @Override protected void updateItem(AnnotationItem item, boolean empty) {
                super.updateItem(item, empty);
                setText(empty || item == null ? null : item.displayName());
                setStyle("-fx-text-fill: black;");
            }
        });
        Label roiLabel = new Label("ROI Filter:");
        roiLabel.setStyle("-fx-text-fill: black; -fx-font-size: 10;");

        qualityFilterPane = new QualityFilterPane(gateTree.getQualityFilter());
        qualityFilterPane.setOnFilterChanged(filter -> onQualityFilterChanged());

        VBox leftPane = new VBox(4, treeView, addRootBtn, roiLabel, roiComboBox, qualityFilterPane);
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
        saveBtn.setTooltip(new Tooltip("Save gate tree to JSON file (Ctrl+S)"));
        Button loadBtn = new Button("Load JSON");
        loadBtn.setOnAction(e -> loadTree());
        loadBtn.setTooltip(new Tooltip("Load gate tree from JSON file (Ctrl+O)"));
        Button exportBtn = new Button("Export CSV");
        exportBtn.setOnAction(e -> exportCsv());
        exportBtn.setTooltip(new Tooltip("Export phenotype assignments to CSV (Ctrl+E)"));

        HBox toolbar = new HBox(8, saveBtn, loadBtn, new Separator(Orientation.VERTICAL), exportBtn);
        toolbar.setPadding(new Insets(6));

        VBox bottomBox = new VBox(statusBar, toolbar);
        setBottom(bottomBox);

        // --- Keyboard shortcuts ---
        setOnKeyPressed(e -> {
            if (new KeyCodeCombination(KeyCode.S, KeyCombination.SHORTCUT_DOWN).match(e)) {
                saveTree(); e.consume();
            } else if (new KeyCodeCombination(KeyCode.O, KeyCombination.SHORTCUT_DOWN).match(e)) {
                loadTree(); e.consume();
            } else if (new KeyCodeCombination(KeyCode.E, KeyCombination.SHORTCUT_DOWN).match(e)) {
                exportCsv(); e.consume();
            }
        });

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
        // Remove old hierarchy listener from previous image
        if (hierarchyListener != null && listenerImageData != null) {
            listenerImageData.getHierarchy().removeListener(hierarchyListener);
            hierarchyListener = null;
            listenerImageData = null;
        }

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

        // Refresh annotation list and compute ROI mask
        refreshAnnotationList();
        recomputeRoiMask();

        // Compute quality mask and stats (using combined mask)
        recomputeQualityMask();
        markerStats = MarkerStats.compute(cellIndex, getCombinedMask());

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

        // Setup preview service
        previewService.setCellIndex(cellIndex);
        previewService.setMarkerStats(markerStats);
        previewService.setRoiMask(cachedRoiMask);
        previewService.setGateTree(gateTree);
        previewService.setImageData(imageData);
        previewService.setUseZScore(editorPane.isUseZScore());
        previewService.setOnUpdateComplete(this::onPreviewUpdated);
        previewService.setOnStatsRecomputed(() -> {
            editorPane.setMarkerStats(previewService.getMarkerStats());
            recomputeQualityMask();
        });

        // Listen for annotation changes (add/remove) to refresh ROI list
        hierarchyListener = event -> {
            if (!event.isChanging()) {
                Platform.runLater(this::refreshAnnotationList);
            }
        };
        listenerImageData = imageData;
        imageData.getHierarchy().addListener(hierarchyListener);

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

    // --- ROI filtering ---

    private void refreshAnnotationList() {
        ImageData<?> imageData = qupath.getImageData();
        AnnotationItem currentSelection = roiComboBox.getValue();
        String currentId = (currentSelection != null) ? currentSelection.id() : null;
        roiComboBox.getItems().clear();
        roiComboBox.getItems().add(ALL_CELLS);

        if (imageData != null) {
            for (PathObject ann : imageData.getHierarchy().getAnnotationObjects()) {
                String name = ann.getDisplayedName();
                if (name != null && !name.isEmpty()) {
                    roiComboBox.getItems().add(new AnnotationItem(ann.getID().toString(), name));
                }
            }
        }

        // Restore previous selection by UUID if still available
        if (currentId != null) {
            for (AnnotationItem item : roiComboBox.getItems()) {
                if (currentId.equals(item.id())) {
                    roiComboBox.setValue(item);
                    return;
                }
            }
        }
        roiComboBox.getSelectionModel().select(0);
    }

    private void recomputeRoiMask() {
        if (cellIndex == null) {
            cachedRoiMask = null;
            return;
        }

        String roiId = gateTree.getRoiFilterId();
        if (roiId == null || roiId.isEmpty()) {
            cachedRoiMask = null;
            previewService.setRoiMask(null);
            return;
        }

        ImageData<?> imageData = qupath.getImageData();
        if (imageData == null) {
            cachedRoiMask = null;
            previewService.setRoiMask(null);
            return;
        }

        ROI roi = null;
        for (PathObject ann : imageData.getHierarchy().getAnnotationObjects()) {
            if (roiId.equals(ann.getID().toString())) {
                roi = ann.getROI();
                break;
            }
        }
        cachedRoiMask = GatingEngine.computeRoiMask(cellIndex, roi);
        previewService.setRoiMask(cachedRoiMask);
    }

    private boolean[] getCombinedMask() {
        if (cachedQualityMask == null) return cachedRoiMask;
        if (cachedRoiMask == null) return cachedQualityMask;
        return GatingEngine.combineMasks(cachedQualityMask, cachedRoiMask);
    }

    private void onRoiSelectionChanged() {
        AnnotationItem selected = roiComboBox.getValue();
        gateTree.setRoiFilterId(selected != null ? selected.id() : null);
        recomputeRoiMask();
        // Recompute stats with new combined mask, then trigger preview
        previewService.recomputeStats();
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
    }

    private void recomputeQualityMask() {
        if (cellIndex == null) {
            cachedQualityMask = null;
            return;
        }
        cachedQualityMask = GatingEngine.computeQualityMask(cellIndex, gateTree.getQualityFilter());
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
        int excluded = previewService.getLastExcludedCount();
        int gateCount = countGates(gateTree.getRoots());
        AnnotationItem roiItem = roiComboBox.getValue();
        String roiInfo = (roiItem != null && roiItem.id() != null)
            ? " | ROI: " + roiItem.displayName()
            : "";
        statusBar.setText(String.format("Total: %,d cells | Excluded: %,d | Gates: %d%s",
            total, excluded, gateCount, roiInfo));
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

            // Restore ROI selection by UUID
            String roiId = gateTree.getRoiFilterId();
            boolean restored = false;
            if (roiId != null) {
                for (AnnotationItem item : roiComboBox.getItems()) {
                    if (roiId.equals(item.id())) {
                        roiComboBox.setValue(item);
                        restored = true;
                        break;
                    }
                }
            }
            if (!restored) {
                roiComboBox.getSelectionModel().select(0);
                gateTree.setRoiFilterId(null);
            }
            recomputeRoiMask();

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
                gateTree, cellIndex, markerStats, editorPane.isUseZScore(), cachedRoiMask);
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
