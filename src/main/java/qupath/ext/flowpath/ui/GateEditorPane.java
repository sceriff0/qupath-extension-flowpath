package qupath.ext.flowpath.ui;

import javafx.geometry.Insets;
import javafx.scene.Node;
import javafx.scene.control.*;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import qupath.ext.flowpath.model.BooleanGate;
import qupath.ext.flowpath.model.Branch;
import qupath.ext.flowpath.model.CellIndex;
import qupath.ext.flowpath.model.ColorUtils;
import qupath.ext.flowpath.model.EllipseGate;
import qupath.ext.flowpath.model.GateNode;
import qupath.ext.flowpath.model.MarkerStats;
import qupath.ext.flowpath.model.PolygonGate;
import qupath.ext.flowpath.model.QuadrantGate;
import qupath.ext.flowpath.model.RectangleGate;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;
import java.util.function.IntConsumer;

/**
 * Right-side editor panel for configuring a single gate node.
 * Swaps controls based on gate type: threshold shows histogram + slider,
 * quadrant shows 2-channel controls, boolean shows operation picker, etc.
 */
public class GateEditorPane extends VBox {

    // --- Shared controls ---
    private final Label gateTypeLabel;
    private final Spinner<Double> clipLowSpinner;
    private final Spinner<Double> clipHighSpinner;
    private final CheckBox excludeOutliersBox;
    private final VBox gateSpecificArea;
    private final VBox branchNamesArea;
    private final VBox actionButtonArea;

    // --- Threshold-specific controls (lazily shown) ---
    private final ComboBox<String> channelCombo;
    private final ToggleGroup modeGroup;
    private final RadioButton rawModeBtn;
    private final RadioButton zscoreModeBtn;
    private final HistogramCanvas histogram;
    private final Slider thresholdSlider;
    private final TextField thresholdValueField;
    private final Label populationCountsLabel;
    private final Label hoverLabel;

    private GateNode currentNode;
    private CellIndex cellIndex;
    private MarkerStats markerStats;
    private boolean suppressEvents = false;

    private Consumer<GateNode> onNodeChanged;
    private Runnable onAddToPositive;
    private Runnable onAddToNegative;
    private IntConsumer onAddToBranch;
    private Runnable onRemoveGate;

    public GateEditorPane() {
        setSpacing(8);
        setPadding(new Insets(10));
        setStyle("-fx-background-color: #2a2a2a;");

        // Gate type indicator
        gateTypeLabel = new Label("No gate selected");
        gateTypeLabel.setStyle("-fx-text-fill: #80b0d0; -fx-font-size: 11; -fx-font-weight: bold;");

        // --- Threshold-specific controls (always created, shown/hidden as needed) ---
        channelCombo = new ComboBox<>();
        channelCombo.setPrefWidth(200);
        channelCombo.setTooltip(new Tooltip("Select the marker channel for this gate"));

        modeGroup = new ToggleGroup();
        rawModeBtn = new RadioButton("Raw");
        rawModeBtn.setToggleGroup(modeGroup);
        rawModeBtn.setStyle("-fx-text-fill: white;");
        rawModeBtn.setTooltip(new Tooltip("Compare marker raw intensity values against threshold"));
        zscoreModeBtn = new RadioButton("Z-score");
        zscoreModeBtn.setToggleGroup(modeGroup);
        zscoreModeBtn.setSelected(true);
        zscoreModeBtn.setStyle("-fx-text-fill: white;");
        zscoreModeBtn.setTooltip(new Tooltip("Compare z-score normalized values against threshold (recommended)"));
        modeGroup.selectedToggleProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setThresholdIsZScore(zscoreModeBtn.isSelected());
                updateHistogram();
                fireNodeChanged();
            }
        });

        histogram = new HistogramCanvas();
        hoverLabel = new Label(" ");
        hoverLabel.setStyle("-fx-text-fill: #aaaaaa; -fx-font-size: 9;");
        histogram.setOnMouseHover(val -> hoverLabel.setText(String.format("Value: %.4f", val)));

        thresholdSlider = new Slider(-5, 5, 0);
        thresholdSlider.setPrefWidth(300);
        thresholdSlider.setBlockIncrement(0.01);
        thresholdValueField = new TextField("0.0000");
        thresholdValueField.setPrefWidth(80);
        thresholdValueField.setStyle("-fx-text-fill: white; -fx-font-family: monospace; -fx-background-color: #3a3a3a;");
        thresholdSlider.valueProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setThreshold(val.doubleValue());
                thresholdValueField.setText(String.format("%.4f", val.doubleValue()));
                histogram.setThreshold(val.doubleValue());
                fireNodeChanged();
                updatePopulationCounts();
            }
        });
        thresholdValueField.setOnAction(e -> applyThresholdFromField());
        thresholdValueField.focusedProperty().addListener((obs, old, focused) -> {
            if (!focused) applyThresholdFromField();
        });

        populationCountsLabel = new Label("Positive: -- | Negative: --");
        populationCountsLabel.setStyle("-fx-text-fill: #aaaaaa; -fx-font-size: 10;");

        // --- Shared: Outlier Clipping ---
        clipLowSpinner = new Spinner<>(0.0, 50.0, 1.0, 0.5);
        clipLowSpinner.setPrefWidth(75);
        clipLowSpinner.setEditable(true);
        clipHighSpinner = new Spinner<>(50.0, 100.0, 99.0, 0.5);
        clipHighSpinner.setPrefWidth(75);
        clipHighSpinner.setEditable(true);
        excludeOutliersBox = new CheckBox("Exclude outliers");
        excludeOutliersBox.setStyle("-fx-text-fill: white;");
        excludeOutliersBox.setTooltip(new Tooltip(
            "When enabled, cells with marker values outside the clip percentile range\n" +
            "are classified as 'Excluded' and removed from gating and CSV export."));

        clipLowSpinner.valueProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                double clamped = Math.min(val, clipHighSpinner.getValue() - 0.5);
                if (clamped != val) { clipLowSpinner.getValueFactory().setValue(clamped); return; }
                currentNode.setClipPercentileLow(val);
                updateHistogram();
                fireNodeChanged();
            }
        });
        clipHighSpinner.valueProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                double clamped = Math.max(val, clipLowSpinner.getValue() + 0.5);
                if (clamped != val) { clipHighSpinner.getValueFactory().setValue(clamped); return; }
                currentNode.setClipPercentileHigh(val);
                updateHistogram();
                fireNodeChanged();
            }
        });
        excludeOutliersBox.selectedProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setExcludeOutliers(val);
                fireNodeChanged();
            }
        });

        HBox clipRow = new HBox(6,
            new Label("Clip:") {{ setStyle("-fx-text-fill: white;"); }},
            clipLowSpinner, new Label("% to") {{ setStyle("-fx-text-fill: white;"); }},
            clipHighSpinner, new Label("%") {{ setStyle("-fx-text-fill: white;"); }},
            excludeOutliersBox);

        // Swappable areas
        gateSpecificArea = new VBox(4);
        branchNamesArea = new VBox(4);
        actionButtonArea = new VBox(4);

        // Wire deferred callbacks
        channelCombo.setOnAction(e -> {
            if (!suppressEvents && currentNode != null) {
                String oldChannel = currentNode.getChannel();
                String newChannel = channelCombo.getValue();
                currentNode.setChannel(newChannel);
                // Only auto-rename branches if they still match the old default pattern
                if (newChannel != null && currentNode.getBranches().size() >= 2) {
                    Branch pos = currentNode.getBranches().get(0);
                    Branch neg = currentNode.getBranches().get(1);
                    if (oldChannel != null && pos.getName().equals(oldChannel + "+")) {
                        pos.setName(newChannel + "+");
                    }
                    if (oldChannel != null && neg.getName().equals(oldChannel + "-")) {
                        neg.setName(newChannel + "-");
                    }
                }
                updateHistogram();
                fireNodeChanged();
                // Refresh branch names in editor
                buildBranchNamesEditor(currentNode);
            }
        });
        histogram.setOnThresholdChanged(val -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setThreshold(val);
                thresholdSlider.setValue(val);
                thresholdValueField.setText(String.format("%.4f", val));
                fireNodeChanged();
                updatePopulationCounts();
            }
        });

        // Assemble
        getChildren().addAll(
            gateTypeLabel,
            gateSpecificArea,
            createSectionHeader("Outlier Clipping"), clipRow,
            new Separator(),
            branchNamesArea,
            new Separator(),
            actionButtonArea
        );

        setDisabled(true);
    }

    /**
     * Populate the editor with a gate node's current values.
     * Rebuilds the gate-specific UI section based on gate type.
     */
    public void setGateNode(GateNode node) {
        this.currentNode = node;
        if (node == null) {
            withSuppressedEvents(() -> setDisabled(true));
            gateTypeLabel.setText("No gate selected");
            gateSpecificArea.getChildren().clear();
            Label hint = new Label("Select a gate from the tree to edit it,\nor click '+ Add Root Gate' to create one.");
            hint.setStyle("-fx-text-fill: #888888; -fx-font-size: 11;");
            hint.setWrapText(true);
            gateSpecificArea.getChildren().add(hint);
            branchNamesArea.getChildren().clear();
            actionButtonArea.getChildren().clear();
            return;
        }
        withSuppressedEvents(() -> {
            setDisabled(false);

            // Shared controls
            clipLowSpinner.getValueFactory().setValue(node.getClipPercentileLow());
            clipHighSpinner.getValueFactory().setValue(node.getClipPercentileHigh());
            excludeOutliersBox.setSelected(node.isExcludeOutliers());

            // Gate type label
            String typeDisplay = switch (node.getGateType()) {
                case "threshold" -> "Threshold Gate";
                case "quadrant" -> "Quadrant Gate";
                case "boolean" -> "Boolean Gate";
                case "polygon" -> "Polygon Gate";
                case "rectangle" -> "Rectangle Gate";
                case "ellipse" -> "Ellipse Gate";
                default -> "Gate";
            };
            gateTypeLabel.setText(typeDisplay);

            // Rebuild gate-specific area
            gateSpecificArea.getChildren().clear();
            if (node instanceof QuadrantGate qg) {
                buildQuadrantEditor(qg);
            } else if (node instanceof BooleanGate bg) {
                buildBooleanEditor(bg);
            } else if (node instanceof PolygonGate || node instanceof RectangleGate || node instanceof EllipseGate) {
                build2DEditor(node);
            } else {
                buildThresholdEditor(node);
            }

            // Rebuild branch names/colors
            buildBranchNamesEditor(node);

            // Rebuild action buttons
            buildActionButtons(node);
        });

        if (!(node instanceof BooleanGate)) {
            updateHistogram();
        }
    }

    // ---- Gate-type-specific editor builders ----

    private void buildThresholdEditor(GateNode node) {
        Label chLabel = new Label("Channel:");
        chLabel.setStyle("-fx-text-fill: white;");
        HBox channelRow = new HBox(8, chLabel, channelCombo);
        channelCombo.setValue(node.getChannel());

        HBox modeRow = new HBox(12, new Label("Mode:") {{ setStyle("-fx-text-fill: white;"); }}, rawModeBtn, zscoreModeBtn);
        if (node.isThresholdIsZScore()) zscoreModeBtn.setSelected(true);
        else rawModeBtn.setSelected(true);

        thresholdSlider.setValue(node.getThreshold());
        thresholdValueField.setText(String.format("%.4f", node.getThreshold()));
        histogram.setPosColor(ColorUtils.intToColor(node.getPositiveColor()));
        histogram.setNegColor(ColorUtils.intToColor(node.getNegativeColor()));

        HBox threshRow = new HBox(8,
            new Label("Threshold:") {{ setStyle("-fx-text-fill: white;"); }},
            thresholdSlider, thresholdValueField);
        HBox.setHgrow(thresholdSlider, Priority.ALWAYS);

        gateSpecificArea.getChildren().addAll(
            channelRow, modeRow,
            createSectionHeader("Histogram"), histogram, hoverLabel,
            createSectionHeader("Threshold"), threshRow, populationCountsLabel
        );
    }

    private void buildQuadrantEditor(QuadrantGate gate) {
        Label chXLabel = new Label("Channel X:");
        chXLabel.setStyle("-fx-text-fill: white;");
        ComboBox<String> chXCombo = new ComboBox<>(channelCombo.getItems());
        chXCombo.setValue(gate.getChannelX());
        chXCombo.setPrefWidth(150);
        chXCombo.setOnAction(e -> { if (!suppressEvents) { gate.setChannelX(chXCombo.getValue()); fireNodeChanged(); }});

        Label chYLabel = new Label("Channel Y:");
        chYLabel.setStyle("-fx-text-fill: white;");
        ComboBox<String> chYCombo = new ComboBox<>(channelCombo.getItems());
        chYCombo.setValue(gate.getChannelY());
        chYCombo.setPrefWidth(150);
        chYCombo.setOnAction(e -> { if (!suppressEvents) { gate.setChannelY(chYCombo.getValue()); fireNodeChanged(); }});

        Slider sliderX = new Slider(-5, 5, gate.getThresholdX());
        sliderX.setBlockIncrement(0.01);
        Label valX = new Label(String.format("%.3f", gate.getThresholdX()));
        valX.setStyle("-fx-text-fill: white; -fx-font-family: monospace;");
        sliderX.valueProperty().addListener((obs, old, val) -> {
            if (!suppressEvents) { gate.setThresholdX(val.doubleValue()); valX.setText(String.format("%.3f", val.doubleValue())); fireNodeChanged(); }
        });

        Slider sliderY = new Slider(-5, 5, gate.getThresholdY());
        sliderY.setBlockIncrement(0.01);
        Label valY = new Label(String.format("%.3f", gate.getThresholdY()));
        valY.setStyle("-fx-text-fill: white; -fx-font-family: monospace;");
        sliderY.valueProperty().addListener((obs, old, val) -> {
            if (!suppressEvents) { gate.setThresholdY(val.doubleValue()); valY.setText(String.format("%.3f", val.doubleValue())); fireNodeChanged(); }
        });

        HBox modeRow = new HBox(12, new Label("Mode:") {{ setStyle("-fx-text-fill: white;"); }}, rawModeBtn, zscoreModeBtn);
        if (gate.isThresholdIsZScore()) zscoreModeBtn.setSelected(true);
        else rawModeBtn.setSelected(true);

        gateSpecificArea.getChildren().addAll(
            new HBox(8, chXLabel, chXCombo), new HBox(8, chYLabel, chYCombo),
            modeRow,
            createSectionHeader("Threshold X"), new HBox(8, sliderX, valX),
            createSectionHeader("Threshold Y"), new HBox(8, sliderY, valY)
        );
    }

    private void buildBooleanEditor(BooleanGate gate) {
        Label opLabel = new Label("Operation:");
        opLabel.setStyle("-fx-text-fill: white;");
        ComboBox<String> opCombo = new ComboBox<>();
        opCombo.getItems().addAll("AND", "OR", "NOT");
        opCombo.setValue(gate.getOperation().name());
        opCombo.setOnAction(e -> {
            if (!suppressEvents) {
                gate.setOperation(BooleanGate.Op.valueOf(opCombo.getValue()));
                fireNodeChanged();
            }
        });

        Label info = new Label("Boolean gates combine other gates.\nAdd operand gates as children of the 'Match' branch.");
        info.setStyle("-fx-text-fill: #aaaaaa; -fx-font-size: 10;");
        info.setWrapText(true);

        Label operandCount = new Label("Operands: " + gate.getOperands().size());
        operandCount.setStyle("-fx-text-fill: white;");

        gateSpecificArea.getChildren().addAll(
            new HBox(8, opLabel, opCombo), info, operandCount
        );
    }

    private void build2DEditor(GateNode node) {
        String chX = node.getChannels().isEmpty() ? "" : node.getChannels().get(0);
        String chY = node.getChannels().size() > 1 ? node.getChannels().get(1) : "";

        Label info = new Label("2D " + node.getGateType() + " gate on " + chX + " vs " + chY);
        info.setStyle("-fx-text-fill: white;");

        // Show scatter plot if we have data
        if (cellIndex != null && !chX.isEmpty() && !chY.isEmpty()) {
            int mxIdx = cellIndex.getMarkerIndex(chX);
            int myIdx = cellIndex.getMarkerIndex(chY);
            if (mxIdx >= 0 && myIdx >= 0) {
                ScatterPlotCanvas scatter = new ScatterPlotCanvas();
                scatter.setData(cellIndex.getMarkerValues(mxIdx), cellIndex.getMarkerValues(myIdx), chX, chY);

                if (node instanceof PolygonGate pg && pg.getVertices().size() >= 3) {
                    scatter.setPolygonOverlay(pg.getVertices());
                } else if (node instanceof RectangleGate rg) {
                    scatter.setRectangleOverlay(rg.getMinX(), rg.getMaxX(), rg.getMinY(), rg.getMaxY());
                } else if (node instanceof EllipseGate eg) {
                    scatter.setEllipseOverlay(eg.getCenterX(), eg.getCenterY(), eg.getRadiusX(), eg.getRadiusY());
                }

                // Drawing toolbar
                ToggleGroup toolGroup = new ToggleGroup();
                ToggleButton polygonBtn = new ToggleButton("Polygon");
                polygonBtn.setToggleGroup(toolGroup);
                ToggleButton rectBtn = new ToggleButton("Rectangle");
                rectBtn.setToggleGroup(toolGroup);
                ToggleButton ellipseBtn = new ToggleButton("Ellipse");
                ellipseBtn.setToggleGroup(toolGroup);

                HBox drawToolbar = new HBox(4, polygonBtn, rectBtn, ellipseBtn);

                // Auto-select based on gate type
                if (node instanceof PolygonGate) polygonBtn.setSelected(true);
                else if (node instanceof RectangleGate) rectBtn.setSelected(true);
                else if (node instanceof EllipseGate) ellipseBtn.setSelected(true);

                // Wire toolbar to drawing mode
                toolGroup.selectedToggleProperty().addListener((obs, old, val) -> {
                    if (val == polygonBtn) scatter.setDrawingMode(ScatterPlotCanvas.DrawingMode.POLYGON);
                    else if (val == rectBtn) scatter.setDrawingMode(ScatterPlotCanvas.DrawingMode.RECTANGLE);
                    else if (val == ellipseBtn) scatter.setDrawingMode(ScatterPlotCanvas.DrawingMode.ELLIPSE);
                    else scatter.setDrawingMode(ScatterPlotCanvas.DrawingMode.NONE);
                });

                // Wire callbacks to update gate model
                scatter.setOnPolygonDrawn(vertices -> {
                    if (node instanceof PolygonGate pg) {
                        pg.setVertices(new ArrayList<>(vertices));
                        scatter.setPolygonOverlay(pg.getVertices());
                        fireNodeChanged();
                    }
                });
                scatter.setOnRectangleDrawn(bounds -> {
                    if (node instanceof RectangleGate rg) {
                        rg.setMinX(bounds[0]); rg.setMaxX(bounds[1]);
                        rg.setMinY(bounds[2]); rg.setMaxY(bounds[3]);
                        scatter.setRectangleOverlay(rg.getMinX(), rg.getMaxX(), rg.getMinY(), rg.getMaxY());
                        fireNodeChanged();
                    }
                });
                scatter.setOnEllipseDrawn(params -> {
                    if (node instanceof EllipseGate eg) {
                        eg.setCenterX(params[0]); eg.setCenterY(params[1]);
                        eg.setRadiusX(params[2]); eg.setRadiusY(params[3]);
                        scatter.setEllipseOverlay(eg.getCenterX(), eg.getCenterY(), eg.getRadiusX(), eg.getRadiusY());
                        fireNodeChanged();
                    }
                });

                // Set initial drawing mode based on gate type
                if (node instanceof PolygonGate) scatter.setDrawingMode(ScatterPlotCanvas.DrawingMode.POLYGON);
                else if (node instanceof RectangleGate) scatter.setDrawingMode(ScatterPlotCanvas.DrawingMode.RECTANGLE);
                else if (node instanceof EllipseGate) scatter.setDrawingMode(ScatterPlotCanvas.DrawingMode.ELLIPSE);

                gateSpecificArea.getChildren().addAll(info, drawToolbar, scatter);
                return;
            }
        }

        Label noData = new Label("Load an image to see the scatter plot");
        noData.setStyle("-fx-text-fill: #888888;");
        gateSpecificArea.getChildren().addAll(info, noData);
    }

    // ---- Branch names/colors editor (generic for any gate type) ----

    private void buildBranchNamesEditor(GateNode node) {
        branchNamesArea.getChildren().clear();
        branchNamesArea.getChildren().add(createSectionHeader("Branch Names & Colors"));

        GridPane grid = new GridPane();
        grid.setHgap(6);
        grid.setVgap(4);

        List<Branch> branches = node.getBranches();
        for (int i = 0; i < branches.size(); i++) {
            Branch branch = branches.get(i);
            int idx = i;

            // Contextual branch label based on gate type and index
            String labelText;
            if (node instanceof QuadrantGate) {
                labelText = new String[]{"Q1 (++):", "Q2 (-+):", "Q3 (+-):", "Q4 (--):"} [Math.min(i, 3)];
            } else if (node instanceof BooleanGate) {
                labelText = i == 0 ? "Match:" : "No Match:";
            } else if (node instanceof PolygonGate || node instanceof RectangleGate || node instanceof EllipseGate) {
                labelText = i == 0 ? "Inside:" : "Outside:";
            } else {
                labelText = i == 0 ? "Positive:" : "Negative:";
            }
            Color labelColor = ColorUtils.intToColor(branch.getColor());
            Label label = new Label(labelText);
            label.setStyle("-fx-text-fill: " + toWebColor(labelColor) + ";");

            TextField nameField = new TextField(branch.getName());
            nameField.setPrefWidth(120);
            nameField.textProperty().addListener((obs, old, val) -> {
                if (!suppressEvents && val != null && !val.isBlank()) {
                    branch.setName(val);
                    fireNodeChanged();
                }
            });

            ColorPicker colorPicker = new ColorPicker(ColorUtils.intToColor(branch.getColor()));
            colorPicker.setPrefWidth(80);
            colorPicker.valueProperty().addListener((obs, old, val) -> {
                if (!suppressEvents) {
                    branch.setColor(ColorUtils.colorToInt(val));
                    fireNodeChanged();
                }
            });

            Label countLabel = new Label(String.format("%,d", branch.getCount()));
            countLabel.setStyle("-fx-text-fill: #aaaaaa; -fx-font-size: 10;");

            grid.add(label, 0, i);
            grid.add(nameField, 1, i);
            grid.add(colorPicker, 2, i);
            grid.add(countLabel, 3, i);
        }

        branchNamesArea.getChildren().add(grid);
    }

    // ---- Action buttons (generic for any gate type) ----

    private void buildActionButtons(GateNode node) {
        actionButtonArea.getChildren().clear();

        List<Branch> branches = node.getBranches();
        HBox buttonRow = new HBox(8);

        // Add "Add child gate to [branch]" button for each branch
        for (int i = 0; i < branches.size(); i++) {
            Branch branch = branches.get(i);
            int branchIdx = i;
            Button addBtn = new Button("+ " + branch.getName());
            addBtn.setStyle("-fx-base: #003300;");
            addBtn.setTooltip(new Tooltip("Add a child gate to '" + branch.getName() + "'"));
            addBtn.setOnAction(e -> {
                if (onAddToBranch != null) onAddToBranch.accept(branchIdx);
            });
            buttonRow.getChildren().add(addBtn);
        }

        Button removeBtn = new Button("Remove Gate");
        removeBtn.setStyle("-fx-base: #440000;");
        removeBtn.setTooltip(new Tooltip("Remove this gate and all its children (Del)"));
        removeBtn.setOnAction(e -> { if (onRemoveGate != null) onRemoveGate.run(); });
        buttonRow.getChildren().add(removeBtn);

        actionButtonArea.getChildren().add(buttonRow);
    }

    // ---- Public API ----

    public void setChannelNames(List<String> names) { channelCombo.getItems().setAll(names); }
    public void setCellIndex(CellIndex index) { this.cellIndex = index; }
    public void setMarkerStats(MarkerStats stats) { this.markerStats = stats; }
    public void setOnNodeChanged(Consumer<GateNode> callback) { this.onNodeChanged = callback; }
    public void setOnAddToPositive(Runnable callback) { this.onAddToPositive = callback; }
    public void setOnAddToNegative(Runnable callback) { this.onAddToNegative = callback; }
    public void setOnAddToBranch(IntConsumer callback) { this.onAddToBranch = callback; }
    public void setOnRemoveGate(Runnable callback) { this.onRemoveGate = callback; }

    public boolean isUseZScore() { return zscoreModeBtn.isSelected(); }

    public void updatePopulationCounts() {
        if (currentNode == null) {
            populationCountsLabel.setText("Positive: -- | Negative: --");
            return;
        }
        List<Branch> branches = currentNode.getBranches();
        if (branches.size() == 2) {
            int pos = branches.get(0).getCount();
            int neg = branches.get(1).getCount();
            int total = pos + neg;
            if (total > 0) {
                populationCountsLabel.setText(String.format(
                    "%s: %,d (%.1f%%) | %s: %,d (%.1f%%)",
                    branches.get(0).getName(), pos, 100.0 * pos / total,
                    branches.get(1).getName(), neg, 100.0 * neg / total));
            } else {
                populationCountsLabel.setText(branches.get(0).getName() + ": 0 | " + branches.get(1).getName() + ": 0");
            }
        }
    }

    // ---- Internal ----

    private void updateHistogram() {
        if (currentNode == null || cellIndex == null || markerStats == null) return;
        String channel = currentNode.getChannel();
        if (channel == null) return;
        int markerIdx = cellIndex.getMarkerIndex(channel);
        if (markerIdx < 0) return;

        double[] rawValues = cellIndex.getMarkerValues(markerIdx);
        boolean useZ = currentNode.isThresholdIsZScore();

        double[] displayValues;
        if (useZ && markerStats.getStd(channel) > 1e-10) {
            displayValues = new double[rawValues.length];
            double mean = markerStats.getMean(channel);
            double std = markerStats.getStd(channel);
            for (int i = 0; i < rawValues.length; i++) {
                displayValues[i] = (rawValues[i] - mean) / std;
            }
        } else {
            displayValues = rawValues;
        }

        double[] sorted = displayValues.clone();
        java.util.Arrays.sort(sorted);
        double clipMin = percentile(sorted, currentNode.getClipPercentileLow());
        double clipMax = percentile(sorted, currentNode.getClipPercentileHigh());

        histogram.setData(displayValues, clipMin, clipMax);
        histogram.setThreshold(currentNode.getThreshold());
        thresholdSlider.setMin(clipMin);
        thresholdSlider.setMax(clipMax);
        updatePopulationCounts();
    }

    private double percentile(double[] sorted, double pct) {
        if (sorted.length == 0) return 0;
        double idx = (pct / 100.0) * (sorted.length - 1);
        int lo = (int) Math.floor(idx);
        int hi = (int) Math.ceil(idx);
        if (lo == hi || hi >= sorted.length) return sorted[Math.min(lo, sorted.length - 1)];
        double frac = idx - lo;
        return sorted[lo] * (1 - frac) + sorted[hi] * frac;
    }

    private void withSuppressedEvents(Runnable action) {
        suppressEvents = true;
        try { action.run(); } finally { suppressEvents = false; }
    }

    private void fireNodeChanged() {
        if (onNodeChanged != null && currentNode != null) onNodeChanged.accept(currentNode);
    }

    private void applyThresholdFromField() {
        if (suppressEvents || currentNode == null) return;
        try {
            double val = Double.parseDouble(thresholdValueField.getText().trim());
            withSuppressedEvents(() -> {
                currentNode.setThreshold(val);
                thresholdSlider.setValue(val);
                histogram.setThreshold(val);
            });
            fireNodeChanged();
            updatePopulationCounts();
        } catch (NumberFormatException ex) {
            thresholdValueField.setText(String.format("%.4f", currentNode.getThreshold()));
        }
    }

    private static String toWebColor(Color c) {
        return String.format("#%02x%02x%02x",
            (int)(c.getRed() * 255), (int)(c.getGreen() * 255), (int)(c.getBlue() * 255));
    }

    private static Label createSectionHeader(String text) {
        Label header = new Label(text);
        header.setStyle("-fx-text-fill: #888888; -fx-font-size: 10; -fx-font-weight: bold;");
        header.setPadding(new Insets(4, 0, 0, 0));
        return header;
    }
}
