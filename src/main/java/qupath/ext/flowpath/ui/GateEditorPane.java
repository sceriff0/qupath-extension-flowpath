package qupath.ext.flowpath.ui;

import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import qupath.ext.flowpath.model.CellIndex;
import qupath.ext.flowpath.model.ColorUtils;
import qupath.ext.flowpath.model.GateNode;
import qupath.ext.flowpath.model.MarkerStats;

import java.util.List;
import java.util.function.Consumer;

/**
 * Right-side editor panel for configuring a single gate node.
 * Shows channel selector, value mode toggle, histogram, threshold slider,
 * positive/negative naming and colors, clip percentiles, and action buttons.
 */
public class GateEditorPane extends VBox {

    private final ComboBox<String> channelCombo;
    private final ToggleGroup modeGroup;
    private final RadioButton rawModeBtn;
    private final RadioButton zscoreModeBtn;
    private final HistogramCanvas histogram;
    private final Slider thresholdSlider;
    private final TextField thresholdValueField;
    private final Label populationCountsLabel;
    private final TextField posNameField;
    private final TextField negNameField;
    private final ColorPicker posColorPicker;
    private final ColorPicker negColorPicker;
    private final Spinner<Double> clipLowSpinner;
    private final Spinner<Double> clipHighSpinner;
    private final CheckBox hideOutliersBox;
    private final Label hoverLabel;

    private final Button addToPosBtn;
    private final Button addToNegBtn;
    private final Button removeGateBtn;

    private GateNode currentNode;
    private CellIndex cellIndex;
    private MarkerStats markerStats;
    private boolean suppressEvents = false;

    private Consumer<GateNode> onNodeChanged;
    private Runnable onAddToPositive;
    private Runnable onAddToNegative;
    private Runnable onRemoveGate;

    public GateEditorPane() {
        setSpacing(8);
        setPadding(new Insets(10));
        setStyle("-fx-background-color: #2a2a2a;");

        // Channel selector
        Label channelLabel = new Label("Channel:");
        channelLabel.setStyle("-fx-text-fill: white;");
        channelCombo = new ComboBox<>();
        channelCombo.setPrefWidth(200);
        // Action listener wired after all fields are initialized (see bottom of constructor)

        HBox channelRow = new HBox(8, channelLabel, channelCombo);

        // Value mode toggle
        modeGroup = new ToggleGroup();
        rawModeBtn = new RadioButton("Raw");
        rawModeBtn.setToggleGroup(modeGroup);
        rawModeBtn.setStyle("-fx-text-fill: white;");
        zscoreModeBtn = new RadioButton("Z-score");
        zscoreModeBtn.setToggleGroup(modeGroup);
        zscoreModeBtn.setSelected(true);
        zscoreModeBtn.setStyle("-fx-text-fill: white;");
        modeGroup.selectedToggleProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setThresholdIsZScore(zscoreModeBtn.isSelected());
                updateHistogram();
                fireNodeChanged();
            }
        });
        HBox modeRow = new HBox(12, new Label("Mode:") {{ setStyle("-fx-text-fill: white;"); }}, rawModeBtn, zscoreModeBtn);

        // Histogram section header
        Label histogramHeader = createSectionHeader("Histogram");

        // Histogram
        histogram = new HistogramCanvas();
        // Threshold callback wired after all fields are initialized (see bottom of constructor)
        hoverLabel = new Label(" ");
        hoverLabel.setStyle("-fx-text-fill: #aaaaaa; -fx-font-size: 9;");
        histogram.setOnMouseHover(val -> hoverLabel.setText(String.format("Value: %.4f", val)));

        // Threshold section header
        Label thresholdHeader = createSectionHeader("Threshold");

        // Threshold slider
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
        // Allow typing exact threshold values
        thresholdValueField.setOnAction(e -> applyThresholdFromField());
        thresholdValueField.focusedProperty().addListener((obs, old, focused) -> {
            if (!focused) applyThresholdFromField();
        });

        HBox threshRow = new HBox(8,
            new Label("Threshold:") {{ setStyle("-fx-text-fill: white;"); }},
            thresholdSlider, thresholdValueField);
        HBox.setHgrow(thresholdSlider, Priority.ALWAYS);

        // Population counts label
        populationCountsLabel = new Label("Positive: -- | Negative: --");
        populationCountsLabel.setStyle("-fx-text-fill: #aaaaaa; -fx-font-size: 10;");

        // Outlier Clipping section header
        Label clipHeader = createSectionHeader("Outlier Clipping");

        // Clip percentiles
        clipLowSpinner = new Spinner<>(0.0, 50.0, 1.0, 0.5);
        clipLowSpinner.setPrefWidth(75);
        clipLowSpinner.setEditable(true);
        clipHighSpinner = new Spinner<>(50.0, 100.0, 99.0, 0.5);
        clipHighSpinner.setPrefWidth(75);
        clipHighSpinner.setEditable(true);
        hideOutliersBox = new CheckBox("Hide outliers");
        hideOutliersBox.setStyle("-fx-text-fill: white;");

        clipLowSpinner.valueProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setClipPercentileLow(val);
                updateHistogram();
                fireNodeChanged();
            }
        });
        clipHighSpinner.valueProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setClipPercentileHigh(val);
                updateHistogram();
                fireNodeChanged();
            }
        });
        hideOutliersBox.selectedProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setHideOutliers(val);
                fireNodeChanged();
            }
        });

        HBox clipRow = new HBox(6,
            new Label("Clip:") {{ setStyle("-fx-text-fill: white;"); }},
            clipLowSpinner, new Label("% to") {{ setStyle("-fx-text-fill: white;"); }},
            clipHighSpinner, new Label("%") {{ setStyle("-fx-text-fill: white;"); }},
            hideOutliersBox);

        // Population Names section header
        Label namesHeader = createSectionHeader("Population Names");

        // Positive / Negative names and colors
        GridPane namesGrid = new GridPane();
        namesGrid.setHgap(6);
        namesGrid.setVgap(4);

        Label posLabel = new Label("Positive:");
        posLabel.setStyle("-fx-text-fill: #00cc00;");
        posNameField = new TextField("Pos");
        posNameField.setPrefWidth(120);
        posColorPicker = new ColorPicker(Color.rgb(0, 200, 0));
        posColorPicker.setPrefWidth(80);

        Label negLabel = new Label("Negative:");
        negLabel.setStyle("-fx-text-fill: #999999;");
        negNameField = new TextField("Neg");
        negNameField.setPrefWidth(120);
        negColorPicker = new ColorPicker(Color.rgb(160, 160, 160));
        negColorPicker.setPrefWidth(80);

        namesGrid.add(posLabel, 0, 0);
        namesGrid.add(posNameField, 1, 0);
        namesGrid.add(posColorPicker, 2, 0);
        namesGrid.add(negLabel, 0, 1);
        namesGrid.add(negNameField, 1, 1);
        namesGrid.add(negColorPicker, 2, 1);

        // Wire name/color changes
        posNameField.textProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                String name = (val == null || val.isBlank()) ? currentNode.getChannel() + "+" : val;
                currentNode.setPositiveName(name);
                fireNodeChanged();
            }
        });
        negNameField.textProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                String name = (val == null || val.isBlank()) ? currentNode.getChannel() + "-" : val;
                currentNode.setNegativeName(name);
                fireNodeChanged();
            }
        });
        posColorPicker.valueProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setPositiveColor(ColorUtils.colorToInt(val));
                histogram.setPosColor(val);
                fireNodeChanged();
            }
        });
        negColorPicker.valueProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setNegativeColor(ColorUtils.colorToInt(val));
                histogram.setNegColor(val);
                fireNodeChanged();
            }
        });

        // Action buttons
        addToPosBtn = new Button("Add Gate to +");
        addToPosBtn.setStyle("-fx-base: #004400;");
        addToNegBtn = new Button("Add Gate to -");
        addToNegBtn.setStyle("-fx-base: #444444;");
        removeGateBtn = new Button("Remove Gate");
        removeGateBtn.setStyle("-fx-base: #440000;");

        addToPosBtn.setOnAction(e -> { if (onAddToPositive != null) onAddToPositive.run(); });
        addToNegBtn.setOnAction(e -> { if (onAddToNegative != null) onAddToNegative.run(); });
        removeGateBtn.setOnAction(e -> { if (onRemoveGate != null) onRemoveGate.run(); });

        HBox buttonRow = new HBox(8, addToPosBtn, addToNegBtn, removeGateBtn);

        // Assemble
        getChildren().addAll(
            channelRow, modeRow,
            histogramHeader, histogram, hoverLabel,
            thresholdHeader, threshRow, populationCountsLabel,
            clipHeader, clipRow,
            new Separator(),
            namesHeader, namesGrid,
            new Separator(),
            buttonRow
        );

        // Wire deferred callbacks (after all fields are initialized)
        channelCombo.setOnAction(e -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setChannel(channelCombo.getValue());
                String ch = channelCombo.getValue();
                if (ch != null) {
                    posNameField.setText(ch + "+");
                    negNameField.setText(ch + "-");
                    currentNode.setPositiveName(ch + "+");
                    currentNode.setNegativeName(ch + "-");
                }
                updateHistogram();
                fireNodeChanged();
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

        setDisabled(true);
    }

    /**
     * Populate the editor with a gate node's current values.
     */
    public void setGateNode(GateNode node) {
        this.currentNode = node;
        if (node == null) {
            withSuppressedEvents(() -> setDisabled(true));
            return;
        }
        withSuppressedEvents(() -> {
            setDisabled(false);
            channelCombo.setValue(node.getChannel());
            if (node.isThresholdIsZScore()) {
                zscoreModeBtn.setSelected(true);
            } else {
                rawModeBtn.setSelected(true);
            }
            thresholdSlider.setValue(node.getThreshold());
            thresholdValueField.setText(String.format("%.4f", node.getThreshold()));
            posNameField.setText(node.getPositiveName());
            negNameField.setText(node.getNegativeName());
            posColorPicker.setValue(ColorUtils.intToColor(node.getPositiveColor()));
            negColorPicker.setValue(ColorUtils.intToColor(node.getNegativeColor()));
            clipLowSpinner.getValueFactory().setValue(node.getClipPercentileLow());
            clipHighSpinner.getValueFactory().setValue(node.getClipPercentileHigh());
            hideOutliersBox.setSelected(node.isHideOutliers());
            histogram.setPosColor(ColorUtils.intToColor(node.getPositiveColor()));
            histogram.setNegColor(ColorUtils.intToColor(node.getNegativeColor()));
        });
        updateHistogram();
    }

    public void setChannelNames(List<String> names) {
        channelCombo.getItems().setAll(names);
    }

    public void setCellIndex(CellIndex index) {
        this.cellIndex = index;
    }

    public void setMarkerStats(MarkerStats stats) {
        this.markerStats = stats;
    }

    public void setOnNodeChanged(Consumer<GateNode> callback) {
        this.onNodeChanged = callback;
    }

    public void setOnAddToPositive(Runnable callback) {
        this.onAddToPositive = callback;
    }

    public void setOnAddToNegative(Runnable callback) {
        this.onAddToNegative = callback;
    }

    public void setOnRemoveGate(Runnable callback) {
        this.onRemoveGate = callback;
    }

    public boolean isUseZScore() {
        return zscoreModeBtn.isSelected();
    }

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

        // Clip range
        double clipLow = currentNode.getClipPercentileLow();
        double clipHigh = currentNode.getClipPercentileHigh();

        // Compute clip values from percentiles
        double[] sorted = displayValues.clone();
        java.util.Arrays.sort(sorted);
        double clipMin = percentile(sorted, clipLow);
        double clipMax = percentile(sorted, clipHigh);

        histogram.setData(displayValues, clipMin, clipMax);
        histogram.setThreshold(currentNode.getThreshold());

        // Update slider range to match display
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
        try {
            action.run();
        } finally {
            suppressEvents = false;
        }
    }

    private void fireNodeChanged() {
        if (onNodeChanged != null && currentNode != null) {
            onNodeChanged.accept(currentNode);
        }
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

    /**
     * Refresh the population counts label from the current GateNode.
     */
    public void updatePopulationCounts() {
        if (currentNode == null) {
            populationCountsLabel.setText("Positive: -- | Negative: --");
            return;
        }
        int pos = currentNode.getPosCount();
        int neg = currentNode.getNegCount();
        int total = pos + neg;
        if (total > 0) {
            double posPct = 100.0 * pos / total;
            double negPct = 100.0 * neg / total;
            populationCountsLabel.setText(String.format(
                "Positive: %,d (%.1f%%) | Negative: %,d (%.1f%%)",
                pos, posPct, neg, negPct));
        } else {
            populationCountsLabel.setText("Positive: 0 (0.0%) | Negative: 0 (0.0%)");
        }
    }

    private static Label createSectionHeader(String text) {
        Label header = new Label(text);
        header.setStyle("-fx-text-fill: #888888; -fx-font-size: 10; -fx-font-weight: bold;");
        header.setPadding(new Insets(4, 0, 0, 0));
        return header;
    }

}
