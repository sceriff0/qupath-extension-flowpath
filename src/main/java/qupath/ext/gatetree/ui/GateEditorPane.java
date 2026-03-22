package qupath.ext.gatetree.ui;

import javafx.geometry.Insets;
import javafx.scene.control.*;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import qupath.ext.gatetree.model.CellIndex;
import qupath.ext.gatetree.model.GateNode;
import qupath.ext.gatetree.model.MarkerStats;

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
    private final Label thresholdValueLabel;
    private final TextField posNameField;
    private final TextField negNameField;
    private final ColorPicker posColorPicker;
    private final ColorPicker negColorPicker;
    private final Spinner<Double> clipLowSpinner;
    private final Spinner<Double> clipHighSpinner;
    private final CheckBox excludeOutliersBox;
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

        // Histogram
        histogram = new HistogramCanvas();
        // Threshold callback wired after all fields are initialized (see bottom of constructor)
        hoverLabel = new Label(" ");
        hoverLabel.setStyle("-fx-text-fill: #aaaaaa; -fx-font-size: 9;");
        histogram.setOnMouseHover(val -> hoverLabel.setText(String.format("Value: %.4f", val)));

        // Threshold slider
        thresholdSlider = new Slider(-5, 5, 0);
        thresholdSlider.setPrefWidth(300);
        thresholdSlider.setBlockIncrement(0.01);
        thresholdValueLabel = new Label("0.0000");
        thresholdValueLabel.setStyle("-fx-text-fill: white; -fx-font-family: monospace;");
        thresholdSlider.valueProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setThreshold(val.doubleValue());
                thresholdValueLabel.setText(String.format("%.4f", val.doubleValue()));
                histogram.setThreshold(val.doubleValue());
                fireNodeChanged();
            }
        });
        HBox threshRow = new HBox(8,
            new Label("Threshold:") {{ setStyle("-fx-text-fill: white;"); }},
            thresholdSlider, thresholdValueLabel);
        HBox.setHgrow(thresholdSlider, Priority.ALWAYS);

        // Clip percentiles
        clipLowSpinner = new Spinner<>(0.0, 50.0, 1.0, 0.5);
        clipLowSpinner.setPrefWidth(75);
        clipLowSpinner.setEditable(true);
        clipHighSpinner = new Spinner<>(50.0, 100.0, 99.0, 0.5);
        clipHighSpinner.setPrefWidth(75);
        clipHighSpinner.setEditable(true);
        excludeOutliersBox = new CheckBox("Exclude outliers");
        excludeOutliersBox.setStyle("-fx-text-fill: white;");

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

        // Positive / Negative names and colors
        GridPane namesGrid = new GridPane();
        namesGrid.setHgap(6);
        namesGrid.setVgap(4);

        Label posLabel = new Label("Positive:");
        posLabel.setStyle("-fx-text-fill: #00cc00;");
        posNameField = new TextField("Pos");
        posNameField.setPrefWidth(120);
        posColorPicker = new ColorPicker(Color.rgb(0, 200, 0));
        posColorPicker.setPrefWidth(60);

        Label negLabel = new Label("Negative:");
        negLabel.setStyle("-fx-text-fill: #999999;");
        negNameField = new TextField("Neg");
        negNameField.setPrefWidth(120);
        negColorPicker = new ColorPicker(Color.rgb(160, 160, 160));
        negColorPicker.setPrefWidth(60);

        namesGrid.add(posLabel, 0, 0);
        namesGrid.add(posNameField, 1, 0);
        namesGrid.add(posColorPicker, 2, 0);
        namesGrid.add(negLabel, 0, 1);
        namesGrid.add(negNameField, 1, 1);
        namesGrid.add(negColorPicker, 2, 1);

        // Wire name/color changes
        posNameField.textProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setPositiveName(val);
                fireNodeChanged();
            }
        });
        negNameField.textProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setNegativeName(val);
                fireNodeChanged();
            }
        });
        posColorPicker.valueProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setPositiveColor(colorToInt(val));
                histogram.setPosColor(val);
                fireNodeChanged();
            }
        });
        negColorPicker.valueProperty().addListener((obs, old, val) -> {
            if (!suppressEvents && currentNode != null) {
                currentNode.setNegativeColor(colorToInt(val));
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
            histogram, hoverLabel,
            threshRow, clipRow,
            new Separator(),
            namesGrid,
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
                thresholdValueLabel.setText(String.format("%.4f", val));
                fireNodeChanged();
            }
        });

        setDisabled(true);
    }

    /**
     * Populate the editor with a gate node's current values.
     */
    public void setGateNode(GateNode node) {
        this.currentNode = node;
        suppressEvents = true;

        if (node == null) {
            setDisabled(true);
            suppressEvents = false;
            return;
        }

        setDisabled(false);

        channelCombo.setValue(node.getChannel());
        if (node.isThresholdIsZScore()) {
            zscoreModeBtn.setSelected(true);
        } else {
            rawModeBtn.setSelected(true);
        }

        thresholdSlider.setValue(node.getThreshold());
        thresholdValueLabel.setText(String.format("%.4f", node.getThreshold()));

        posNameField.setText(node.getPositiveName());
        negNameField.setText(node.getNegativeName());
        posColorPicker.setValue(intToColor(node.getPositiveColor()));
        negColorPicker.setValue(intToColor(node.getNegativeColor()));

        clipLowSpinner.getValueFactory().setValue(node.getClipPercentileLow());
        clipHighSpinner.getValueFactory().setValue(node.getClipPercentileHigh());
        excludeOutliersBox.setSelected(node.isExcludeOutliers());

        histogram.setPosColor(intToColor(node.getPositiveColor()));
        histogram.setNegColor(intToColor(node.getNegativeColor()));

        suppressEvents = false;
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

    private void fireNodeChanged() {
        if (onNodeChanged != null && currentNode != null) {
            onNodeChanged.accept(currentNode);
        }
    }

    static int colorToInt(Color c) {
        int r = (int) (c.getRed() * 255);
        int g = (int) (c.getGreen() * 255);
        int b = (int) (c.getBlue() * 255);
        return (r << 16) | (g << 8) | b;
    }

    static Color intToColor(int packed) {
        int r = (packed >> 16) & 0xFF;
        int g = (packed >> 8) & 0xFF;
        int b = packed & 0xFF;
        return Color.rgb(r, g, b);
    }
}
