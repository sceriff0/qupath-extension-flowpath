package qupath.ext.flowpath.ui;

import javafx.geometry.Insets;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.Slider;
import javafx.scene.control.TitledPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.VBox;
import qupath.ext.flowpath.model.QualityFilter;

import java.util.function.Consumer;

/**
 * Quality filter panel with sliders for area, eccentricity, solidity, and total intensity.
 * Sits below the TreeView and provides global pre-gating quality control.
 */
public class QualityFilterPane extends TitledPane {

    private QualityFilter filter;
    private boolean suppressEvents = false;

    private final Slider areaMinSlider;
    private final Slider areaMaxSlider;
    private final Slider totalIntensitySlider;
    private final Slider eccentricitySlider;
    private final Slider soliditySlider;

    private final Label areaMinLabel;
    private final Label areaMaxLabel;
    private final Label totalIntLabel;
    private final Label eccentLabel;
    private final Label solidLabel;

    private Consumer<QualityFilter> onFilterChanged;

    public QualityFilterPane(QualityFilter filter) {
        this.filter = filter;
        setText("Quality Filter");
        setCollapsible(true);
        setExpanded(true);

        GridPane grid = new GridPane();
        grid.setHgap(8);
        grid.setVgap(4);
        grid.setPadding(new Insets(6));

        int row = 0;

        // Area min
        grid.add(new Label("Area min:"), 0, row);
        areaMinSlider = createSlider(0, 1000, filter.getMinArea());
        areaMinLabel = new Label(fmt(filter.getMinArea()));
        grid.add(areaMinSlider, 1, row);
        grid.add(areaMinLabel, 2, row);
        row++;

        // Area max
        grid.add(new Label("Area max:"), 0, row);
        areaMaxSlider = createSlider(0, 50000, filter.getMaxArea() == Double.MAX_VALUE ? 50000 : filter.getMaxArea());
        areaMaxLabel = new Label(filter.getMaxArea() == Double.MAX_VALUE ? "off" : fmt(filter.getMaxArea()));
        grid.add(areaMaxSlider, 1, row);
        grid.add(areaMaxLabel, 2, row);
        row++;

        // Total intensity
        grid.add(new Label("Total int:"), 0, row);
        totalIntensitySlider = createSlider(0, 5000, filter.getMinTotalIntensity());
        totalIntLabel = new Label(fmt(filter.getMinTotalIntensity()));
        grid.add(totalIntensitySlider, 1, row);
        grid.add(totalIntLabel, 2, row);
        row++;

        // Eccentricity max
        grid.add(new Label("Eccent max:"), 0, row);
        eccentricitySlider = createSlider(0, 1, filter.getMaxEccentricity());
        eccentLabel = new Label(fmt(filter.getMaxEccentricity()));
        grid.add(eccentricitySlider, 1, row);
        grid.add(eccentLabel, 2, row);
        row++;

        // Solidity min
        grid.add(new Label("Solidity min:"), 0, row);
        soliditySlider = createSlider(0, 1, filter.getMinSolidity());
        solidLabel = new Label(fmt(filter.getMinSolidity()));
        grid.add(soliditySlider, 1, row);
        grid.add(solidLabel, 2, row);
        row++;

        // Reset button
        Button resetBtn = new Button("Reset Filters");
        resetBtn.setMaxWidth(Double.MAX_VALUE);
        resetBtn.setOnAction(e -> resetToDefaults());
        grid.add(resetBtn, 0, row, 3, 1);

        setContent(grid);

        // Wire listeners
        areaMinSlider.valueProperty().addListener((obs, old, val) -> {
            if (suppressEvents) return;
            filter.setMinArea(val.doubleValue());
            areaMinLabel.setText(fmt(val.doubleValue()));
            fireChanged();
        });
        areaMaxSlider.valueProperty().addListener((obs, old, val) -> {
            if (suppressEvents) return;
            double v = val.doubleValue();
            if (v >= areaMaxSlider.getMax() - 1) {
                filter.setMaxArea(Double.MAX_VALUE);
                areaMaxLabel.setText("off");
            } else {
                filter.setMaxArea(v);
                areaMaxLabel.setText(fmt(v));
            }
            fireChanged();
        });
        totalIntensitySlider.valueProperty().addListener((obs, old, val) -> {
            if (suppressEvents) return;
            filter.setMinTotalIntensity(val.doubleValue());
            totalIntLabel.setText(fmt(val.doubleValue()));
            fireChanged();
        });
        eccentricitySlider.valueProperty().addListener((obs, old, val) -> {
            if (suppressEvents) return;
            filter.setMaxEccentricity(val.doubleValue());
            eccentLabel.setText(fmt(val.doubleValue()));
            fireChanged();
        });
        soliditySlider.valueProperty().addListener((obs, old, val) -> {
            if (suppressEvents) return;
            filter.setMinSolidity(val.doubleValue());
            solidLabel.setText(fmt(val.doubleValue()));
            fireChanged();
        });
    }

    /**
     * Replace the filter reference and sync all UI controls to its values.
     */
    public void setFilter(QualityFilter newFilter) {
        suppressEvents = true;
        this.filter = newFilter;
        areaMinSlider.setValue(newFilter.getMinArea());
        areaMinLabel.setText(fmt(newFilter.getMinArea()));
        areaMaxSlider.setValue(newFilter.getMaxArea() == Double.MAX_VALUE ? areaMaxSlider.getMax() : newFilter.getMaxArea());
        areaMaxLabel.setText(newFilter.getMaxArea() == Double.MAX_VALUE ? "off" : fmt(newFilter.getMaxArea()));
        totalIntensitySlider.setValue(newFilter.getMinTotalIntensity());
        totalIntLabel.setText(fmt(newFilter.getMinTotalIntensity()));
        eccentricitySlider.setValue(newFilter.getMaxEccentricity());
        eccentLabel.setText(fmt(newFilter.getMaxEccentricity()));
        soliditySlider.setValue(newFilter.getMinSolidity());
        solidLabel.setText(fmt(newFilter.getMinSolidity()));
        suppressEvents = false;
    }

    /**
     * Update slider ranges based on actual data statistics.
     * Suppresses events to avoid silently corrupting filter values.
     */
    public void updateRanges(double maxArea, double maxTotalIntensity) {
        suppressEvents = true;
        areaMinSlider.setMax(maxArea);
        areaMaxSlider.setMax(maxArea * 1.1);
        totalIntensitySlider.setMax(maxTotalIntensity);
        suppressEvents = false;
    }

    /**
     * Enable/disable QC metric sliders based on data availability.
     */
    public void setAvailableMetrics(boolean hasEccentricity, boolean hasSolidity) {
        eccentricitySlider.setDisable(!hasEccentricity);
        soliditySlider.setDisable(!hasSolidity);
        eccentLabel.setText(hasEccentricity ? fmt(filter.getMaxEccentricity()) : "\u2014");
        solidLabel.setText(hasSolidity ? fmt(filter.getMinSolidity()) : "\u2014");
    }

    public void setOnFilterChanged(Consumer<QualityFilter> callback) {
        this.onFilterChanged = callback;
    }

    public QualityFilter getFilter() {
        return filter;
    }

    /**
     * Reset all filter sliders to their default (no-op) values.
     */
    public void resetToDefaults() {
        suppressEvents = true;
        areaMinSlider.setValue(0);
        areaMinLabel.setText(fmt(0));
        filter.setMinArea(0);

        areaMaxSlider.setValue(areaMaxSlider.getMax());
        areaMaxLabel.setText("off");
        filter.setMaxArea(Double.MAX_VALUE);

        totalIntensitySlider.setValue(0);
        totalIntLabel.setText(fmt(0));
        filter.setMinTotalIntensity(0);

        eccentricitySlider.setValue(1.0);
        eccentLabel.setText(fmt(1.0));
        filter.setMaxEccentricity(1.0);

        soliditySlider.setValue(0.0);
        solidLabel.setText(fmt(0.0));
        filter.setMinSolidity(0.0);

        suppressEvents = false;
        updateActiveIndicator();
        fireChanged();
    }

    /**
     * Update the TitledPane text and color depending on whether any filter is active.
     */
    private void updateActiveIndicator() {
        boolean active = filter.getMinArea() > 0
                || filter.getMaxArea() != Double.MAX_VALUE
                || filter.getMinTotalIntensity() > 0
                || filter.getMaxEccentricity() < 1.0
                || filter.getMinSolidity() > 0.0;
        if (active) {
            setText("Quality Filter (active)");
            setStyle("-fx-text-fill: #ff9900;");
        } else {
            setText("Quality Filter");
            setStyle("");
        }
    }

    private void fireChanged() {
        updateActiveIndicator();
        if (onFilterChanged != null) {
            onFilterChanged.accept(filter);
        }
    }

    private Slider createSlider(double min, double max, double value) {
        Slider slider = new Slider(min, max, value);
        slider.setPrefWidth(150);
        slider.setBlockIncrement((max - min) / 100);
        return slider;
    }

    private String fmt(double v) {
        if (v == (long) v) return String.format("%.0f", v);
        return String.format("%.2f", v);
    }
}
