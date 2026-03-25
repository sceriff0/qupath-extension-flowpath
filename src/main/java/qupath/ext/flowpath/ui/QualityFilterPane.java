package qupath.ext.flowpath.ui;

import javafx.geometry.Insets;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.control.Slider;
import javafx.scene.control.TitledPane;
import javafx.scene.layout.GridPane;
import qupath.ext.flowpath.model.QualityFilter;

import java.util.function.Consumer;

/**
 * Quality filter panel with min+max sliders for area, eccentricity, solidity,
 * total intensity, and perimeter. Sits below the TreeView and provides global
 * pre-gating quality control.
 */
public class QualityFilterPane extends TitledPane {

    private QualityFilter filter;
    private boolean suppressEvents = false;

    private final Slider areaMinSlider, areaMaxSlider;
    private final Slider eccentMinSlider, eccentMaxSlider;
    private final Slider solidMinSlider, solidMaxSlider;
    private final Slider totalIntMinSlider, totalIntMaxSlider;
    private final Slider perimMinSlider, perimMaxSlider;

    private final Label areaMinLabel, areaMaxLabel;
    private final Label eccentMinLabel, eccentMaxLabel;
    private final Label solidMinLabel, solidMaxLabel;
    private final Label totalIntMinLabel, totalIntMaxLabel;
    private final Label perimMinLabel, perimMaxLabel;

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

        // Area
        grid.add(new Label("Area:"), 0, row);
        areaMinSlider = createSlider(0, 1000, filter.getMinArea());
        areaMinLabel = new Label(fmt(filter.getMinArea()));
        areaMaxSlider = createSlider(0, 50000, filter.getMaxArea() == Double.MAX_VALUE ? 50000 : filter.getMaxArea());
        areaMaxLabel = new Label(filter.getMaxArea() == Double.MAX_VALUE ? "off" : fmt(filter.getMaxArea()));
        grid.add(areaMinSlider, 1, row);
        grid.add(areaMinLabel, 2, row);
        grid.add(areaMaxSlider, 3, row);
        grid.add(areaMaxLabel, 4, row);
        row++;

        // Eccentricity
        grid.add(new Label("Eccent:"), 0, row);
        eccentMinSlider = createSlider(0, 1, filter.getMinEccentricity());
        eccentMinLabel = new Label(fmt(filter.getMinEccentricity()));
        eccentMaxSlider = createSlider(0, 1, filter.getMaxEccentricity());
        eccentMaxLabel = new Label(filter.getMaxEccentricity() >= 1.0 ? "off" : fmt(filter.getMaxEccentricity()));
        grid.add(eccentMinSlider, 1, row);
        grid.add(eccentMinLabel, 2, row);
        grid.add(eccentMaxSlider, 3, row);
        grid.add(eccentMaxLabel, 4, row);
        row++;

        // Solidity
        grid.add(new Label("Solidity:"), 0, row);
        solidMinSlider = createSlider(0, 1, filter.getMinSolidity());
        solidMinLabel = new Label(fmt(filter.getMinSolidity()));
        solidMaxSlider = createSlider(0, 1, filter.getMaxSolidity());
        solidMaxLabel = new Label(filter.getMaxSolidity() >= 1.0 ? "off" : fmt(filter.getMaxSolidity()));
        grid.add(solidMinSlider, 1, row);
        grid.add(solidMinLabel, 2, row);
        grid.add(solidMaxSlider, 3, row);
        grid.add(solidMaxLabel, 4, row);
        row++;

        // Total intensity
        grid.add(new Label("Total int:"), 0, row);
        totalIntMinSlider = createSlider(0, 5000, filter.getMinTotalIntensity());
        totalIntMinLabel = new Label(fmt(filter.getMinTotalIntensity()));
        totalIntMaxSlider = createSlider(0, 50000, filter.getMaxTotalIntensity() == Double.MAX_VALUE ? 50000 : filter.getMaxTotalIntensity());
        totalIntMaxLabel = new Label(filter.getMaxTotalIntensity() == Double.MAX_VALUE ? "off" : fmt(filter.getMaxTotalIntensity()));
        grid.add(totalIntMinSlider, 1, row);
        grid.add(totalIntMinLabel, 2, row);
        grid.add(totalIntMaxSlider, 3, row);
        grid.add(totalIntMaxLabel, 4, row);
        row++;

        // Perimeter
        grid.add(new Label("Perim:"), 0, row);
        perimMinSlider = createSlider(0, 1000, filter.getMinPerimeter());
        perimMinLabel = new Label(fmt(filter.getMinPerimeter()));
        perimMaxSlider = createSlider(0, 5000, filter.getMaxPerimeter() == Double.MAX_VALUE ? 5000 : filter.getMaxPerimeter());
        perimMaxLabel = new Label(filter.getMaxPerimeter() == Double.MAX_VALUE ? "off" : fmt(filter.getMaxPerimeter()));
        grid.add(perimMinSlider, 1, row);
        grid.add(perimMinLabel, 2, row);
        grid.add(perimMaxSlider, 3, row);
        grid.add(perimMaxLabel, 4, row);
        row++;

        // Reset button spanning all 5 columns
        Button resetBtn = new Button("Reset Filters");
        resetBtn.setMaxWidth(Double.MAX_VALUE);
        resetBtn.setOnAction(e -> resetToDefaults());
        grid.add(resetBtn, 0, row, 5, 1);

        setContent(grid);

        // Wire all 10 listeners
        wireMinSlider(areaMinSlider, areaMinLabel, filter::setMinArea);
        wireMinSlider(eccentMinSlider, eccentMinLabel, filter::setMinEccentricity);
        wireMinSlider(solidMinSlider, solidMinLabel, filter::setMinSolidity);
        wireMinSlider(totalIntMinSlider, totalIntMinLabel, filter::setMinTotalIntensity);
        wireMinSlider(perimMinSlider, perimMinLabel, filter::setMinPerimeter);
        wireMaxSlider(areaMaxSlider, areaMaxLabel, filter::setMaxArea);
        wireMaxSlider(eccentMaxSlider, eccentMaxLabel, filter::setMaxEccentricity);
        wireMaxSlider(solidMaxSlider, solidMaxLabel, filter::setMaxSolidity);
        wireMaxSlider(totalIntMaxSlider, totalIntMaxLabel, filter::setMaxTotalIntensity);
        wireMaxSlider(perimMaxSlider, perimMaxLabel, filter::setMaxPerimeter);
    }

    private void wireMinSlider(Slider slider, Label label, java.util.function.DoubleConsumer setter) {
        slider.valueProperty().addListener((obs, old, val) -> {
            if (suppressEvents) return;
            setter.accept(val.doubleValue());
            label.setText(fmt(val.doubleValue()));
            fireChanged();
        });
    }

    private void wireMaxSlider(Slider slider, Label label, java.util.function.DoubleConsumer setter) {
        slider.valueProperty().addListener((obs, old, val) -> {
            if (suppressEvents) return;
            double v = val.doubleValue();
            double range = slider.getMax() - slider.getMin();
            // "off" when slider is within 0.1% of its max (works for both [0,1] and [0,50000] ranges)
            if (v >= slider.getMax() - range * 0.001) {
                setter.accept(Double.MAX_VALUE);
                label.setText("off");
            } else {
                setter.accept(v);
                label.setText(fmt(v));
            }
            fireChanged();
        });
    }

    /**
     * Replace the filter reference and sync all UI controls to its values.
     */
    public void setFilter(QualityFilter newFilter) {
        suppressEvents = true;
        this.filter = newFilter;
        syncSlider(areaMinSlider, areaMinLabel, newFilter.getMinArea(), false);
        syncSlider(areaMaxSlider, areaMaxLabel, newFilter.getMaxArea(), true);
        syncSlider(eccentMinSlider, eccentMinLabel, newFilter.getMinEccentricity(), false);
        syncSlider(eccentMaxSlider, eccentMaxLabel, newFilter.getMaxEccentricity(), true);
        syncSlider(solidMinSlider, solidMinLabel, newFilter.getMinSolidity(), false);
        syncSlider(solidMaxSlider, solidMaxLabel, newFilter.getMaxSolidity(), true);
        syncSlider(totalIntMinSlider, totalIntMinLabel, newFilter.getMinTotalIntensity(), false);
        syncSlider(totalIntMaxSlider, totalIntMaxLabel, newFilter.getMaxTotalIntensity(), true);
        syncSlider(perimMinSlider, perimMinLabel, newFilter.getMinPerimeter(), false);
        syncSlider(perimMaxSlider, perimMaxLabel, newFilter.getMaxPerimeter(), true);
        suppressEvents = false;
    }

    private void syncSlider(Slider slider, Label label, double value, boolean isMax) {
        if (isMax && value == Double.MAX_VALUE) {
            slider.setValue(slider.getMax());
            label.setText("off");
        } else if (isMax && value >= 1.0 && slider.getMax() == 1.0) {
            slider.setValue(slider.getMax());
            label.setText("off");
        } else {
            slider.setValue(value);
            label.setText(fmt(value));
        }
    }

    /**
     * Update slider ranges based on actual data statistics.
     * Suppresses events to avoid silently corrupting filter values.
     */
    public void updateRanges(double maxArea, double maxTotalIntensity, double maxPerimeter) {
        suppressEvents = true;
        areaMinSlider.setMax(maxArea);
        areaMaxSlider.setMax(maxArea * 1.1);
        totalIntMinSlider.setMax(maxTotalIntensity);
        totalIntMaxSlider.setMax(maxTotalIntensity * 1.1);
        perimMinSlider.setMax(maxPerimeter);
        perimMaxSlider.setMax(maxPerimeter * 1.1);
        suppressEvents = false;
    }

    /**
     * Enable/disable QC metric sliders based on data availability.
     */
    public void setAvailableMetrics(boolean hasEccentricity, boolean hasSolidity, boolean hasPerimeter) {
        eccentMinSlider.setDisable(!hasEccentricity);
        eccentMaxSlider.setDisable(!hasEccentricity);
        eccentMinLabel.setText(hasEccentricity ? fmt(filter.getMinEccentricity()) : "\u2014");
        eccentMaxLabel.setText(hasEccentricity ? fmt(filter.getMaxEccentricity()) : "\u2014");
        solidMinSlider.setDisable(!hasSolidity);
        solidMaxSlider.setDisable(!hasSolidity);
        solidMinLabel.setText(hasSolidity ? fmt(filter.getMinSolidity()) : "\u2014");
        solidMaxLabel.setText(hasSolidity ? fmt(filter.getMaxSolidity()) : "\u2014");
        perimMinSlider.setDisable(!hasPerimeter);
        perimMaxSlider.setDisable(!hasPerimeter);
        perimMinLabel.setText(hasPerimeter ? fmt(filter.getMinPerimeter()) : "\u2014");
        perimMaxLabel.setText(hasPerimeter ? fmt(filter.getMaxPerimeter()) : "\u2014");
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
        resetSlider(areaMinSlider, areaMinLabel, 0, false);
        filter.setMinArea(0);
        resetSlider(areaMaxSlider, areaMaxLabel, areaMaxSlider.getMax(), true);
        filter.setMaxArea(Double.MAX_VALUE);

        resetSlider(eccentMinSlider, eccentMinLabel, 0, false);
        filter.setMinEccentricity(0);
        resetSlider(eccentMaxSlider, eccentMaxLabel, 1.0, true);
        filter.setMaxEccentricity(1.0);

        resetSlider(solidMinSlider, solidMinLabel, 0, false);
        filter.setMinSolidity(0);
        resetSlider(solidMaxSlider, solidMaxLabel, 1.0, true);
        filter.setMaxSolidity(1.0);

        resetSlider(totalIntMinSlider, totalIntMinLabel, 0, false);
        filter.setMinTotalIntensity(0);
        resetSlider(totalIntMaxSlider, totalIntMaxLabel, totalIntMaxSlider.getMax(), true);
        filter.setMaxTotalIntensity(Double.MAX_VALUE);

        resetSlider(perimMinSlider, perimMinLabel, 0, false);
        filter.setMinPerimeter(0);
        resetSlider(perimMaxSlider, perimMaxLabel, perimMaxSlider.getMax(), true);
        filter.setMaxPerimeter(Double.MAX_VALUE);

        suppressEvents = false;
        updateActiveIndicator();
        fireChanged();
    }

    private void resetSlider(Slider slider, Label label, double value, boolean showOff) {
        slider.setValue(value);
        label.setText(showOff ? "off" : fmt(value));
    }

    /**
     * Update the TitledPane text and color depending on whether any filter is active.
     */
    private void updateActiveIndicator() {
        boolean active = filter.getMinArea() > 0
                || filter.getMaxArea() != Double.MAX_VALUE
                || filter.getMinEccentricity() > 0.0
                || filter.getMaxEccentricity() < 1.0
                || filter.getMinSolidity() > 0.0
                || filter.getMaxSolidity() < 1.0
                || filter.getMinTotalIntensity() > 0
                || filter.getMaxTotalIntensity() != Double.MAX_VALUE
                || filter.getMinPerimeter() > 0
                || filter.getMaxPerimeter() != Double.MAX_VALUE;
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
