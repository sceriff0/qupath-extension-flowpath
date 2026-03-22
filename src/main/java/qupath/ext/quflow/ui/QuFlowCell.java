package qupath.ext.quflow.ui;

import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.Label;
import javafx.scene.control.TreeCell;
import javafx.scene.layout.Background;
import javafx.scene.layout.BackgroundFill;
import javafx.scene.layout.CornerRadii;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import qupath.ext.quflow.model.ColorUtils;
import qupath.ext.quflow.model.GateNode;

/**
 * Custom TreeCell for rendering gate tree items with a polished dark-theme design (v0.2.0).
 * <p>
 * Gate nodes render as full-width colored bars (muted blue-gray) with bold channel name
 * and lighter threshold text. Branch items render as colored pill badges whose background
 * matches the user-set positive/negative color, with right-aligned counts and a star
 * marker for leaf phenotypes.
 */
public class QuFlowCell extends TreeCell<Object> {

    private static final Color GATE_BAR_COLOR = Color.web("#3a4a5a");
    private static final CornerRadii BAR_RADII = new CornerRadii(6);
    private static final CornerRadii PILL_RADII = new CornerRadii(10);
    private static final Insets BAR_PADDING = new Insets(5, 10, 5, 10);
    private static final Insets PILL_PADDING = new Insets(2, 10, 2, 10);
    private static final Insets CELL_PADDING = new Insets(2, 0, 2, 0);
    private static final String STAR = "\u2605";

    /**
     * Wrapper for a branch (positive or negative side of a gate).
     */
    public static class BranchItem {
        public final GateNode parentGate;
        public final boolean isPositive;

        public BranchItem(GateNode parentGate, boolean isPositive) {
            this.parentGate = parentGate;
            this.isPositive = isPositive;
        }
    }

    @Override
    protected void updateItem(Object item, boolean empty) {
        super.updateItem(item, empty);

        if (empty || item == null) {
            setText(null);
            setGraphic(null);
            setBackground(Background.EMPTY);
            setPadding(Insets.EMPTY);
            return;
        }

        setPadding(CELL_PADDING);

        if (item instanceof GateNode node) {
            setGraphic(buildGateNodeGraphic(node));
            setText(null);
        } else if (item instanceof BranchItem branch) {
            setGraphic(buildBranchGraphic(branch));
            setText(null);
        } else {
            setText(item.toString());
            setGraphic(null);
        }
    }

    // ---- Gate node: full-width colored bar ------------------------------------------------

    private HBox buildGateNodeGraphic(GateNode node) {
        HBox bar = new HBox(8);
        bar.setAlignment(Pos.CENTER_LEFT);
        bar.setPadding(BAR_PADDING);
        bar.setBackground(new Background(new BackgroundFill(GATE_BAR_COLOR, BAR_RADII, Insets.EMPTY)));
        // Stretch to fill available width
        bar.setMaxWidth(Double.MAX_VALUE);
        HBox.setHgrow(bar, Priority.ALWAYS);

        // Channel name - bold white
        Label channelLabel = new Label(node.getChannel());
        channelLabel.setFont(Font.font(null, FontWeight.BOLD, 13));
        channelLabel.setTextFill(Color.WHITE);

        // Threshold value - lighter text
        String thresholdText = node.isThresholdIsZScore()
                ? String.format("z = %.3f", node.getThreshold())
                : String.format("t = %.3f", node.getThreshold());
        Label thresholdLabel = new Label(thresholdText);
        thresholdLabel.setFont(Font.font(null, FontWeight.NORMAL, 11));
        thresholdLabel.setTextFill(Color.web("#a0b0c0"));

        bar.getChildren().addAll(channelLabel, thresholdLabel);
        return bar;
    }

    // ---- Branch item: colored pill/badge --------------------------------------------------

    private HBox buildBranchGraphic(BranchItem branch) {
        GateNode gate = branch.parentGate;
        boolean isPos = branch.isPositive;

        int colorInt = isPos ? gate.getPositiveColor() : gate.getNegativeColor();
        Color pillColor = ColorUtils.intToColor(colorInt);
        String name = isPos ? gate.getPositiveName() : gate.getNegativeName();
        int count = isPos ? gate.getPosCount() : gate.getNegCount();
        boolean isLeaf = isPos ? gate.getPositiveChildren().isEmpty() : gate.getNegativeChildren().isEmpty();

        // Outer row spanning full width
        HBox row = new HBox(8);
        row.setAlignment(Pos.CENTER_LEFT);
        row.setMaxWidth(Double.MAX_VALUE);
        HBox.setHgrow(row, Priority.ALWAYS);

        // Pill badge
        HBox pill = new HBox(4);
        pill.setAlignment(Pos.CENTER_LEFT);
        pill.setPadding(PILL_PADDING);
        pill.setBackground(new Background(new BackgroundFill(pillColor, PILL_RADII, Insets.EMPTY)));

        String displayName = isLeaf ? (STAR + " " + name) : name;
        Label nameLabel = new Label(displayName);
        nameLabel.setTextFill(Color.WHITE);
        nameLabel.setFont(Font.font(null, isLeaf ? FontWeight.BOLD : FontWeight.NORMAL, 12));

        pill.getChildren().add(nameLabel);

        // Spacer pushes count to the right
        Region spacer = new Region();
        HBox.setHgrow(spacer, Priority.ALWAYS);

        // Count label - right-aligned gray text
        Label countLabel = new Label(String.format("%,d", count));
        countLabel.setFont(Font.font(null, FontWeight.NORMAL, 11));
        countLabel.setTextFill(Color.web("#888888"));

        row.getChildren().addAll(pill, spacer, countLabel);
        return row;
    }
}
