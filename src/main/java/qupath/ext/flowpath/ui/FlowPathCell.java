package qupath.ext.flowpath.ui;

import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.Label;
import javafx.scene.control.TreeCell;
import javafx.scene.control.TreeItem;
import javafx.scene.layout.Background;
import javafx.scene.layout.BackgroundFill;
import javafx.scene.layout.CornerRadii;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Priority;
import javafx.scene.layout.Region;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.scene.text.FontWeight;
import qupath.ext.flowpath.model.Branch;
import qupath.ext.flowpath.model.ColorUtils;
import qupath.ext.flowpath.model.GateNode;
import qupath.ext.flowpath.model.QuadrantGate;

/**
 * Custom TreeCell for rendering gate tree items with a polished dark-theme design.
 * <p>
 * Gate nodes render as full-width colored bars with bold channel name and threshold text.
 * Branch items render as colored pill badges with right-aligned counts and a star marker
 * for leaf phenotypes. Supports both threshold (2 branches) and quadrant (4 branches) gates.
 */
public class FlowPathCell extends TreeCell<Object> {

    private static final Color GATE_BAR_COLOR = Color.web("#3a4a5a");
    private static final Color GATE_BAR_DISABLED_COLOR = Color.web("#2a2a2a");
    private static final CornerRadii BAR_RADII = new CornerRadii(6);
    private static final CornerRadii PILL_RADII = new CornerRadii(10);
    private static final Insets BAR_PADDING = new Insets(5, 10, 5, 10);
    private static final Insets PILL_PADDING = new Insets(2, 10, 2, 10);
    private static final Insets CELL_PADDING = new Insets(2, 0, 2, 0);
    private static final String STAR = "\u2605";

    /**
     * Wrapper for a branch of a gate (generic — works for any gate type).
     */
    public static class BranchItem {
        public final GateNode parentGate;
        public final Branch branch;
        public final int branchIndex;
        /** @deprecated Use the Branch-based constructor instead. */
        public final boolean isPositive;

        /**
         * New generic constructor using Branch reference.
         */
        public BranchItem(GateNode parentGate, Branch branch, int branchIndex) {
            this.parentGate = parentGate;
            this.branch = branch;
            this.branchIndex = branchIndex;
            this.isPositive = (branchIndex == 0);
        }

        /**
         * Backward-compatible constructor for threshold gates.
         */
        public BranchItem(GateNode parentGate, boolean isPositive) {
            this.parentGate = parentGate;
            this.branchIndex = isPositive ? 0 : 1;
            this.branch = parentGate.getBranches().get(this.branchIndex);
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
        } else if (item instanceof BranchItem bi) {
            setGraphic(buildBranchGraphic(bi));
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
        Color barColor = node.isEnabled() ? GATE_BAR_COLOR : GATE_BAR_DISABLED_COLOR;
        bar.setBackground(new Background(new BackgroundFill(barColor, BAR_RADII, Insets.EMPTY)));
        bar.setMaxWidth(Double.MAX_VALUE);
        bar.setOpacity(node.isEnabled() ? 1.0 : 0.5);
        HBox.setHgrow(bar, Priority.ALWAYS);

        if (node instanceof QuadrantGate qg) {
            // Quadrant gate: show both channels
            Label channelLabel = new Label(qg.getChannelX() + " / " + qg.getChannelY());
            channelLabel.setFont(Font.font(null, FontWeight.BOLD, 13));
            channelLabel.setTextFill(Color.WHITE);

            Label typeLabel = new Label("[Quadrant]");
            typeLabel.setFont(Font.font(null, FontWeight.NORMAL, 10));
            typeLabel.setTextFill(Color.web("#80b0d0"));

            bar.getChildren().addAll(channelLabel, typeLabel);
        } else {
            // Threshold gate
            Label channelLabel = new Label(node.getChannel());
            channelLabel.setFont(Font.font(null, FontWeight.BOLD, 13));
            channelLabel.setTextFill(Color.WHITE);

            String thresholdText = node.isThresholdIsZScore()
                    ? String.format("z = %.3f", node.getThreshold())
                    : String.format("t = %.3f", node.getThreshold());
            Label thresholdLabel = new Label(thresholdText);
            thresholdLabel.setFont(Font.font(null, FontWeight.NORMAL, 11));
            thresholdLabel.setTextFill(Color.web("#a0b0c0"));

            bar.getChildren().addAll(channelLabel, thresholdLabel);
        }

        return bar;
    }

    // ---- Branch item: colored pill/badge --------------------------------------------------

    private HBox buildBranchGraphic(BranchItem bi) {
        Branch branch = bi.branch;
        Color pillColor = ColorUtils.intToColor(branch.getColor());
        String name = branch.getName();
        int count = branch.getCount();
        boolean isLeaf = branch.isLeaf();

        HBox row = new HBox(8);
        row.setAlignment(Pos.CENTER_LEFT);
        row.setMaxWidth(Double.MAX_VALUE);
        HBox.setHgrow(row, Priority.ALWAYS);

        HBox pill = new HBox(4);
        pill.setAlignment(Pos.CENTER_LEFT);
        pill.setPadding(PILL_PADDING);
        pill.setBackground(new Background(new BackgroundFill(pillColor, PILL_RADII, Insets.EMPTY)));

        String displayName = isLeaf ? (STAR + " " + name) : name;
        Label nameLabel = new Label(displayName);
        nameLabel.setTextFill(Color.WHITE);
        nameLabel.setFont(Font.font(null, isLeaf ? FontWeight.BOLD : FontWeight.NORMAL, 12));

        pill.getChildren().add(nameLabel);

        Region spacer = new Region();
        HBox.setHgrow(spacer, Priority.ALWAYS);

        // Show count with percentage relative to parent gate total
        int totalParent = 0;
        for (Branch b : bi.parentGate.getBranches()) {
            totalParent += b.getCount();
        }
        String countText = totalParent > 0
            ? String.format("%,d (%.1f%%)", count, 100.0 * count / totalParent)
            : String.format("%,d", count);
        Label countLabel = new Label(countText);
        countLabel.setFont(Font.font(null, FontWeight.NORMAL, 11));
        countLabel.setTextFill(Color.web("#888888"));

        row.getChildren().addAll(pill, spacer, countLabel);
        return row;
    }
}
