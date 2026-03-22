package qupath.ext.quflow;

import javafx.scene.Scene;
import javafx.scene.control.MenuItem;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyCodeCombination;
import javafx.scene.input.KeyCombination;
import javafx.stage.Stage;
import qupath.ext.quflow.ui.QuFlowPane;
import qupath.lib.gui.QuPathGUI;
import qupath.lib.gui.extensions.QuPathExtension;

/**
 * QuPath extension entry point for the interactive tree-based gating tool.
 * Registers a menu item under Extensions and opens a floating Stage.
 */
public class QuFlowExtension implements QuPathExtension {

    private static final String NAME = "QuFlow";
    private static final String DESCRIPTION = "Interactive tree-based cell phenotyping gating";

    private Stage stage;
    private QuFlowPane gateTreePane;

    @Override
    public void installExtension(QuPathGUI qupath) {
        var menuItem = new MenuItem(NAME);
        menuItem.setOnAction(e -> showGateTreeWindow(qupath));
        menuItem.setAccelerator(new KeyCodeCombination(KeyCode.G, KeyCombination.CONTROL_DOWN));
        qupath.getMenu("Extensions", true).getItems().add(menuItem);
    }

    private void showGateTreeWindow(QuPathGUI qupath) {
        if (stage != null && stage.isShowing()) {
            stage.toFront();
            stage.requestFocus();
            return;
        }

        gateTreePane = new QuFlowPane(qupath);

        stage = new Stage();
        stage.setTitle("QuFlow \u2014 Cell Phenotyping");
        stage.initOwner(qupath.getStage());
        stage.setScene(new Scene(gateTreePane, 900, 700));
        stage.setMinWidth(700);
        stage.setMinHeight(500);

        stage.setOnCloseRequest(e -> {
            if (gateTreePane != null) {
                gateTreePane.shutdown();
            }
        });

        stage.show();
    }

    @Override
    public String getName() {
        return NAME;
    }

    @Override
    public String getDescription() {
        return DESCRIPTION;
    }
}
