package qupath.ext.flowpath.ui;

import javafx.scene.canvas.Canvas;
import javafx.scene.canvas.GraphicsContext;
import javafx.scene.input.MouseEvent;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

/**
 * Canvas-based 2D scatter plot for visualizing cells on two marker axes.
 * Supports overlay of gate boundaries (polygon, rectangle, ellipse).
 * Subsamples for display if cell count exceeds MAX_DISPLAY_POINTS.
 */
public class ScatterPlotCanvas extends Canvas {

    private static final int MAX_DISPLAY_POINTS = 20000;
    private static final double PADDING_LEFT = 45;
    private static final double PADDING_RIGHT = 10;
    private static final double PADDING_TOP = 10;
    private static final double PADDING_BOTTOM = 30;
    private static final double DOT_SIZE = 2.0;

    private double[] xValues;
    private double[] yValues;
    private double minX, maxX, minY, maxY;
    private String labelX = "X";
    private String labelY = "Y";

    // Gate overlays
    private List<double[]> polygonVertices;
    private double[] rectBounds;  // [minX, maxX, minY, maxY]
    private double[] ellipseParams; // [centerX, centerY, radiusX, radiusY]

    private Color insideColor = Color.rgb(0, 200, 0, 0.6);
    private Color outsideColor = Color.rgb(128, 128, 128, 0.3);

    // Drawing interaction
    public enum DrawingMode { NONE, POLYGON, RECTANGLE, ELLIPSE }

    private DrawingMode drawingMode = DrawingMode.NONE;
    private final List<double[]> drawingVertices = new ArrayList<>();
    private double[] dragStart;
    private double[] dragCurrent;
    private Consumer<List<double[]>> onPolygonDrawn;
    private Consumer<double[]> onRectangleDrawn;
    private Consumer<double[]> onEllipseDrawn;

    // Handle editing existing overlays
    private int dragHandleIndex = -1;
    private static final double HANDLE_RADIUS = 6.0;
    private static final double HANDLE_HIT_RADIUS = 8.0;

    public ScatterPlotCanvas() {
        super(380, 300);
        widthProperty().addListener((obs, o, n) -> repaint());
        heightProperty().addListener((obs, o, n) -> repaint());

        setOnMousePressed(this::handleMousePressed);
        setOnMouseDragged(this::handleMouseDragged);
        setOnMouseReleased(this::handleMouseReleased);
    }

    @Override public boolean isResizable() { return true; }
    @Override public double prefWidth(double h) { return 380; }
    @Override public double prefHeight(double w) { return 300; }
    @Override public double minWidth(double h) { return 200; }
    @Override public double minHeight(double w) { return 150; }
    @Override public double maxWidth(double h) { return Double.MAX_VALUE; }
    @Override public double maxHeight(double w) { return Double.MAX_VALUE; }

    @Override
    public void resize(double width, double height) {
        setWidth(width);
        setHeight(height);
        repaint();
    }

    public void setData(double[] xValues, double[] yValues, String labelX, String labelY) {
        this.xValues = xValues;
        this.yValues = yValues;
        this.labelX = labelX;
        this.labelY = labelY;

        if (xValues == null || yValues == null || xValues.length == 0) {
            repaint();
            return;
        }

        // Compute ranges
        minX = Double.MAX_VALUE; maxX = -Double.MAX_VALUE;
        minY = Double.MAX_VALUE; maxY = -Double.MAX_VALUE;
        for (int i = 0; i < xValues.length; i++) {
            if (!Double.isNaN(xValues[i]) && !Double.isNaN(yValues[i])) {
                minX = Math.min(minX, xValues[i]);
                maxX = Math.max(maxX, xValues[i]);
                minY = Math.min(minY, yValues[i]);
                maxY = Math.max(maxY, yValues[i]);
            }
        }
        if (maxX <= minX) maxX = minX + 1;
        if (maxY <= minY) maxY = minY + 1;

        // Add 5% padding
        double padX = (maxX - minX) * 0.05;
        double padY = (maxY - minY) * 0.05;
        minX -= padX; maxX += padX;
        minY -= padY; maxY += padY;

        repaint();
    }

    public void setPolygonOverlay(List<double[]> vertices) {
        this.polygonVertices = vertices;
        this.rectBounds = null;
        this.ellipseParams = null;
        repaint();
    }

    public void setRectangleOverlay(double minX, double maxX, double minY, double maxY) {
        this.rectBounds = new double[]{minX, maxX, minY, maxY};
        this.polygonVertices = null;
        this.ellipseParams = null;
        repaint();
    }

    public void setEllipseOverlay(double cx, double cy, double rx, double ry) {
        this.ellipseParams = new double[]{cx, cy, rx, ry};
        this.polygonVertices = null;
        this.rectBounds = null;
        repaint();
    }

    public void clearOverlay() {
        this.polygonVertices = null;
        this.rectBounds = null;
        this.ellipseParams = null;
        repaint();
    }

    public void setInsideColor(Color c) { this.insideColor = c; repaint(); }
    public void setOutsideColor(Color c) { this.outsideColor = c; repaint(); }

    // ---- Coordinate conversion helpers ----

    public static double valueToPixel(double value, double min, double max, double plotSize) {
        if (max <= min) return 0;
        return ((value - min) / (max - min)) * plotSize;
    }

    public static double pixelToValue(double pixel, double min, double max, double plotSize) {
        if (plotSize <= 0) return min;
        return min + (pixel / plotSize) * (max - min);
    }

    public double dataXToScreenX(double dataX) {
        double plotW = getWidth() - PADDING_LEFT - PADDING_RIGHT;
        return PADDING_LEFT + valueToPixel(dataX, minX, maxX, plotW);
    }

    public double dataYToScreenY(double dataY) {
        double plotH = getHeight() - PADDING_TOP - PADDING_BOTTOM;
        return PADDING_TOP + plotH - valueToPixel(dataY, minY, maxY, plotH);
    }

    public double screenXToDataX(double screenX) {
        double plotW = getWidth() - PADDING_LEFT - PADDING_RIGHT;
        return pixelToValue(screenX - PADDING_LEFT, minX, maxX, plotW);
    }

    public double screenYToDataY(double screenY) {
        double plotH = getHeight() - PADDING_TOP - PADDING_BOTTOM;
        return pixelToValue(PADDING_TOP + plotH - screenY, minY, maxY, plotH);
    }

    // ---- Drawing mode API ----

    public void setDrawingMode(DrawingMode mode) {
        this.drawingMode = mode;
        drawingVertices.clear();
        dragStart = null;
        dragCurrent = null;
        dragHandleIndex = -1;
        repaint();
    }

    public DrawingMode getDrawingMode() { return drawingMode; }

    public void setOnPolygonDrawn(Consumer<List<double[]>> cb) { this.onPolygonDrawn = cb; }
    public void setOnRectangleDrawn(Consumer<double[]> cb) { this.onRectangleDrawn = cb; }
    public void setOnEllipseDrawn(Consumer<double[]> cb) { this.onEllipseDrawn = cb; }

    // ---- Mouse interaction ----

    private void handleMousePressed(MouseEvent e) {
        if (drawingMode == DrawingMode.NONE) return;

        double sx = e.getX();
        double sy = e.getY();

        // Check for handle hit first (editing existing overlay), but not while drawing a new polygon
        int handle = (drawingMode == DrawingMode.POLYGON && !drawingVertices.isEmpty()) ? -1 : findHandle(sx, sy);
        if (handle >= 0) {
            dragHandleIndex = handle;
            dragStart = new double[]{sx, sy};
            dragCurrent = new double[]{sx, sy};
            e.consume();
            return;
        }

        if (drawingMode == DrawingMode.POLYGON) {
            if (e.getClickCount() == 2) {
                // Double-click: close polygon if enough vertices
                if (drawingVertices.size() >= 3 && onPolygonDrawn != null) {
                    onPolygonDrawn.accept(new ArrayList<>(drawingVertices));
                }
                drawingVertices.clear();
                repaint();
            } else if (e.getClickCount() == 1) {
                // Single-click: add vertex
                double dx = screenXToDataX(sx);
                double dy = screenYToDataY(sy);
                drawingVertices.add(new double[]{dx, dy});
                repaint();
            }
        } else {
            // RECTANGLE or ELLIPSE: start drag
            dragStart = new double[]{sx, sy};
            dragCurrent = new double[]{sx, sy};
        }
        e.consume();
    }

    private void handleMouseDragged(MouseEvent e) {
        if (drawingMode == DrawingMode.NONE) return;

        double sx = e.getX();
        double sy = e.getY();

        if (dragHandleIndex >= 0) {
            // Dragging an existing handle — update overlay
            updateHandleDrag(sx, sy);
            repaint();
            e.consume();
            return;
        }

        if (drawingMode == DrawingMode.RECTANGLE || drawingMode == DrawingMode.ELLIPSE) {
            dragCurrent = new double[]{sx, sy};
            repaint();
        }
        e.consume();
    }

    private void handleMouseReleased(MouseEvent e) {
        if (drawingMode == DrawingMode.NONE) return;

        double sx = e.getX();
        double sy = e.getY();

        if (dragHandleIndex >= 0) {
            updateHandleDrag(sx, sy);
            fireHandleDragComplete();
            dragHandleIndex = -1;
            dragStart = null;
            dragCurrent = null;
            repaint();
            e.consume();
            return;
        }

        if (drawingMode == DrawingMode.RECTANGLE && dragStart != null) {
            double x1 = screenXToDataX(dragStart[0]);
            double y1 = screenYToDataY(dragStart[1]);
            double x2 = screenXToDataX(sx);
            double y2 = screenYToDataY(sy);
            double rMinX = Math.min(x1, x2), rMaxX = Math.max(x1, x2);
            double rMinY = Math.min(y1, y2), rMaxY = Math.max(y1, y2);
            if (rMaxX - rMinX > 1e-10 && rMaxY - rMinY > 1e-10 && onRectangleDrawn != null) {
                onRectangleDrawn.accept(new double[]{rMinX, rMaxX, rMinY, rMaxY});
            }
            dragStart = null;
            dragCurrent = null;
            repaint();
        } else if (drawingMode == DrawingMode.ELLIPSE && dragStart != null) {
            double x1 = screenXToDataX(dragStart[0]);
            double y1 = screenYToDataY(dragStart[1]);
            double x2 = screenXToDataX(sx);
            double y2 = screenYToDataY(sy);
            double cx = (x1 + x2) / 2.0;
            double cy = (y1 + y2) / 2.0;
            double rx = Math.abs(x2 - x1) / 2.0;
            double ry = Math.abs(y2 - y1) / 2.0;
            if (rx > 1e-10 && ry > 1e-10 && onEllipseDrawn != null) {
                onEllipseDrawn.accept(new double[]{cx, cy, rx, ry});
            }
            dragStart = null;
            dragCurrent = null;
            repaint();
        }
        e.consume();
    }

    // ---- Handle hit-testing and dragging ----

    private int findHandle(double sx, double sy) {
        if (polygonVertices != null && polygonVertices.size() >= 3) {
            for (int i = 0; i < polygonVertices.size(); i++) {
                double hx = dataXToScreenX(polygonVertices.get(i)[0]);
                double hy = dataYToScreenY(polygonVertices.get(i)[1]);
                if (Math.hypot(sx - hx, sy - hy) <= HANDLE_HIT_RADIUS) return i;
            }
        }
        if (rectBounds != null) {
            // Handles: 0=topLeft, 1=topRight, 2=bottomRight, 3=bottomLeft
            double[][] corners = getRectHandleScreenPositions();
            for (int i = 0; i < corners.length; i++) {
                if (Math.hypot(sx - corners[i][0], sy - corners[i][1]) <= HANDLE_HIT_RADIUS) return i;
            }
        }
        if (ellipseParams != null) {
            // Handles: 0=top, 1=right, 2=bottom, 3=left (cardinal points)
            double[][] cardinals = getEllipseHandleScreenPositions();
            for (int i = 0; i < cardinals.length; i++) {
                if (Math.hypot(sx - cardinals[i][0], sy - cardinals[i][1]) <= HANDLE_HIT_RADIUS) return i;
            }
        }
        return -1;
    }

    private double[][] getRectHandleScreenPositions() {
        if (rectBounds == null) return new double[0][];
        double x0 = dataXToScreenX(rectBounds[0]);
        double x1 = dataXToScreenX(rectBounds[1]);
        double y0 = dataYToScreenY(rectBounds[3]); // maxY -> top of screen
        double y1 = dataYToScreenY(rectBounds[2]); // minY -> bottom of screen
        return new double[][]{{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
    }

    private double[][] getEllipseHandleScreenPositions() {
        if (ellipseParams == null) return new double[0][];
        double cx = dataXToScreenX(ellipseParams[0]);
        double cy = dataYToScreenY(ellipseParams[1]);
        double plotW = getWidth() - PADDING_LEFT - PADDING_RIGHT;
        double plotH = getHeight() - PADDING_TOP - PADDING_BOTTOM;
        double erx = valueToPixel(ellipseParams[2], 0, maxX - minX, plotW);
        double ery = valueToPixel(ellipseParams[3], 0, maxY - minY, plotH);
        // top, right, bottom, left
        return new double[][]{{cx, cy - ery}, {cx + erx, cy}, {cx, cy + ery}, {cx - erx, cy}};
    }

    private void updateHandleDrag(double sx, double sy) {
        double dx = screenXToDataX(sx);
        double dy = screenYToDataY(sy);

        if (polygonVertices != null && dragHandleIndex >= 0 && dragHandleIndex < polygonVertices.size()) {
            polygonVertices.get(dragHandleIndex)[0] = dx;
            polygonVertices.get(dragHandleIndex)[1] = dy;
        } else if (rectBounds != null && dragHandleIndex >= 0 && dragHandleIndex < 4) {
            // Move corner: 0=topLeft, 1=topRight, 2=bottomRight, 3=bottomLeft
            switch (dragHandleIndex) {
                case 0 -> { rectBounds[0] = dx; rectBounds[3] = dy; }
                case 1 -> { rectBounds[1] = dx; rectBounds[3] = dy; }
                case 2 -> { rectBounds[1] = dx; rectBounds[2] = dy; }
                case 3 -> { rectBounds[0] = dx; rectBounds[2] = dy; }
            }
            // Ensure min < max
            if (rectBounds[0] > rectBounds[1]) { double t = rectBounds[0]; rectBounds[0] = rectBounds[1]; rectBounds[1] = t; }
            if (rectBounds[2] > rectBounds[3]) { double t = rectBounds[2]; rectBounds[2] = rectBounds[3]; rectBounds[3] = t; }
        } else if (ellipseParams != null && dragHandleIndex >= 0 && dragHandleIndex < 4) {
            // Cardinal points: 0=top, 1=right, 2=bottom, 3=left
            switch (dragHandleIndex) {
                case 0 -> ellipseParams[3] = Math.abs(ellipseParams[1] - dy); // top: adjust ry
                case 1 -> ellipseParams[2] = Math.abs(dx - ellipseParams[0]); // right: adjust rx
                case 2 -> ellipseParams[3] = Math.abs(dy - ellipseParams[1]); // bottom: adjust ry
                case 3 -> ellipseParams[2] = Math.abs(ellipseParams[0] - dx); // left: adjust rx
            }
        }
    }

    private void fireHandleDragComplete() {
        if (polygonVertices != null && onPolygonDrawn != null) {
            onPolygonDrawn.accept(new ArrayList<>(polygonVertices));
        } else if (rectBounds != null && onRectangleDrawn != null) {
            onRectangleDrawn.accept(new double[]{rectBounds[0], rectBounds[1], rectBounds[2], rectBounds[3]});
        } else if (ellipseParams != null && onEllipseDrawn != null) {
            onEllipseDrawn.accept(new double[]{ellipseParams[0], ellipseParams[1], ellipseParams[2], ellipseParams[3]});
        }
    }

    private void repaint() {
        GraphicsContext gc = getGraphicsContext2D();
        double w = getWidth();
        double h = getHeight();

        gc.setFill(Color.rgb(30, 30, 30));
        gc.fillRect(0, 0, w, h);

        if (xValues == null || yValues == null || xValues.length == 0) {
            gc.setFill(Color.gray(0.5));
            gc.setFont(Font.font(12));
            gc.fillText("No data", w / 2 - 20, h / 2);
            return;
        }

        double plotW = w - PADDING_LEFT - PADDING_RIGHT;
        double plotH = h - PADDING_TOP - PADDING_BOTTOM;

        // Draw dots
        int step = Math.max(1, xValues.length / MAX_DISPLAY_POINTS);
        for (int i = 0; i < xValues.length; i += step) {
            if (Double.isNaN(xValues[i]) || Double.isNaN(yValues[i])) continue;
            double px = PADDING_LEFT + valueToPixel(xValues[i], minX, maxX, plotW);
            double py = PADDING_TOP + plotH - valueToPixel(yValues[i], minY, maxY, plotH);

            boolean inside = isInsideOverlay(xValues[i], yValues[i]);
            gc.setFill(inside ? insideColor : outsideColor);
            gc.fillOval(px - DOT_SIZE / 2, py - DOT_SIZE / 2, DOT_SIZE, DOT_SIZE);
        }

        // Draw gate overlay
        drawOverlay(gc, plotW, plotH);

        // Draw in-progress drawing preview
        drawDrawingPreview(gc, plotW, plotH);

        // Border
        gc.setStroke(Color.gray(0.3));
        gc.setLineWidth(1);
        gc.strokeRect(PADDING_LEFT, PADDING_TOP, plotW, plotH);

        // Axis labels
        gc.setFill(Color.gray(0.7));
        gc.setFont(Font.font(9));
        gc.fillText(String.format("%.2f", minX), PADDING_LEFT, h - 3);
        gc.fillText(String.format("%.2f", maxX), PADDING_LEFT + plotW - 30, h - 3);
        gc.fillText(labelX, PADDING_LEFT + plotW / 2 - 10, h - 3);

        gc.save();
        gc.translate(10, PADDING_TOP + plotH / 2);
        gc.rotate(-90);
        gc.fillText(labelY, 0, 0);
        gc.restore();
    }

    private boolean isInsideOverlay(double x, double y) {
        if (polygonVertices != null && polygonVertices.size() >= 3) {
            return pointInPolygon(x, y, polygonVertices);
        }
        if (rectBounds != null) {
            return x >= rectBounds[0] && x <= rectBounds[1] && y >= rectBounds[2] && y <= rectBounds[3];
        }
        if (ellipseParams != null) {
            double dx = (x - ellipseParams[0]) / ellipseParams[2];
            double dy = (y - ellipseParams[1]) / ellipseParams[3];
            return dx * dx + dy * dy <= 1.0;
        }
        return true; // No overlay = all inside
    }

    private void drawOverlay(GraphicsContext gc, double plotW, double plotH) {
        gc.setStroke(Color.YELLOW);
        gc.setLineWidth(1.5);

        if (polygonVertices != null && polygonVertices.size() >= 3) {
            double[] xp = new double[polygonVertices.size()];
            double[] yp = new double[polygonVertices.size()];
            for (int i = 0; i < polygonVertices.size(); i++) {
                xp[i] = dataXToScreenX(polygonVertices.get(i)[0]);
                yp[i] = dataYToScreenY(polygonVertices.get(i)[1]);
            }
            gc.strokePolygon(xp, yp, xp.length);

            // Draw handles at vertices
            drawHandles(gc, xp, yp);
        }

        if (rectBounds != null) {
            double rx = dataXToScreenX(rectBounds[0]);
            double ry = dataYToScreenY(rectBounds[3]);
            double rw = valueToPixel(rectBounds[1], minX, maxX, plotW) - valueToPixel(rectBounds[0], minX, maxX, plotW);
            double rh = valueToPixel(rectBounds[3], minY, maxY, plotH) - valueToPixel(rectBounds[2], minY, maxY, plotH);
            gc.strokeRect(rx, ry, rw, rh);

            // Draw handles at corners
            double[][] corners = getRectHandleScreenPositions();
            drawHandles(gc, arrayCol(corners, 0), arrayCol(corners, 1));
        }

        if (ellipseParams != null) {
            double cx = dataXToScreenX(ellipseParams[0]);
            double cy = dataYToScreenY(ellipseParams[1]);
            double erx = valueToPixel(ellipseParams[2], 0, maxX - minX, plotW);
            double ery = valueToPixel(ellipseParams[3], 0, maxY - minY, plotH);
            gc.strokeOval(cx - erx, cy - ery, erx * 2, ery * 2);

            // Draw handles at cardinal points
            double[][] cardinals = getEllipseHandleScreenPositions();
            drawHandles(gc, arrayCol(cardinals, 0), arrayCol(cardinals, 1));
        }
    }

    private void drawHandles(GraphicsContext gc, double[] hx, double[] hy) {
        gc.setFill(Color.WHITE);
        gc.setStroke(Color.YELLOW);
        gc.setLineWidth(1.0);
        for (int i = 0; i < hx.length; i++) {
            gc.fillOval(hx[i] - HANDLE_RADIUS, hy[i] - HANDLE_RADIUS, HANDLE_RADIUS * 2, HANDLE_RADIUS * 2);
            gc.strokeOval(hx[i] - HANDLE_RADIUS, hy[i] - HANDLE_RADIUS, HANDLE_RADIUS * 2, HANDLE_RADIUS * 2);
        }
        gc.setLineWidth(1.5);
    }

    private static double[] arrayCol(double[][] arr, int col) {
        double[] result = new double[arr.length];
        for (int i = 0; i < arr.length; i++) result[i] = arr[i][col];
        return result;
    }

    private void drawDrawingPreview(GraphicsContext gc, double plotW, double plotH) {
        gc.setStroke(Color.CYAN);
        gc.setLineWidth(1.5);
        gc.setLineDashes(6, 4);

        if (drawingMode == DrawingMode.POLYGON && !drawingVertices.isEmpty()) {
            // Draw lines connecting vertices
            for (int i = 0; i < drawingVertices.size() - 1; i++) {
                double x1 = dataXToScreenX(drawingVertices.get(i)[0]);
                double y1 = dataYToScreenY(drawingVertices.get(i)[1]);
                double x2 = dataXToScreenX(drawingVertices.get(i + 1)[0]);
                double y2 = dataYToScreenY(drawingVertices.get(i + 1)[1]);
                gc.strokeLine(x1, y1, x2, y2);
            }
            // Draw small circles at each vertex
            gc.setFill(Color.CYAN);
            for (double[] v : drawingVertices) {
                double sx = dataXToScreenX(v[0]);
                double sy = dataYToScreenY(v[1]);
                gc.fillOval(sx - 4, sy - 4, 8, 8);
            }
        }

        if ((drawingMode == DrawingMode.RECTANGLE) && dragStart != null && dragCurrent != null) {
            double x = Math.min(dragStart[0], dragCurrent[0]);
            double y = Math.min(dragStart[1], dragCurrent[1]);
            double rw = Math.abs(dragCurrent[0] - dragStart[0]);
            double rh = Math.abs(dragCurrent[1] - dragStart[1]);
            gc.strokeRect(x, y, rw, rh);
        }

        if ((drawingMode == DrawingMode.ELLIPSE) && dragStart != null && dragCurrent != null) {
            double x = Math.min(dragStart[0], dragCurrent[0]);
            double y = Math.min(dragStart[1], dragCurrent[1]);
            double ew = Math.abs(dragCurrent[0] - dragStart[0]);
            double eh = Math.abs(dragCurrent[1] - dragStart[1]);
            gc.strokeOval(x, y, ew, eh);
        }

        gc.setLineDashes(null);
    }

    private static boolean pointInPolygon(double x, double y, List<double[]> vertices) {
        boolean inside = false;
        int n = vertices.size();
        for (int i = 0, j = n - 1; i < n; j = i++) {
            double xi = vertices.get(i)[0], yi = vertices.get(i)[1];
            double xj = vertices.get(j)[0], yj = vertices.get(j)[1];
            if (((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)) {
                inside = !inside;
            }
        }
        return inside;
    }
}
