package qupath.ext.flowpath.model;

import java.util.ArrayList;
import java.util.List;

/**
 * A 2D ellipse gate that classifies cells based on whether their (channelX, channelY)
 * marker values fall inside an elliptical region. Produces 2 branches: inside/outside.
 * <p>
 * Containment test: ((x-cx)/rx)^2 + ((y-cy)/ry)^2 &lt;= 1
 */
public class EllipseGate extends GateNode {

    private String channelX;
    private String channelY;
    private double centerX, centerY, radiusX, radiusY;

    private final Branch insideBranch;
    private final Branch outsideBranch;

    public EllipseGate() {
        int green = (0 << 16) | (200 << 8) | 0;
        int gray = (128 << 16) | (128 << 8) | 128;
        this.insideBranch = new Branch("Inside", green);
        this.outsideBranch = new Branch("Outside", gray);
    }

    public EllipseGate(String channelX, String channelY,
                        double centerX, double centerY, double radiusX, double radiusY) {
        this.channelX = channelX;
        this.channelY = channelY;
        this.centerX = centerX;
        this.centerY = centerY;
        this.radiusX = radiusX;
        this.radiusY = radiusY;
        int green = (0 << 16) | (200 << 8) | 0;
        int gray = (128 << 16) | (128 << 8) | 128;
        this.insideBranch = new Branch(channelX + "/" + channelY + " (in)", green);
        this.outsideBranch = new Branch(channelX + "/" + channelY + " (out)", gray);
    }

    @Override public List<Branch> getBranches() { return List.of(insideBranch, outsideBranch); }
    @Override public List<String> getChannels() {
        List<String> ch = new ArrayList<>();
        if (channelX != null) ch.add(channelX);
        if (channelY != null) ch.add(channelY);
        return ch;
    }
    @Override public String getGateType() { return "ellipse"; }
    @Override public String getChannel() { return channelX; }
    @Override public void setChannel(String channel) { this.channelX = channel; }

    public String getChannelX() { return channelX; }
    public void setChannelX(String v) { this.channelX = v; }
    public String getChannelY() { return channelY; }
    public void setChannelY(String v) { this.channelY = v; }
    public double getCenterX() { return centerX; }
    public void setCenterX(double v) { this.centerX = v; }
    public double getCenterY() { return centerY; }
    public void setCenterY(double v) { this.centerY = v; }
    public double getRadiusX() { return radiusX; }
    public void setRadiusX(double v) { this.radiusX = v; }
    public double getRadiusY() { return radiusY; }
    public void setRadiusY(double v) { this.radiusY = v; }
    public Branch getInsideBranch() { return insideBranch; }
    public Branch getOutsideBranch() { return outsideBranch; }

    public boolean contains(double x, double y) {
        if (radiusX <= 0 || radiusY <= 0) return false;
        double dx = (x - centerX) / radiusX;
        double dy = (y - centerY) / radiusY;
        return (dx * dx + dy * dy) <= 1.0;
    }

    @Override
    public GateNode deepCopy() {
        EllipseGate copy = new EllipseGate();
        copy.channelX = this.channelX;
        copy.channelY = this.channelY;
        copy.centerX = this.centerX; copy.centerY = this.centerY;
        copy.radiusX = this.radiusX; copy.radiusY = this.radiusY;
        copySharedFieldsTo(copy);
        copyBranchesTo(copy);
        return copy;
    }
}
