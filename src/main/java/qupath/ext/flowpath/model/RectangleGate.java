package qupath.ext.flowpath.model;

import java.util.ArrayList;
import java.util.List;

/**
 * A 2D rectangle gate that classifies cells based on whether their (channelX, channelY)
 * marker values fall inside a rectangular region. Produces 2 branches: inside/outside.
 */
public class RectangleGate extends GateNode {

    private String channelX;
    private String channelY;
    private double minX, maxX, minY, maxY;

    private final Branch insideBranch;
    private final Branch outsideBranch;

    public RectangleGate() {
        int green = (0 << 16) | (200 << 8) | 0;
        int gray = (128 << 16) | (128 << 8) | 128;
        this.insideBranch = new Branch("Inside", green);
        this.outsideBranch = new Branch("Outside", gray);
    }

    public RectangleGate(String channelX, String channelY, double minX, double maxX, double minY, double maxY) {
        this.channelX = channelX;
        this.channelY = channelY;
        this.minX = minX;
        this.maxX = maxX;
        this.minY = minY;
        this.maxY = maxY;
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
    @Override public String getGateType() { return "rectangle"; }
    @Override public String getChannel() { return channelX; }
    @Override public void setChannel(String channel) { this.channelX = channel; }

    public String getChannelX() { return channelX; }
    public void setChannelX(String v) { this.channelX = v; }
    public String getChannelY() { return channelY; }
    public void setChannelY(String v) { this.channelY = v; }
    public double getMinX() { return minX; }
    public void setMinX(double v) { this.minX = v; }
    public double getMaxX() { return maxX; }
    public void setMaxX(double v) { this.maxX = v; }
    public double getMinY() { return minY; }
    public void setMinY(double v) { this.minY = v; }
    public double getMaxY() { return maxY; }
    public void setMaxY(double v) { this.maxY = v; }
    public Branch getInsideBranch() { return insideBranch; }
    public Branch getOutsideBranch() { return outsideBranch; }

    public boolean contains(double x, double y) {
        return x >= minX && x <= maxX && y >= minY && y <= maxY;
    }

    @Override
    public GateNode deepCopy() {
        RectangleGate copy = new RectangleGate();
        copy.channelX = this.channelX;
        copy.channelY = this.channelY;
        copy.minX = this.minX; copy.maxX = this.maxX;
        copy.minY = this.minY; copy.maxY = this.maxY;
        copySharedFieldsTo(copy);
        copyBranchesTo(copy);
        return copy;
    }
}
