package qupath.ext.flowpath.model;

import java.util.ArrayList;
import java.util.List;

/**
 * A quadrant gate that splits cells into 4 populations based on thresholds
 * on two marker channels (X and Y axes).
 * <p>
 * The four quadrants are:
 * <ul>
 *   <li>Q1 (++): X >= thresholdX AND Y >= thresholdY (top-right)</li>
 *   <li>Q2 (-+): X < thresholdX AND Y >= thresholdY (top-left)</li>
 *   <li>Q3 (+-): X >= thresholdX AND Y < thresholdY (bottom-right)</li>
 *   <li>Q4 (--): X < thresholdX AND Y < thresholdY (bottom-left)</li>
 * </ul>
 */
public class QuadrantGate extends GateNode {

    private String channelX;
    private String channelY;
    private double thresholdX;
    private double thresholdY;
    private boolean thresholdIsZScore = true;

    // 4 branches: PP(++), NP(-+), PN(+-), NN(--)
    private final Branch branchPP;
    private final Branch branchNP;
    private final Branch branchPN;
    private final Branch branchNN;

    public QuadrantGate() {
        int green = (0 << 16) | (200 << 8) | 0;
        int blue = (0 << 16) | (100 << 8) | 200;
        int orange = (200 << 16) | (150 << 8) | 0;
        int gray = (128 << 16) | (128 << 8) | 128;
        this.branchPP = new Branch("", green);
        this.branchNP = new Branch("", blue);
        this.branchPN = new Branch("", orange);
        this.branchNN = new Branch("", gray);
    }

    public QuadrantGate(String channelX, String channelY) {
        this(channelX, channelY, 0.0, 0.0);
    }

    public QuadrantGate(String channelX, String channelY, double thresholdX, double thresholdY) {
        this.channelX = channelX;
        this.channelY = channelY;
        this.thresholdX = thresholdX;
        this.thresholdY = thresholdY;
        this.thresholdIsZScore = true;

        int green = (0 << 16) | (200 << 8) | 0;
        int blue = (0 << 16) | (100 << 8) | 200;
        int orange = (200 << 16) | (150 << 8) | 0;
        int gray = (128 << 16) | (128 << 8) | 128;

        this.branchPP = new Branch(channelX + "+/" + channelY + "+", green);
        this.branchNP = new Branch(channelX + "-/" + channelY + "+", blue);
        this.branchPN = new Branch(channelX + "+/" + channelY + "-", orange);
        this.branchNN = new Branch(channelX + "-/" + channelY + "-", gray);
    }

    @Override
    public List<Branch> getBranches() {
        return List.of(branchPP, branchNP, branchPN, branchNN);
    }

    @Override
    public List<String> getChannels() {
        List<String> channels = new ArrayList<>();
        if (channelX != null) channels.add(channelX);
        if (channelY != null) channels.add(channelY);
        return channels;
    }

    @Override
    public String getGateType() {
        return "quadrant";
    }

    @Override
    public String getChannel() {
        return channelX;
    }

    @Override
    public void setChannel(String channel) {
        this.channelX = channel;
    }

    // ---- Quadrant-specific getters/setters ----

    public String getChannelX() { return channelX; }
    public void setChannelX(String channelX) { this.channelX = channelX; }

    public String getChannelY() { return channelY; }
    public void setChannelY(String channelY) { this.channelY = channelY; }

    public double getThresholdX() { return thresholdX; }
    public void setThresholdX(double thresholdX) { this.thresholdX = thresholdX; }

    public double getThresholdY() { return thresholdY; }
    public void setThresholdY(double thresholdY) { this.thresholdY = thresholdY; }

    @Override
    public double getThreshold() { return thresholdX; }

    @Override
    public void setThreshold(double threshold) { this.thresholdX = threshold; }

    @Override
    public boolean isThresholdIsZScore() { return thresholdIsZScore; }

    @Override
    public void setThresholdIsZScore(boolean v) { this.thresholdIsZScore = v; }

    // Branch accessors
    public Branch getBranchPP() { return branchPP; }
    public Branch getBranchNP() { return branchNP; }
    public Branch getBranchPN() { return branchPN; }
    public Branch getBranchNN() { return branchNN; }

    /**
     * Evaluate which quadrant a cell falls into based on its X and Y values.
     * Returns the branch index: 0=PP, 1=NP, 2=PN, 3=NN.
     */
    public int evaluateQuadrant(double valueX, double valueY) {
        boolean xPos = valueX >= thresholdX;
        boolean yPos = valueY >= thresholdY;
        if (xPos && yPos) return 0;   // PP
        if (!xPos && yPos) return 1;  // NP
        if (xPos) return 2;           // PN
        return 3;                      // NN
    }

    @Override
    public GateNode deepCopy() {
        QuadrantGate copy = new QuadrantGate();
        copy.channelX = this.channelX;
        copy.channelY = this.channelY;
        copy.thresholdX = this.thresholdX;
        copy.thresholdY = this.thresholdY;
        copy.thresholdIsZScore = this.thresholdIsZScore;
        copySharedFieldsTo(copy);
        copyBranchesTo(copy);
        return copy;
    }
}
