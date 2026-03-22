package qupath.ext.flowpath.model;

import java.util.ArrayList;
import java.util.List;

public class GateNode {

    private String channel;
    private double threshold;
    private boolean thresholdIsZScore;
    private String positiveName;
    private String negativeName;
    private int positiveColor = (0 << 16) | (200 << 8) | 0; // green
    private int negativeColor = (128 << 16) | (128 << 8) | 128; // gray
    private double clipPercentileLow = 1.0;
    private double clipPercentileHigh = 99.0;
    private boolean hideOutliers = false;
    private List<GateNode> positiveChildren = new ArrayList<>();
    private List<GateNode> negativeChildren = new ArrayList<>();

    private transient int posCount;
    private transient int negCount;

    public GateNode() {
    }

    public GateNode(String channel) {
        this(channel, 0.0);
    }

    public GateNode(String channel, double threshold) {
        this.channel = channel;
        this.threshold = threshold;
        this.thresholdIsZScore = true;
        this.positiveName = channel + "+";
        this.negativeName = channel + "-";
    }

    public String getChannel() {
        return channel;
    }

    public void setChannel(String channel) {
        this.channel = channel;
    }

    public double getThreshold() {
        return threshold;
    }

    public void setThreshold(double threshold) {
        this.threshold = threshold;
    }

    public boolean isThresholdIsZScore() {
        return thresholdIsZScore;
    }

    public void setThresholdIsZScore(boolean thresholdIsZScore) {
        this.thresholdIsZScore = thresholdIsZScore;
    }

    public String getPositiveName() {
        return positiveName;
    }

    public void setPositiveName(String positiveName) {
        this.positiveName = positiveName;
    }

    public String getNegativeName() {
        return negativeName;
    }

    public void setNegativeName(String negativeName) {
        this.negativeName = negativeName;
    }

    public int getPositiveColor() {
        return positiveColor;
    }

    public void setPositiveColor(int positiveColor) {
        this.positiveColor = positiveColor;
    }

    public int getNegativeColor() {
        return negativeColor;
    }

    public void setNegativeColor(int negativeColor) {
        this.negativeColor = negativeColor;
    }

    public double getClipPercentileLow() {
        return clipPercentileLow;
    }

    public void setClipPercentileLow(double clipPercentileLow) {
        this.clipPercentileLow = clipPercentileLow;
    }

    public double getClipPercentileHigh() {
        return clipPercentileHigh;
    }

    public void setClipPercentileHigh(double clipPercentileHigh) {
        this.clipPercentileHigh = clipPercentileHigh;
    }

    public boolean isHideOutliers() {
        return hideOutliers;
    }

    public void setHideOutliers(boolean hideOutliers) {
        this.hideOutliers = hideOutliers;
    }

    public List<GateNode> getPositiveChildren() {
        return positiveChildren;
    }

    public void setPositiveChildren(List<GateNode> positiveChildren) {
        this.positiveChildren = positiveChildren;
    }

    public List<GateNode> getNegativeChildren() {
        return negativeChildren;
    }

    public void setNegativeChildren(List<GateNode> negativeChildren) {
        this.negativeChildren = negativeChildren;
    }

    public int getPosCount() {
        return posCount;
    }

    public void setPosCount(int posCount) {
        this.posCount = posCount;
    }

    public int getNegCount() {
        return negCount;
    }

    public void setNegCount(int negCount) {
        this.negCount = negCount;
    }

    public boolean isLeaf() {
        return positiveChildren.isEmpty() && negativeChildren.isEmpty();
    }

    /**
     * Create a deep copy of this node and all descendants.
     * Transient counts are not copied (they belong to the engine run).
     */
    public GateNode deepCopy() {
        GateNode copy = new GateNode();
        copy.channel = this.channel;
        copy.threshold = this.threshold;
        copy.thresholdIsZScore = this.thresholdIsZScore;
        copy.positiveName = this.positiveName;
        copy.negativeName = this.negativeName;
        copy.positiveColor = this.positiveColor;
        copy.negativeColor = this.negativeColor;
        copy.clipPercentileLow = this.clipPercentileLow;
        copy.clipPercentileHigh = this.clipPercentileHigh;
        copy.hideOutliers = this.hideOutliers;
        copy.positiveChildren = new ArrayList<>();
        for (GateNode child : this.positiveChildren) {
            copy.positiveChildren.add(child.deepCopy());
        }
        copy.negativeChildren = new ArrayList<>();
        for (GateNode child : this.negativeChildren) {
            copy.negativeChildren.add(child.deepCopy());
        }
        return copy;
    }

    /**
     * Copy transient counts from another node into this one.
     */
    public void transferCountsFrom(GateNode source) {
        this.posCount = source.posCount;
        this.negCount = source.negCount;
    }

    public void collectLeafNames(List<String> out) {
        if (positiveChildren.isEmpty()) {
            out.add(positiveName);
        } else {
            for (GateNode child : positiveChildren) {
                child.collectLeafNames(out);
            }
        }
        if (negativeChildren.isEmpty()) {
            out.add(negativeName);
        } else {
            for (GateNode child : negativeChildren) {
                child.collectLeafNames(out);
            }
        }
    }
}
