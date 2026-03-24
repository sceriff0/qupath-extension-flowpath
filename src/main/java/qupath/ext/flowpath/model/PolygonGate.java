package qupath.ext.flowpath.model;

import java.util.ArrayList;
import java.util.List;

/**
 * A 2D polygon gate that classifies cells based on whether their (channelX, channelY)
 * marker values fall inside a user-drawn polygon. Produces 2 branches: inside/outside.
 */
public class PolygonGate extends GateNode {

    private String channelX;
    private String channelY;
    private List<double[]> vertices = new ArrayList<>(); // [[x0,y0], [x1,y1], ...]

    private final Branch insideBranch;
    private final Branch outsideBranch;

    public PolygonGate() {
        int green = (0 << 16) | (200 << 8) | 0;
        int gray = (128 << 16) | (128 << 8) | 128;
        this.insideBranch = new Branch("Inside", green);
        this.outsideBranch = new Branch("Outside", gray);
    }

    public PolygonGate(String channelX, String channelY) {
        this.channelX = channelX;
        this.channelY = channelY;
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
    @Override public String getGateType() { return "polygon"; }
    @Override public String getChannel() { return channelX; }
    @Override public void setChannel(String channel) { this.channelX = channel; }

    public String getChannelX() { return channelX; }
    public void setChannelX(String v) { this.channelX = v; }
    public String getChannelY() { return channelY; }
    public void setChannelY(String v) { this.channelY = v; }
    public List<double[]> getVertices() { return vertices; }
    public void setVertices(List<double[]> v) { this.vertices = v; }
    public Branch getInsideBranch() { return insideBranch; }
    public Branch getOutsideBranch() { return outsideBranch; }

    /**
     * Point-in-polygon test using ray casting algorithm.
     */
    public boolean contains(double x, double y) {
        if (vertices.size() < 3) return false;
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

    @Override
    public GateNode deepCopy() {
        PolygonGate copy = new PolygonGate();
        copy.channelX = this.channelX;
        copy.channelY = this.channelY;
        copy.vertices = new ArrayList<>();
        for (double[] v : this.vertices) {
            copy.vertices.add(new double[]{v[0], v[1]});
        }
        copySharedFieldsTo(copy);
        copyBranchesTo(copy);
        return copy;
    }
}
