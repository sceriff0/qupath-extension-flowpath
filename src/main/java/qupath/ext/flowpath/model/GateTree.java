package qupath.ext.flowpath.model;

import java.util.ArrayList;
import java.util.List;

public class GateTree {

    private List<GateNode> roots = new ArrayList<>();
    private QualityFilter qualityFilter = new QualityFilter();
    private String roiFilterId;

    public GateTree() {
    }

    public void addRoot(GateNode node) {
        roots.add(node);
    }

    public void removeRoot(GateNode node) {
        roots.remove(node);
    }

    public List<GateNode> getRoots() {
        return roots;
    }

    public void setRoots(List<GateNode> roots) {
        this.roots = roots;
    }

    public QualityFilter getQualityFilter() {
        return qualityFilter;
    }

    public void setQualityFilter(QualityFilter qualityFilter) {
        this.qualityFilter = qualityFilter;
    }

    public String getRoiFilterId() {
        return roiFilterId;
    }

    public void setRoiFilterId(String roiFilterId) {
        this.roiFilterId = roiFilterId;
    }

    /**
     * Create a deep copy of this tree (all nodes and the quality filter are cloned).
     */
    public GateTree deepCopy() {
        GateTree copy = new GateTree();
        copy.qualityFilter = this.qualityFilter.deepCopy();
        copy.roiFilterId = this.roiFilterId;
        copy.roots = new ArrayList<>();
        for (GateNode root : this.roots) {
            copy.roots.add(root.deepCopy());
        }
        return copy;
    }

    /**
     * Transfer transient counts from a copy's nodes back to the originals.
     * Walks both trees in parallel; stops gracefully if structures differ.
     */
    public static void transferCounts(List<GateNode> originals, List<GateNode> copies) {
        for (int i = 0; i < originals.size() && i < copies.size(); i++) {
            originals.get(i).transferCountsFrom(copies.get(i));
            List<Branch> origBranches = originals.get(i).getBranches();
            List<Branch> copyBranches = copies.get(i).getBranches();
            for (int b = 0; b < origBranches.size() && b < copyBranches.size(); b++) {
                transferCounts(origBranches.get(b).getChildren(), copyBranches.get(b).getChildren());
            }
        }
    }

    public List<String> collectLeafNames() {
        List<String> names = new ArrayList<>();
        for (GateNode root : roots) {
            root.collectLeafNames(names);
        }
        return names;
    }
}
