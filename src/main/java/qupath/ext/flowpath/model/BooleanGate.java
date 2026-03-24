package qupath.ext.flowpath.model;

import java.util.ArrayList;
import java.util.List;

/**
 * A boolean combination gate that combines other gates using AND, OR, or NOT logic.
 * <p>
 * Boolean gates don't evaluate marker thresholds directly. Instead, they reference
 * operand gates (stored as children of the match branch) and combine their results:
 * <ul>
 *   <li><b>AND</b>: Cell must pass ALL operand gates</li>
 *   <li><b>OR</b>: Cell must pass ANY operand gate</li>
 *   <li><b>NOT</b>: Cell must NOT pass the first operand gate</li>
 * </ul>
 * <p>
 * The boolean gate itself produces 2 branches: match (cells satisfying the condition)
 * and no-match (cells not satisfying it).
 */
public class BooleanGate extends GateNode {

    public enum Op { AND, OR, NOT }

    private Op operation = Op.AND;
    private List<GateNode> operands = new ArrayList<>();

    private final Branch matchBranch;
    private final Branch noMatchBranch;

    public BooleanGate() {
        int green = (0 << 16) | (200 << 8) | 0;
        int gray = (128 << 16) | (128 << 8) | 128;
        this.matchBranch = new Branch("Match", green);
        this.noMatchBranch = new Branch("No Match", gray);
    }

    public BooleanGate(Op operation, String name) {
        this.operation = operation;
        int green = (0 << 16) | (200 << 8) | 0;
        int gray = (128 << 16) | (128 << 8) | 128;
        this.matchBranch = new Branch(name + " (match)", green);
        this.noMatchBranch = new Branch(name + " (no match)", gray);
    }

    @Override
    public List<Branch> getBranches() {
        return List.of(matchBranch, noMatchBranch);
    }

    @Override
    public List<String> getChannels() {
        // Boolean gates derive channels from their operands
        List<String> channels = new ArrayList<>();
        for (GateNode operand : operands) {
            for (String ch : operand.getChannels()) {
                if (!channels.contains(ch)) channels.add(ch);
            }
        }
        return channels;
    }

    @Override
    public String getGateType() {
        return "boolean";
    }

    @Override
    public String getChannel() {
        return operation.name();
    }

    public Op getOperation() { return operation; }
    public void setOperation(Op operation) { this.operation = operation; }

    public List<GateNode> getOperands() { return operands; }
    public void setOperands(List<GateNode> operands) { this.operands = operands; }

    public Branch getMatchBranch() { return matchBranch; }
    public Branch getNoMatchBranch() { return noMatchBranch; }

    @Override
    public GateNode deepCopy() {
        BooleanGate copy = new BooleanGate();
        copy.operation = this.operation;
        copy.operands = new ArrayList<>();
        for (GateNode op : this.operands) {
            copy.operands.add(op.deepCopy());
        }
        copySharedFieldsTo(copy);
        copyBranchesTo(copy);
        return copy;
    }
}
