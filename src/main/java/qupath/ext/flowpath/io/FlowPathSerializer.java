package qupath.ext.flowpath.io;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import qupath.ext.flowpath.model.ColorUtils;
import qupath.ext.flowpath.model.GateNode;
import qupath.ext.flowpath.model.GateTree;
import qupath.ext.flowpath.model.QualityFilter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Serializes and deserializes {@link GateTree} instances to/from JSON files.
 */
public class FlowPathSerializer {

    private static final int CURRENT_VERSION = 1;

    private FlowPathSerializer() {
        // static utility class
    }

    /**
     * Save a gate tree to a JSON file.
     *
     * @param tree the gate tree to save
     * @param file the destination file
     * @throws IOException if writing fails
     */
    public static void save(GateTree tree, File file) throws IOException {
        JsonObject root = new JsonObject();
        root.addProperty("version", CURRENT_VERSION);
        root.add("qualityFilter", serializeQualityFilter(tree.getQualityFilter()));
        if (tree.getRoiFilterName() != null) {
            root.addProperty("roiFilterName", tree.getRoiFilterName());
        }
        root.add("gates", serializeNodeList(tree.getRoots()));

        Gson gson = new GsonBuilder().setPrettyPrinting().create();
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(file))) {
            writer.write(gson.toJson(root));
        }
    }

    /**
     * Load a gate tree from a JSON file.
     *
     * @param file the source file
     * @return the deserialized gate tree
     * @throws IOException if reading fails or the format is invalid
     */
    public static GateTree load(File file) throws IOException {
        JsonObject root;
        try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
            root = JsonParser.parseReader(reader).getAsJsonObject();
        }

        // Version check (currently only version 1 is supported)
        int version = root.has("version") ? root.get("version").getAsInt() : 1;
        if (version > CURRENT_VERSION) {
            throw new IOException("Unsupported gate tree version: " + version
                    + " (maximum supported: " + CURRENT_VERSION + ")");
        }

        GateTree tree = new GateTree();

        if (root.has("qualityFilter")) {
            tree.setQualityFilter(deserializeQualityFilter(root.getAsJsonObject("qualityFilter")));
        }

        if (root.has("roiFilterName")) {
            tree.setRoiFilterName(root.get("roiFilterName").getAsString());
        }

        if (root.has("gates")) {
            tree.setRoots(deserializeNodeList(root.getAsJsonArray("gates")));
        }

        return tree;
    }

    // -----------------------------------------------------------------------
    //  Quality filter
    // -----------------------------------------------------------------------

    private static JsonObject serializeQualityFilter(QualityFilter qf) {
        JsonObject obj = new JsonObject();
        obj.addProperty("minArea", qf.getMinArea());
        obj.addProperty("maxArea", qf.getMaxArea());
        obj.addProperty("minTotalIntensity", qf.getMinTotalIntensity());
        obj.addProperty("maxEccentricity", qf.getMaxEccentricity());
        obj.addProperty("minSolidity", qf.getMinSolidity());
        obj.addProperty("hideFiltered", qf.isHideFiltered());
        return obj;
    }

    private static QualityFilter deserializeQualityFilter(JsonObject obj) {
        QualityFilter qf = new QualityFilter();
        if (obj.has("minArea"))
            qf.setMinArea(obj.get("minArea").getAsDouble());
        if (obj.has("maxArea"))
            qf.setMaxArea(obj.get("maxArea").getAsDouble());
        if (obj.has("minTotalIntensity"))
            qf.setMinTotalIntensity(obj.get("minTotalIntensity").getAsDouble());
        if (obj.has("maxEccentricity"))
            qf.setMaxEccentricity(obj.get("maxEccentricity").getAsDouble());
        if (obj.has("minSolidity"))
            qf.setMinSolidity(obj.get("minSolidity").getAsDouble());
        if (obj.has("hideFiltered"))
            qf.setHideFiltered(obj.get("hideFiltered").getAsBoolean());
        return qf;
    }

    // -----------------------------------------------------------------------
    //  Gate nodes
    // -----------------------------------------------------------------------

    private static JsonArray serializeNodeList(List<GateNode> nodes) {
        JsonArray array = new JsonArray();
        for (GateNode node : nodes) {
            array.add(serializeNode(node));
        }
        return array;
    }

    private static JsonObject serializeNode(GateNode node) {
        JsonObject obj = new JsonObject();
        obj.addProperty("channel", node.getChannel());
        obj.addProperty("threshold", node.getThreshold());
        obj.addProperty("thresholdIsZScore", node.isThresholdIsZScore());
        obj.addProperty("positiveName", node.getPositiveName());
        obj.addProperty("negativeName", node.getNegativeName());
        obj.add("positiveColor", ColorUtils.toJsonArray(node.getPositiveColor()));
        obj.add("negativeColor", ColorUtils.toJsonArray(node.getNegativeColor()));
        obj.addProperty("clipPercentileLow", node.getClipPercentileLow());
        obj.addProperty("clipPercentileHigh", node.getClipPercentileHigh());
        obj.addProperty("excludeOutliers", node.isExcludeOutliers());
        obj.add("positiveChildren", serializeNodeList(node.getPositiveChildren()));
        obj.add("negativeChildren", serializeNodeList(node.getNegativeChildren()));
        return obj;
    }

    private static List<GateNode> deserializeNodeList(JsonArray array) {
        List<GateNode> nodes = new ArrayList<>();
        for (JsonElement elem : array) {
            nodes.add(deserializeNode(elem.getAsJsonObject()));
        }
        return nodes;
    }

    private static GateNode deserializeNode(JsonObject obj) {
        GateNode node = new GateNode();

        if (obj.has("channel"))
            node.setChannel(obj.get("channel").getAsString());
        if (obj.has("threshold"))
            node.setThreshold(obj.get("threshold").getAsDouble());
        if (obj.has("thresholdIsZScore"))
            node.setThresholdIsZScore(obj.get("thresholdIsZScore").getAsBoolean());
        if (obj.has("positiveName"))
            node.setPositiveName(obj.get("positiveName").getAsString());
        if (obj.has("negativeName"))
            node.setNegativeName(obj.get("negativeName").getAsString());
        if (obj.has("positiveColor"))
            node.setPositiveColor(ColorUtils.fromJsonArray(obj.getAsJsonArray("positiveColor")));
        if (obj.has("negativeColor"))
            node.setNegativeColor(ColorUtils.fromJsonArray(obj.getAsJsonArray("negativeColor")));
        if (obj.has("clipPercentileLow"))
            node.setClipPercentileLow(obj.get("clipPercentileLow").getAsDouble());
        if (obj.has("clipPercentileHigh"))
            node.setClipPercentileHigh(obj.get("clipPercentileHigh").getAsDouble());
        if (obj.has("excludeOutliers"))
            node.setExcludeOutliers(obj.get("excludeOutliers").getAsBoolean());
        else if (obj.has("hideOutliers"))
            node.setExcludeOutliers(obj.get("hideOutliers").getAsBoolean());
        if (obj.has("positiveChildren"))
            node.setPositiveChildren(deserializeNodeList(obj.getAsJsonArray("positiveChildren")));
        if (obj.has("negativeChildren"))
            node.setNegativeChildren(deserializeNodeList(obj.getAsJsonArray("negativeChildren")));

        return node;
    }

}
