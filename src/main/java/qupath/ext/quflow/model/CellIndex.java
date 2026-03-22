package qupath.ext.quflow.model;

import qupath.lib.objects.PathObject;

import java.util.Collection;
import java.util.List;
import java.util.Map;

public class CellIndex {

    private final PathObject[] objects;
    private final String[] markerNames;
    private final double[][] values; // [markerIndex][cellIndex]
    private final double[] areas;
    private final double[] perimeters;
    private final double[] eccentricities;
    private final double[] solidities;
    private final double[] totalIntensities;
    private final double[] centroidX;
    private final double[] centroidY;
    private final int size;

    private CellIndex(PathObject[] objects, String[] markerNames, double[][] values,
                      double[] areas, double[] perimeters, double[] eccentricities,
                      double[] solidities, double[] totalIntensities,
                      double[] centroidX, double[] centroidY) {
        this.objects = objects;
        this.markerNames = markerNames;
        this.values = values;
        this.areas = areas;
        this.perimeters = perimeters;
        this.eccentricities = eccentricities;
        this.solidities = solidities;
        this.totalIntensities = totalIntensities;
        this.centroidX = centroidX;
        this.centroidY = centroidY;
        this.size = objects.length;
    }

    public static CellIndex build(Collection<PathObject> detections, List<String> markerNames) {
        int n = detections.size();
        int m = markerNames.size();

        PathObject[] objects = detections.toArray(new PathObject[0]);
        String[] markers = markerNames.toArray(new String[0]);
        double[][] values = new double[m][n];
        double[] areas = new double[n];
        double[] perimeters = new double[n];
        double[] eccentricities = new double[n];
        double[] solidities = new double[n];
        double[] totalIntensities = new double[n];
        double[] centroidX = new double[n];
        double[] centroidY = new double[n];

        int i = 0;
        for (PathObject obj : objects) {
            Map<String, Number> measurements = getMeasurements(obj);

            double area = findMeasurement(measurements, "area");
            double convexArea = findMeasurement(measurements, "convex_area");
            double eccentricity = findMeasurement(measurements, "eccentricity");
            double perimeter = findMeasurement(measurements, "perimeter");

            areas[i] = area;
            perimeters[i] = perimeter;
            eccentricities[i] = eccentricity;
            // Solidity = area / convex_area; NaN if either is missing
            if (!Double.isNaN(area) && !Double.isNaN(convexArea) && convexArea > 0) {
                solidities[i] = area / convexArea;
            } else {
                solidities[i] = Double.NaN;
            }

            // Spatial coordinates
            centroidX[i] = findMeasurement(measurements, "Centroid X");
            centroidY[i] = findMeasurement(measurements, "Centroid Y");

            double totalIntensity = 0;
            for (int j = 0; j < m; j++) {
                double v = findMarkerValue(measurements, markers[j]);
                values[j][i] = v;
                totalIntensity += v;
            }
            totalIntensities[i] = totalIntensity;

            i++;
        }

        return new CellIndex(objects, markers, values, areas, perimeters, eccentricities,
                solidities, totalIntensities, centroidX, centroidY);
    }

    private static Map<String, Number> getMeasurements(PathObject obj) {
        try {
            var m = obj.getMeasurements();
            if (m != null) return m;
        } catch (Exception ignored) {
        }
        return Map.of();
    }

    /**
     * Find a marker intensity value by channel name.
     * Tries exact match first, then looks for "[layer] channel" patterns
     * (from import_phenotype.groovy layer-prefixed measurements).
     */
    private static double findMarkerValue(Map<String, Number> measurements, String channel) {
        // Exact match
        Number val = measurements.get(channel);
        if (val != null) return val.doubleValue();

        // Layer-prefixed match: "[something] channel"
        String suffix = "] " + channel;
        for (Map.Entry<String, Number> entry : measurements.entrySet()) {
            if (entry.getKey().endsWith(suffix) && entry.getValue() != null) {
                return entry.getValue().doubleValue();
            }
        }
        return 0.0;
    }

    /**
     * Find a morphological measurement by key name.
     * Tries exact match, then layer-prefixed "[layer] key", then prefix match
     * (e.g., "area" matches "area µm²"). Returns NaN if not found.
     */
    private static double findMeasurement(Map<String, Number> measurements, String key) {
        Number val = measurements.get(key);
        if (val != null) return val.doubleValue();

        // Layer-prefixed match: "[layer] key"
        String suffix = "] " + key;
        for (Map.Entry<String, Number> entry : measurements.entrySet()) {
            if (entry.getKey().endsWith(suffix) && entry.getValue() != null) {
                return entry.getValue().doubleValue();
            }
        }

        // Prefix match: "area" matches "area µm²"
        for (Map.Entry<String, Number> entry : measurements.entrySet()) {
            if (entry.getKey().startsWith(key) && entry.getValue() != null) {
                return entry.getValue().doubleValue();
            }
        }

        return Double.NaN;
    }

    public double[] getMarkerValues(int markerIndex) {
        return values[markerIndex];
    }

    public int getMarkerIndex(String name) {
        for (int i = 0; i < markerNames.length; i++) {
            if (markerNames[i].equals(name))
                return i;
        }
        return -1;
    }

    public PathObject getObject(int cellIndex) {
        return objects[cellIndex];
    }

    public PathObject[] getObjects() {
        return objects;
    }

    public String[] getMarkerNames() {
        return markerNames;
    }

    public double[] getAreas() {
        return areas;
    }

    public double[] getPerimeters() {
        return perimeters;
    }

    public double[] getEccentricities() {
        return eccentricities;
    }

    public double[] getSolidities() {
        return solidities;
    }

    public double[] getTotalIntensities() {
        return totalIntensities;
    }

    public int getSize() {
        return size;
    }

    public int size() {
        return size;
    }

    public double getArea(int i) {
        return areas[i];
    }

    public double getPerimeter(int i) {
        return perimeters[i];
    }

    public double getEccentricity(int i) {
        return eccentricities[i];
    }

    public double getSolidity(int i) {
        return solidities[i];
    }

    public double getTotalIntensity(int i) {
        return totalIntensities[i];
    }

    public double getCentroidX(int i) {
        return centroidX[i];
    }

    public double getCentroidY(int i) {
        return centroidY[i];
    }
}
