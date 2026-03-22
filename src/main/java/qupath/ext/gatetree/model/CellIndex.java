package qupath.ext.gatetree.model;

import qupath.lib.objects.PathObject;

import java.util.Collection;
import java.util.List;
import java.util.Map;

public class CellIndex {

    private final PathObject[] objects;
    private final String[] markerNames;
    private final double[][] values; // [markerIndex][cellIndex]
    private final double[] areas;
    private final double[] eccentricities;
    private final double[] solidities;
    private final double[] totalIntensities;
    private final int size;

    private CellIndex(PathObject[] objects, String[] markerNames, double[][] values,
                      double[] areas, double[] eccentricities, double[] solidities,
                      double[] totalIntensities) {
        this.objects = objects;
        this.markerNames = markerNames;
        this.values = values;
        this.areas = areas;
        this.eccentricities = eccentricities;
        this.solidities = solidities;
        this.totalIntensities = totalIntensities;
        this.size = objects.length;
    }

    public static CellIndex build(Collection<PathObject> detections, List<String> markerNames) {
        int n = detections.size();
        int m = markerNames.size();

        PathObject[] objects = detections.toArray(new PathObject[0]);
        String[] markers = markerNames.toArray(new String[0]);
        double[][] values = new double[m][n];
        double[] areas = new double[n];
        double[] eccentricities = new double[n];
        double[] solidities = new double[n];
        double[] totalIntensities = new double[n];

        int i = 0;
        for (PathObject obj : objects) {
            Map<String, Number> measurements = getMeasurements(obj);

            double area = findMeasurement(measurements, "area");
            double convexArea = findMeasurement(measurements, "convex_area");
            double eccentricity = findMeasurement(measurements, "eccentricity");

            areas[i] = area;
            eccentricities[i] = eccentricity;
            solidities[i] = (convexArea > 0) ? area / convexArea : 1.0;

            double totalIntensity = 0;
            for (int j = 0; j < m; j++) {
                Number val = measurements.get(markers[j]);
                double v = (val != null) ? val.doubleValue() : 0.0;
                values[j][i] = v;
                totalIntensity += v;
            }
            totalIntensities[i] = totalIntensity;

            i++;
        }

        return new CellIndex(objects, markers, values, areas, eccentricities, solidities, totalIntensities);
    }

    private static Map<String, Number> getMeasurements(PathObject obj) {
        try {
            var m = obj.getMeasurements();
            if (m != null) return m;
        } catch (Exception ignored) {
        }
        return Map.of();
    }

    private static double findMeasurement(Map<String, Number> measurements, String key) {
        // Try exact match first
        Number val = measurements.get(key);
        if (val != null)
            return val.doubleValue();

        // Case-insensitive partial match
        String keyLower = key.toLowerCase();
        for (Map.Entry<String, Number> entry : measurements.entrySet()) {
            if (entry.getKey().toLowerCase().contains(keyLower) && entry.getValue() != null) {
                return entry.getValue().doubleValue();
            }
        }
        return 0.0;
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

    public double getEccentricity(int i) {
        return eccentricities[i];
    }

    public double getSolidity(int i) {
        return solidities[i];
    }

    public double getTotalIntensity(int i) {
        return totalIntensities[i];
    }
}
