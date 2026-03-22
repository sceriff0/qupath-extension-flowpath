package qupath.ext.flowpath.model;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class MarkerStats {

    private static final int HISTOGRAM_BINS = 200;

    private final Map<String, double[]> sortedValues;
    private final Map<String, Double> means;
    private final Map<String, Double> stds;
    private final Map<String, Double> mins;
    private final Map<String, Double> maxs;
    private final Map<String, double[]> histogramBins;
    private final Map<String, double[]> histogramCounts;

    private MarkerStats(Map<String, double[]> sortedValues,
                        Map<String, Double> means, Map<String, Double> stds,
                        Map<String, Double> mins, Map<String, Double> maxs,
                        Map<String, double[]> histogramBins,
                        Map<String, double[]> histogramCounts) {
        this.sortedValues = sortedValues;
        this.means = means;
        this.stds = stds;
        this.mins = mins;
        this.maxs = maxs;
        this.histogramBins = histogramBins;
        this.histogramCounts = histogramCounts;
    }

    public static MarkerStats compute(CellIndex index, boolean[] qualityMask) {
        String[] markers = index.getMarkerNames();
        int n = index.getSize();

        // Count passing cells
        int passCount = 0;
        for (int i = 0; i < n; i++) {
            if (qualityMask[i]) passCount++;
        }

        Map<String, double[]> sortedValues = new HashMap<>();
        Map<String, Double> means = new HashMap<>();
        Map<String, Double> stds = new HashMap<>();
        Map<String, Double> mins = new HashMap<>();
        Map<String, Double> maxs = new HashMap<>();
        Map<String, double[]> histogramBins = new HashMap<>();
        Map<String, double[]> histogramCounts = new HashMap<>();

        for (int m = 0; m < markers.length; m++) {
            String name = markers[m];
            double[] raw = index.getMarkerValues(m);

            // Extract passing values
            double[] passing = new double[passCount];
            int idx = 0;
            for (int i = 0; i < n; i++) {
                if (qualityMask[i]) {
                    passing[idx++] = raw[i];
                }
            }

            Arrays.sort(passing);
            sortedValues.put(name, passing);

            if (passCount == 0) {
                means.put(name, 0.0);
                stds.put(name, 0.0);
                mins.put(name, 0.0);
                maxs.put(name, 0.0);
                histogramBins.put(name, new double[HISTOGRAM_BINS + 1]);
                histogramCounts.put(name, new double[HISTOGRAM_BINS]);
                continue;
            }

            double min = passing[0];
            double max = passing[passCount - 1];
            mins.put(name, min);
            maxs.put(name, max);

            double sum = 0;
            for (double v : passing) sum += v;
            double mean = sum / passCount;
            means.put(name, mean);

            double sumSq = 0;
            for (double v : passing) sumSq += (v - mean) * (v - mean);
            double std = Math.sqrt(sumSq / passCount);
            stds.put(name, std);

            // Histogram
            double[] bins = new double[HISTOGRAM_BINS + 1];
            double[] counts = new double[HISTOGRAM_BINS];
            double range = max - min;
            if (range < 1e-10) range = 1.0;
            double binWidth = range / HISTOGRAM_BINS;

            for (int b = 0; b <= HISTOGRAM_BINS; b++) {
                bins[b] = min + b * binWidth;
            }
            for (double v : passing) {
                int bin = (int) ((v - min) / binWidth);
                if (bin >= HISTOGRAM_BINS) bin = HISTOGRAM_BINS - 1;
                if (bin < 0) bin = 0;
                counts[bin]++;
            }

            histogramBins.put(name, bins);
            histogramCounts.put(name, counts);
        }

        return new MarkerStats(sortedValues, means, stds, mins, maxs, histogramBins, histogramCounts);
    }

    public double toZScore(String channel, double rawValue) {
        double mean = means.getOrDefault(channel, 0.0);
        double std = stds.getOrDefault(channel, 0.0);
        if (std < 1e-10) return 0.0;
        return (rawValue - mean) / std;
    }

    public double fromZScore(String channel, double zScore) {
        double mean = means.getOrDefault(channel, 0.0);
        double std = stds.getOrDefault(channel, 0.0);
        return zScore * std + mean;
    }

    public double getPercentileValue(String channel, double percentile) {
        double[] sorted = sortedValues.get(channel);
        if (sorted == null || sorted.length == 0) return 0.0;
        double idx = (percentile / 100.0) * (sorted.length - 1);
        int lo = (int) Math.floor(idx);
        int hi = (int) Math.ceil(idx);
        if (lo < 0) lo = 0;
        if (hi >= sorted.length) hi = sorted.length - 1;
        double frac = idx - lo;
        return sorted[lo] + frac * (sorted[hi] - sorted[lo]);
    }

    public double[] getHistogramBins(String channel) {
        return histogramBins.get(channel);
    }

    public double[] getHistogramCounts(String channel) {
        return histogramCounts.get(channel);
    }

    public double getMean(String channel) {
        return means.getOrDefault(channel, 0.0);
    }

    public double getStd(String channel) {
        return stds.getOrDefault(channel, 0.0);
    }

    public double getMin(String channel) {
        return mins.getOrDefault(channel, 0.0);
    }

    public double getMax(String channel) {
        return maxs.getOrDefault(channel, 0.0);
    }

    public Set<String> getMarkerNames() {
        return means.keySet();
    }
}
