package qupath.ext.quflow.model;

public class QualityFilter {

    private double minArea = 0;
    private double maxArea = Double.MAX_VALUE;
    private double minTotalIntensity = 0;
    private double maxEccentricity = 1.0;
    private double minSolidity = 0.0;
    private boolean hideFiltered = true;
    private boolean excludeFromCsv = true;

    public QualityFilter() {
    }

    public boolean passes(double area, double eccentricity, double solidity, double totalIntensity) {
        if (!Double.isNaN(area) && (area < minArea || area > maxArea)) return false;
        if (!Double.isNaN(eccentricity) && eccentricity > maxEccentricity) return false;
        if (!Double.isNaN(solidity) && solidity < minSolidity) return false;
        if (!Double.isNaN(totalIntensity) && totalIntensity < minTotalIntensity) return false;
        return true;
    }

    public double getMinArea() {
        return minArea;
    }

    public void setMinArea(double minArea) {
        this.minArea = minArea;
    }

    public double getMaxArea() {
        return maxArea;
    }

    public void setMaxArea(double maxArea) {
        this.maxArea = maxArea;
    }

    public double getMinTotalIntensity() {
        return minTotalIntensity;
    }

    public void setMinTotalIntensity(double minTotalIntensity) {
        this.minTotalIntensity = minTotalIntensity;
    }

    public double getMaxEccentricity() {
        return maxEccentricity;
    }

    public void setMaxEccentricity(double maxEccentricity) {
        this.maxEccentricity = maxEccentricity;
    }

    public double getMinSolidity() {
        return minSolidity;
    }

    public void setMinSolidity(double minSolidity) {
        this.minSolidity = minSolidity;
    }

    public boolean isHideFiltered() {
        return hideFiltered;
    }

    public void setHideFiltered(boolean hideFiltered) {
        this.hideFiltered = hideFiltered;
    }

    public boolean isExcludeFromCsv() {
        return excludeFromCsv;
    }

    public void setExcludeFromCsv(boolean excludeFromCsv) {
        this.excludeFromCsv = excludeFromCsv;
    }
}
