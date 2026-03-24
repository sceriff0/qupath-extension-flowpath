package qupath.ext.flowpath.model;

public class QualityFilter {

    private double minArea = 0;
    private double maxArea = Double.MAX_VALUE;
    private double minTotalIntensity = 0;
    private double maxEccentricity = 1.0;
    private double minSolidity = 0.0;

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

    /**
     * Create a deep copy of this filter with all fields copied.
     */
    public QualityFilter deepCopy() {
        QualityFilter copy = new QualityFilter();
        copy.minArea = this.minArea;
        copy.maxArea = this.maxArea;
        copy.minTotalIntensity = this.minTotalIntensity;
        copy.maxEccentricity = this.maxEccentricity;
        copy.minSolidity = this.minSolidity;
        return copy;
    }

}
