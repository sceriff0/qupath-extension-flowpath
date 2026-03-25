package qupath.ext.flowpath.model;

public class QualityFilter {

    private double minArea = 0;
    private double maxArea = Double.MAX_VALUE;
    private double minTotalIntensity = 0;
    private double maxEccentricity = 1.0;
    private double minSolidity = 0.0;
    private double minEccentricity = 0.0;
    private double maxSolidity = 1.0;
    private double maxTotalIntensity = Double.MAX_VALUE;
    private double minPerimeter = 0;
    private double maxPerimeter = Double.MAX_VALUE;

    public QualityFilter() {
    }

    public boolean passes(double area, double eccentricity, double solidity, double totalIntensity, double perimeter) {
        if (!Double.isNaN(area) && (area < minArea || area > maxArea)) return false;
        if (!Double.isNaN(eccentricity) && (eccentricity < minEccentricity || eccentricity > maxEccentricity)) return false;
        if (!Double.isNaN(solidity) && (solidity < minSolidity || solidity > maxSolidity)) return false;
        if (!Double.isNaN(totalIntensity) && (totalIntensity < minTotalIntensity || totalIntensity > maxTotalIntensity)) return false;
        if (!Double.isNaN(perimeter) && (perimeter < minPerimeter || perimeter > maxPerimeter)) return false;
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

    public double getMinEccentricity() { return minEccentricity; }
    public void setMinEccentricity(double minEccentricity) { this.minEccentricity = minEccentricity; }

    public double getMaxSolidity() { return maxSolidity; }
    public void setMaxSolidity(double maxSolidity) { this.maxSolidity = maxSolidity; }

    public double getMaxTotalIntensity() { return maxTotalIntensity; }
    public void setMaxTotalIntensity(double maxTotalIntensity) { this.maxTotalIntensity = maxTotalIntensity; }

    public double getMinPerimeter() { return minPerimeter; }
    public void setMinPerimeter(double minPerimeter) { this.minPerimeter = minPerimeter; }

    public double getMaxPerimeter() { return maxPerimeter; }
    public void setMaxPerimeter(double maxPerimeter) { this.maxPerimeter = maxPerimeter; }

    /**
     * Create a deep copy of this filter with all fields copied.
     */
    public QualityFilter deepCopy() {
        QualityFilter copy = new QualityFilter();
        copy.minArea = this.minArea;
        copy.maxArea = this.maxArea;
        copy.minTotalIntensity = this.minTotalIntensity;
        copy.maxTotalIntensity = this.maxTotalIntensity;
        copy.minEccentricity = this.minEccentricity;
        copy.maxEccentricity = this.maxEccentricity;
        copy.minSolidity = this.minSolidity;
        copy.maxSolidity = this.maxSolidity;
        copy.minPerimeter = this.minPerimeter;
        copy.maxPerimeter = this.maxPerimeter;
        return copy;
    }

}
