package qupath.ext.flowpath.model;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class QualityFilterTest {

    @Test
    void defaultFilterPassesEverything() {
        var qf = new QualityFilter();
        assertTrue(qf.passes(100, 0.5, 0.9, 5000, 0.0));
    }

    @Test
    void rejectsAreaBelowMin() {
        var qf = new QualityFilter();
        qf.setMinArea(50);
        assertFalse(qf.passes(10, 0.5, 0.9, 5000, 0.0));
    }

    @Test
    void rejectsAreaAboveMax() {
        var qf = new QualityFilter();
        qf.setMaxArea(200);
        assertFalse(qf.passes(300, 0.5, 0.9, 5000, 0.0));
    }

    @Test
    void rejectsEccentricityAboveMax() {
        var qf = new QualityFilter();
        qf.setMaxEccentricity(0.8);
        assertFalse(qf.passes(100, 0.95, 0.9, 5000, 0.0));
    }

    @Test
    void rejectsSolidityBelowMin() {
        var qf = new QualityFilter();
        qf.setMinSolidity(0.7);
        assertFalse(qf.passes(100, 0.5, 0.3, 5000, 0.0));
    }

    @Test
    void rejectsTotalIntensityBelowMin() {
        var qf = new QualityFilter();
        qf.setMinTotalIntensity(1000);
        assertFalse(qf.passes(100, 0.5, 0.9, 500, 0.0));
    }

    @Test
    void nanValuesAreSkipped() {
        var qf = new QualityFilter();
        qf.setMinArea(50);
        qf.setMaxEccentricity(0.8);
        qf.setMinSolidity(0.7);
        qf.setMinTotalIntensity(1000);
        // NaN values should not trigger rejection
        assertTrue(qf.passes(Double.NaN, Double.NaN, Double.NaN, Double.NaN, Double.NaN));
    }

    @Test
    void boundaryValuesPass() {
        var qf = new QualityFilter();
        qf.setMinArea(50);
        qf.setMaxArea(200);
        qf.setMaxEccentricity(0.8);
        qf.setMinSolidity(0.7);
        qf.setMinTotalIntensity(1000);
        // Exact boundary values should pass
        assertTrue(qf.passes(50, 0.8, 0.7, 1000, 0.0));
        assertTrue(qf.passes(200, 0.0, 1.0, 9999, 0.0));
    }

    @Test
    void rejectsEccentricityBelowMin() {
        var qf = new QualityFilter();
        qf.setMinEccentricity(0.3);
        assertFalse(qf.passes(100, 0.1, 0.9, 5000, 0.0));
        assertTrue(qf.passes(100, 0.5, 0.9, 5000, 0.0));
    }

    @Test
    void rejectsSolidityAboveMax() {
        var qf = new QualityFilter();
        qf.setMaxSolidity(0.8);
        assertFalse(qf.passes(100, 0.5, 0.95, 5000, 0.0));
        assertTrue(qf.passes(100, 0.5, 0.7, 5000, 0.0));
    }

    @Test
    void rejectsTotalIntensityAboveMax() {
        var qf = new QualityFilter();
        qf.setMaxTotalIntensity(3000);
        assertFalse(qf.passes(100, 0.5, 0.9, 5000, 0.0));
        assertTrue(qf.passes(100, 0.5, 0.9, 2000, 0.0));
    }

    @Test
    void rejectsPerimeterBelowMin() {
        var qf = new QualityFilter();
        qf.setMinPerimeter(10);
        assertFalse(qf.passes(100, 0.5, 0.9, 5000, 5));
        assertTrue(qf.passes(100, 0.5, 0.9, 5000, 15));
    }

    @Test
    void rejectsPerimeterAboveMax() {
        var qf = new QualityFilter();
        qf.setMaxPerimeter(100);
        assertFalse(qf.passes(100, 0.5, 0.9, 5000, 150));
        assertTrue(qf.passes(100, 0.5, 0.9, 5000, 50));
    }

    @Test
    void boundaryValuesIncludeNewFields() {
        var qf = new QualityFilter();
        qf.setMinEccentricity(0.2);
        qf.setMaxSolidity(0.9);
        qf.setMaxTotalIntensity(5000);
        qf.setMinPerimeter(10);
        qf.setMaxPerimeter(200);
        assertTrue(qf.passes(100, 0.2, 0.9, 5000, 10));
        assertTrue(qf.passes(100, 1.0, 0.0, 0, 200));
    }

    @Test
    void deepCopyIncludesNewFields() {
        var qf = new QualityFilter();
        qf.setMinEccentricity(0.3);
        qf.setMaxSolidity(0.8);
        qf.setMaxTotalIntensity(3000);
        qf.setMinPerimeter(10);
        qf.setMaxPerimeter(200);
        var copy = qf.deepCopy();
        assertEquals(0.3, copy.getMinEccentricity());
        assertEquals(0.8, copy.getMaxSolidity());
        assertEquals(3000, copy.getMaxTotalIntensity());
        assertEquals(10, copy.getMinPerimeter());
        assertEquals(200, copy.getMaxPerimeter());
    }
}
