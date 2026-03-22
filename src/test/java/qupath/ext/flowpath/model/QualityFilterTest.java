package qupath.ext.flowpath.model;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class QualityFilterTest {

    @Test
    void defaultFilterPassesEverything() {
        var qf = new QualityFilter();
        assertTrue(qf.passes(100, 0.5, 0.9, 5000));
    }

    @Test
    void rejectsAreaBelowMin() {
        var qf = new QualityFilter();
        qf.setMinArea(50);
        assertFalse(qf.passes(10, 0.5, 0.9, 5000));
    }

    @Test
    void rejectsAreaAboveMax() {
        var qf = new QualityFilter();
        qf.setMaxArea(200);
        assertFalse(qf.passes(300, 0.5, 0.9, 5000));
    }

    @Test
    void rejectsEccentricityAboveMax() {
        var qf = new QualityFilter();
        qf.setMaxEccentricity(0.8);
        assertFalse(qf.passes(100, 0.95, 0.9, 5000));
    }

    @Test
    void rejectsSolidityBelowMin() {
        var qf = new QualityFilter();
        qf.setMinSolidity(0.7);
        assertFalse(qf.passes(100, 0.5, 0.3, 5000));
    }

    @Test
    void rejectsTotalIntensityBelowMin() {
        var qf = new QualityFilter();
        qf.setMinTotalIntensity(1000);
        assertFalse(qf.passes(100, 0.5, 0.9, 500));
    }

    @Test
    void nanValuesAreSkipped() {
        var qf = new QualityFilter();
        qf.setMinArea(50);
        qf.setMaxEccentricity(0.8);
        qf.setMinSolidity(0.7);
        qf.setMinTotalIntensity(1000);
        // NaN values should not trigger rejection
        assertTrue(qf.passes(Double.NaN, Double.NaN, Double.NaN, Double.NaN));
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
        assertTrue(qf.passes(50, 0.8, 0.7, 1000));
        assertTrue(qf.passes(200, 0.0, 1.0, 9999));
    }
}
