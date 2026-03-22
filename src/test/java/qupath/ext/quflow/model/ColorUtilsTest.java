package qupath.ext.quflow.model;

import com.google.gson.JsonArray;
import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class ColorUtilsTest {

    @Test
    void toJsonArrayAndBack() {
        int color = (200 << 16) | (100 << 8) | 50; // R=200, G=100, B=50
        JsonArray arr = ColorUtils.toJsonArray(color);

        assertEquals(3, arr.size());
        assertEquals(200, arr.get(0).getAsInt());
        assertEquals(100, arr.get(1).getAsInt());
        assertEquals(50, arr.get(2).getAsInt());

        int roundTripped = ColorUtils.fromJsonArray(arr);
        assertEquals(color, roundTripped);
    }

    @Test
    void toJsonArrayBlack() {
        JsonArray arr = ColorUtils.toJsonArray(0);
        assertEquals(0, arr.get(0).getAsInt());
        assertEquals(0, arr.get(1).getAsInt());
        assertEquals(0, arr.get(2).getAsInt());
    }

    @Test
    void toJsonArrayWhite() {
        int white = (255 << 16) | (255 << 8) | 255;
        JsonArray arr = ColorUtils.toJsonArray(white);
        assertEquals(255, arr.get(0).getAsInt());
        assertEquals(255, arr.get(1).getAsInt());
        assertEquals(255, arr.get(2).getAsInt());
    }

    @Test
    void fromJsonArrayNullReturnsDefaultGray() {
        int result = ColorUtils.fromJsonArray(null);
        assertEquals(0xFF808080, result);
    }

    @Test
    void fromJsonArrayTooShortReturnsDefaultGray() {
        var arr = new JsonArray();
        arr.add(255);
        arr.add(0);
        // Only 2 elements
        assertEquals(0xFF808080, ColorUtils.fromJsonArray(arr));
    }

    @Test
    void fromJsonArrayExtractsRgb() {
        var arr = new JsonArray();
        arr.add(10);
        arr.add(20);
        arr.add(30);
        int expected = (10 << 16) | (20 << 8) | 30;
        assertEquals(expected, ColorUtils.fromJsonArray(arr));
    }
}
