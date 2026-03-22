package qupath.ext.quflow.model;

import com.google.gson.JsonArray;
import javafx.scene.paint.Color;
import qupath.lib.common.ColorTools;

/**
 * Shared color conversion utilities.
 * Packed int format: (R << 16) | (G << 8) | B.
 */
public final class ColorUtils {

    private ColorUtils() {
    }

    /** JavaFX Color → packed int. */
    public static int colorToInt(Color c) {
        int r = (int) (c.getRed() * 255);
        int g = (int) (c.getGreen() * 255);
        int b = (int) (c.getBlue() * 255);
        return (r << 16) | (g << 8) | b;
    }

    /** Packed int → JavaFX Color. */
    public static Color intToColor(int packed) {
        int r = (packed >> 16) & 0xFF;
        int g = (packed >> 8) & 0xFF;
        int b = packed & 0xFF;
        return Color.rgb(r, g, b);
    }

    /** Packed int → [R, G, B] JSON array. */
    public static JsonArray toJsonArray(int color) {
        int r = (color >> 16) & 0xFF;
        int g = (color >> 8) & 0xFF;
        int b = color & 0xFF;
        JsonArray arr = new JsonArray();
        arr.add(r);
        arr.add(g);
        arr.add(b);
        return arr;
    }

    /** [R, G, B] JSON array → packed int. */
    public static int fromJsonArray(JsonArray arr) {
        if (arr == null || arr.size() < 3)
            return 0xFF808080; // default gray
        int r = arr.get(0).getAsInt();
        int g = arr.get(1).getAsInt();
        int b = arr.get(2).getAsInt();
        return (r << 16) | (g << 8) | b;
    }

    /** Packed int → QuPath color (via ColorTools). */
    public static int toQuPathColor(int packed) {
        return ColorTools.packRGB(
                (packed >> 16) & 0xFF,
                (packed >> 8) & 0xFF,
                packed & 0xFF);
    }
}
