package adduct;

public class Adduct {

    /**
     * Calculate the mass to search depending on the adduct hypothesis
     *
     * @param mz mz
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return the monoisotopic mass of the experimental mass mz with the adduct @param adduct
     */
    public static Double getMonoisotopicMassFromMZ(Double mz, String adduct) {
        if (mz == null || adduct == null || adduct.isEmpty()) return null;

        // 1. Try to get the adduct mass
        Double shift = AdductList.MAPMZPOSITIVEADDUCTS.get(adduct);
        if (shift == null) {
            shift = AdductList.MAPMZNEGATIVEADDUCTS.get(adduct);
        }
        if (shift == null) {
            System.out.println("Unknown adduct: " + adduct);
            return null;
        }
        // 2. Determine multimer (e.g., 2M = 2 molecules)
        int multimer = 1;
        if (adduct.contains("3M")) {
            multimer = 3;
        } else if (adduct.contains("2M")) {
            multimer = 2;
        }

        // 3. Determine charge (default = 1)
        int charge = 1;
        if (adduct.matches(".*3[+−-]")) {
            //.* — anything before followed by the number 3 (or 2) followed by one of these characters: +, − or -
            charge = 3;
        } else if (adduct.matches(".*2[+−-]")) {
            charge = 2;
        }

        // 4. Calculate monoisotopic mass
        double mass = (mz * charge + shift) / multimer;

        return mass;
    }

    /**
     * Calculate the mz of a monoisotopic mass with the corresponding adduct
     *
     * @param monoMass
     * @param adduct adduct name ([M+H]+, [2M+H]+, [M+2H]2+, etc..)
     *
     * @return
     */
    public static Double getMZFromMonoisotopicMass(Double monoMass, String adduct) {
        if (monoMass == null || adduct == null || adduct.isEmpty()) return null;

        // 1. Try to retrieve the adduct mass
        Double shift = AdductList.MAPMZPOSITIVEADDUCTS.get(adduct);
        if (shift == null) {
            shift = AdductList.MAPMZNEGATIVEADDUCTS.get(adduct);
        }
        if (shift == null) {
            System.out.println("Unknown adduct: " + adduct);
            return null;
        }

        // 2. Detect multimer count (default = 1)
        int multimer = 1;
        if (adduct.contains("3M")) {
            multimer = 3;
        } else if (adduct.contains("2M")) {
            multimer = 2;
        }

        // 3. Detect charge state (default = 1)
        int charge = 1;
        if (adduct.matches(".*3[+−-]")) {
            charge = 3;
        } else if (adduct.matches(".*2[+−-]")) {
            charge = 2;
        }

        // 4. Compute the theoretical m/z
        double mz = ((monoMass * multimer) + shift) / charge;
        return mz;
    }

    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * @param experimentalMass    Mass measured by MS
     * @param theoreticalMass Theoretical mass of the compound
     */
    public static int calculatePPMIncrement(Double experimentalMass, Double theoreticalMass) {
        int ppmIncrement;
        ppmIncrement = (int) Math.round(Math.abs((experimentalMass - theoreticalMass) * 1000000
                / theoreticalMass));
        return ppmIncrement;
    }

    /**
     * Returns the ppm difference between measured mass and theoretical mass
     *
     * @param experimentalMass    Mass measured by MS
     * @param ppm ppm of tolerance
     */
    public static double calculateDeltaPPM(Double experimentalMass, int ppm) {
        double deltaPPM;
        deltaPPM =  Math.round(Math.abs((experimentalMass * ppm) / 1000000));
        return deltaPPM;

    }
}
