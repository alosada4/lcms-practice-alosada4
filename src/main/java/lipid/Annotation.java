package lipid;

import adduct.Adduct;
import adduct.AdductList;

import java.util.*;

/**
 * Class to represent the annotation over a lipid
 */
public class Annotation {

    private final int PPM_TOLERANCE = 10;
    private final Lipid lipid;
    private final double mz;
    private final double intensity; // intensity of the most abundant peak in the groupedPeaks
    private final double rtMin;
    private final IonizationMode ionizationMode;
    private String adduct;
    private final Set<Peak> groupedSignals;
    private int score;
    private int totalScoresApplied;


    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IonizationMode ionizationMode) {
        this(lipid, mz, intensity, retentionTime, ionizationMode, Collections.emptySet());
    }

    /**
     * @param lipid
     * @param mz
     * @param intensity
     * @param retentionTime
     * @param ionizationMode
     * @param groupedSignals
     */
    public Annotation(Lipid lipid, double mz, double intensity, double retentionTime, IonizationMode ionizationMode, Set<Peak> groupedSignals) {
        this.lipid = lipid;
        this.mz = mz;
        this.rtMin = retentionTime;
        this.intensity = intensity;
        this.ionizationMode = ionizationMode;
        this.groupedSignals = new TreeSet<>(groupedSignals);
        this.score = 0;
        this.totalScoresApplied = 0;
        detectAdduct();
    }

    public Lipid getLipid() {
        return lipid;
    }

    public double getMz() {
        return mz;
    }

    public double getRtMin() {
        return rtMin;
    }

    public String getAdduct() {
        return adduct;
    }

    public void setAdduct(String adduct) {
        this.adduct = adduct;
    }

    public double getIntensity() {
        return intensity;
    }

    public IonizationMode getIonizationMode() {
        return ionizationMode;
    }

    public Set<Peak> getGroupedSignals() {
        return Collections.unmodifiableSet(groupedSignals);
    }


    public int getScore() {
        return score;
    }

    public void setScore(int score) {
        this.score = score;
    }

    // !CHECK Take into account that the score should be normalized between -1 and 1
    public void addScore(int delta) {
        this.score += delta;
        this.totalScoresApplied++;
    }

    /**
     * @return The normalized score between 0 and 1 that consists on the final number divided into the times that the rule
     * has been applied.
     */
    public double getNormalizedScore() {
        return (double) this.score / this.totalScoresApplied;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Annotation)) return false;
        Annotation that = (Annotation) o;
        return Double.compare(that.mz, mz) == 0 &&
                Double.compare(that.rtMin, rtMin) == 0 &&
                Objects.equals(lipid, that.lipid);
    }

    @Override
    public int hashCode() {
        return Objects.hash(lipid, mz, rtMin);
    }

    @Override
    public String toString() {
        return String.format("Annotation(%s, mz=%.4f, RT=%.2f, adduct=%s, intensity=%.1f, score=%d)",
                lipid.getName(), mz, rtMin, adduct, intensity, score);
    }

    /**
     * Infers the adduct by comparing grouped peak m/z values and matching their monoisotopic masses.
     */
    public void detectAdduct() {
        if (this.groupedSignals == null || this.groupedSignals.size() < 2) return;

        List<Peak> peaks = new ArrayList<>(this.groupedSignals);
        Map<String, Double> adducts = (this.ionizationMode == IonizationMode.POSITIVE)
                ? AdductList.MAPMZPOSITIVEADDUCTS
                : AdductList.MAPMZNEGATIVEADDUCTS;

        for (String adduct1 : adducts.keySet()) {
            for (String adduct2 : adducts.keySet()) {
                if (adduct1.equals(adduct2)) continue;

                for (Peak peak1 : peaks) {
                    for (Peak peak2 : peaks) {
                        if (peak1.equals(peak2)) continue;

                        double mz1 = peak1.getMz();
                        double mz2 = peak2.getMz();

                        Double mass1 = Adduct.getMonoisotopicMassFromMZ(mz1, adduct1);
                        Double mass2 = Adduct.getMonoisotopicMassFromMZ(mz2, adduct2);

                        if (mass1 == null || mass2 == null) continue;

                        int ppmMassDiff = Adduct.calculatePPMIncrement(mass1, mass2);
                        if (ppmMassDiff > PPM_TOLERANCE) continue;

                        int ppm1 = Adduct.calculatePPMIncrement(this.mz, mz1);
                        int ppm2 = Adduct.calculatePPMIncrement(this.mz, mz2);

                        if (ppm1 <= PPM_TOLERANCE) {
                            this.adduct = adduct1;
                            //System.out.println("Adduct detected from pair: " + adduct1);
                            return;
                        } else if (ppm2 <= PPM_TOLERANCE) {
                            this.adduct = adduct2;
                            //System.out.println("Adduct detected from pair: " + adduct2);
                            return;
                        }
                    }
                }
            }
        }

        // If no adduct was detected from peak pairs, try detecting from this.mz alone
        for (String adduct : adducts.keySet()) {
            Double monoMass = Adduct.getMonoisotopicMassFromMZ(this.mz, adduct);
            if (monoMass != null) {
                Double expectedMz = Adduct.getMZFromMonoisotopicMass(monoMass, adduct);
                if (expectedMz != null) {
                    int ppm = Adduct.calculatePPMIncrement(this.mz, expectedMz);
                    if (ppm <= PPM_TOLERANCE) {
                        this.adduct = adduct;
                        //System.out.println("Adduct detected from direct match: " + adduct);
                        return;
                    }
                }
            }
        }
    }
}


