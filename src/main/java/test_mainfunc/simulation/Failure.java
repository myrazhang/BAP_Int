package test_mainfunc.simulation;

public class Failure {
    public StochNum downTime;
    public StochNum upTime;

    // @Erica
    // todo: specify the parameters of uptime and downtime.
    // todo: write the methods getOneSample(){...} for desired distributions of uptime and downtime in the class 'StochNum', see 'private class Normal{}' as example.


    public Failure() {
        downTime = new StochNum();
        upTime = new StochNum();
    }

    // This method generates ONLY repair time
    // Repair time can be zero, means no repair activity; or positive, means there failure occurance.
    public void repairTimeGeneration(double[] processingTimeSamples, double[] repairTimeSamples){
        int sampleSize = processingTimeSamples.length;
        repairTimeSamples = new double [sampleSize];

        double ttf = upTime.getOneSample();
        double P = 0;

        for (int i = 1; i< sampleSize;i++){
            P = P + processingTimeSamples[i];
            if (P > ttf) {
                P = P - ttf;
                repairTimeSamples[i] = downTime.getOneSample();
                ttf = upTime.getOneSample();
            }
            else
                repairTimeSamples[i] = 0;
        }
    }

}
