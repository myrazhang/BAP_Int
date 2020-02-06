package test_mainfunc.simulation;

public class Failure {
    public StochNum downTime;
    public StochNum upTime;
    public double[] repairTimeSamples;

    // @Erica
    // todo: specify the parameters of uptime and downtime.
    // DONE: write the methods getOneSample(){...} for desired distributions of uptime and downtime in the class 'StochNum', see 'private class Normal{}' as example.


    public Failure() {
        downTime = new StochNum();
        upTime = new StochNum();
    }

    // This method generates ONLY repair time
    // Repair time can be zero, means no repair activity; or positive, means there failure occurence.
    public void repairTimeGeneration(double[] processingTimeSamples, double upRate, double downRate){
        int sampleSize = processingTimeSamples.length;
        repairTimeSamples = new double [sampleSize];

        double ttf = upTime.getOneSample(upRate, 0);
        double P = 0;

        for (int i = 1; i< sampleSize;i++){
            P = P + processingTimeSamples[i];
            if (P > ttf) {
                P = P - ttf;
                repairTimeSamples[i] = downTime.getOneSample(downRate, 0);
                ttf = upTime.getOneSample(upRate, 0);
            }
            else
                repairTimeSamples[i] = 0;
        }
    }

    public void ProctimeUpdateWithRep(double[] processingTimeSamples, double[] repairTimeSamples)
    {
        int sampleSize = processingTimeSamples.length;
        System.out.println("samplesize: "+ sampleSize);
        if(sampleSize != repairTimeSamples.length)
            throw new UnsupportedOperationException("vectors of different length");

        for (int i = 1; i < sampleSize; i ++)
        {
            processingTimeSamples[i] = processingTimeSamples[i] + repairTimeSamples[i];
            if( i < 20 )
                System.out.println("pt: " +processingTimeSamples[i]  );

        }
    }

}
