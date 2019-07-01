package test_mainfunc.simulation;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.random.RandomGeneratorFactory;

import java.util.Date;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import java.util.Random;

public class StochNum {
    public String distribution;
    public double para1, para2, para3, para4;


    void iidGeneration(int sampleSize, double[] generatedSamples, int seed){

        RandomGenerator generator = RandomGeneratorFactory.createRandomGenerator(new Random());
        generator.setSeed(seed);

        double p;
        switch (this.distribution) {
            case "Norm":
                NormalDistribution normalSample = new NormalDistribution(generator, this.para1, this.para1 * this.para2);
                p = normalSample.sample();
                for (int i = 0; i <= sampleSize; i++) {
                    while (p < 0 || p > 2 * this.para1) {
                        p = normalSample.sample();
                    }
                    generatedSamples[i] = p;
                    p = normalSample.sample();
            }
            break;

            case "Beta":
                BetaDistribution betaSample=new BetaDistribution(generator, this.para1, this.para2);
                p = betaSample.sample();
                for (int i = 0; i <= sampleSize; i++) {
                    generatedSamples[i] =this.para3+ p * (this.para4 - this.para3);
                    p = betaSample.sample();
                }
                break;

            case "Exp":
                ExponentialDistribution expSample= new ExponentialDistribution(generator,this.para1);
                for (int i = 0; i <= sampleSize; i++)
                    generatedSamples[i] =expSample.sample();
                break;

            case "LogNorm":
                //para1 is the mean, para2 is the coefficient of variation
                LogNormalDistribution lognormalSample = new LogNormalDistribution(generator, this.para1, this.para1*this.para2 );
                p = lognormalSample.sample();
                for (int i = 0; i <= sampleSize; i++) {
                    generatedSamples[i] = p;
                    p = lognormalSample.sample();
                }
                break;
        }
    }
}


