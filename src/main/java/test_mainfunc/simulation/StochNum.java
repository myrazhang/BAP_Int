package test_mainfunc.simulation;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.random.RandomGeneratorFactory;

import java.util.Date;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import java.util.Random;

import static java.lang.Math.exp;
import static java.lang.Math.pow;

public class StochNum {
    public String distribution;
    public double para1, para2, para3, para4;
    public String para1_label, para2_label, para3_label, para4_label;


    void iidGeneration(int sampleSize, double[] generatedSamples, int seed){

        switch (this.distribution) {
            case "Norm":
                Normal normalSample=new Normal();
                normalSample.iidGeneration(sampleSize, generatedSamples, seed);
            break;

            case "Beta":
                Beta betaSample=new Beta();
                betaSample.iidGeneration(sampleSize, generatedSamples, seed);
                break;

            case "Exp":
                Exponential expSample=new Exponential();
                expSample.iidGeneration(sampleSize, generatedSamples, seed);
                break;

            case "LogNorm":
                LogNormal logNormSample=new LogNormal();
                logNormSample.iidGeneration(sampleSize, generatedSamples, seed);
                break;
            case "Deterministic":
                Deterministic deterministicSample= new Deterministic();
                deterministicSample.iidGeneration(sampleSize, generatedSamples, seed);
                break;
            default:
                throw new UnsupportedOperationException("Distribution is not supported!");

        }
    }

    public double getMean(){
        double mu=0;
        switch (this.distribution){
            case "Norm":
                Normal normalSample=new Normal();
                return normalSample.getMean();

            case "Beta":
                Beta betaSample=new Beta();
                return betaSample.getMean();

            case "Exp":
                Exponential expSample=new Exponential();
                return expSample.getMean();

            case "LogNorm":
                LogNormal logNormSample=new LogNormal();
                return logNormSample.getMean();
            case "Deterministic":
                Deterministic deterministicSample=new Deterministic();
                return  deterministicSample.getMean();

            default:
                throw new UnsupportedOperationException("Distribution is not supported!");
        }

    }

    double getOneSample(double rate, int seed){

            RandomGenerator generator = RandomGeneratorFactory.createRandomGenerator(new Random());
          //  generator.setSeed(seed);

            double p;
            ExponentialDistribution expSample= new ExponentialDistribution(generator,rate);
            double generatedSamples;
            generatedSamples = expSample.sample();
            return generatedSamples;
    }

    private class Normal{
        void iidGeneration(int sampleSize, double[] generatedSamples, int seed){
            RandomGenerator generator = RandomGeneratorFactory.createRandomGenerator(new Random());
            generator.setSeed(seed);

            double p;
            NormalDistribution normalSample = new NormalDistribution(generator, para1, para1 * para2);
            p = normalSample.sample();
            for (int i = 0; i <= sampleSize; i++) {
                while (p < 0 || p > 2 * para1) {
                    p = normalSample.sample();
                }
                generatedSamples[i] = p;
                p = normalSample.sample();
            }
        }

        double getMean(){
            return para1;
        }

        double getOneSample(){
            RandomGenerator generator = RandomGeneratorFactory.createRandomGenerator(new Random());
            NormalDistribution normalSample = new NormalDistribution(generator, para1, para1 * para2);
            return normalSample.sample();
        }
    }

    private class Exponential{
        void iidGeneration(int sampleSize, double[] generatedSamples, int seed){
            RandomGenerator generator = RandomGeneratorFactory.createRandomGenerator(new Random());
            generator.setSeed(seed);

            double p;
            ExponentialDistribution expSample= new ExponentialDistribution(generator,para1);
            for (int i = 0; i <= sampleSize; i++)
                generatedSamples[i] =expSample.sample();
        }

        double getMean(){
            return para1;
        }
    }

    private class Deterministic{
        void iidGeneration(int sampleSize, double[] generatedSamples, int seed){
            for (int i = 0; i <= sampleSize; i++)
                generatedSamples[i] = para1;
        }

        double getMean(){
            return para1;
        }
    }

    private class Beta{
        void iidGeneration(int sampleSize, double[] generatedSamples, int seed){
            RandomGenerator generator = RandomGeneratorFactory.createRandomGenerator(new Random());
            generator.setSeed(seed);

            double p;
            BetaDistribution betaSample=new BetaDistribution(generator, para1, para2);
            p = betaSample.sample();
            for (int i = 0; i <= sampleSize; i++) {
                generatedSamples[i] =para3+ p * (para4 - para3);
                p = betaSample.sample();
            }
        }

        double getMean(){
            return para1/(para1+para2);
        }
    }

    private class LogNormal{
        //para1 is the mean, para2 is the coefficient of variation of the NORMAL VARIATE

        void iidGeneration(int sampleSize, double[] generatedSamples, int seed){
            RandomGenerator generator = RandomGeneratorFactory.createRandomGenerator(new Random());
            generator.setSeed(seed);

            double p;
            LogNormalDistribution lognormalSample = new LogNormalDistribution(generator, para1, para1*para2 );
            p = lognormalSample.sample();
            for (int i = 0; i <= sampleSize; i++) {
                generatedSamples[i] = p;
                p = lognormalSample.sample();
            }
        }

        double getMean(){
            return exp(para1+pow(para1*para2,2)/2);
        }

    }
}


