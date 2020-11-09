package test_mainfunc.util;

public class Stopwatch {
    public float elapseTimeSeconds;
    private long startTime;
    private boolean isRunning=false;

    public void start(){
        if (isRunning){
            System.err.println("Cannot start Stopwatch while it is running!");

        }
        else{
            elapseTimeSeconds=0;
            this.startTime=System.currentTimeMillis();
            this.isRunning=true;
        }
    }

    public void restart(){
        if (isRunning){
            System.err.println("Cannot restart Stopwatch while it is running!");
        }
        else{
            this.startTime=System.currentTimeMillis();
            this.isRunning=true;
        }
    }

    public void pause(){
        if (!isRunning)
            System.err.println("Cannot pause Stopwatch while it is not running!");
        else{
            this.elapseTimeSeconds+=(System.currentTimeMillis()-this.startTime)/1000F;
            this.isRunning=false;
        }

    }

    public void stop(){
        if (!isRunning) {
            System.err.println("Cannot stop Stopwatch while it is not running!");
            new Exception("").printStackTrace();
        }
        else{
            this.elapseTimeSeconds+=(System.currentTimeMillis()-this.startTime)/1000F;
            this.isRunning=false;
        }

    }

}
