package test_mainfunc;

public class StopWatch {
    private long startTime;
    public float elapseTimeSeconds=0;
    private boolean isRunning=false;

    public void start(){
        if (isRunning)
            System.err.println("Cannot start StopWatch while it is running!");
        else{
            this.startTime=System.currentTimeMillis();
            this.isRunning=true;
        }
    }

    public void restart(){
        if (isRunning)
            System.err.println("Cannot restart StopWatch while it is running!");
        else{
            this.startTime=System.currentTimeMillis();
            this.isRunning=true;
        }
    }

    public void pause(){
        if (!isRunning)
            System.err.println("Cannot pause StopWatch while it is not running!");
        else{
            this.elapseTimeSeconds=(System.currentTimeMillis()-this.startTime)/1000F;
            this.isRunning=false;
        }

    }

    public void stop(){
        if (isRunning)
            System.err.println("Cannot pause StopWatch while it is not running!");
        else{
            this.elapseTimeSeconds=(System.currentTimeMillis()-this.startTime)/1000F;
            this.isRunning=false;
        }

    }

}
