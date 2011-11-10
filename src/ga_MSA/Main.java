package ga_MSA;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Date;

public class Main {

	
	public static int  maxIter = 2000;
	public static int runNo = 100;
	//public static double [] average = new double [maxIter];
	public static double [] endRun = new double [runNo];
	public static double [] endRunT = new double [runNo];
	public static double [] endRunG = new double [runNo];
	public static double successRate = 0;
	
	public static double sdFast ( double[] data )
    {
	    // sd is sqrt of sum of (values-mean) squared divided by n - 1
	    // Calculate the mean
	    double mean = 0;
	    final int n = data.length;
	    if ( n < 2 )
	       {
	       return Double.NaN;
	       }
	    for ( int i=0; i<n; i++ )
	       {
	       mean += data[i];
	       }
	    mean /= n;
	    // calculate the sum of squares
	    double sum = 0;
	    for ( int i=0; i<n; i++ )
	       {
	       final double v = data[i] - mean;
	       sum += v * v;
	       }
	    // Change to ( n - 1 ) to n if you have complete data instead of a sample.
	    return Math.sqrt( sum / ( n - 1 ) );
    }
	public static void main(String [] argv){
		try {
			//FileOutputStream ouA = new FileOutputStream(new File("./result/end-of-run-AGA.txt"));
			FileOutputStream ouB = new FileOutputStream(new File("./result/all-AGA.txt"));
			
			PrintStream all = new PrintStream(ouB);
			all.println("Average:\tSD:\tMax:\tMin:\tTime:\tSuccessRate");

			

			//Knapsack knap = new Knapsack(250, 10);
			MSA msa = new MSA();
			for (int k = 10; k < 11; k++) {
				FileOutputStream ouA = new FileOutputStream(new File("./result/end-of-run-AGA-P=" +k +".txt"));
				PrintStream  end = new PrintStream(ouA);
				successRate = 0;
				end.println("\nN = " + k + "\t *************");
				System.out.println("\nN = " + k + "\t *************");
				for (int i = 0; i < runNo; i++) {
					System.out.println("Run " + String.valueOf(i + 1));
					//GA ga = new GA(k);
					MPGA mpga = new MPGA(8,3, msa);
					Date startT = new Date();
					//ga.runGA();
					mpga.runGA();
					Date endT = new Date();
					//endRun[i] = ga.getBestF();
					endRun[i] = mpga.getBestF();
					endRunT[i]= endT.getTime()-startT.getTime();
					endRunG[i]= mpga.generation;
					end.println(endRun[i] + "\t" + endRunT[i]);
					//System.out.println(endRunT[i]);
				} 
				double aveBest = 0;
				double aveTime = 0;
				double aveGen = 0;
				double maxBest = -Double.MAX_VALUE;
				double minBest = Double.MAX_VALUE;
				
				for (int i = 0; i < endRun.length; i++) {
					if(endRun[i] > maxBest)maxBest = endRun[i];
					if(endRun[i] < minBest)minBest = endRun[i];
					aveBest += endRun[i];
					aveTime += endRunT[i];
					aveGen += endRunG[i];
				}
				end.println("\nAverage:\t" + (aveBest/runNo) + "\tSD:\t" + sdFast(endRun) + 
						"\tMax:\t" + maxBest + "\tMin:\t" + minBest  +"\tTime:\t" + (aveTime/runNo) + "\tGen:\t" + (aveGen/runNo)+ "\tSuccess Rate:" + (1.0*successRate/runNo));
				
				all.println((aveBest/runNo) + "\t" + sdFast(endRun) + "\t" + maxBest + "\t" + minBest + "\t" + (aveTime/runNo) + "\t" + (aveGen/runNo)+ "\t" + (1.0*successRate/runNo));
				System.out.println(String.valueOf((aveBest/runNo)) + "\t" + (aveTime/runNo) + "\t" + (aveGen/runNo) +"\t" + (1.0*successRate/runNo));
				System.out.println("Time SD: " + sdFast(endRunT) + "\tGeneration SD: " + sdFast(endRunG));

		}//end k
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
