package ga_MSA;

public class MSA_Test {
	
	public static void main(String[] args) {
		MSA msa = new MSA();
		
		System.out.println("Run Started");
		System.out.println("maxL = " + msa.maxL);
		MPGA mpga = new MPGA(msa.maxL,msa.maxW, msa);
		System.out.println("Best Individuals:");
		for (int i = 0; i < mpga.elites.length; i++)
			System.out.println((i+1) + ": " + mpga.elites[i]);
		System.out.println("Best Individual:");
		printElite(mpga, msa);
		System.out.println("Running Adaptive MPGA....");
		double lastBest = mpga.eliteF;
		for (int k = 1; k < 501; k++) {
			mpga.singleGARun();
			//if (k%10 != 0 && k > 11) continue;
			/*System.out.println("Best Individuals in run " + k);
			for (int i = 0; i < mpga.elites.length; i++) {
				System.out.println((i) + ": " + mpga.elites[i] + "   Avg: " + 
						mpga.avgFit[i] + "    Min: " + mpga.minFit[i] + "   Max: " + mpga.maxFit[i] + 
						"     Pop: " + mpga.fitness.get(i).size());
			}*/
			//System.out.println("Best Individual:");
			//printElite(mpga, msa);
			//if (mpga.eliteF > lastBest)
				System.out.println("Best Score so far in run " + k + " = " + mpga.eliteF);
			lastBest = mpga.eliteF;
		}
		System.out.println("Best Individual:");
		printElite(mpga, msa);
		
	}

	private static void printElite(MPGA mpga, MSA msa) {
		Integer[][] indiv = mpga.eliteIndiv;
		for (int l = 0; l < indiv.length; l++) {
			int size = -1;
			for (int c = 0; c < indiv[0].length; c++) {
				char c1;
				if (indiv[l][c] == 1) {
					size++;
					c1 = msa.getChar(l, size);
				} else
					c1 = '-';
				System.out.print(c1);
			}
			System.out.println();
		}
		System.out.println("Score = " + mpga.eliteF);
	}
}

