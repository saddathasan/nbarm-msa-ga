package ga;

public class knapsackTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		Knapsack knap = new Knapsack(250, 10);
		//knap.printBest();
		System.out.println("Run Started");
		MPGA mpga = new MPGA(250, knap);
		System.out.println("Best Individuals:");
		for (int i = 0; i < mpga.elites.length; i++)
			System.out.println((i+1) + ": " + mpga.elites[i][0]);
		//System.out.println("Best Individual = " + mpga.eliteF);
		//System.out.println("Applying repair...");
		//mpga.repairAll();
		//System.out.println("Best Individuals:");
		//for (int i = 0; i < mpga.elites.length; i++)
			//System.out.println((i+1) + ": " + mpga.elites[i][0]);
		System.out.println("Best Individual = " + mpga.eliteF + " with weight = " + mpga.eliteW);
		System.out.println("Running Adaptive MPGA....");
		for (int k = 1; k < 1001; k++) {
			mpga.singleGARun();
			//if (k%10 != 0 && k > 11) continue;
			System.out.println("Best Individuals in run " + k);
			for (int i = 0; i < mpga.elites.length; i++) {
				System.out.println((i) + ": " + mpga.elites[i][0] + "   Avg: " + 
						mpga.avgFit[i] + "    Min: " + mpga.minFit[i] + "   Max: " + mpga.maxFit[i] + 
						"     Pop: " + mpga.fitness.get(i).size());
			}
			System.out.println("Best Individual = " + mpga.eliteF + " with weight = " + mpga.eliteW);
		}
		/*Integer[] number = new Integer [20];
		Integer[] order = new Integer [20];
		for (int i = 0; i < 20; i++) {
			order[i] = i;
			number[i] = i * i;
		}
		mpga.sort(number, order);
		for (int i = 0; i < 20; i++) {
			System.out.println("order " + i + " = " + number[order[i]]);
		}*/
	}

}
