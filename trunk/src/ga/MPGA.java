package ga;

import java.util.Vector;

public class MPGA {
	
	private double mutation = 0.005;
	private double cross = 0.025;
	private double[] migrationRate;
	private int migrationInterval = 1;
	private double[] new_FG; //new fitness gain
	private double[] old_FG;
	public int generation;
	public static int l = 250;					
	public static int n = 500;				
	public static int p = 5;
	//public int [][] populationGA;
	public Vector<Vector<Integer[]>> GA;
	public Vector<Vector<Double>> fitness;
	public Vector<Vector<Integer[]>> GA_New;
	public Vector<Vector<Double>> fitness_New;
	public double eliteF;
	public double eliteW;
	public double [][] elites;
	private double [][] old_elites;
	public int optimal = 200;	// Set the optimum here	
	private static Knapsack knap; //knapsack problem
	private Integer[][][] populationNew;
	private Double[][] fitnessNew;
	private Integer[][] order;
	public double[] avgFit;
	public double[] minFit;
	public double[] maxFit;
	
	public MPGA(int N, Knapsack K) {
		//populationGA = new int [n][l];
		l = N;
		GA = new Vector<Vector<Integer[]>>();
		fitness = new Vector<Vector<Double>>();
		GA_New = new Vector<Vector<Integer[]>>();
		fitness_New = new Vector<Vector<Double>>();
		for (int i = 0; i < p; i++) {
			GA.add(new Vector<Integer[]>());
			fitness.add(new Vector<Double>());
			GA_New.add(new Vector<Integer[]>());
			fitness_New.add(new Vector<Double>());
		}
		knap = K;
		generation = 0;
		avgFit = new double[p];
		minFit = new double[p];
		maxFit = new double[p];
		migrationRate = new double[p];
		new_FG = new double[p];
		old_FG = new double[p];
		randomGAPop();
		repairAll();
		eliteF = -1;
		getBestF();
		old_elites = elites;
		for (int i = 0; i < p; i++) {
			migrationRate[i] = 0.15;
			old_FG[i] = 0;
		}
	}
	
	public void randomGAPop(){
		for (int i = 0; i < p; i++) {
			for (int j = 0; j < n/p; j++) {
				Integer [] temp = new Integer[l];
				for (int k = 0; k < l; k++) {
					temp[k] = (Math.random()>0.5)? 1 : 0;
				}
				GA.get(i).add(temp);
				fitness.get(i).add(getFitness(temp));
			}
		}
	}
	public static double getFitness(Integer [] indiv){
		double sum = 0;
		for (int i = 0; i < l; i++) {
			if (indiv[i] == 1) {
				sum += knap.profits[i];
			}
		}		
		return sum;	
	}
	
	public double getWeight(Integer [] indiv) {
		double weight = 0;
		for (int i = 0; i < l; i++) {
			if (indiv[i] == 1) {
				weight += knap.weights[i];
			}
		}		
		return weight;
	}
	
	public void repairAll() {
		for (int k = 0; k < p; k++) {
			for (int n = 0; n < GA.get(k).size(); n++) {
				Integer[] indiv = GA.get(k).get(n);
				double weight = getWeight(indiv);
				if (weight > knap.capacity) {
					fitness.get(k).set(n, repair(indiv));
				}
			}
		}
		eliteF = -1;
		getBestF();
	}
	
	public double repair(Integer[] indiv) {
		double weight = 0;
		double sum = 0;
		for (int i = 0; i < l; i++) {
			if (indiv[knap.order[i]] == 1) {
				if (weight + knap.weights[knap.order[i]] <= knap.capacity) {
					weight += knap.weights[knap.order[i]];
					sum += knap.profits[knap.order[i]];
				} else {
					indiv[knap.order[i]] = 0;
				}
			}
		}
		return sum;
	}
	
	public void singleGARun(){
		generation++;
		crossover();
		mutation();
		setElite();
		replacement();
		if (generation % migrationInterval == 0) 
			//migration_ring_fixed();
			//migration_ring_adaptive();
			//migration_mix_fixed();
			migration_mix_adaptive();
	}
	public void runGA(){
		for (int i = 0; i < Main.maxIter; i++) {
			if(eliteF == optimal){
				Main.successRate++;
				break;
			}
			singleGARun();
		}
	}
	
	public void migration_ring_adaptive() {
		//System.out.println("========Migration at time " + generation + "============");
		//testing ring topology
		//p1->p2, p2->p3, ...., p5->p1
		int size = 0, dist1, dist2;
		//setting migration rates
		for (int k = 0; k < p; k++) {
			new_FG[k] = (elites[k][0] - old_elites[k][0]);
			if (new_FG[k] > old_FG[k] && migrationRate[k] <= 0.70)
				migrationRate[k] += 0.15;
			else if (new_FG[k] < old_FG[k] && migrationRate[k] >= 0.3)
				migrationRate[k] -= 0.15;
			old_FG[k] = new_FG[k];
		}
		old_elites = elites;
		//Adding emigrants to the GA_New population so they will be added during selection
		for (int k = 0; k < p; k++) {
			size = (int) (GA.get(k).size() * migrationRate[k]);
			if (k == 0) {
				dist1 = p - 1;
				dist2 = k + 1;
			} else if (k == p - 1) {
				dist1 = k - 1;
				dist2 = 0;
			} else {
				dist1 = k - 1;
				dist2 = k + 1;
			}
			for (int i = 0; i < size; i++) {
				GA_New.get(dist1).add(populationNew[k][order[k][i]]);
				fitness_New.get(dist1).add(fitnessNew[k][order[k][i]]);
				GA_New.get(dist2).add(populationNew[k][order[k][i]]);
				fitness_New.get(dist2).add(fitnessNew[k][order[k][i]]);
			}
		}
	}
	
	public void migration_ring_fixed() {
		//testing ring topology
		//p1->p2, p2->p3, ...., p5->p1
		int size = 0, dist1, dist2; 
		//Adding emigrants to the GA_New population so they will be added during selection
		for (int k = 0; k < p; k++) {
			size = (int) (GA.get(k).size() * migrationRate[k]);
			if (k == 0) {
				dist1 = p - 1;
				dist2 = k + 1;
			} else if (k == p - 1) {
				dist1 = k - 1;
				dist2 = 0;
			} else {
				dist1 = k - 1;
				dist2 = k + 1;
			}
			for (int i = 0; i < size; i++) {
				GA_New.get(dist1).add(populationNew[k][order[k][i]]);
				fitness_New.get(dist1).add(fitnessNew[k][order[k][i]]);
				GA_New.get(dist2).add(populationNew[k][order[k][i]]);
				fitness_New.get(dist2).add(fitnessNew[k][order[k][i]]);
			}
		}
	}
	
	public void migration_mix_adaptive() {
		//testing dynamic topology
		int size = 0;
		//setting migration rates
		for (int k = 0; k < p; k++) {
			new_FG[k] = (elites[k][0] - old_elites[k][0]);
			if (new_FG[k] > old_FG[k] && migrationRate[k] <= 0.7)
				migrationRate[k] += 0.15;
			else if (new_FG[k] < old_FG[k] && migrationRate[k] >= 0.3)
				migrationRate[k] -= 0.15;
			old_FG[k] = new_FG[k];
		}
		old_elites = elites;
		
		//now choosing the population with the top indiv to migrate from
		double max = optimal / 10;
		for (int i = 0; i < p; i++) {
			int top = i;
			double part[] = new double[p];
			
			for (int k = 0; k < p; k++) {
				if (k == top)
					continue; //skipping sender
				part[k] = (optimal - maxFit[k]) / max;
			}
			for (int k = 0; k < p; k++) {
				if (k == top)
					continue; //skipping sender
				size = (int) (GA.get(top).size() * migrationRate[top] * part[k]);
				if (size > GA.get(top).size())
					size = GA.get(top).size();
				for (int n = 0; n < size; n++) {
					GA_New.get(k).add(populationNew[top][order[top][n]]);
					fitness_New.get(k).add(fitnessNew[top][order[top][n]]);
				}
			}
		}

	}
	
	public void migration_mix_fixed() {
		int size = 0;
		
		//now choosing the population with the top indiv to migrate from
		for (int i = 0; i < p; i++) {
			int top = i;
			for (int k = 0; k < p; k++) {
				if (k == top)
					continue; //skipping sender
				size = (int) (GA.get(top).size() * .05);
				if (size > GA.get(top).size())
					size = GA.get(top).size();
				for (int n = 0; n < size; n++) {
					GA_New.get(k).add(populationNew[top][order[top][n]]);
					fitness_New.get(k).add(fitnessNew[top][order[top][n]]);
				}
			}
		}

	}
	
	
	public void mutation(){
		for (int k = 0; k < p; k++) {
			for (int i = 0; i < GA.get(k).size(); i++) {
				Integer[] indiv = GA.get(k).get(i);
				for (int j = 0; j < l; j++) {
					// Flip the bit if get mutated (using bitwise XOR with 1)
					if(Math.random() < mutation) {
						//populationGA[i][j] = populationGA[i][j]^1;
						indiv[j] = indiv[j]^1;
					}
				}
				double weight = getWeight(indiv);
				if (weight > knap.capacity) {
					repair(indiv);
				}
				GA.get(k).set(i, indiv);
				/*int Age = age.get(k).get(i);
				double LT = getLifetime(i, k);
				if (Age > LT) {
					GA.get(k).remove(i);
					fitness.get(k).remove(i);
					age.get(k).remove(i);
				} else
					i++; //else move on*/
					
			}
			for (int i = 0; i < GA_New.get(k).size(); i++) {
				Integer[] indiv = GA_New.get(k).get(i);
				for (int j = 0; j < l; j++) {
					// Flip the bit if get mutated (using bitwise XOR with 1)
					if(Math.random() < mutation) {
						//populationGA[i][j] = populationGA[i][j]^1;
						indiv[j] = indiv[j]^1;
					}
				}
				double weight = getWeight(indiv);
				if (weight > knap.capacity) {
					repair(indiv);
				}
				GA_New.get(k).set(i, indiv);
			}
		}
	}

	public void crossover(){
		Integer [] child1 = new Integer[l];
		Integer [] child2 = new Integer[l];
		for (int k = 0; k < p; k++) {
			int size = GA.get(k).size();
			for (int i = 0; i < size; i++) {
				if(Math.random() < cross){
					// Select a mate based on tournament selection and a random crossover point
					int randIndex = (int)Math.floor(Math.random()*size);
					Integer [] mate = GA.get(k).get(randIndex);
					Integer [] indiv = GA.get(k).get(i);
					Integer[] order = new Integer[(int)(size*.3)];
					Double[] fit = new Double[(int)(size*0.3)];
					Integer[][] genePool = new Integer[(int)(size*0.3)][l];
					int count = 0;
					while (count < genePool.length) {
						randIndex = (int)Math.floor(Math.random()*size);
						genePool[count]=GA.get(k).get(randIndex);
						fit[count] = fitness.get(k).get(randIndex);
						order[count] = count++;
					}
					//order the set
					sort(fit, order);
					//pick 2
					Integer [] mate1, mate2;
					mate1 = genePool[order[0]];
					mate2 = genePool[order[1]];
					/*int found = 0;
					for (int j = 0; j < genePool.length; j++)
						if(Math.random() < cross){
							if (found == 0) {
								mate1 = genePool[order[j]];
								found++;
							} else if (found == 1) {
								mate2 = genePool[order[j]];
								found++;
							} else
								break;
						}*/
					//if (Math.random() < 0.5)
						mate = mate1;
					//else
						//mate = mate2;
					int crossPoint = (int)Math.floor((Math.random()*(l-1)) + 1);
					
					// Perform crossover to generate offsprings
					for (int j = 0; j < crossPoint; j++) {
						child1[j] = indiv[j];
						child2[j] = mate[j];
					}
					for (int j = crossPoint; j < l; j++) {
						child1[j] = mate[j];
						child2[j] = indiv[j];
					}
					double val1 = getFitness(child1);
					double val2 = getFitness(child2);
					GA_New.get(k).add(child1);
					GA_New.get(k).add(child2);
					fitness_New.get(k).add(val1);
					fitness_New.get(k).add(val2);
					
					/*int bestF = repair(child1);
					int temp = repair(child2);
					if (temp > bestF)
						bestF = temp;
					if(Math.random()< cross){
						if(bestF == getFitness(child1))
							System.arraycopy(child1, 0, indiv, 0, l);
						else
							System.arraycopy(child2, 0, indiv, 0, l);
					}*/
					
					//Select one of children randomly and replace it with parent
					//if(Math.random()>0.5)
						//System.arraycopy(child1, 0, indiv, 0, l);
					//else
						//System.arraycopy(child2, 0, indiv, 0, l);
				}
			}
		}
	}
	
	/**
	 * Selecting the top to continue from current population and 
	 * new population (children and immigrants) 
	 */
	public void replacement(){
		populationNew = new Integer[p][][];
		fitnessNew = new Double[p][];
		order = new Integer[p][];
		for (int k = 0; k < p; k++) { 
			int n = GA.get(k).size();
			int n2 = GA_New.get(k).size();
			populationNew[k] = new Integer [n+n2][l];
			fitnessNew[k] = new Double[n+n2];
			order[k] = new Integer[n+n2];
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < l; j++)
					populationNew[k][i][j] = GA.get(k).get(i)[j];
				fitnessNew[k][i] = fitness.get(k).get(i);
				order[k][i] = i;
			}
			for (int i = 0; i < n2; i++) {
				for (int j = 0; j < l; j++)
					populationNew[k][n+i][j] = GA_New.get(k).get(i)[j];
				fitnessNew[k][n+i] = fitness_New.get(k).get(i);
				order[k][n+i] = n+i;
			}
			//sorting based on fitness
			sort(fitnessNew[k], order[k]);
			//System.out.println("==K="+k+"=G="+generation+"==");
			//for (int i = 0; i < 5; i++)
				//System.out.println("Order " + i + " = " + fitnessNew[k][order[k][i]]);
			//System.out.println("============");
			//now coping n individuals to continue
			GA.get(k).clear();
			fitness.get(k).clear();
			GA_New.get(k).clear();
			fitness_New.get(k).clear();
			avgFit[k] = 0; 
			for (int i = 0; i < n; i++) {
				int index = (int)Math.random()*i;
				GA.get(k).add(index, populationNew[k][order[k][i]]);
				fitness.get(k).add(index, fitnessNew[k][order[k][i]]);
				avgFit[k] += fitnessNew[k][order[k][i]];
			}
			
			//calculating avgFitness
			avgFit[k] /= n;
			minFit[k] = fitnessNew[k][order[k][n-1]];
			maxFit[k] = fitnessNew[k][order[k][0]];
			
			/*double [] expectedCount = new double [n];
			int [] actualCount = new int [n];
			//int [][] populationNew = new int [n][l];
			populationNew[k] = new Integer [n][l];
			double sumF = getSumFitness(fitness.get(k));
			int count = 0;
			
			//Set the expected counts of each string based upon fitness proportional selection
			for (int i = 0; i < n; i++){
				expectedCount[i] = n*fitness.get(k).get(i)/sumF;
				actualCount[i] = (int) expectedCount[i];
				count += actualCount[i];
			}
			
			do
			  {
			    for (int i=0; i<n; i++)
			    {
			      if (Math.random() > (expectedCount[i]-actualCount[i])) continue;
			      if (count == n) break;
			      actualCount[i]++;
			      count++;
			    }
			 } while (count < n);
			
			// Construct the gene pool
			count = 0;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < actualCount[i]; j++) {
					System.arraycopy(GA.get(k).get(i),0,populationNew[k][count],0,l);
					count++;
				}
			}
			
			// Reset the population back
			for (int i = 0; i < n; i++)
				System.arraycopy(populationNew[k][i],0,GA.get(k).get(i),0,l);*/
			
		}//for each sub population
		
	}
	
	/*
	 * 1. Update fitness values of each individual
	 * 2. Construct Elite
	 */
	public void setElite(){
		for (int k = 0; k < p; k++) {
			// Update fitness
			int size = fitness.get(k).size();
			for (int i = 0; i < size; i++)
				fitness.get(k).set(i, getFitness(GA.get(k).get(i)));
			size = fitness_New.get(k).size();
			for (int i = 0; i < size; i++)
				fitness_New.get(k).set(i, getFitness(GA_New.get(k).get(i)));
		}	
		// Construct Elite
		getBestF();
	}
	public double getBestF(){
		double[][] best = new double[p][2];
		for (int k = 0; k < p; k++) {
			best[k][0] = -1.0;
			for (int i = 0; i < fitness.get(k).size(); i++) {
				if(best[k][0] < fitness.get(k).get(i)) {
					best[k][0] = fitness.get(k).get(i);
					best[k][1] = getWeight(GA.get(k).get(i));
				}
			}
			for (int i = 0; i < fitness_New.get(k).size(); i++) {
				if(best[k][0] < fitness_New.get(k).get(i)) {
					best[k][0] = fitness_New.get(k).get(i);
					best[k][1] = getWeight(GA_New.get(k).get(i));
				}
			}
		}
		elites = best;
		for (int i = 0; i < p; i++) {
			if (elites[i][0] > eliteF) {
				eliteF = elites[i][0];
				eliteW = elites[i][1];
			}
		}
		return eliteF;
	}
	public double getSumFitness(Vector<Double> pop){
		double sum = 0;
		for (int i = 0; i < pop.size(); i++) {
			sum += pop.get(i); 
		}
		return sum;
	}
	
	/*
	 * Bubble sorting array values with indices in arg2
	 */
	public void sort(Double[] arg, Integer[] arg2) {
		boolean done = true;
		for (int i = 0; i < arg.length - 1; i++) {
			Double a = arg[arg2[i]];
			Double b = arg[arg2[i+1]];
			if (b > a) {
				int temp = arg2[i];
				arg2[i] = arg2[i+1];
				arg2[i+1] = temp;
				done = false;
			}
		}
		if (done)
			return;
		else
			sort(arg, arg2);
	}

}
