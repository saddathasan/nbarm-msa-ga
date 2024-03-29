package ga_MSA;

import java.util.Vector;

public class MPGA {
	//##############################################
	// RUNNING TEST ON GAP SCORES for 1000 GEN's
	// using 0.0005 mut and 0.25 cross
	// Open: Extend: Columns matching: Score
	//  -1:  0:   5 : +909
	//	-1: -0.01:6 : +864 
	//	-1: -0.05:3 : +781
	//	-1: -0.1: 7 : +747
	//  -1: -0.2: 6 : +850
	//	-1:	-0.3: 3 : +620
	//  -2:  0:   5 : +659
	//	-2: -0.01:6 : +580
	//  -2: -0.1: 5 : +591
	//  -2: -0.2: 4 : +709
	// 	-3:  0	: 5 : +544
	//	-3: -0.01:5 : +434
	//	-3: -0.1: 6	: +620
	//  -3: -0.15:3 : +427
	//	-3:	-0.2: 7	: +433
	//  -3: -0.25:6 : +391
	//	-3:	-0.3: 6 : +590
	//	-3: -0.4: 3 : +450
	//-2.9: -0.01:4 : +410
	//	-4: -0.01:6 : +474
	//	-4: -0.05:5 : +306
	//	-4: -0.1: 4 : +330
	//	-4: -0.2: 5 : +407
	//	-4: -0.3: 4 : +341
	//  -5:-0.01: 6 : +443
	//  -5: -0.1: 6 : +296
	//	-5: -0.2: 4 : +348
	private double mutation = 0.0005;
	private double cross = 0.40;
	private double[] migrationRate;
	private int migrationInterval = 5;
	private double[] new_FG; //new fitness gain
	private double[] old_FG;
	public int generation;
	public static int l, w;					
	public static int n = 500;				
	public static int p = 5;
	//public int [][] populationGA;
	public Vector<Vector<Integer[][]>> GA;
	public Vector<Vector<Double>> fitness;
	public Vector<Vector<Integer[][]>> GA_New;
	public Vector<Vector<Double>> fitness_New;
	public double eliteF;
	public int[] eliteI;
	public Integer[][] eliteIndiv;
	public double [] elites;
	private double [] old_elites;
	public int optimal = 600;	// Set the optimum here	
	//private static Knapsack knap; //knapsack problem
	private static MSA msa; //multi sequence alignment problem
	private Integer[][][][] populationNew;
	private Double[][] fitnessNew;
	private Integer[][] order;
	public double[] avgFit;
	public double[] minFit;
	public double[] maxFit;
	
	public MPGA(int L, int W, MSA M) {
		//populationGA = new int [n][l];
		l = L;
		w = W;
		GA = new Vector<Vector<Integer[][]>>();
		fitness = new Vector<Vector<Double>>();
		GA_New = new Vector<Vector<Integer[][]>>();
		fitness_New = new Vector<Vector<Double>>();
		for (int i = 0; i < p; i++) {
			GA.add(new Vector<Integer[][]>());
			fitness.add(new Vector<Double>());
			GA_New.add(new Vector<Integer[][]>());
			fitness_New.add(new Vector<Double>());
		}
		//knap = K;
		msa = M;
		generation = 0;
		avgFit = new double[p];
		minFit = new double[p];
		maxFit = new double[p];
		migrationRate = new double[p];
		new_FG = new double[p];
		old_FG = new double[p];
		randomGAPop();
		eliteF = -100000;
		eliteI = new int[3];
		eliteI[0] = eliteI[1] = eliteI[2] = 0;
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
				Integer [][] temp = new Integer[w][l];
				for (int k = 0; k < w; k++) {
					/*int limit = msa.length(k);
					//adding 1's first
					for (int n = 0; n < limit; n++)
						temp[k][n] = 1;
					//adding zero's
					for (int n = limit; n < l; n++)
						temp[k][n] = 0;
					//random shuffling
					for (int n = 0; n < l; n++) {
						int r = (int)Math.floor(Math.random()*l);
						int t = temp[k][r];
						temp[k][r] = temp[k][n];
						temp[k][n] = t;
					}*/
					for (int n = 0; n < l; n++) {
						temp[k][n] = 1;
					}
				}
				repair(temp, null);
				GA.get(i).add(temp);
				fitness.get(i).add(getFitness(temp));
			}
		}
	}
	
	//Using the sum of pairs for calculating
	public static double getFitness(Integer [][] indiv){
		double sumTotal = 0;
		double[] sumsOfPairs = new double[l];
		int[] size = new int[w];
		char[]	column = new char[indiv.length];
		int matches = 0;
		for (int n = 0; n < w; n++) {
			size[n] = -1;
		}
		for (int c = 0; c < l; c++) {
			boolean lastAllGaps = true;
			for (int n = 0; n < w; n++) {
				if (indiv[n][c] == 1) {
					size[n]++;
					column[n] = msa.getChar(n, size[n]);
				} else
					column[n] = '-';
				if (c > 0 && indiv[n][c-1] == 1)
					lastAllGaps = false;
			}
			sumsOfPairs[c] = 0;
			for (int s = 0; s < w - 1; s++) {
				char c1, c2;
				if (indiv[s][c] == 1)
					c1 = msa.getChar(s, size[s]);
				else
					c1 = '-';
				for (int t = s + 1; t < w; t++) {
					if (indiv[t][c] == 1)
						c2 = msa.getChar(t, size[t]);
					else
						c2 = '-';
					double dist = 0;
					if (c1 == '-' && c2 == '-')
						dist = 0;
					else if (c1 == '-') {
						if (c > 0 && indiv[s][c-1] == 0 && !lastAllGaps)
							dist = -0.2;
						else
							dist = -3;
					} else if (c2 == '-') {
						if (c > 0 && indiv[t][c-1] == 0 && !lastAllGaps)
							dist = -0.2;
						else
							dist = -3;
					} else if (c1 != '-' && c2 != '-')
						try {
							/*if (c1 == 'u' || c2 == 'u') {
								System.out.println(c1 + " " + c2);
								System.out.println(s + " " + t + " " + c);
							}*/
							dist = Blosum.getDistance(c1, c2);
						} catch (Exception e) {
							e.printStackTrace();
						}

				    //System.out.println(c1 + " : " + c2 + " = " + dist);
					sumsOfPairs[c] += dist; 
				}
			}

			boolean match = true;
			char e = column[0];
			for (int i = 1; i < indiv.length; i++) {
				if (column[i] != e)
					match = false;
			}
			if (match && e != '-') {
				matches++;
				//multiplying score
				sumsOfPairs[c] *= 2;
			}
		}
		//summing up the column score		
		for (int c = 0; c < l; c++)
			sumTotal += sumsOfPairs[c];
				
		return sumTotal;
	}
	
	public void singleGARun(){
		generation++;
		crossover();
		mutation();
		setElite();
		replacement();
		if (generation % migrationInterval == 0)
			migration_ring_fixed();
			//migration_ring_adaptive();
			//migration_mix_fixed();
			//migration_mix_adaptive();
			//replacement();
	}
	public void runGA(){
		for (int i = 0; i < Main.maxIter; i++) {
			if(eliteF >= optimal){
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
			new_FG[k] = (elites[k] - old_elites[k]);
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
			new_FG[k] = (elites[k] - old_elites[k]);
			if (new_FG[k] > old_FG[k] && migrationRate[k] <= 0.95)
				migrationRate[k] += 0.05;
			else if (new_FG[k] < old_FG[k] && migrationRate[k] >= 0.1)
				migrationRate[k] -= 0.05;
			old_FG[k] = new_FG[k];
		}
		old_elites = elites;
		
		//now choosing the population with the top indiv to migrate from
		for (int i = 0; i < p; i++) {
			int top = i;
			double part[] = new double[p];
			
			for (int k = 0; k < p; k++) {
				if (k == top)
					continue; //skipping sender
				//part[k] = Math.abs(850 - maxFit[k]) / max;
				//if (avgFit[top] > avgFit[k])
					part[k] = (1 + ((avgFit[top] - avgFit[k]) / avgFit[k]));
					if (part[k] < 0)
						part[k] = 0;
			}
			for (int k = 0; k < p; k++) {
				if (k == top)
					continue; //skipping sender
				size = (int) (GA.get(top).size() * migrationRate[top] * part[k]);
				//size = (int) (GA.get(top).size() * part[k]);
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
				for (int n = 0; n < size; n++) {
					GA_New.get(k).add(populationNew[top][order[top][n]]);
					fitness_New.get(k).add(fitnessNew[top][order[top][n]]);
				}
			}
		}

	}
	
	
	public void mutation(){
		int type = 0;
		if (Math.random() >= 2.0/3.0 ) {
			type = 1;
		} else if (Math.random() >= 1.0/3.0 ) {
			type = 2;
		}
		int max = 1;
		for (int k = 0; k < p; k++) {
			for (int i = 0; i < GA.get(k).size(); i++) {
				Integer[][] indiv = GA.get(k).get(i);
				int T = type;
				for (int m = 0; m < max; m++) {
					Integer[][] indivTemp = new Integer[w][l];
					for (int a = 0; a < w; a++) {
						indivTemp[a] = copy(indiv[a]);
					}
					for (int n = 0; n < w; n++) {
						for (int j = 0; j < l; j++) {
							//if (indiv[n][j] == 1 && type != 2)
								//continue;
							if(Math.random() < mutation) {
								//populationGA[i][j] = populationGA[i][j]^1;
								int r = (int)Math.floor(Math.random()*l);
								int t;
								// Swapping gap
								if (T == 1) {
									t = indivTemp[n][r];
									indivTemp[n][r] = indivTemp[n][j];
									indivTemp[n][j] = t;
									continue;
								} else if (T == 2) {
									indivTemp[n][j] = indivTemp[n][j]^1;
									continue;
								} else if (T == 3) {
									//adding a gap
								}
								//else, shifting sequence
								//if (indiv[n][j] == 1)
								//	continue;
								int start = j;
								int stop = j;
								//find 1st gap
								while (start > 0) {
									if (indiv[n][start-1] == 0)
										start--;
									else
										break;
								}
								//find last gap
								while (stop + 1 < l) {
									if (indiv[n][stop+1] == 0)
										stop++;
									else
										break;
								}
								//choose random new position
								int size = stop - start + 1;
								int rand = (int)Math.floor(Math.random()*(l));
								int[] temp = new int[size];
								for (int e = 0; e < size; e++) {
									temp[e] = indivTemp[n][start+e];
								}
								//copy 
								if (rand > start && rand <= stop && stop + rand - start < l) {
									//shifting elements after gaps
									for (int e = 0; e < rand - start; e++) {
										indivTemp[n][start+e] = indivTemp[n][stop + 1 + e];
									}
									for (int e = 0; e < size; e++) {
										indivTemp[n][rand + e] = temp[e];
									}
								} else if (rand > stop) {
									//shifting elements after gaps
									for (int e = 0; e < rand - stop; e++) {
										indivTemp[n][start+e] = indivTemp[n][stop + 1 + e];
									}
									for (int e = 0; e < size; e++) {
										indivTemp[n][rand - size + e] = temp[e];
									}
								} else if (rand < start) {
									for (int e = 0; e < start - rand ; e++) {
										indivTemp[n][stop - e] = indivTemp[n][start-1-e];
									}
									for (int e = 0; e < size; e++) {
										indivTemp[n][rand + e] = temp[e];
									}
								}
								
								/*int stop = 0;
								if (r > j) {
									stop = r - j - 1; 
									t = indivTemp[n][j];
									for (int e = 0; e < stop; e++) {
										indivTemp[n][j+e] = indivTemp[n][j+e+1];
									}
									indivTemp[n][r-1] = t;
								} else if (j > r) {
									stop = j - r - 1; 
									t = indivTemp[n][j];
									for (int e = stop + 1; e > 0; e--) {
										indivTemp[n][r+e] = indivTemp[n][r+e-1];
									}
									indivTemp[n][r] = t;
								}*/
							}//mut
						}//l
					}//w
					if (T == 2)
						repair(indivTemp, null);
					//if (getFitness(indivTemp) > fitness.get(k).get(i)) {
						indiv = indivTemp;
						GA.get(k).set(i, indiv);
					/*	break;
					}  else {
						//new type
						if (Math.random() >= 0.5)
							T = (T+1)%3;
						else
							T = (T+2)%3;
					}*/
				}//m times
			}
			
			for (int i = 0; i < GA_New.get(k).size(); i++) {
				Integer[][] indiv = GA_New.get(k).get(i);
				int T = type;
				for (int m = 0; m < max; m++) {
					Integer[][] indivTemp = new Integer[w][l];
					for (int a = 0; a < w; a++) {
						indivTemp[a] = copy(indiv[a]);
					}
					for (int n = 0; n < w; n++) {
						for (int j = 0; j < l; j++) {
							//if (indiv[n][j] == 1 && type != 2)
								//continue;
							if(Math.random() < mutation) {
								//populationGA[i][j] = populationGA[i][j]^1;
								int r = (int)Math.floor(Math.random()*l);
								int t;
								// Swapping gap
								if (T == 1) {
									t = indivTemp[n][r];
									indivTemp[n][r] = indivTemp[n][j];
									indivTemp[n][j] = t;
									continue;
								} else if (T == 2) {
									indivTemp[n][j] = indivTemp[n][j]^1;
									continue;
								}
								//else, shifting sequence
								//if (indiv[n][j] == 1)
								//	continue;
								int start = j;
								int stop = j;
								//find 1st gap
								while (start - 1 >= 0) {
									if (indiv[n][start-1] == 0)
										start--;
									else
										break;
								}
								//find last gap
								while (stop + 1 < l) {
									if (indiv[n][stop+1] == 0)
										stop++;
									else
										break;
								}
								//choose random new position
								int size = stop - start + 1;
								int rand = (int)Math.floor(Math.random()*(l));
								int[] temp = new int[size];
								for (int e = 0; e < size; e++) {
									temp[e] = indivTemp[n][start+e];
								}
								//copy 
								if (rand > start && rand <= stop && stop + rand - start < l) {
									//shifting elements after gaps
									for (int e = 0; e < rand - start; e++) {
										indivTemp[n][start+e] = indivTemp[n][stop + 1 + e];
									}
									for (int e = 0; e < size; e++) {
										indivTemp[n][rand + e] = temp[e];
									}
								} else if (rand > stop) {
									//shifting elements after gaps
									for (int e = 0; e < rand - stop; e++) {
										indivTemp[n][start+e] = indivTemp[n][stop + 1 + e];
									}
									for (int e = 0; e < size; e++) {
										indivTemp[n][rand - size + e] = temp[e];
									}
								} else if (rand < start) {
									for (int e = 0; e < start - rand ; e++) {
										indivTemp[n][stop - e] = indivTemp[n][start-1-e];
									}
									for (int e = 0; e < size; e++) {
										indivTemp[n][rand + e] = temp[e];
									}
								}
								//else, shifting sequence
								/*int stop = 0;
								if (r > j) {
									stop = r - j - 1; 
									t = indivTemp[n][j];
									for (int e = 0; e < stop; e++) {
										indivTemp[n][j+e] = indivTemp[n][j+e+1];
									}
									indivTemp[n][r-1] = t;
								} else if (j > r) {
									stop = j - r - 1; 
									t = indivTemp[n][j];
									for (int e = stop + 1; e > 0; e--) {
										indivTemp[n][r+e] = indivTemp[n][r+e-1];
									}
									indivTemp[n][r] = t;
								}*/
							}//mut
						}//l
					}//w
					if (T == 2)
						repair(indivTemp, null);
					//if (getFitness(indivTemp) > fitness_New.get(k).get(i)) {
						indiv = indivTemp;
						GA_New.get(k).set(i, indiv);
					/*	break;
					} else {
						//new type
						if (Math.random() >= 0.5)
							T = (T+1)%3;
						else
							T = (T+2)%3;
					}*/
				}//m times
			}//indiv
		}//pop
	}

	public void crossover(){
		Integer [][] child1;
		Integer [][] child2;
		int type = 0;
		/*if (Math.random() >= 2.0/3.0 ) {
			type = 1;
		} else if (Math.random() >= 1.0/3.0 ) {
			type = 2;
		}*/
		for (int k = 0; k < p; k++) {
			int size = GA.get(k).size();
			for (int i = 0; i < size; i++) {
				if(Math.random() < cross){
					// Select a mate based on tournament selection and a random crossover point
					int randIndex = (int)Math.floor(Math.random()*size);
					Integer [][] mate = GA.get(k).get(randIndex);
					Integer [][] indiv = GA.get(k).get(i);
					double ratio = 0.30;
					Integer[] order = new Integer[(int)(size*ratio)];
					Double[] fit = new Double[(int)(size*ratio)];
					Integer[][][] genePool = new Integer[(int)(size*ratio)][w][l];
					int count = 0;
					while (count < genePool.length) {
						randIndex = (int)Math.floor(Math.random()*size);
						genePool[count]= GA.get(k).get(randIndex);
						fit[count] = fitness.get(k).get(randIndex);
						order[count] = count++;
					}
					//order the set
					sort(fit, order);
					//pick 2
					Integer [][] mate1;//, mate2;
					mate1 = genePool[order[0]];
					//mate2 = genePool[order[1]];
					mate = mate1;

					child1 = new Integer[w][l];
					child2 = new Integer[w][l];
					// Perform random row crossover to generate offsprings
					if (type == 1) {
						int crossPoint = (int)Math.floor(Math.random()*(w));
						for (int j = 0; j <= crossPoint; j++) {
							child1[j] = copy(indiv[j]);
							child2[j] = copy(mate[j]);
						}
						for (int j = crossPoint + 1; j < w; j++) {
							child1[j] = copy(mate[j]);
							child2[j] = copy(indiv[j]);
						}
					} else if (type == 2) {
						int[] cp = new int[w];
						for (int n = 0; n < w; n++) {
							int crossPoint = (int)Math.floor(Math.random()*(l));
							cp[n] = crossPoint;
							for (int j = 0; j < crossPoint; j++) {
								child1[n][j] = indiv[n][j];
								child2[n][j] = mate[n][j];
							}
							for (int j = crossPoint; j < l; j++) {
								child1[n][j] = mate[n][j];
								child2[n][j] = indiv[n][j];
							}
						}
						repair(child1, cp);
						repair(child2, cp);
					} else {
						//else do random column crossover
						int[] crossPoints = new int[w];
						int[] gaps = new int[w];
						int[] gapsSoFar = new int[w];
						int[] aminoSoFar = new int[w];
						int[] gapsSoFarMate = new int[w];
						int[] aminoSoFarMate = new int[w];
						for (int j = 0; j < w; j++) {
							crossPoints[j] = (int)Math.floor((Math.random()*(l-1)) + 1);
							gaps[j] = l - msa.length(j);
							gapsSoFar[j] = 0;
							aminoSoFar[j] = 0;
							for (int c = 0; c <= crossPoints[j]; c++) {
								if (indiv[j][c] == 1)
									aminoSoFar[j]++;
								else
									gapsSoFar[j]++;
								child1[j][c] = indiv[j][c];
							}
							gapsSoFarMate[j] = 0;
							aminoSoFarMate[j] = 0;
							int cross = 0;
							while (aminoSoFarMate[j] < aminoSoFar[j]) {
								if (mate[j][cross] == 1)
									aminoSoFarMate[j]++;
								else
									gapsSoFarMate[j]++;
								child2[j][cross] = mate[j][cross];
								cross++;
							}
							//System.out.println("cross+1" +cross);
							//start coping 2nd part
							int gapsExtra = Math.abs(gapsSoFarMate[j] - gapsSoFar[j]);
							int gapsNeeded = gapsExtra;
							if (gapsSoFar[j] > gapsSoFarMate[j]) {
								for (int c = cross, d = crossPoints[j]+1; c < l; c++) {
									if (gapsExtra > 0 && mate[j][c] == 0) {
										gapsExtra--;
										continue;
									}
									child1[j][d++] = mate[j][c];
								}
								for (int e = 0; e < gapsNeeded; e++) {
									child2[j][cross++] = 0;
								}
								for (int c = cross, d = crossPoints[j]+1; c < l; c++) {
									child2[j][c] = indiv[j][d++];
								}
							} else if (gapsSoFarMate[j] > gapsSoFar[j]) {
								for (int c = cross, d = crossPoints[j]+1; d < l; d++) {
									if (gapsExtra > 0 && indiv[j][d] == 0) {
										gapsExtra--;
										continue;
									}
									child2[j][c++] = indiv[j][d];
								}
								for (int e = 0 ; e < gapsNeeded; e++) {
									child1[j][++crossPoints[j]] = 0;
								}
								for (int c = cross, d = crossPoints[j]+1; c < l; c++) {
									child1[j][d++] = mate[j][c];
								}
							} else {
								for (int c = cross, d = crossPoints[j]+1; c < l; c++) {
									child1[j][d++] = mate[j][c];
								}
								for (int c = cross, d = crossPoints[j]+1; c < l; c++) {
									child2[j][c] = indiv[j][d++];
								}							
							}
						}
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
		populationNew = new Integer[p][][][];
		fitnessNew = new Double[p][];
		order = new Integer[p][];
		for (int k = 0; k < p; k++) { 
			int n = GA.get(k).size();
			int n2 = GA_New.get(k).size();
			populationNew[k] = new Integer [n+n2][w][l];
			fitnessNew[k] = new Double[n+n2];
			order[k] = new Integer[n+n2];
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < w; j++)
					populationNew[k][i][j] = copy(GA.get(k).get(i)[j]);
				fitnessNew[k][i] = fitness.get(k).get(i);
				order[k][i] = i;
			}
			for (int i = 0; i < n2; i++) {
				for (int j = 0; j < w; j++)
					populationNew[k][n+i][j] = copy(GA_New.get(k).get(i)[j]);
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
		double[] best = new double[p];
		double bestAll = -100000;
		for (int k = 0; k < p; k++) {
			best[k] = -100000;
			for (int i = 0; i < fitness.get(k).size(); i++) {
				if(best[k] < fitness.get(k).get(i)) {
					best[k] = fitness.get(k).get(i);
					if(eliteF < fitness.get(k).get(i)) {
						eliteI[0] = k; 
						eliteI[1] = i;
						eliteI[2] = 0;
						eliteIndiv = new Integer[w][l];
						for (int j = 0; j < w; j++) {
							eliteIndiv[j] = copy(GA.get(eliteI[0]).get(eliteI[1])[j]);
						}
					}
				}
			}
			for (int i = 0; i < fitness_New.get(k).size(); i++) {
				if(best[k] < fitness_New.get(k).get(i)) {
					best[k] = fitness_New.get(k).get(i);
					if(eliteF < fitness_New.get(k).get(i)) {
						eliteI[0] = k; 
						eliteI[1] = i;
						eliteI[2] = 1;
						eliteIndiv = new Integer[w][l];
						for (int j = 0; j < w; j++) {
							eliteIndiv[j] = copy(GA_New.get(eliteI[0]).get(eliteI[1])[j]);
						}
					}
				}
			}
		}
		elites = best;
		for (int i = 0; i < p; i++) {
			if (elites[i] > eliteF) {
				eliteF = elites[i];
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
	
	//returns a copy of an integer array
	public Integer[] copy(Integer[] arg1) {
		Integer[] toReturn = new Integer[arg1.length];
		for (int i = 0; i < arg1.length; i++) {
			toReturn[i] = arg1[i];
		}
		return toReturn;
	}
	
	private void repair(Integer[][] arg1, int[] arg2) {

		int[] count = new int[w];
		for (int n = 0; n < w; n++) {
			count[n] = 0;
			for (int j = 0; j < l; j++) {
				if (arg1[n][j] == 1)
					count[n]++;
			}
		}
		//check if fix is needed
		for (int n = 0; n < w; n++) {
			//System.out.println("count c = " + count[n] + ", length c = " +  msa.length(n));
			if (count[n] < msa.length(n)) {
				//need to add ones
				while (count[n] < msa.length(n)){
					int start = (arg2==null)? 0: arg2[n];
					for (int j = start; j < l; j++) {
						if (Math.random() < 1.0 / l) {
							if (arg1[n][j] == 0) {
								arg1[n][j] = 1;
								count[n]++;
								if (count[n] >= msa.length(n))
									break;
							}
						}
					}
				}
			} else if (count[n] > msa.length(n)) {
				//need to remove ones
				while (count[n] > msa.length(n)){
					int start = (arg2==null)? 0: arg2[n];
					for (int j = start; j < l; j++) {
						if (Math.random() < 1.0 / l) {
							if (arg1[n][j] == 1) {
								arg1[n][j] = 0;
								count[n]--;
								if (count[n] <= msa.length(n))
									break;
							}
						}
					}
				}
			} 
		}
	}//repair

}
