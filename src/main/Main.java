package main;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.util.Precision;

import java.math.BigDecimal;

public class Main {
    private static MersenneTwister r = new MersenneTwister();
    public static void main(String[] args) {
        Space s = new Space(10);
        long time1 = System.nanoTime();
    	s.readDB();
    	s.readCFG();
    	if (s.useInput) {
    		//Read molecules from Input.xyz, do not propagate
			s.readInput();
		}
    	else {
			//If molecules can't be placed in space of given size within 20 tries, increase size by 10% and retry
			while (!s.propagate()) {
				s.size = s.size * 1.1;
			}
		}
        s.makeDirectory();
        long time2 = System.nanoTime();
        String stamp = timestamp(time1, time2);
        String str = s.write("Output0.xyz");
        String initText = "Initialization ";
        if (!s.useInput){
        	initText += "and propagation ";
		}
        s.log(initText + "done in " + stamp);
        s.log("Cluster movement constrained within cube with side length " + s.size * 1.5);
        double startingEnergy = s.calcEnergy();
		s.log("Starting Energy: " + startingEnergy);
        time1 = System.nanoTime();
        if (s.staticTemp){
        	s.numTeeth = 1;
        	s.writeEnergy(startingEnergy);
		}
        sawtoothAnneal(s, s.maxTemperature, s.movePerPoint, s.pointsPerTooth, s.pointIncrement, s.numTeeth, s.tempDecreasePerTooth, s.maxTransDist, s.magwalkFactorTrans, s.magwalkProbTrans, s.maxRotDegree, s.magwalkProbRot, startingEnergy);
        time2 = System.nanoTime();
        stamp = timestamp(time1, time2);
		s.log("Annealing done in " + stamp + ".");
    }
    
    public static void sawtoothAnneal(Space s, double maxTemp, double numMovesPerPoint, int ptsPerTooth, int ptsIncrement, int numTeeth, double toothScale, double maxD, double magwalkFactorTrans, double magwalkProbTrans, double maxRot, double magwalkProbRot, double startingEnergy) {
    	double t = maxTemp;
    	double saveT = t;
    	//Boolean used to check if final cycle
    	boolean isDoubleCycle = false;
    	for (int x = 0; x < numTeeth; x++) {
    		double delT = t / (ptsPerTooth - 1);
    		for (int y = 0; y < ptsPerTooth; y++) {
    			for (int z = 0; z < numMovesPerPoint; z++) {
    				double energyChange = 0;
    				Molecule m = s.randMolecule(); //Pick a random molecule
    				if (r.nextDouble() >= 0.5 && m.atoms.size() > 1) { //If the rotation of m matters (more than 1 atom), 50% for either rotate or translate; otherwise just translate
                        if (r.nextDouble() >= magwalkProbRot) { //Chance to magwalk from config
							energyChange = s.rotate(m, maxRot, t);
                        }
                        else {
							energyChange = s.rotate(m, 2 * Math.PI, t); //Magwalking sets rotation maximum to 2PI
                        }
                    }
                    else{
                        if (r.nextDouble() >= magwalkProbTrans) { //Chance to magwalk from config
							energyChange = s.move(m, maxD, t);
                        }
                        else{
                            energyChange = s.move(m, maxD * magwalkFactorTrans, t); //Magwalking multiplies distance maximum by magwalk factor (specified in config file)
                        }
                    }
					startingEnergy += energyChange;
        			s.writeEnergy(startingEnergy);
    			}
    			if (!s.staticTemp){
					t -= delT; //Decrease temperature by decrement factor
				}
    			if (t < 0) {
    				t = 0; //Prevents temperature from becoming negative, which causes issues
    			}

    			/**double percentAccepted = c / numMovesPerPoint;
				if (percentAccepted < 0.3) {
    				maxDTemp *= 0.8;
    				maxRotTemp *= 0.8;
    				//System.out.println(percentAccepted);
    			}
    			else if (percentAccepted > 0.7) {
    				maxDTemp *= 1.2;
    				maxRotTemp *= 1.2;
    				if (maxRotTemp > Math.PI * 2) {
    					maxRotTemp = Math.PI * 2;
    				}
    				//System.out.println(percentAccepted);
    			}**/
    			//c = 0;

    			s.writeMovie("Output" + x + "_" + (x + 1) + "_Movie.xyz");
    		}
    		saveT *= toothScale;
    		t = saveT;
    		ptsPerTooth += ptsIncrement;
			if (x == numTeeth - 1 && !isDoubleCycle && s.extraCycle){
				t = 0;
				isDoubleCycle = true;
				x--;
			}
			else{
				s.write("Output" + (x + 1) + ".xyz");
				s.log("Energy at end of tooth " + (x + 1) + ": " + s.calcEnergy());
			}
    	}
    }
    public static String timestamp(long time1, long time2){
    	String ret = "";
		float runtime = (time2 - time1) / 1000000000f;
		double seconds = Precision.round(runtime % 60.0, 3, BigDecimal.ROUND_HALF_UP);
		int minutes = (int) runtime / 60;
		int hours = (int) minutes / 60;
		if (hours > 0) ret += hours + "h ";
		if (minutes > 0) ret += minutes + "m ";
		return ret + seconds + "s";
	}
}
