package main;
import org.apache.commons.math3.util.Pair;
import org.apache.commons.math3.util.Precision;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.net.URLDecoder;
import java.util.HashMap;
import java.util.Map;

public class Main {
    public static void main(String[] args) {
    	long procStart = System.nanoTime();
		Space s = new Space();
    	try {
			long time1 = System.nanoTime();
			Map<String, String> parsedArgs = getArgs(args);
			if (parsedArgs.containsKey("seed")) {
				try {
					s.setSeed(Long.parseLong(parsedArgs.get("seed")));
				} catch (NumberFormatException e) {
					throw new RuntimeException("Error: --set-seed argument must be a valid long value.");
				}
			}
			s.makeDirectoryName(parsedArgs);
			s.makeDirectory(parsedArgs);
			System.out.println("Writing output to: " + s.getDir());
			s.readDB(parsedArgs.get("dbase"));
			Config.parseConfig(parsedArgs.get("config"), s);
//			s.readCFG(parsedArgs.get("config"));
			if (!Config.useInput && parsedArgs.get("inputIncluded").equals("yes")) {
				throw new RuntimeException("Error: You cannot provide a parameter for --input if \"Use Input.xyz\" is not selected in your config.");
			}
			if (!Config.chooseParams && parsedArgs.get("paramsIncluded").equals("yes")) {
				throw new RuntimeException("Error: You cannot provide a parameter for --params if \"Choose All Interaction Parameters\" is not selected in your config.");
			}
			if (Config.chooseParams) s.readParams(parsedArgs.get("interactionParams"));
			else s.setupPairVals();
			s.writeParams();
			if (Config.useInput) {
				//Read molecules from Input.xyz, do not propagate
				s.readInput(parsedArgs.get("input"));
			} else {
				//If molecules can't be placed in space of given size within 20 tries, increase size by 10% and retry
				while (!s.propagate()) {
					Config.spaceLen *= 1.1;
				}
			}
			long time2 = System.nanoTime();
			String stamp = timestamp(time1, time2);
			s.write(0);
			String initText = "Initialization ";
			if (!Config.useInput) {
				initText += "and propagation ";
			}
			s.log(initText + "done in " + stamp);
			s.log("Cluster movement constrained within cube with side length " + Config.spaceLen * 1.5);
			double startingEnergy = s.calcEnergy();
			s.log("Starting Energy: " + startingEnergy);
			time1 = System.nanoTime();
			if (Config.staticTemp) {
				Config.numTeeth = 1;
				Config.ptsPerTooth = 1;
				if (Config.writeEnergiesEnabled) {
					s.writeEnergy(startingEnergy);
				}
			}
			sawtoothAnneal(s, startingEnergy);
			time2 = System.nanoTime();
			stamp = timestamp(time1, time2);
			s.log("Annealing done in " + stamp + ".");
		} catch (Exception exc) {
    		System.err.println(exc.getMessage());
		} finally {
    		s.writeExecTime(procStart);
		}
    }

    public static Map<String, String> getArgs(String[] args) throws IOException {
		//Construct path to file, rather complicated but has to work with .jar or project files for testing
		String pathDir = Main.class.getProtectionDomain().getCodeSource().getLocation().getPath();
		pathDir = URLDecoder.decode(pathDir, "utf-8");
		pathDir = "/" + pathDir.substring(1, pathDir.lastIndexOf("/"));

		Map<String, String> parsed = new HashMap<>();

		parsed.put("output", pathDir);
		parsed.put("config", pathDir + "/config.txt");
		parsed.put("dbase", pathDir + "/dbase.txt");
		parsed.put("input", pathDir + "/Input.xyz");
		parsed.put("inputIncluded", "no");
		parsed.put("interactionParams", pathDir + "/interaction_params.txt");
		parsed.put("paramsIncluded", "no");

		for (int i = 0, argsLength = args.length; i < argsLength; i++) {
			String arg = args[i];
			String argName = null;
			String fileType = null;
			switch (arg) {
				case "-i":
				case "--input":
					argName = "input";
					fileType = ".xyz";
					parsed.put("inputIncluded", "yes");
				case "-p":
				case "--params":
					if (argName == null) {
						argName = "interactionParams";
						parsed.put("paramsIncluded", "yes");
					}
					if (fileType == null) fileType = ".txt";
				case "-d":
				case "--dbase":
					if (argName == null) argName = "dbase";
					if (fileType == null) fileType = ".txt";
				case "-c":
				case "--config":
					if (argName == null) argName = "config";
					if (fileType == null) fileType = ".txt";
					if (i == argsLength - 1) throw new RuntimeException(String.format("Error: Expected file path for \"%s\" parameter", argName));
					i++;
					String value = args[i];
					File file = new File(value);
					if (!file.exists()) throw new RuntimeException(String.format("Error: File not found: %s", file.getCanonicalPath()));
					if (!file.canRead()) throw new RuntimeException(String.format("Error: Cannot read file %s", file.getCanonicalPath()));
					if (!file.getPath().endsWith(fileType)) throw new RuntimeException(String.format("Error: Bad filetype for \"%s\" parameter; expected %s file", argName, fileType));
					parsed.put(argName, value);
					break;
				case "-o":
				case "--output":
					if (i == argsLength - 1) throw new RuntimeException("Error: Expected output directory path for \"output\" parameter");
					i++;
					value = args[i];
					file = new File(value);
					if (!file.exists()) throw new RuntimeException(String.format("Error: Directory not found: %s", file.getCanonicalPath()));
					if (!file.isDirectory()) throw new RuntimeException(String.format("Error: %s is not a directory", file.getCanonicalPath()));
					if (!file.canWrite()) throw new RuntimeException(String.format("Error: Cannot write to directory %s", file.getCanonicalPath()));
					parsed.put("output", value);
					break;
				case "-s":
				case "--set-seed":
					i++;
					parsed.put("seed", args[i]);
					break;
				default:
					if (arg.startsWith("-")) {
						throw new RuntimeException(String.format("Error: Unknown flag: %s", arg));
					} else {
						throw new RuntimeException(String.format("Error: Unexpected argument: %s", arg));
					}
			}
		}

		// test existence of required files
		File cFile = new File(parsed.get("config"));
		if (!cFile.exists()) throw new RuntimeException(String.format("Error: File not found: %s", cFile.getCanonicalPath()));
		File dFile = new File(parsed.get("dbase"));
		if (!dFile.exists()) throw new RuntimeException(String.format("Error: File not found: %s", dFile.getCanonicalPath()));

		return parsed;
	}

    public static void sawtoothAnneal(Space s, double startingEnergy) {
    	double t = Config.maxTemp;
    	double saveT = t;
    	//Boolean used to check if final cycle
    	boolean isDoubleCycle = false;
		int ptsPerTooth = Config.ptsPerTooth;
		double magwalkProbRot = Config.magwalkProbRot;
		double magwalkProbTrans = Config.magwalkProbTrans;
		double maxD = Config.maxTrans;
		double maxRot = Config.maxRot;
    	for (int x = 0; x < Config.numTeeth; x++) {
    		double delT = t / (ptsPerTooth - 1);
    		int accepted = 0;
    		int total = 0;
			for (int y = 0; y < ptsPerTooth; y++) {
				//Energy counters used for average energy
				double totalEnergy = 0;
				double totalSquaredEnergy = 0;
    			for (int z = 0; z < Config.movePerPt; z++) {
    				total++;
    				Molecule m = s.randMolecule(); //Pick a random molecule
					Pair<Double, Integer> values;
    				if (Space.getR().nextDouble() >= 0.5 && m.atoms.size() > 1) { //If the rotation of m matters (more than 1 atom), 50% for either rotate or translate; otherwise just translate
                        if (Space.getR().nextDouble() >= magwalkProbRot) { //Chance to magwalk from config
							values = s.rotate(m, maxRot, t);
                        }
                        else {
							values = s.rotate(m, 2 * Math.PI, t); //Magwalking sets rotation maximum to 2PI
                        }
                    }
                    else{
                        if (Space.getR().nextDouble() >= magwalkProbTrans) { //Chance to magwalk from config
							values = s.move(m, maxD, t);
                        }
                        else{
							values = s.move(m, maxD * Config.magwalkFactorTrans, t); //Magwalking multiplies distance maximum by magwalk factor (specified in config file)
                        }
                    }
					startingEnergy += values.getFirst();
					if (!Config.staticTemp || z > Config.eqConfigs) {
						totalEnergy += startingEnergy;
						totalSquaredEnergy += Math.pow(startingEnergy, 2);
					}
                    if (Config.writeEnergiesEnabled) {
						s.writeEnergy(startingEnergy);
					}
                    accepted += values.getSecond();
    			}
    			s.writeAcceptance(t, accepted, total);
				if (Config.writeConfHeatCapacitiesEnabled) {
					double roundedTemperature = new BigDecimal(t).setScale(2, RoundingMode.HALF_UP).doubleValue();
					if (roundedTemperature == 0) continue;
					double avgEnergy = totalEnergy / Config.movePerPt;
					double avgSquaredEnergy = totalSquaredEnergy / Config.movePerPt;
					double configurationalHeatCapacity = (avgSquaredEnergy - Math.pow(avgEnergy, 2)) / (Space.BOLTZMANN_CONSTANT * Math.pow(t, 2));
					s.writeConfigurationalHeatCapacity(roundedTemperature, avgEnergy, configurationalHeatCapacity);
				}
    			if (!Config.staticTemp){
					t -= delT; //Decrease temperature by decrement factor
				}
    			if (t < 0) {
    				t = 0; //Prevents temperature from becoming negative, which causes issues
    			}
    			s.writeMovie(x, x+1);
    		}
    		saveT *= Config.tempDecreasePerTooth;
    		t = saveT;
    		ptsPerTooth += Config.ptsAddedPerTooth;
			if (x == Config.numTeeth - 1 && !isDoubleCycle && Config.finale){
				t = 0;
				isDoubleCycle = true;
				x--;
				maxD /= 100;
				maxRot /= 100;
				magwalkProbRot = 0;
				magwalkProbTrans = 0;
			}
			else{
				s.write(x + 1);
				s.log("Energy at end of tooth " + (x + 1) + ": " + s.calcEnergy());
			}
    	}
    	//After all teeth, if numTeeth > 1, write minimum energy
		if (Config.numTeeth > 1){
			s.writeMinEnergy();
		}
    }
    public static String timestamp(long time1, long time2){
    	String ret = "";
		float runtime = (time2 - time1) / 1000000000f;
		double seconds = Precision.round(runtime % 60.0, 3, BigDecimal.ROUND_HALF_UP);
		int minutes = (int) (runtime / 60);
		int hours = minutes / 60;
		if (hours > 0){
			ret += hours + "h ";
			minutes %= 60;
		}
		if (minutes > 0) ret += minutes + "m ";
		return ret + seconds + "s";
	}
}
