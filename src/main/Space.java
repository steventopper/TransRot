package main;
import java.io.*;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.net.URLDecoder;
import java.nio.file.Files;
import java.util.*;

import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.util.Pair;
import org.apache.commons.math3.util.Precision;
import java.time.LocalDateTime;

public class Space {
    //All set in config.txt
    double size;
    private int maxPropFailures;
    double maxTemperature;
    int movePerPoint;
    int pointsPerTooth;
    int pointIncrement;
    int numTeeth;
    double tempDecreasePerTooth;
    double maxTransDist;
    double maxRotDegree;
    double magwalkFactorTrans;
    double magwalkProbTrans;
    double magwalkProbRot;

    boolean useInput;
    boolean extraCycle;
    boolean staticTemp;
    boolean writeEnergiesEnabled;
    boolean writeConfHeatCapacitiesEnabled;
    int eqConfigs;
    boolean writeAcceptanceRatios;

    private String dir; //Directory of .xyz files to be saved in, created at runtime

    //For saving lowest energy structure separately
    private double minEnergy = Double.MAX_VALUE;
    private int minEnergyToothNum;

    private final ArrayList<Molecule> list;
    private final ArrayList<Molecule> space;
    private final ArrayList<Molecule> moveableMolecules;
    private final ArrayList<Molecule> dbase;
    private final HashMap<Pair<UUID, UUID>, double[]> pairwiseDbase;
    private static final MersenneTwister r = new MersenneTwister(); //Used instead of Java.Random for greater accuracy in randomization
    static final double BOLTZMANN_CONSTANT = 0.0019872;
    public Space(double s){
        size = s;
        list = new ArrayList<>();
        space = new ArrayList<>();
        moveableMolecules = new ArrayList<>();
        dbase = new ArrayList<>();
        pairwiseDbase = new HashMap<>();
    }
    //Adds molecule to space from config file string
    public boolean add(String s){ //Adds larger molecules to be placed first, to lower chance of overlap later
        int c = 0;
        for (Molecule a : dbase){
            if (a.name.equals(s)){
                Molecule m = new Molecule(a);
                while (c < list.size()){
                    if (m.radius > list.get(c).radius){ //If radius of m larger than radius of checked molecule, place m before checked molecule
                        list.add(c, new Molecule(m));
                        return true;
                    }
                    c++;
                }
                if (c == list.size()){ //If m smaller than all other molecules, place m last in list
                    list.add(list.size(), m);
                    return true;
                }
            }
        }
        return false;
    }
    //Places molecules in list into space with random positions and rotations
    public boolean propagate(){
        for (Molecule m : list){
            int c = 0;
            while (true) {
                //Determine random position
                double x = (r.nextDouble() * size) - (size / 2);
                double y = (r.nextDouble() * size) - (size / 2);
                double z = (r.nextDouble() * size) - (size / 2);
                if (space.isEmpty()){ //If placing first molecule, no reason to do comparisons, just place the molecule
                    m.put(x, y, z);
                    space.add(m);
                    moveableMolecules.add(m);
                    break;
                }
                boolean ok = true;
                for (Molecule n : space) {
                    //Check for overlap
                    double min = m.radius + n.radius;
                    double d = Math.sqrt(Math.pow(n.x - x, 2) + Math.pow(n.y - y, 2) + Math.pow(n.z - z, 2)); //Distance between molecules
                    if (d < min){ //If distance lower than sum of molecule radiuses, invalid placement
                        c++;
                        ok = false;
                        break;
                    }
                }
                if (ok) {
                    m.put(x, y, z);
                    space.add(m);
                    moveableMolecules.add(m);
                    break;
                }
                if (c >= maxPropFailures){ //Need larger space to place molecules
                	space.clear();
                	moveableMolecules.clear();
                    return false;
                }
            }
        }
        for (Molecule m : space){
            //Rotate each molecule randomly
            m.rotateTemp(2 * Math.PI);
            m.setTemps();
        }
        list.clear(); //No reason for this to take up space anymore
        return true;
    }
    //Reads molecules from dbase.txt
    public void readDB(){
    	try {
    	    //Construct path to file, rather complicated but has to work with .jar or project files for testing
            String pathName = Main.class.getProtectionDomain().getCodeSource().getLocation().getPath();
            pathName = URLDecoder.decode(pathName, "utf-8");
            pathName = "/" + pathName.substring(1,pathName.lastIndexOf("/")) + "/dbase.txt";
            Scanner scanner = new Scanner(new File(pathName));
            int currLine = 0; //Used to print line during error detection
            while (scanner.hasNextLine()){
                //Parse file based on expected file format; if an issue is encountered, assume the file is improperly set up and throw an error
                try {
                    String[] words = scanner.nextLine().split(" " + " +");
                    currLine++;
                    if (words.length == 1) {
                        int n;
                        try {
                            n = Integer.parseInt(words[0]);
                        } catch (Exception exc) {
                            throw new Exception();
                        }
                        ArrayList<Atom> atoms = new ArrayList<>();
                        if (!scanner.hasNextLine()) {
                            throw new Exception();
                        }
                        String[] s = scanner.nextLine().split(" " + " +");
                        currLine++;
                        if (s.length != 2) {
                            throw new Exception();
                        }
                        double radius = Double.parseDouble(s[1]);
                        String name = s[0];
                        int numGhosts = 0;
                        for (int x = 0; x < n; x++) {
                            String[] atom = scanner.nextLine().split(" " + " +");
                            currLine++;
                            if (atom.length != 10) {
                                throw new Exception();
                            }
                            Atom a = new Atom(atom[0], Double.parseDouble(atom[1]), Double.parseDouble(atom[2]), Double.parseDouble(atom[3]), Double.parseDouble(atom[4]), Double.parseDouble(atom[5]), Double.parseDouble(atom[6]), Double.parseDouble(atom[7]), Double.parseDouble(atom[8]), Double.parseDouble(atom[9]));
                            atoms.add(a);
                            if (a.symbol.contains("*")){
                                numGhosts++;
                            }
                        }
                        if (numGhosts == n){
                            System.out.println("Error on line" + currLine + " in dbase.txt: Molecule " + name + " cannot be comprised of only ghost atoms.");
                            System.exit(0);
                        }
                        Molecule m = new Molecule(name, radius, atoms);
                        dbase.add(m);
                    }
                }
                catch (Exception exc){
                    System.out.println("Error on line " + currLine + " in dbase.txt: File incorrectly formatted.");
                    scanner.close();
                    System.exit(0);
                }
            }
            scanner.close();
    	}
    	catch (Exception exc) {
    		System.out.println("Error: File " + dir + "/dbase.txt not found.");
            System.exit(0);
    	}
    }
    //Reads annealing config from config.txt
    public void readCFG(){
        try {
            //Construct path to file, rather complicated but has to work with .jar or project files for testing
            String pathName = Main.class.getProtectionDomain().getCodeSource().getLocation().getPath();
            pathName = URLDecoder.decode(pathName, "utf-8");
            pathName = "/" + pathName.substring(1, pathName.lastIndexOf("/")) + "/config.txt";
            Scanner scanner = new Scanner(new File(pathName));
            int currLine = 3; //Tracks current line for error printing
            //Skip 2 lines for comments
            scanner.nextLine();
            scanner.nextLine();
            try {
                //Read configs; splits when 2 or more spaces are found
                maxTemperature = Double.parseDouble(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                movePerPoint = Integer.parseInt(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                pointsPerTooth = Integer.parseInt(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                pointIncrement = Integer.parseInt(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                numTeeth = Integer.parseInt(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                tempDecreasePerTooth = Double.parseDouble(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                maxTransDist = Double.parseDouble(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                magwalkFactorTrans = Double.parseDouble(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                magwalkProbTrans = Double.parseDouble(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                maxRotDegree = Double.parseDouble(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                magwalkProbRot = Double.parseDouble(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                size = Double.parseDouble(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                maxPropFailures = Integer.parseInt(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                useInput = Boolean.parseBoolean(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                extraCycle = Boolean.parseBoolean(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                staticTemp = Boolean.parseBoolean(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                writeEnergiesEnabled = Boolean.parseBoolean(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                writeConfHeatCapacitiesEnabled = Boolean.parseBoolean(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                eqConfigs = Integer.parseInt(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                writeAcceptanceRatios = Boolean.parseBoolean(scanner.nextLine().split(" " + " +")[1]);
                currLine++;
                //Skip 3 lines for comments
                scanner.nextLine();
                scanner.nextLine();
                scanner.nextLine();
                currLine += 2;
                //Read molecules to place
                while (scanner.hasNextLine()) {
                    String line = scanner.nextLine();
                    currLine++;
                    String name = line.split(" ")[0];
                    int num = Integer.parseInt(line.split(" ")[1]);
                    for (int x = 0; x < num; x++) {
                        if (!add(name)) {
                            System.out.println("Error on line " + currLine + " in config.txt: Unable to find molecule '" + name + "' in dbase.txt.");
                            System.exit(0);
                        }
                    }
                }
            }
            catch (Exception exc){
                System.out.println("Error on line " + currLine + " in config.txt: File incorrectly formatted.");
                System.exit(0);
            }
        }
        catch (Exception exc){
            System.out.println("Error: File config.txt not found.");
            System.exit(0);
        }
    }
    public void readInput(){
        try{
            //Construct path to file, rather complicated but has to work with .jar or project files for testing
            String pathName = Main.class.getProtectionDomain().getCodeSource().getLocation().getPath();
            pathName = URLDecoder.decode(pathName, "utf-8");
            pathName = "/" + pathName.substring(1, pathName.lastIndexOf("/")) + "/Input.xyz";
            Scanner scanner = new Scanner(new File(pathName));
            int currLine = 1; //Tracks current line for error printing
            //Skip next line, since we don't care about the total number of atoms, we get it all from the comment on the next line
            scanner.nextLine();
            try {
                String[] atomNums = scanner.nextLine().split(" ");
                currLine++;
                for (String num : atomNums) {
                    int atomNum;
                    boolean unmoving = false;
                    try {
                        atomNum = Integer.parseInt(num);
                    } catch (NumberFormatException exc) {
                        atomNum = Integer.parseInt(num.substring(0, num.length() - 1));
                        unmoving = true;
                    }
                    //Four arrays used for temporary storage, never get too large since they're overwritten for each molecule
                    String[] atomSyms = new String[atomNum];
                    double[] atomXs = new double[atomNum];
                    double[] atomYs = new double[atomNum];
                    double[] atomZs = new double[atomNum];
                    for (int y = 0; y < atomNum; y++) {
                        String[] atomVals = scanner.nextLine().split(" +"); //splits by all sets of spaces > 1
                        currLine++;
                        //If line is preceded by space, everything breaks since we get atomVals[0] == "". So, if that's true, start from 1.
                        int c = 0;
                        if (atomVals[0].isEmpty()) {
                            c = 1;
                        }
                        atomSyms[y] = atomVals[c];
                        atomXs[y] = Double.parseDouble(atomVals[c + 1]);
                        atomYs[y] = Double.parseDouble(atomVals[c + 2]);
                        atomZs[y] = Double.parseDouble(atomVals[c + 3]);
                    }
                    //Set up hashmap to compare to atomSymbols in Molecule
                    HashMap<String, Integer> matcher = new HashMap<>();
                    for (String s : atomSyms) {
                        if (matcher.containsKey(s)) {
                            matcher.replace(s, matcher.get(s) + 1);
                        } else {
                            matcher.put(s, 1);
                        }
                    }
                    boolean dbContainsMolecule = false;
                    for (Molecule molecule : dbase) {
                        if (molecule.atomSymbols.equals(matcher)) {
                            //Copy molecule from dbase
                            dbContainsMolecule = true;
                            Molecule m = new Molecule(molecule);
                            //Move each atom in m to positions defined in Input.xyz
                            ArrayList<Atom> oldAtoms = m.atoms;
                            ArrayList<Atom> newAtoms = new ArrayList<>();
                            for (int z = 0; z < atomSyms.length; z++) {
                                for (Atom atom : oldAtoms) {
                                    if (atom.symbol.equals(atomSyms[z])) {
                                        atom.x = atomXs[z];
                                        atom.y = atomYs[z];
                                        atom.z = atomZs[z];
                                        oldAtoms.remove(atom);
                                        newAtoms.add(atom);
                                        break;
                                    }
                                }
                            }
                            m.atoms = newAtoms;
                            m.setCenterOfMass();
                            space.add(m);
                            if (!unmoving) {
                                moveableMolecules.add(m);
                            }
                            break;
                        }
                    }
                    if (!dbContainsMolecule) {
                        System.out.println("Error: Unable to find one or more molecules from Input.xyz in dbase.txt.");
                        System.exit(0);
                    }
                }
            }
            catch (Exception exc){
                System.out.println("Error on line" + currLine + " in Input.xyz: File incorrectly formatted.");
                System.exit(0);
            }
        }
        catch (Exception exc){
            System.out.println("Error: File Input.xyz not found.");
            System.exit(0);
        }
    }

    //Set up pair
    public void setupPairVals(){

        //Compile unique atoms from molecules in dbase
        ArrayList<Atom> allAtoms = new ArrayList<>();
        for (Molecule molecule : dbase) {
            allAtoms.addAll(molecule.atoms);
        }
        for (Atom atom1 : allAtoms) {
            for (Atom atom2 : allAtoms) {
                double A = Math.sqrt(atom1.a * atom2.a);
                double B = (atom1.b + atom2.b) / 2;
                double C = Math.sqrt(atom1.c * atom2.c);
                double D = Math.sqrt(atom1.d * atom2.d);
                Pair<UUID, UUID> key1 = new Pair<>(atom1.uuid, atom2.uuid);
                Pair<UUID, UUID> key2 = new Pair<>(atom2.uuid, atom1.uuid);
                double[] value = {A, B, C, D};
                pairwiseDbase.put(key1, value);
                pairwiseDbase.put(key2, value);
            }
        }
    }
    //Creates new directory to place output files into
    public void makeDirectory(){
        //Get current datetime and create directory name
        LocalDateTime time = LocalDateTime.now();
        String name = time.getYear() + "_" + time.getMonthValue() + "_" + time.getDayOfMonth() + "_" + time.getHour() + "_" + time.getMinute() + "_" + time.getSecond();
        //Get path of object
        try {
            String dirPath = Main.class.getProtectionDomain().getCodeSource().getLocation().getPath();
            dirPath = URLDecoder.decode(dirPath, "utf-8");
            dir = "/" + dirPath.substring(1, dirPath.lastIndexOf("/")) + "/" + name;
            File f = new File(dir);
            if (!f.mkdir()){
                throw new Exception();
            }
            //Copy config.txt to new directory
            String cfgPath = Main.class.getProtectionDomain().getCodeSource().getLocation().getPath();
            cfgPath = URLDecoder.decode(cfgPath, "utf-8");
            cfgPath = "/" + cfgPath.substring(1, cfgPath.lastIndexOf("/")) + "/config.txt";
            Scanner scanner = new Scanner(new File(cfgPath));
            StringBuilder out = new StringBuilder();
            while (scanner.hasNextLine()){
                out.append(scanner.nextLine()).append("\n");
            }
            String copyPath = dir + "/config.txt";
            FileWriter writer = new FileWriter(copyPath);
            writer.write(out.toString());
            writer.close();
            if (useInput) { //If enabled, copy Input.xyz to new directory
                String inputPath = Main.class.getProtectionDomain().getCodeSource().getLocation().getPath();
                inputPath = URLDecoder.decode(inputPath, "utf-8");
                inputPath = "/" + inputPath.substring(1, inputPath.lastIndexOf("/")) + "/Input.xyz";
                scanner = new Scanner(new File(inputPath));
                out = new StringBuilder();
                while (scanner.hasNextLine()) {
                    out.append(scanner.nextLine()).append("\n");
                }
                copyPath = dir + "/Input.xyz";
                writer = new FileWriter(copyPath);
                writer.write(out.toString());
                writer.close();
            }
        }
        catch (Exception exc){
            System.out.println("Error: Failed to create new directory.");
            System.exit(0);
        }
    }
    //Writes atom placements to .xyz file. Programs that read .xyz files will figure out what atoms go to what molecules, so that information is unnecessary
    public void write(int outputFileNumber){
        try{
        	String pathName = dir + "/Output" + outputFileNumber + ".xyz"; //dir specified in makeDirectory()
            FileWriter writer = new FileWriter(pathName);
            double toothEnergy = calcEnergy();
            //Write to file in correct .xyz output format
            StringBuilder content = new StringBuilder("          " + numAtoms() + "\nEnergy: " + toothEnergy + " Kcal/mole");
            for (Molecule m : space){
                for (Atom a : m.atoms){
                    if (a.symbol.contains("*")){
                        continue;
                    }
                    content.append("\n ").append(a.symbol);
                    if (a.symbol.length() == 1) {
                    	content.append(" ");
                    }
                    //Round all doubles to 10 decimal places to keep file neat
                    String[] x = Double.toString(Precision.round(a.x, 10, BigDecimal.ROUND_HALF_UP)).split("\\.");
                    content.append("  ");
                    for (int b = 0; b + x[0].length() < 4; b++) {
                		content.append(" ");
                	}
                    content.append(x[0]).append(".").append(x[1]);
                	for (int b = 0; b + x[1].length() < 10; b++) {
                		content.append("0");
                	}

                	String[] y = Double.toString(Precision.round(a.y, 10, BigDecimal.ROUND_HALF_UP)).split("\\.");
                    content.append("  ");
                    for (int b = 0; b + y[0].length() < 4; b++) {
                		content.append(" ");
                	}
                    content.append(y[0]).append(".").append(y[1]);
                	for (int b = 0; b + y[1].length() < 10; b++) {
                		content.append("0");
                	}

                	String[] z = Double.toString(Precision.round(a.z, 10, BigDecimal.ROUND_HALF_UP)).split("\\.");
                    content.append("  ");
                    for (int b = 0; b + z[0].length() < 4; b++) {
                		content.append(" ");
                	}
                    content.append(z[0]).append(".").append(z[1]);
                	for (int b = 0; b + z[1].length() < 10; b++) {
                		content.append("0");
                	}
                }
            }
            writer.write(content + "\n");
            if (toothEnergy < minEnergy){
                minEnergy = toothEnergy;
                minEnergyToothNum = outputFileNumber;
            }
            writer.close();
        }
        catch(Exception exc){
            System.out.println("Error: Failed to write to " + dir + "/Output" + outputFileNumber + ".xyz");
            System.exit(0);
        }
    }
    //Return number of atoms in space, used for write() as part of .xyz file format
    private int numAtoms(){
        int ret = 0;
        for (Molecule molecule : space) {
            for (int y = 0; y < molecule.atoms.size(); y++) {
                if (!molecule.atoms.get(y).symbol.contains("*")) { //Exclude ghosts from count
                    ret++;
                }
            }
        }
        return ret;
    }
    //Prints database List as string, for debugging only
    public String printDbase() {
    	StringBuilder ret = new StringBuilder();
    	for (int x = 0; x < dbase.size(); x++) {
    		if ((x != 0) && (x % 3 == 0)) {
    			ret.append("\n");
    		}
    		ret.append(dbase.get(x).name).append("\t");
    	}
    	return ret.toString();
    }
    //Calculates total energy of space, used for console output at end of propagation and each sawtooth
    public double calcEnergy() {
    	double totEnergy = 0;
    	for (int x = 0; x < space.size(); x++) {
    		for (int y = 0; y < space.size(); y++) {
    			if (y > x) {
    				totEnergy += space.get(x).calcEnergy(space.get(y), pairwiseDbase);
    			}
    		}
    	}
    	return totEnergy;
    }
    //Chooses whether to rotate molecule or not based on change in energy
    public Pair<Double, Integer> rotate(Molecule m, double maxRot, double temp) {
        m.rotateTemp(maxRot); //Temporarily rotate molecule
        for (Atom a : m.atoms){
            if (a.tempx > size * 0.75 || a.tempx < size * -0.75 || a.tempy > size * 0.75 || a.tempy < size * -0.75 ||a.tempz > size * 0.75 || a.tempz < size * -0.75){
                m.resetTemps();
                return new Pair<>(0.0, 0);
            }
        }
        double eStart = 0;
        double eEnd = 0;
        for (Molecule n : space) { //Calculate energy based on both current position and temporary position
            if (m != n) {
                eStart += m.calcEnergy(n, pairwiseDbase);
                eEnd += m.calcTempEnergy(n, pairwiseDbase);
            }
        }
        if (eEnd - eStart > 0) { //If energy increases, certain chance to accept move based on boltzmann constant and temperature
            double rand = r.nextDouble();
            if (temp == 0) { //If temperature is zero, formula fails but only decreases in energy should be accepted, so reject manually
                m.resetTemps();
                return new Pair<>(0.0, 0);
            }
            double test = Math.pow(Math.E, -1 * (eEnd - eStart)/(BOLTZMANN_CONSTANT * temp));
            if (rand > test) {
                m.resetTemps();
                return new Pair<>(0.0, 0);
            }
        }
        m.setTemps();
        return new Pair<>(eEnd - eStart, 1);
    }
    //Chooses whether to move molecule or not based on change in energy
    public Pair<Double, Integer> move(Molecule m, double maxD, double temp) {
    	//Temporarily move molecule
        double x = (r.nextDouble() * 2 * maxD) - maxD;
    	double y = (r.nextDouble() * 2 * maxD) - maxD;
    	double z = (r.nextDouble() * 2 * maxD) - maxD;
    	m.tempx += x;
    	m.tempy += y;
    	m.tempz += z;
    	for (Atom a : m.atoms) {
    		a.tempx += x;
    		a.tempy += y;
    		a.tempz += z;
    	}
    	for (Atom a : m.atoms){
    	    if (a.tempx > size * 0.75 || a.tempx < size * -0.75 || a.tempy > size * 0.75 || a.tempy < size * -0.75 ||a.tempz > size * 0.75 || a.tempz < size * -0.75){
                m.resetTemps();
                return new Pair<>(0.0, 0);
            }
        }
        double eStart = 0;
    	double eEnd = 0;
		for (Molecule n : space) {
			if (m != n) {
				eStart += m.calcEnergy(n, pairwiseDbase);
				eEnd += m.calcTempEnergy(n, pairwiseDbase);
			}
		}
		//double delEnergy = eEnd - eStart;
		if (eEnd - eStart > 0) { //If energy increases, certain chance to accept move based on boltzmann constant and temperature
			double rand = r.nextDouble();
			if (temp == 0) { //If temperature is zero, formula fails but only decreases in energy should be accepted, so reject manually
				m.resetTemps();
                return new Pair<>(0.0, 0);
			}
			double test = Math.pow(Math.E, -1 * (eEnd - eStart)/(BOLTZMANN_CONSTANT * temp));
			if (rand > test) {
				m.resetTemps();
                return new Pair<>(0.0, 0);
			}
		}
		m.setTemps();
		return new Pair<>(eEnd - eStart, 1);
    }
    public Molecule randMolecule() {
        int x = r.nextInt(moveableMolecules.size());
        return moveableMolecules.get(x);
    }
    //Writes by appending to file in order to create a movie .xyz file that can be viewed in Avogadro
    public void writeMovie(int outputStartNumber, int outputEndNumber){
        try{
            String pathName = dir + "/Output" + outputStartNumber + "_" + outputEndNumber + "_Movie.xyz"; //dir specified in makeDirectory()
            FileWriter writer = new FileWriter(pathName, true);
            //Write to file in correct .xyz output format
            StringBuilder content = new StringBuilder("          " + numAtoms() + "\nEnergy: " + calcEnergy() + " Kcal/mole");
            for (Molecule m : space){
                for (Atom a : m.atoms){
                    if (a.symbol.contains("*")){
                        continue;
                    }
                    content.append("\n ").append(a.symbol);
                    if (a.symbol.length() == 1) {
                    	content.append(" ");
                    }
                    //Round all doubles to 10 decimal places to keep file neat
                    String[] x = Double.toString(Precision.round(a.x, 10, BigDecimal.ROUND_HALF_UP)).split("\\.");
                    content.append("  ");
                    for (int b = 0; b + x[0].length() < 4; b++) {
                		content.append(" ");
                	}
                    content.append(x[0]).append(".").append(x[1]);
                	for (int b = 0; b + x[1].length() < 10; b++) {
                		content.append("0");
                	}

                	String[] y = Double.toString(Precision.round(a.y, 10, BigDecimal.ROUND_HALF_UP)).split("\\.");
                    content.append("  ");
                    for (int b = 0; b + y[0].length() < 4; b++) {
                		content.append(" ");
                	}
                    content.append(y[0]).append(".").append(y[1]);
                	for (int b = 0; b + y[1].length() < 10; b++) {
                		content.append("0");
                	}

                	String[] z = Double.toString(Precision.round(a.z, 10, BigDecimal.ROUND_HALF_UP)).split("\\.");
                    content.append("  ");
                    for (int b = 0; b + z[0].length() < 4; b++) {
                		content.append(" ");
                	}
                    content.append(z[0]).append(".").append(z[1]);
                	for (int b = 0; b + z[1].length() < 10; b++) {
                		content.append("0");
                	}
                }
            }
            writer.write(content + "\n");
            writer.close();
        }
        catch(Exception exc){
            System.out.println("Error: Failed to write to " + dir + "/Output" + outputStartNumber + "_" + outputEndNumber + "_Movie.xyz");
            System.exit(0);
        }
    }
    //Writes output and logs to appropriate file in dir, then prints to terminal
    public void log(String text){
        try {
            String pathName = dir + "/log.txt";
            FileWriter writer = new FileWriter(pathName, true);
            writer.write(text + "\n");
            writer.close();
            System.out.println(text);
        }
        catch (Exception exc){
            System.out.println("Error: Failed to write to " + dir + "/log.txt.");
            System.exit(0);
        }
    }
    //Appends energy values to energies.txt in dir if staticTemp = true
    public void writeEnergy(double energy){
        if (staticTemp){
            try {
                String pathName = dir + "/energies.txt";
                FileWriter writer = new FileWriter(pathName, true);
                writer.write(energy + "\n");
                writer.close();
            }
            catch (Exception exc){
                System.out.println("Error: Failed to write to " + dir + "/energies.txt.");
                System.exit(0);
            }
        }
    }
    public void writeAcceptance(double temperature, int accepted, int total){
        if(writeAcceptanceRatios){ //TODO: fix this
            try{
                double ratio = (double) accepted / (double) total;
                String pathName = dir + "/acceptance_ratios.txt";
                FileWriter writer = new FileWriter(pathName, true);
                writer.write(new BigDecimal(temperature).setScale(2, RoundingMode.HALF_UP) + " " + new BigDecimal(ratio).setScale(5, RoundingMode.HALF_UP) + "\n");
                writer.close();
            }
            catch (Exception exc){
                System.out.println("Error: Failed to write to " + dir + "/acceptance_ratios.txt");
                System.exit(0);
            }
        }
    }
    //Writes min energy structure to min_structure.xyz
    public void writeMinEnergy(){
        String writePath = dir + "/Min_Energy_Structure_" + minEnergyToothNum + ".xyz";
        try {
            Files.copy(new File(dir + "/Output" + minEnergyToothNum + ".xyz").toPath(), new File(writePath).toPath());
        }
        catch (Exception exc){
            System.out.println("Error: Failed to write to " + dir + "/min_energy_structure" + minEnergyToothNum + ".xyz");
            System.exit(0);
        }
    }

    public void writeConfigurationalHeatCapacity(double temperature, double avgEnergy, double configurationalHeatCapacity){
        String writePath = dir + "/configurational_heat_capacities.txt";
        String output = temperature + "\t" + avgEnergy + "\t" + configurationalHeatCapacity + "\n";

        //Write configurational heat capacity
        try {
            FileWriter writer = new FileWriter(writePath, true);
            writer.write(output);
            writer.close();
        }
        catch (Exception exc){
            System.out.println("Error: Failed to write to " + dir + "/configurational_heat_capacities.txt");
            System.exit(1);
        }
    }
}
