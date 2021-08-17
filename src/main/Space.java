package main;
import java.io.*;
import java.math.BigDecimal;
import java.net.URLDecoder;
import java.util.ArrayList;
import java.util.Scanner;
import org.apache.commons.math3.random.MersenneTwister;
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

    private String dir; //Directory of .xyz files to be saved in, created at runtime

    private ArrayList<Molecule> list;
    private ArrayList<Molecule> space;
    private ArrayList<Molecule> dbase;
    private MersenneTwister r = new MersenneTwister(); //Used instead of Java.Random for greater accuracy in randomization
    private static final double BOLTZMANN_CONSTANT = 0.0019872;
    public Space(double s){
        size = s;
        list = new ArrayList<>();
        space = new ArrayList<>();
        dbase = new ArrayList<>();
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
                double x = r.nextDouble() * size;
                double y = r.nextDouble() * size;
                double z = r.nextDouble() * size;
                if (space.size() == 0){ //If placing first molecule, no reason to do comparisons, just place the molecule
                    m.put(x, y, z);
                    space.add(m);
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
                    break;
                }
                if (c >= maxPropFailures){ //Need larger space to place molecules
                	space.clear();
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
            while (scanner.hasNextLine()){
                //Parse file based on expected file format; if an issue is encountered, assume the file is improperly set up and throw an error
                try {
                    String[] words = scanner.nextLine().split("\t");
                    if (words.length == 1) {
                        int n;
                        try {
                            n = Integer.parseInt(words[0]);
                        } catch (Exception exc) {
                            continue;
                        }
                        ArrayList<Atom> atoms = new ArrayList<>();
                        if (!scanner.hasNextLine()) {
                            throw new Exception();
                        }
                        String[] s = scanner.nextLine().split("\t");
                        double radius = 0;
                        if (s.length != 2) {
                            throw new Exception();
                        }
                        radius = Double.parseDouble(s[1]);
                        String name = s[0];
                        for (int x = 0; x < n; x++) {
                            String[] atom = scanner.nextLine().split("\t");
                            if (atom.length != 10) {
                                throw new Exception();
                            }
                            Atom a = new Atom(atom[0], Double.parseDouble(atom[1]), Double.parseDouble(atom[2]), Double.parseDouble(atom[3]), Double.parseDouble(atom[4]), Double.parseDouble(atom[5]), Double.parseDouble(atom[6]), Double.parseDouble(atom[7]), Double.parseDouble(atom[8]), Double.parseDouble(atom[9]));
                            atoms.add(a);
                        }
                        Molecule m = new Molecule(name, radius, atoms);
                        dbase.add(m);
                    }
                }
                catch (Exception exc){
                    System.out.println("Error: File dbase.txt incorrectly formatted.");
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
            //Skip 2 lines for comments
            scanner.nextLine();
            scanner.nextLine();
            try{
                //Read configs
                maxTemperature = Double.parseDouble(scanner.nextLine().split("\t")[1]);
                movePerPoint = Integer.parseInt(scanner.nextLine().split("\t")[1]);
                pointsPerTooth = Integer.parseInt(scanner.nextLine().split("\t")[1]);
                pointIncrement = Integer.parseInt(scanner.nextLine().split("\t")[1]);
                numTeeth = Integer.parseInt(scanner.nextLine().split("\t")[1]);
                tempDecreasePerTooth = Double.parseDouble(scanner.nextLine().split("\t")[1]);
                maxTransDist = Double.parseDouble(scanner.nextLine().split("\t")[1]);
                magwalkFactorTrans = Double.parseDouble(scanner.nextLine().split("\t")[1]);
                magwalkProbTrans = Double.parseDouble(scanner.nextLine().split("\t")[1]);
                maxRotDegree = Double.parseDouble(scanner.nextLine().split("\t")[1]);
                magwalkProbRot = Double.parseDouble(scanner.nextLine().split("\t")[1]);
                size = Double.parseDouble(scanner.nextLine().split("\t")[1]);
                maxPropFailures = Integer.parseInt(scanner.nextLine().split("\t")[1]);
                //Skip 3 lines for comments
                scanner.nextLine();
                scanner.nextLine();
                scanner.nextLine();
                //Read molecules to place
                while (scanner.hasNextLine()){
                    String line = scanner.nextLine();
                    String name = line.split(" ")[0];
                    int num = Integer.parseInt(line.split(" ")[1]);
                    for (int x = 0; x < num; x++){
                        if (!add(name)){
                            System.out.println("Error: unable to find molecule '" + name + "' in database.");
                            System.exit(0);
                        }
                    }
                }
            }
            catch (Exception exc){
                System.out.println("Error: File config.txt incorrectly formatted.");
                System.exit(0);
            }
        }
        catch (Exception exc){
            System.out.println("Error: File  " + dir + "/config.txt not found.");
            System.exit(0);
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
            String cfgPath = Main.class.getProtectionDomain().getCodeSource().getLocation().getPath();
            cfgPath = URLDecoder.decode(cfgPath, "utf-8");
            cfgPath = "/" + cfgPath.substring(1, cfgPath.lastIndexOf("/")) + "/config.txt";
            Scanner scanner = new Scanner(new File(cfgPath));
            String out = "";
            while (scanner.hasNextLine()){
                out += scanner.nextLine() + "\n";
            }
            String copyPath = dir + "/config.txt";
            FileWriter writer = new FileWriter(new File(copyPath));
            writer.write(out);
            writer.close();
        }
        catch (Exception exc){
            System.out.println("Error: Failed to create directory.");
            System.exit(0);
        }
    }
    //Writes atom placements to .xyz file. Programs that read .xyz files will figure out what atoms go to what molecules, so that information is unnecessary
    public String write(String fileName){
        try{
        	String pathName = dir + "/" + fileName; //dir specified in makeDirectory()
            FileWriter writer = new FileWriter(new File(pathName));
            //Write to file in correct .xyz output format
            String content = "          " + numAtoms() + "\nEnergy: " + calcEnergy() + " Kcal/mole";
            for (Molecule m : space){
                for (Atom a : m.atoms){
                    content += "\n " + a.symbol;
                    if (a.symbol.length() == 1) {
                    	content += " ";
                    }
                    //Round all doubles to 10 decimal places to keep file neat
                    String[] x = Double.toString(Precision.round(a.x, 10, BigDecimal.ROUND_HALF_UP)).split("\\.");
                    content += "  ";
                    for (int b = 0; b + x[0].length() < 4; b++) {
                		content += " ";
                	}
                    content += x[0] + "." + x[1];
                	for (int b = 0; b + x[1].length() < 10; b++) {
                		content += "0";
                	}
                    
                	String[] y = Double.toString(Precision.round(a.y, 10, BigDecimal.ROUND_HALF_UP)).split("\\.");
                    content += "  ";
                    for (int b = 0; b + y[0].length() < 4; b++) {
                		content += " ";
                	}
                    content += y[0] + "." + y[1];
                	for (int b = 0; b + y[1].length() < 10; b++) {
                		content += "0";
                	}
                	
                	String[] z = Double.toString(Precision.round(a.z, 10, BigDecimal.ROUND_HALF_UP)).split("\\.");
                    content += "  ";
                    for (int b = 0; b + z[0].length() < 4; b++) {
                		content += " ";
                	}
                    content += z[0] + "." + z[1];
                	for (int b = 0; b + z[1].length() < 10; b++) {
                		content += "0";
                	}
                }
            }
            writer.write(content + "\n");
            writer.close();
            return pathName;
        }
        catch(Exception exc){
            return "Error";
        }
    }
    //Return number of atoms in space, used for write() as part of .xyz file format
    private int numAtoms(){
        int ret = 0;
        for (int x = 0; x < space.size(); x++){
            for (int y = 0; y < space.get(x).atoms.size(); y++){
                ret++;
            }
        }
        return ret;
    }
    //Prints database List as string, for debugging only
    public String printDbase() {
    	String ret = "";
    	for (int x = 0; x < dbase.size(); x++) {
    		if ((x != 0) && (x % 3 == 0)) {
    			ret += "\n";
    		}
    		ret += dbase.get(x).name + "\t";
    	}
    	return ret;
    }
    //Calculates total energy of space, used for console output at end of propagation and each sawtooth
    public double calcEnergy() {
    	double totEnergy = 0;
    	for (int x = 0; x < space.size(); x++) {
    		for (int y = 0; y < space.size(); y++) {
    			if (y > x) {
    				totEnergy += space.get(x).calcEnergy(space.get(y));
    			}
    		}
    	}
    	return totEnergy;
    }
    //Chooses whether to rotate molecule or not based on change in energy
    public boolean rotate(Molecule m, double maxRot, double temp) {
        m.rotateTemp(maxRot); //Temporarily rotate molecule
        for (Atom a : m.atoms){
            if (a.tempx > size * 1.5 || a.tempx < 0 || a.tempy > size * 1.5 || a.tempy < 0 ||a.tempz > size * 1.5 || a.tempz < 0){
                m.resetTemps();
                return false;
            }
        }
        double eStart = 0;
        double eEnd = 0;
        for (Molecule n : space) { //Calculate energy based on both current position and temporary position
            if (m != n) {
                eStart += m.calcEnergy(n);
                eEnd += m.calcTempEnergy(n);
            }
        }
        if (eEnd - eStart > 0) { //If energy increases, certain chance to accept move based on boltzmann constant and temperature
            double rand = r.nextDouble();
            if (temp == 0) { //If temperature is zero, formula fails but only decreases in energy should be accepted, so reject manually
                m.resetTemps();
                return false;
            }
            double test = Math.pow(Math.E, -1 * (eEnd - eStart)/(BOLTZMANN_CONSTANT * temp));
            if (rand > test) {
                m.resetTemps();
                return false;
            }
        }
        m.setTemps();
        return true;
    }
    //Chooses whether to move molecule or not based on change in energy
    public boolean move(Molecule m, double maxD, double temp) {
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
    	    if (a.tempx > size * 1.5 || a.tempx < 0 || a.tempy > size * 1.5 || a.tempy < 0 ||a.tempz > size * 1.5 || a.tempz < 0){
                m.resetTemps();
                return false;
            }
        }
    	/**if (m.tempx > size * 1.5 || m.tempy > size * 1.5 || m.tempz > size * 1.5 || m.tempx < 0 || m.tempy < 0 || m.tempz < 0){ //If molecule would be moved out of valid space, reject the move
            m.resetTemps();
            return false;
        }**/
    	double eStart = 0;
    	double eEnd = 0;
		for (Molecule n : space) {
			if (m != n) {
				eStart += m.calcEnergy(n);
				eEnd += m.calcTempEnergy(n);
			}
		}
		//double delEnergy = eEnd - eStart;
		if (eEnd - eStart > 0) { //If energy incrases, certain chance to accept move based on boltzmann constant and temperature
			double rand = r.nextDouble();
			if (temp == 0) { //If temperature is zero, formula fails but only decreases in energy should be accepted, so reject manually
				m.resetTemps();
				return false;
			}
			double test = Math.pow(Math.E, -1 * (eEnd - eStart)/(BOLTZMANN_CONSTANT * temp));
			if (rand > test) {
				m.resetTemps();
				return false;
			}
		}
		m.setTemps();
		return true;
    }
    public Molecule randMolecule() {
    	int x = r.nextInt(space.size());
    	return space.get(x);
    }
    //Writes by appending to file in order to create a movie .xyz file that can be viewed in Avogadro
    public String writeMovie(String fileName){
        try{
            String pathName = dir + "/" + fileName; //dir specified in makeDirectory()
            FileWriter writer = new FileWriter(new File(pathName), true);
            //Write to file in correct .xyz output format
            String content = "          " + numAtoms() + "\nEnergy: " + calcEnergy() + " Kcal/mole";
            for (Molecule m : space){
                for (Atom a : m.atoms){
                    content += "\n " + a.symbol;
                    if (a.symbol.length() == 1) {
                    	content += " ";
                    }
                    //Round all doubles to 10 decimal places to keep file neat
                    String[] x = Double.toString(Precision.round(a.x, 10, BigDecimal.ROUND_HALF_UP)).split("\\.");
                    content += "  ";
                    for (int b = 0; b + x[0].length() < 4; b++) {
                		content += " ";
                	}
                    content += x[0] + "." + x[1];
                	for (int b = 0; b + x[1].length() < 10; b++) {
                		content += "0";
                	}
                    
                	String[] y = Double.toString(Precision.round(a.y, 10, BigDecimal.ROUND_HALF_UP)).split("\\.");
                    content += "  ";
                    for (int b = 0; b + y[0].length() < 4; b++) {
                		content += " ";
                	}
                    content += y[0] + "." + y[1];
                	for (int b = 0; b + y[1].length() < 10; b++) {
                		content += "0";
                	}
                	
                	String[] z = Double.toString(Precision.round(a.z, 10, BigDecimal.ROUND_HALF_UP)).split("\\.");
                    content += "  ";
                    for (int b = 0; b + z[0].length() < 4; b++) {
                		content += " ";
                	}
                    content += z[0] + "." + z[1];
                	for (int b = 0; b + z[1].length() < 10; b++) {
                		content += "0";
                	}
                }
            }
            writer.write(content + "\n");
            writer.close();
            return pathName;
        }
        catch(Exception exc){
            return "Error";
        }
    }
    //Writes output and logs to appropriate file in dir, then prints to terminal
    public void log(String text){
        try {
            String pathName = dir + "/log.txt";
            FileWriter writer = new FileWriter(new File(pathName), true);
            writer.write(text + "\n");
            writer.close();
            System.out.println(text);
        }
        catch (Exception exc){
            System.out.println("Error in writing to log.txt");
            System.exit(0);
        }
    }
}
