package main;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.util.Pair;
import org.apache.commons.math3.util.Precision;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.file.Files;
import java.time.LocalDateTime;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

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
    boolean useParams;
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
    public void readDB(String pathName){
    	try {
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
                            throw new IOException();
                        }
                        ArrayList<Atom> atoms = new ArrayList<>();
                        if (!scanner.hasNextLine()) {
                            throw new IOException();
                        }
                        String[] s = scanner.nextLine().split(" " + " +");
                        currLine++;
                        if (s.length != 2) {
                            throw new IOException();
                        }
                        double radius = Double.parseDouble(s[1]);
                        String name = s[0];
                        int numGhosts = 0;
                        for (int x = 0; x < n; x++) {
                            String[] atom = scanner.nextLine().split(" " + " +");
                            currLine++;
                            if (atom.length != 10) {
                                throw new IOException();
                            }
                            Atom a = new Atom(atom[0], Double.parseDouble(atom[1]), Double.parseDouble(atom[2]), Double.parseDouble(atom[3]), Double.parseDouble(atom[4]), Double.parseDouble(atom[5]), Double.parseDouble(atom[6]), Double.parseDouble(atom[7]), Double.parseDouble(atom[8]), Double.parseDouble(atom[9]));
                            atoms.add(a);
                            if (a.symbol.contains("*")){
                                numGhosts++;
                            }
                        }
                        if (numGhosts == n){
                            throw new RuntimeException("Error on line " + currLine + " in dbase.txt: Molecule " + name + " cannot be comprised of only ghost atoms.");
                        }
                        Molecule m = new Molecule(name, radius, atoms);
                        dbase.add(m);
                    }
                }
                catch (IOException exc){
                    scanner.close();
                    throw new RuntimeException("Error on line " + currLine + " in dbase.txt: File incorrectly formatted.");
                }
            }
            scanner.close();
    	}
    	catch (IOException exc) {
    		throw new RuntimeException("Error: File " + pathName + " not found.");
    	}
    }
    //Reads annealing config from config.txt
    public void readCFG(String pathName){
        try {
            //Construct path to file, rather complicated but has to work with .jar or project files for testing
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
                useParams = Boolean.parseBoolean(scanner.nextLine().split(" " + " +")[1]);
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
                    String name = line.split(" {2,}")[0];
                    int num = Integer.parseInt(line.split(" {2,}")[1]);
                    for (int x = 0; x < num; x++) {
                        if (!add(name)) {
                            scanner.close();
                            throw new RuntimeException("Error on line " + currLine + " in config.txt: Unable to find particle '" + name + "' in dbase.txt.");
                        }
                    }
                }
            }
            catch (NumberFormatException exc){
                scanner.close();
                throw new RuntimeException("Error on line " + currLine + " in config.txt: File incorrectly formatted.");
            }
        }
        catch (IOException exc){
            throw new RuntimeException("Error: File " + pathName + " not found.");
        }
    }
    public void readInput(String pathName){
        try (Scanner scanner = new Scanner(new File(pathName))) {
            int currLine = 1; // Tracks current line for error printing
            // Skip next line, since we don't care about the total number of atoms, we get it all from the comment on the next line
            String firstLine = scanner.nextLine().trim();
            int numAtoms = Integer.parseInt(firstLine);
            int atomsRead = 0;

            try {
                String[] molInfo = scanner.nextLine().replaceFirst("Energy: -?\\d*.?\\d* Kcal/mole", "").split("\\|");
                currLine++;
                for (String info : molInfo) {
                    String trimmed = info.trim();
                    String[] parts = trimmed.split(" ");
                    if (parts.length < 2) throw new RuntimeException("Error: Incorrect format on line 2 of input file");
                    String num = parts[0];
                    boolean frozen = num.toUpperCase().endsWith("F");
                    int numMols = Integer.parseInt(num.substring(0, num.length() - (frozen ? 1 : 0)));
                    String molName = String.join(" ", Arrays.stream(parts).collect(Collectors.toList()).subList(1, parts.length));
                    if (!dbase.stream().map(m -> m.name).collect(Collectors.toList()).contains(molName)) throw new RuntimeException("Error: Unable to find one or more molecules from input file in dbase.txt.");
                    Molecule mol = dbase.stream().filter(m -> m.name.equals(molName)).findFirst().get();

                    for (int i = 0; i < numMols; i++) {
                        Molecule molecule = new Molecule(mol);
                        for (int j = 0; j < molecule.atoms.size(); j++) {
                            String[] atomVals = scanner.nextLine().trim().split(" {2,}"); // splits by all sets of spaces > 1
                            currLine++;
                            String symbol = atomVals[0];
                            if (!symbol.equals(molecule.atoms.get(j).symbol)) throw new RuntimeException("Error: Unexpected atom on line " + currLine + " of input file.");

                            double x = Double.parseDouble(atomVals[1]);
                            double y = Double.parseDouble(atomVals[2]);
                            double z = Double.parseDouble(atomVals[3]);
                            atomsRead++;
                            Atom a = molecule.atoms.get(j);
                            a.x = x;
                            a.y = y;
                            a.z = z;
                        }
                        molecule.setCenterOfMass();
                        space.add(molecule);
                        if (!frozen) moveableMolecules.add(molecule);
                    }
                }
                if (scanner.hasNextLine() && !scanner.nextLine().equals("")) throw new RuntimeException("Error: Molecule counts must match # of lines in input file.");
                if (numAtoms != atomsRead) throw new RuntimeException("Error: Atom count at top of input file does not match number of molecules read.");
            }
            catch (NumberFormatException exc){
                scanner.close();
                throw new RuntimeException("Error on line " + currLine + " in input file: File incorrectly formatted.");
            }
        }
        catch (IOException exc){
            throw new RuntimeException("Error: File " + pathName + " not found.");
        }
    }
    public void readParams(String pathName) {
        String aliasRegex = "^(?<particle>.+) *= *(?<aliases>([A-Za-z0-9-_]+, *)*[A-Za-z0-9-_]+)$";
        Pattern aliasPattern = Pattern.compile(aliasRegex);
        String paramRegex = "^(?<firstID>[A-Za-z0-9-_]+) *- *(?<secondID>[A-Za-z0-9-_]+)(?<values>( {2,}\\d+(.\\d+)?){4})$";
        Pattern paramPattern = Pattern.compile(paramRegex);
        List<String> aliasedParticles = new ArrayList<>();
        Map<String, List<Atom>> aliases = new HashMap<>();
        List<Pair<Pair<String, String>, double[]>> paramsList = new ArrayList<>();

        try (Scanner scanner = new Scanner(new File(pathName))) {
            int currLine = 1; // Tracks current line for error printing
            for (String line = scanner.nextLine().trim();; line = scanner.nextLine(), currLine++) {
                // Skip empty lines
                if (line.length() == 0) continue;
                // Parse alias lines (anywhere in the file)
                Matcher aliasMatcher = aliasPattern.matcher(line);
                Matcher paramMatcher = paramPattern.matcher(line);
                if (aliasMatcher.find()) {
                    String particleName = aliasMatcher.group("particle");
                    if (aliasedParticles.contains(particleName)) {
                        scanner.close();
                        throw new RuntimeException("Error on line " + currLine + " in interaction parameters file: Cannot alias the same particle type twice.");
                    }
                    List<String> aliasList = Arrays.stream(aliasMatcher.group("aliases").split(",")).map(String::trim).collect(Collectors.toList());
                    Molecule aliased = null;
                    for (Molecule particle : dbase) {
                        if (particle.name.equals(particleName)) {
                            aliased = particle;
                            break;
                        }
                    }
                    if (aliased == null) {
                        scanner.close();
                        throw new RuntimeException("Error on line " + currLine + " in interaction parameters file: Aliased particle does not exist.");
                    }
                    if (aliased.atoms.size() != aliasList.size()) {
                        scanner.close();
                        throw new RuntimeException("Error on line " + currLine + " in interaction parameters file: Number of atom aliases does not match number of atoms defined in the input database.");
                    }

                    for (int i = 0; i < aliasList.size(); i++) {
                        String alias = aliasList.get(i);
                        aliases.putIfAbsent(alias, new ArrayList<>());
                        aliases.get(alias).add(aliased.atoms.get(i));
                    }
                    aliasedParticles.add(particleName);
                }
                // Parse parameter definition lines
                else if (paramMatcher.find()) {
                    String firstID = paramMatcher.group("firstID");
                    String secondID = paramMatcher.group("secondID");
                    String[] values = paramMatcher.group("values").trim().split(" {2,}");
                    paramsList.add(new Pair<>(new Pair<>(firstID, secondID), Arrays.stream(values).mapToDouble(Double::parseDouble).toArray()));
                }
                // Otherwise, error
                else {
                    scanner.close();
                    throw new RuntimeException("Error on line " + currLine + " in interaction parameters file: File incorrectly formatted.");
                }
                if (!scanner.hasNextLine()) break;
            }
        }
        catch (IOException exc){
            throw new RuntimeException("Error: File " + pathName + " not found.");
        }

        // Error if not all particles are aliased
        if (aliasedParticles.size() != dbase.size()) throw new RuntimeException("Error in interaction parameters file: Must provide aliases for all defined particles in the database.");

        // Error if not all interactions are defined
        if (paramsList.size() < (aliases.size() * (aliases.size() - 1)) / 2) throw new RuntimeException("Error in interaction parameters file: Not all interactions are defined.");

        HashMap<Pair<UUID, UUID>, double[]> tempDbase = new HashMap<>();
        for (Pair<Pair<String, String>, double[]> interParams : paramsList) {
            String firstID = interParams.getFirst().getFirst();
            String secondID = interParams.getFirst().getSecond();
            double[] params = interParams.getSecond();

            // If either ID does not exist as an alias, error
            if (!aliases.containsKey(firstID)) throw new RuntimeException("Error in interaction parameters file: Atom alias '" + firstID + "' not defined in particle alias lists.");
            if (!aliases.containsKey(secondID)) throw new RuntimeException("Error in interaction parameters file: Atom alias '" + secondID + "' not defined in particle alias lists.");

            for (Atom atom : aliases.get(firstID)) {
                for (Atom atom2 : aliases.get(secondID)) {
                    Pair<UUID, UUID> key1 = new Pair<>(atom.uuid, atom2.uuid);
                    Pair<UUID, UUID> key2 = new Pair<>(atom2.uuid, atom.uuid);
                    if (tempDbase.containsKey(key1) || tempDbase.containsKey(key2))
                        throw new RuntimeException("Error in interaction parameters file: Interaction between '" + firstID + "' and '" + secondID + "' cannot be defined twice.");

                    tempDbase.put(key1, params.clone());
                    tempDbase.put(key2, params.clone());
                }
            }
        }

        // Add interaction parameters for self collisions
        for (String alias : aliases.keySet()) {
            for (Atom atom : aliases.get(alias)) {
                for (Atom atom2 : aliases.get(alias)) {
                    double A = Math.sqrt(atom.a * atom2.a);
                    double B = (atom.b + atom2.b) / 2;
                    double C = Math.sqrt(atom.c * atom2.c);
                    double D = Math.sqrt(atom.d * atom2.d);
                    Pair<UUID, UUID> key1 = new Pair<>(atom.uuid, atom2.uuid);
                    Pair<UUID, UUID> key2 = new Pair<>(atom2.uuid, atom.uuid);
                    double[] value = {A, B, C, D};
                    tempDbase.put(key1, value);
                    tempDbase.put(key2, value);
                }
            }
        }

        pairwiseDbase.putAll(tempDbase);
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
    // creates directory path to be saved to later
    public void makeDirectoryName(Map<String, String> parsedArgs)  {
        // Get current datetime and create directory name
        LocalDateTime time = LocalDateTime.now();
        String name = time.getYear() + "_" + time.getMonthValue() + "_" + time.getDayOfMonth() + "_" + time.getHour() + "_" + time.getMinute() + "_" + time.getSecond();
        String dirName = parsedArgs.get("output") + "/" + name;
        // throws an exception if directory cannot be created
        SecurityManager security = System.getSecurityManager();
        if (security != null) security.checkWrite(dirName);
        dir = dirName;
    }
    // Creates new directory to place output files into
    public void makeDirectory(Map<String, String> parsedArgs){
        // Get path of object
        try {
            File f = new File(dir);
            if (!f.mkdir()){
                throw new Exception();
            }

            // Copy config.txt to new directory;
            Scanner scanner = new Scanner(new File(parsedArgs.get("config")));
            StringBuilder out = new StringBuilder();
            while (scanner.hasNextLine()){
                out.append(scanner.nextLine()).append("\n");
            }
            String copyPath = dir + "/config.txt";
            FileWriter writer = new FileWriter(copyPath);
            writer.write(out.toString());
            writer.close();

            if (useInput) { // If enabled, copy Input.xyz to new directory
                scanner = new Scanner(new File(parsedArgs.get("input")));
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
            throw new RuntimeException("Error: Failed to create new directory.");
        }
    }

    // Writes interaction parameters for a TransRot run to interaction_params.txt
    public void writeParams() throws IOException {
        // First, pairwiseDbase is converted into an adjacency list, which associates
        // each atom with a map containing all of its interactions
        Map<UUID, Map<UUID, double[]>> adjList = new HashMap<>();
        for (Map.Entry<Pair<UUID, UUID>, double[]> pairEntry : pairwiseDbase.entrySet()) {
            UUID first = pairEntry.getKey().getFirst();
            UUID second = pairEntry.getKey().getSecond();
            adjList.putIfAbsent(first, new HashMap<>());
            adjList.putIfAbsent(second, new HashMap<>());
            adjList.get(first).put(second, pairEntry.getValue());
            adjList.get(second).put(first, pairEntry.getValue());
        }

        // Atoms are grouped by identical adjacency lists into idLists
        Map<Map<UUID, double[]>, List<UUID>> idLists = new HashMap<>();
        Map<Map<UUID, double[]>, String> symbols = new HashMap<>();
        for (Map.Entry<UUID, Map<UUID, double[]>> interactionList : adjList.entrySet()) {
            // Find atom by UUID
            Atom atomFromID = null;
            Molecule parent = null;
            for (Molecule molecule : dbase) {
                for (Atom atom : molecule.atoms) {
                    if (atom.uuid == interactionList.getKey()) {
                        parent = molecule;
                        atomFromID = atom;
                        break;
                    }
                }
                if (atomFromID != null) break;
            }
            if (atomFromID == null) return;

            boolean givenID = false;
            for (Map.Entry<Map<UUID, double[]>, List<UUID>> paramList : idLists.entrySet()) {
                // Find paramList atom's parent
                Molecule firstParent = null;
                for (Molecule molecule : dbase) {
                    for (Atom atom : molecule.atoms) {
                        if (atom.uuid == paramList.getValue().get(0)) {
                            firstParent = molecule;
                            break;
                        }
                    }
                    if (firstParent != null) break;
                }

                // The atom is equivalent to other atoms if its interaction parameters are identical
                // it has the same atomic symbol, and it is from the same particle
                if (interactionList.getValue().size() == paramList.getKey().size()
                        && interactionList.getValue().keySet().stream().allMatch(id -> Arrays.equals(interactionList.getValue().get(id), paramList.getKey().get(id)))
                        && firstParent == parent
                        && atomFromID.symbol.equals(symbols.get(paramList.getKey()))) {
                    paramList.getValue().add(interactionList.getKey());
                    givenID = true;
                }
            }
            // If no matching atom group is found, a new group is started
            if (!givenID) {
                idLists.put(interactionList.getValue(), new ArrayList<>(Collections.singletonList(interactionList.getKey())));
                symbols.put(interactionList.getValue(), atomFromID.symbol);
            }
        }

        // Groups are extracted from idLists
        List<List<UUID>> groups = new ArrayList<>(idLists.values());

        // Sorting is done both within groups and between them so output aliases
        // are in a logical and readable order
        List<Atom> atoms1D = dbase.stream().flatMap(mol -> mol.atoms.stream()).collect(Collectors.toList());
        Map<UUID, Integer> orderMap = new HashMap<>();
        for (int i = 0; i < atoms1D.size(); i++) {
            orderMap.put(atoms1D.get(i).uuid, i);
        }
        groups.forEach(group -> group.sort(Comparator.comparingInt(orderMap::get)));
        groups.sort(Comparator.comparingInt(list -> orderMap.get(list.get(0))));

        // Write output string
        StringBuilder returnStr = new StringBuilder();

        Map<String, Integer> symbolOccurences = new HashMap<>();
        Map<UUID, String> definedUUIDs = new HashMap<>();
        // To assign aliases to atoms, the order of definition in dbase.txt is used for consistency
        for (Molecule molecule : dbase) {
            List<String> ids = new ArrayList<>();
            for (Atom atom : molecule.atoms) {
                String symbol = "";
                // For each atom, if its alias is already assigned, simply pull out this alias for further use
                // Otherwise, every atom in its group is assigned an alias id, numerically according to the group's
                // collective atomic symbol.
                if (definedUUIDs.containsKey(atom.uuid)) symbol = definedUUIDs.get(atom.uuid);
                else {
                    for (List<UUID> group : groups) {
                        if (group.contains(atom.uuid)) {
                            String key = atom.symbol.replaceAll("\\*", "");
                            symbolOccurences.putIfAbsent(key, 0);
                            symbolOccurences.put(key, symbolOccurences.get(key) + 1);
                            symbol = key + symbolOccurences.get(key);
                            for (UUID uuid : group) {
                                definedUUIDs.put(uuid, symbol);
                            }
                            break;
                        }
                    }
                }
                ids.add(symbol);
            }
            // For each molecule, a line is added to the output to define its aliases
            returnStr.append(molecule.name).append("=").append(String.join(",", ids)).append("\n");
        }
        returnStr.append("\n");

        // A final loop iterates through group pairs in order, being careful to not have repeats
        // For each interaction, a line is added to the output with retrieved symbols and the original interaction parameters
        for (int i = 0; i < groups.size(); i++) {
            List<UUID> group1 = groups.get(i);
            for (int j = i + 1; j < groups.size(); j++) {
                List<UUID> group2 = groups.get(j);
                returnStr
                        .append(definedUUIDs.get(group1.get(0)))
                        .append("-")
                        .append(definedUUIDs.get(group2.get(0)))
                        .append(Arrays.stream(pairwiseDbase.get(new Pair<>(group1.get(0), group2.get(0)))).mapToObj(d -> "  " + d).collect(Collectors.joining("")))
                        .append("\n");
            }
        }

        // Write interaction_params.txt
        String paramsPath = dir + "/interaction_params.txt";
        FileWriter writer = new FileWriter(paramsPath);
        writer.write(returnStr.toString());
        writer.close();
    }

    // Writes atom placements to .xyz file. Programs that read .xyz files will figure out what atoms go to what molecules, so that information is unnecessary
    public void write(int outputFileNumber){
        try{
        	String pathName = dir + "/Output" + outputFileNumber + ".xyz"; //dir specified in makeDirectory()
            FileWriter writer = new FileWriter(pathName);
            double toothEnergy = calcEnergy();

            StringBuilder molString = new StringBuilder();
            String lastMol = "";
            int lastCount = 0;
            for (Molecule m : space) {
                if (!m.name.equals(lastMol)) {
                    if (lastMol.length() > 0) molString.append(lastCount).append(" ").append(lastMol).append(" | ");
                    lastCount = 0;
                    lastMol = m.name;
                }
                lastCount++;
            }
            molString.append(lastCount).append(" ").append(lastMol).append(" | ");
            molString.delete(molString.length() - 3, molString.length());

            //Write to file in correct .xyz output format
            StringBuilder content = new StringBuilder("          " + numAtoms() + "\nEnergy: " + toothEnergy + " Kcal/mole  " + molString);
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
            throw new RuntimeException("Error: Failed to write to " + dir + "/Output" + outputFileNumber + ".xyz");
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
            throw new RuntimeException("Error: Failed to write to " + dir + "/Output" + outputStartNumber + "_" + outputEndNumber + "_Movie.xyz");
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
            throw new RuntimeException("Error: Failed to write to " + dir + "/log.txt.");
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
                throw new RuntimeException("Error: Failed to write to " + dir + "/energies.txt.");
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
                throw new RuntimeException("Error: Failed to write to " + dir + "/acceptance_ratios.txt");
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
            throw new RuntimeException("Error: Failed to write to " + dir + "/min_energy_structure" + minEnergyToothNum + ".xyz");
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
            throw new RuntimeException("Error: Failed to write to " + dir + "/configurational_heat_capacities.txt");
        }
    }

    public String getDir() {
        return dir;
    }

    public void writeExecTime(long procStart) {
        long procEnd = System.nanoTime();
        try {
            String writePath = dir + "/elapsed_time.log";
            FileWriter writer = new FileWriter(writePath, false);
            writer.write((procEnd - procStart) + "");
            writer.close();
        } catch (Exception exc) {
            System.err.println("Error: Failed to write to " + dir + "/elapsed_time.log");
            System.exit(1);
        }
    }
}
