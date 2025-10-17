package main;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class Config {

    public static double maxTemp;
    public static int movePerPt;
    public static int ptsPerTooth;
    public static int ptsAddedPerTooth;
    public static int numTeeth;
    public static double tempDecreasePerTooth;
    public static double maxTrans;
    public static double magwalkFactorTrans;
    public static double magwalkProbTrans;
    public static double maxRot;
    public static double magwalkProbRot;
    public static double spaceLen;
    public static int maxPropFailures;

    public static boolean useInput;
    public static boolean chooseParams;
    public static boolean finale;
    public static boolean staticTemp;
    public static boolean writeEnergiesEnabled;
    public static boolean writeConfHeatCapacitiesEnabled;
    public static int eqConfigs;
    public static boolean writeAcceptanceRatios;

    // VarType uses RegEx to validate the type of an input value
    public enum VarType {
        INT("^\\d+$"), DOUBLE("^\\d+(\\.\\d+)?$"), BOOLEAN("^(true|false)$");

        private final String regex;

        VarType(String re) {
            regex = re;
        }

        public boolean validate(String var) {
            return var.matches(regex);
        }
    }

    private static final List<String> varNames = Arrays.asList(
        "Max Temperature",
        "Moves per Point",
        "Points per Tooth",
        "Points Added per Tooth",
        "Number of Teeth",
        "Temperature Decrease per Tooth",
        "Max Translation Distance",
        "Magwalk Translation Multiplication Factor",
        "Magwalk Translation Probability",
        "Max Rotation",
        "Magwalk Rotation Probability",
        "Length of Cubic Space",
        "Max Failures During Propagation",
        "Use Input.xyz",
        "Choose All Interaction Parameters",
        "0K Finale",
        "Static Temperature",
        "Write Energies When Static Temperature",
        "Write Configurational Heat Capacities",
        "Number of Equilibration Configurations",
        "Write Acceptance Ratios"
    );
    private static final List<VarType> types = Arrays.asList(
        VarType.DOUBLE, // Max Temperature
        VarType.INT, // Moves per Point
        VarType.INT, // Points per Tooth
        VarType.INT, // Points Added per Tooth
        VarType.INT, // Number of Teeth
        VarType.DOUBLE, // Temperature Decrease per Tooth
        VarType.DOUBLE, // Max Translation Distance
        VarType.DOUBLE, // Magwalk Translation Multiplication Factor
        VarType.DOUBLE, // Magwalk Translation Probability
        VarType.DOUBLE, // Max Rotation
        VarType.DOUBLE, // Magwalk Rotation Probability
        VarType.DOUBLE, // Length of Cubic Space
        VarType.INT, // Max Failures During Propagation
        VarType.BOOLEAN, // Use Input.xyz
        VarType.BOOLEAN, // Choose All Interaction Parameters
        VarType.BOOLEAN, // 0K Finale
        VarType.BOOLEAN, // Static Temperature
        VarType.BOOLEAN, // Write Energies When Static Temperature
        VarType.BOOLEAN, // Write Configurational Heat Capacities
        VarType.INT, // Number of Equilibration Configurations
        VarType.BOOLEAN // Write Acceptance Ratios
    );
    private static final Map<String, VarType> varTypes = varNames.stream().collect(Collectors.toMap(
        varName -> varName,
        varName -> types.get(varNames.indexOf(varName))
    ));

    // Parses a valid configuration .txt file
    public static void parseConfig(String filename, Space s) {
        // Regex is used for line matching for more specific error management
        String settingRegex = "^(?<setting>.+): +(?<value>.+)$";
        Pattern settingPattern = Pattern.compile(settingRegex);
        String particleRegex = "^(?<particle>(\\S+\\s?)+) {2,}(?<count>\\d+)$";
        Pattern particlePattern = Pattern.compile(particleRegex);

        // Config settings are stored first, and then propagated into variables
        Map<String, String> configValues = new HashMap<>();
        List<String> addedParticles = new ArrayList<>();

        // A try-with-resources block automatically closes the scanner before exiting the method
        try (Scanner scanner = new Scanner(new File(filename))) {
            int currLine = 1;
            boolean finishedSettings = false;
            // All lines are trimmed, so leading or trailing space is insignificant
            for (String line = scanner.nextLine().trim();; line = scanner.nextLine().trim(), currLine++) {
                Matcher settingMatcher = settingPattern.matcher(line);
                Matcher particleMatcher = particlePattern.matcher(line);

                // Skip commas and empty lines
                if (line.startsWith("//") || line.length() == 0);
                // Config settings must all be present, in any order in the file
                else if (settingMatcher.find()) {
                    if (finishedSettings)
                        throw new RuntimeException("Error on line " + currLine + " in config file: Particle counts must be defined after all configuration settings in configuration file.");

                    String settingName = settingMatcher.group("setting");
                    String value = settingMatcher.group("value");
                    settingName = settingName.split("\\([a-zA-Z0-9/_-]+\\)")[0].trim();

                    // No new settings are allowed, and all settings must be of the specified type
                    if (!varNames.contains(settingName))
                        throw new RuntimeException("Error on line " + currLine + " in config file: Configuration setting '" + settingName + "' not expected.");

                    if (!varTypes.get(settingName).validate(value))
                        throw new RuntimeException("Error on line " + currLine + " in config file: Unexpected type for '" + settingName + "' configuration");

                    configValues.put(settingName, value);
                // Particle counts are directly added to the Space if valid
                } else if (particleMatcher.find()) {
                    finishedSettings = true;
                    String name = particleMatcher.group("particle");
                    int count = Integer.parseInt(particleMatcher.group("count"));
                    if (addedParticles.contains(name))
                        throw new RuntimeException("Error on line " + currLine + " in config file: Particle counts can only be specified once per particle.");

                    if (!s.add(name, count))
                        throw new RuntimeException("Error on line " + currLine + " in config file: Unable to find particle '" + name + "' in dbase.");

                    addedParticles.add(name);
                } else throw new RuntimeException("Error on line " + currLine + " in config file: Unexpected format.");
                if (!scanner.hasNextLine()) break;
            }
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Error: File " + filename + " not found.");
        }

        // Check all values are defined
        Optional<String> unknownKey = configValues.keySet().stream().filter(key -> !varNames.contains(key)).findFirst();
        if (unknownKey.isPresent())
            throw new RuntimeException("Error in config file: Unexpected setting '" + unknownKey + "'.");

        if (configValues.size() != varNames.size())
            throw new RuntimeException("Error in config file: Missing setting '" + varNames.stream().filter(name -> !configValues.containsKey(name)).findFirst().get() + "'.");

        setVars(configValues);
    }

    // Sets all static configuration variables for ease of use in other classes
    private static void setVars(Map<String, String> values) {
        maxTemp = Double.parseDouble(values.get("Max Temperature"));
        movePerPt = Integer.parseInt(values.get("Moves per Point"));
        ptsPerTooth = Integer.parseInt(values.get("Points per Tooth"));
        ptsAddedPerTooth = Integer.parseInt(values.get("Points Added per Tooth"));
        numTeeth = Integer.parseInt(values.get("Number of Teeth"));
        tempDecreasePerTooth = Double.parseDouble(values.get("Temperature Decrease per Tooth"));
        maxTrans = Double.parseDouble(values.get("Max Translation Distance"));
        magwalkFactorTrans = Double.parseDouble(values.get("Magwalk Translation Multiplication Factor"));
        magwalkProbTrans = Double.parseDouble(values.get("Magwalk Translation Probability"));
        maxRot = Double.parseDouble(values.get("Max Rotation"));
        magwalkProbRot = Double.parseDouble(values.get("Magwalk Rotation Probability"));
        spaceLen = Double.parseDouble(values.get("Length of Cubic Space"));
        maxPropFailures = Integer.parseInt(values.get("Max Failures During Propagation"));
        useInput = Boolean.parseBoolean(values.get("Use Input.xyz"));
        chooseParams = Boolean.parseBoolean(values.get("Choose All Interaction Parameters"));
        finale = Boolean.parseBoolean(values.get("0K Finale"));
        staticTemp = Boolean.parseBoolean(values.get("Static Temperature"));
        writeEnergiesEnabled = Boolean.parseBoolean(values.get("Write Energies When Static Temperature"));
        writeConfHeatCapacitiesEnabled = Boolean.parseBoolean(values.get("Write Configurational Heat Capacities"));
        eqConfigs = Integer.parseInt(values.get("Number of Equilibration Configurations"));
        writeAcceptanceRatios = Boolean.parseBoolean(values.get("Write Acceptance Ratios"));
    }
}
