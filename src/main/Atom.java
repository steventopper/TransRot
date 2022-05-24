package main;

import org.apache.commons.math3.util.Pair;
import java.util.HashMap;

public class Atom {
	private static final double K_VALUE = 332.0637133;
    String symbol;
    double x;
    double y;
    double z;
    //Def (default) variables used to save relative position to center of molecule; ensures that atoms are placed properly compared to center of molecule when propogating space
    double defX;
    double defY;
    double defZ;
    //Temp variables used for checking if moves increase or decrease energy
    double tempx;
    double tempy;
    double tempz;
    //a,b,c,d,q constants used for energy calculations
    double a;
    double b;
    double c;
    double d;
    double q;
    //Composite of K_VALUE and q, to remove one floating point operation per energy calculation
    double kTimesQ;

    double mass;
    public Atom(String s, double x, double y, double z, double A, double B, double C, double D, double Q, double mass){
        symbol = s;
        this.x = x;
        this.y = y;
        this.z = z;
        defX = x;
        defY = y;
        defZ = z;
        resetTemps();
        a = A;
        b = B;
        c = C;
        d = D;
        q = Q;
        this.mass = mass;
        this.kTimesQ = K_VALUE * q;
    }
    //For copying molecules from database List object, called by Molecule.Molecule(Molecule m)
    public Atom (Atom atom){
        symbol = atom.symbol;
        x = atom.x;
        y = atom.y;
        z = atom.z;
        defX = atom.defX;
        defY = atom.defY;
        defZ = atom.defZ;
        resetTemps();
        a = atom.a;
        b = atom.b;
        c = atom.c;
        d = atom.d;
        q = atom.q;
        mass = atom.mass;
        kTimesQ = K_VALUE * q;
    }
    //Calculates energy based on atom full position
    public double calcEnergy(Atom atom, boolean isPairwise, HashMap<Pair<String, String>, double[]> pairwiseDbase) { //Calculates energy compared to rest of system, instead of total energy of system
        double A;
        double B;
        double C;
        double D;
        if (isPairwise){
            Pair<String, String> key = new Pair<>(this.symbol, atom.symbol);
            double[] value = pairwiseDbase.get(key);
            if (value == null){
                System.out.println("Error: pairwise.txt does not contain a pairing for " + this.symbol + " and " + atom.symbol + ".");
                System.exit(0);
            }
            A = value[0];
            B = value[1];
            C = value[2];
            D = value[3];
        }
        else {
            A = Math.sqrt(this.a * atom.a);
            B = (this.b + atom.b) / 2;
            C = Math.sqrt(this.c * atom.c);
            D = Math.sqrt(this.d * atom.d);
        }
    	double r = Math.sqrt(Math.pow(atom.x - x, 2) + Math.pow(atom.y - y, 2) + Math.pow(atom.z - z, 2)); //Distance between atoms
    	return (kTimesQ * atom.q / r) + A * Math.exp(-1 * B * r) - (C / Math.pow(r,  6)) + (D / Math.pow(r,  12));
    }
   // Calculates energy based on atom temp position
    public double calcTempEnergy(Atom atom, boolean isPairwise, HashMap<Pair<String, String>, double[]> pairwiseDbase) { //Calculates energy compared to rest of system, instead of total energy of system
        double A;
        double B;
        double C;
        double D;
        double Q;
        if (isPairwise){
            Pair<String, String> key = new Pair<>(this.symbol, atom.symbol);
            double[] value = pairwiseDbase.get(key);
            if (value == null){
                System.out.println("Error: pairwise.txt does not contain a pairing for " + this.symbol + " and " + atom.symbol + ".");
                System.exit(0);
            }
            A = value[0];
            B = value[1];
            C = value[2];
            D = value[3];
        }
        else {
            A = Math.sqrt(this.a * atom.a);
            B = (this.b + atom.b) / 2;
            C = Math.sqrt(this.c * atom.c);
            D = Math.sqrt(this.d * atom.d);
        }
    	double r = Math.sqrt(Math.pow(atom.x - tempx, 2) + Math.pow(atom.y - tempy, 2) + Math.pow(atom.z - tempz, 2)); //Distance between atoms
    	return (kTimesQ * atom.q / r) + A * Math.exp(-1 * B * r) - (C / Math.pow(r,  6)) + (D / Math.pow(r,  12));
    }
    //Sets temps to atom's current position, called by Molecule.resetTemps()
    public void resetTemps() {
    	tempx = x;
    	tempy = y;
    	tempz = z;
    }
    //Moves molecule based on position of temps, called by Molecule.setTemps()
    public void setTemps() {
    	x = tempx;
    	y = tempy;
    	z = tempz;
    }
}
